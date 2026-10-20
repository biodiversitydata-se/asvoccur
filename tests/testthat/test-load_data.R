tiny <- function(env = parent.frame()) {
  data_path <- withr::local_tempdir(.local_envir = env)
  make_dataset(data_path, "DS-A", n_events = 2, n_asvs = 3)
  data_path
}

test_that("load_data returns the five elements, each a list indexed by datasetID", {
  loaded <- load_data(tiny())

  expect_named(loaded, c("counts", "asvs", "events", "datasets", "emof"))
  for (x in loaded) expect_named(x, "DS-A")
})

test_that("counts is a taxon x event sparse matrix holding the quantities of occurrence.tsv", {
  counts <- load_data(tiny())$counts[["DS-A"]]

  expect_s4_class(counts, "dgCMatrix")
  expect_equal(rownames(counts), asv_ids(2))
  expect_equal(colnames(counts), c("DS-A:S1", "DS-A:S2"))
  expect_equal(as.matrix(counts), matrix(c(7, 14, 21, 28), 2, dimnames = list(asv_ids(2), c("DS-A:S1", "DS-A:S2"))))
})

test_that("asvs is keyed on taxonID, drops dataset_pid and has NA for blank taxonomy cells", {
  asvs <- load_data(tiny())$asvs[["DS-A"]]

  expect_equal(data.table::key(asvs), "taxonID")
  expect_equal(asvs$taxonID, asv_ids(3))
  expect_false("dataset_pid" %in% names(asvs))
  expect_equal(asvs$DNA_sequence, c("ACGT1", "ACGT2", "ACGT3"))
  expect_equal(asvs$class, rep("Gamma", 3))
  expect_true(all(is.na(asvs$phylum)))
  expect_true(all(is.na(asvs$specificEpithet)))
})

test_that("events is keyed on eventID and drops dataset_pid, datasetName and ipt_resource_id", {
  events <- load_data(tiny())$events[["DS-A"]]

  expect_equal(data.table::key(events), "eventID")
  expect_equal(events$eventID, c("DS-A:S1", "DS-A:S2"))
  expect_false(any(c("dataset_pid", "datasetName", "ipt_resource_id") %in% names(events)))
  expect_equal(events$verbatimLocality, c("Arkona Basin", "Bothnian Bay"))
  expect_equal(events$decimalLatitude, c(55.185, 55.185))
  expect_equal(events$remarks, rep("\"Niskin, 5 m", 2))
})

test_that("datasets has the datasetID taken from the eventIDs and the datasetName", {
  datasets <- load_data(tiny())$datasets[["DS-A"]]

  expect_equal(nrow(datasets), 1)
  expect_equal(datasets$datasetID, "DS-A")
  expect_equal(datasets$datasetName, "Name of DS-A")
})

test_that("emof is wide: one column per 'type (unit)' and one row per event", {
  emof <- load_data(tiny())$emof[["DS-A"]]

  expect_equal(data.table::key(emof), "eventID")
  expect_equal(emof$eventID, c("DS-A:S1", "DS-A:S2"))
  expect_equal(emof[["salinity (psu)"]], c(1.5, 2.5))
  expect_equal(emof[["temperature (C)"]], c(3.5, 4.5))
})

test_that("an empty emof gives a table with only eventID and a warning", {
  data_path <- withr::local_tempdir()
  make_dataset(data_path, "DS-A", n_events = 2, n_asvs = 3, empty_emof = TRUE)

  expect_warning(loaded <- load_data(data_path), "Adding empty emof table")
  expect_equal(names(loaded$emof[["DS-A"]]), "eventID")
  expect_equal(loaded$emof[["DS-A"]]$eventID, c("DS-A:S1", "DS-A:S2"))
})

test_that("several datasets are loaded side by side, named by the ZIP name", {
  data_path <- withr::local_tempdir()
  make_dataset(data_path, "DS-A", n_events = 2, n_asvs = 3)
  make_dataset(data_path, "DS-B", n_events = 3, n_asvs = 3)
  loaded <- load_data(data_path)

  expect_setequal(names(loaded$counts), c("DS-A", "DS-B"))
  expect_equal(ncol(loaded$counts[["DS-B"]]), 3)
  expect_equal(loaded$events[["DS-B"]]$eventID, paste0("DS-B:S", 1:3))
})

test_that("a ZIP with a wrapper folder and macOS metadata loads like a plain one", {
  plain <- withr::local_tempdir()
  wrapped <- withr::local_tempdir()
  make_dataset(plain, "DS-A", n_events = 2, n_asvs = 3)
  make_dataset(wrapped, "DS-A", n_events = 2, n_asvs = 3, wrapper = TRUE)

  expect_equal(load_data(wrapped), load_data(plain))
})

test_that("an unzipped folder loads like the ZIP and warns", {
  zipped <- withr::local_tempdir()
  folder <- withr::local_tempdir()
  make_dataset(zipped, "DS-A", n_events = 2, n_asvs = 3)
  make_dataset(folder, "DS-A", n_events = 2, n_asvs = 3, as = "dir")

  expect_warning(from_folder <- load_data(folder), "unzipped dataset folder")
  expect_equal(from_folder, load_data(zipped))
})

test_that("a ZIP wins over a folder with the same ID, with a message", {
  data_path <- withr::local_tempdir()
  make_dataset(data_path, "DS-A", n_events = 2, n_asvs = 3)
  other <- withr::local_tempdir()
  make_dataset(other, "DS-A", as = "dir", n_events = 5, n_asvs = 3)
  file.copy(file.path(other, "DS-A"), data_path, recursive = TRUE)

  expect_message(loaded <- load_data(data_path), "Using ZIP")
  expect_equal(nrow(loaded$events[["DS-A"]]), 2)
})

test_that("load_data stops with a message that names the problem", {
  empty <- withr::local_tempdir()
  bad_id <- withr::local_tempdir()
  make_dataset(bad_id, "DS A")
  no_asv <- withr::local_tempdir()
  make_dataset(no_asv, "DS-A", drop = "asv")

  expect_error(load_data(file.path(empty, "nope")), "was not found")
  expect_error(load_data(empty), "No ZIP files or dataset directories")
  expect_error(load_data(bad_id), "Invalid dataset id")
  expect_error(load_data(no_asv), "missing required files: asv.tsv")
})

test_that("load_data leaves no unzipped files behind", {
  leftovers <- function() list.files(tempdir(), pattern = "^asv_zip_")
  before <- leftovers()
  load_data(tiny())

  expect_equal(leftovers(), before)
})

test_that("a blank cell among filled taxonomy cells is NA, not an empty string", {
  data_path <- withr::local_tempdir()
  root <- file.path(data_path, "DS-A")
  write_dataset_files(root, "DS-A", n_events = 2, n_asvs = 3)
  asv <- utils::read.delim(file.path(root, "asv.tsv"), colClasses = "character", na.strings = NULL)
  asv$genus[3] <- ""
  write_tsv(asv, file.path(root, "asv.tsv"))

  expect_warning(asvs <- load_data(data_path)$asvs[["DS-A"]], "unzipped dataset folder")
  expect_equal(asvs$genus, c("Gen", "Gen", NA))
})
