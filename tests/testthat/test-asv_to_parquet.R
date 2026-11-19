two_datasets <- function(env = parent.frame()) {
  data_path <- withr::local_tempdir(.local_envir = env)
  make_dataset(data_path, "DS-A", n_events = 3, asv_from = 1, n_asvs = 5)
  make_dataset(data_path, "DS-B", n_events = 4, asv_from = 3, n_asvs = 6)
  data_path
}

convert <- function(data_path, env = parent.frame()) {
  out <- file.path(withr::local_tempdir(.local_envir = env), "parquet")
  asv_to_parquet(data_path, out)
  out
}

test_that("occurrence rows equal the loaded counts of every dataset, so no row is lost, doubled or moved", {
  data_path <- two_datasets()
  out <- convert(data_path)
  loaded <- load_data(data_path)
  occ <- read_pq(out, "occurrence")

  for (id in c("DS-A", "DS-B")) {
    got <- sorted(occ[occ$datasetID == id, c("eventID", "taxonID", "organismQuantity")], c("eventID", "taxonID"))
    expect_gt(nrow(got), 0)
    expect_equal(got, triplets(loaded$counts[[id]]), ignore_attr = TRUE)
  }
  expect_equal(nrow(occ), sum(vapply(loaded$counts, Matrix::nnzero, numeric(1))))
})

test_that("organismQuantity is numeric and sums per event to the column sums of the loaded counts", {
  data_path <- two_datasets()
  out <- convert(data_path)
  loaded <- load_data(data_path)
  occ <- read_pq(out, "occurrence")

  expect_true(is.numeric(occ$organismQuantity))
  totals <- tapply(occ$organismQuantity, occ$eventID, sum)
  expected <- unlist(lapply(loaded$counts, function(m) Matrix::colSums(m)), use.names = FALSE)
  names(expected) <- unlist(lapply(loaded$counts, colnames), use.names = FALSE)
  expect_equal(as.numeric(totals[names(expected)]), as.numeric(expected))
})

test_that("the five tables are written as one hive partition per dataset, keyed by datasetID", {
  out <- convert(two_datasets())

  for (tbl in all_tables) {
    expect_setequal(list.files(file.path(out, tbl)), c("datasetID=DS-A", "datasetID=DS-B"))
    expect_gt(length(list.files(file.path(out, tbl, "datasetID=DS-A"), pattern = "\\.parquet$")), 0)
  }
})

test_that("the datasetID of a partition equals the prefix of every eventID in it", {
  out <- convert(two_datasets())

  for (tbl in c("occurrence", "event", "emof")) {
    x <- read_pq(out, tbl)
    expect_equal(sub(":.*$", "", x$eventID), x$datasetID)
  }
})

test_that("dataset_pid is in no table, since it is not unique across datasets", {
  out <- convert(two_datasets())

  for (tbl in all_tables) expect_false("dataset_pid" %in% names(read_pq(out, tbl)))
})

test_that("blank taxonomy cells in asv are NA, as in load_data, and the sequence and name are kept", {
  data_path <- two_datasets()
  out <- convert(data_path)
  loaded <- load_data(data_path)
  asv <- read_pq(out, "asv")
  tax_cols <- c("kingdom", "phylum", "order", "class", "family", "genus", "specificEpithet", "infraspecificEpithet", "otu")

  for (id in c("DS-A", "DS-B")) {
    got <- sorted(asv[asv$datasetID == id, ], "taxonID")
    ref <- sorted(as.data.frame(loaded$asvs[[id]]), "taxonID")
    expect_equal(lapply(got[tax_cols], as.character), lapply(ref[tax_cols], as.character))
    expect_equal(got[c("taxonID", "DNA_sequence", "scientificName")], ref[c("taxonID", "DNA_sequence", "scientificName")])
    expect_setequal(setdiff(names(got), "datasetID"), names(ref))
  }
  expect_true(all(is.na(asv$phylum)))
  expect_false(any(is.na(asv$genus)))
})

test_that("an ASV in two datasets has one row in each, with the same taxonID", {
  asv <- read_pq(convert(two_datasets()), "asv")
  shared <- asv_ids(7)[3:5]

  expect_equal(sort(table(asv$taxonID[asv$taxonID %in% shared])), sort(table(rep(shared, 2))), ignore_attr = TRUE)
  expect_setequal(asv$datasetID[asv$taxonID == shared[1]], c("DS-A", "DS-B"))
})

test_that("event keeps all columns of load_data and a cell that starts with an unclosed quote does not swallow rows", {
  data_path <- two_datasets()
  out <- convert(data_path)
  loaded <- load_data(data_path)
  event <- read_pq(out, "event")

  expect_equal(nrow(event), 3 + 4)
  expect_setequal(setdiff(names(event), "datasetID"), names(loaded$events[["DS-A"]]))
  expect_true(all(event$remarks == "\"Niskin, 5 m"))
  expect_true(all(event$sop == "\"\""))
})

test_that("emof stays long, one row per event and measurement, with the value kept as text", {
  data_path <- two_datasets()
  emof <- read_pq(convert(data_path), "emof")

  expect_equal(nrow(emof), 2 * 3 + 2 * 4)
  expect_setequal(emof$measurementType, c("salinity", "temperature"))
  expect_true(all(c("eventID", "measurementType", "measurementUnit", "measurementValue") %in% names(emof)))
  expect_equal(emof$measurementValue[emof$eventID == "DS-A:S1" & emof$measurementType == "salinity"], "1.5")
})

test_that("a dataset with an empty emof gets an empty partition and the other datasets keep their rows, without a warning", {
  data_path <- withr::local_tempdir()
  make_dataset(data_path, "DS-A")
  make_dataset(data_path, "DS-E", empty_emof = TRUE)
  out <- file.path(withr::local_tempdir(), "parquet")

  expect_no_warning(asv_to_parquet(data_path, out))
  emof <- read_pq(out, "emof")
  expect_equal(unique(emof$datasetID), "DS-A")
  expect_true(dir.exists(file.path(out, "emof", "datasetID=DS-E")))
})

test_that("dataset has one row per dataset with the name from event.tsv, equal to load_data", {
  data_path <- two_datasets()
  out <- convert(data_path)
  loaded <- load_data(data_path)
  ds <- read_pq(out, "dataset")
  ref <- do.call(rbind, lapply(loaded$datasets, as.data.frame))

  expect_equal(nrow(ds), 2)
  expect_equal(sorted(ds[c("datasetID", "datasetName")], "datasetID"), sorted(ref, "datasetID"), ignore_attr = TRUE)
})

test_that("a zip with a wrapper folder and macOS metadata converts like a plain zip", {
  plain <- withr::local_tempdir()
  wrapped <- withr::local_tempdir()
  make_dataset(plain, "DS-A")
  make_dataset(wrapped, "DS-A", wrapper = TRUE)
  a <- convert(plain)
  b <- convert(wrapped)

  for (tbl in all_tables) {
    key <- if (tbl == "dataset") "datasetID" else names(read_pq(a, tbl))[1:3]
    expect_equal(sorted(read_pq(b, tbl), key), sorted(read_pq(a, tbl), key))
  }
})

test_that("an unzipped folder gives the same tables as the zip and warns as load_data does", {
  zipped <- withr::local_tempdir()
  folder <- withr::local_tempdir()
  make_dataset(zipped, "DS-A")
  make_dataset(folder, "DS-A", as = "dir")
  a <- convert(zipped)
  out <- file.path(withr::local_tempdir(), "parquet")

  expect_warning(asv_to_parquet(folder, out), "unzipped dataset folder")
  for (tbl in all_tables) {
    key <- if (tbl == "dataset") "datasetID" else names(read_pq(a, tbl))[1:3]
    expect_equal(sorted(read_pq(out, tbl), key), sorted(read_pq(a, tbl), key))
  }
})

test_that("when a zip and a folder have the same ID the zip is used and the user is told", {
  data_path <- withr::local_tempdir()
  make_dataset(data_path, "DS-A", n_events = 3)
  other <- withr::local_tempdir()
  make_dataset(other, "DS-A", as = "dir", n_events = 6)
  file.copy(file.path(other, "DS-A"), data_path, recursive = TRUE)
  out <- file.path(withr::local_tempdir(), "parquet")

  expect_message(asv_to_parquet(data_path, out), "Using ZIP")
  expect_equal(nrow(read_pq(out, "event")), 3)
})

test_that("running again replaces a changed dataset and keeps the others", {
  data_path <- two_datasets()
  out <- convert(data_path)
  before_b <- sorted(read_pq(out, "occurrence"), c("eventID", "taxonID"))
  before_b <- before_b[before_b$datasetID == "DS-B", ]

  unlink(file.path(data_path, "DS-A.zip"))
  make_dataset(data_path, "DS-A", n_events = 2, asv_from = 1, n_asvs = 5)
  asv_to_parquet(data_path, out)
  occ <- read_pq(out, "occurrence")

  expect_setequal(unique(occ$eventID[occ$datasetID == "DS-A"]), c("DS-A:S1", "DS-A:S2"))
  expect_equal(sorted(occ[occ$datasetID == "DS-B", ], c("eventID", "taxonID")), before_b, ignore_attr = TRUE)
  expect_length(list.files(file.path(out, "occurrence", "datasetID=DS-A")), 1)
})

test_that("a second run on the same input gives the same rows, not twice as many", {
  data_path <- two_datasets()
  out <- convert(data_path)
  n <- vapply(all_tables, function(t) nrow(read_pq(out, t)), numeric(1))

  asv_to_parquet(data_path, out)

  expect_equal(vapply(all_tables, function(t) nrow(read_pq(out, t)), numeric(1)), n)
})

test_that("the function returns out_dir invisibly", {
  data_path <- two_datasets()
  out <- file.path(withr::local_tempdir(), "parquet")

  expect_invisible(res <- asv_to_parquet(data_path, out))
  expect_equal(normalizePath(res), normalizePath(out))
})

test_that("a zip leaves no extracted files behind, also when the conversion fails", {
  leftovers <- function() list.files(tempdir(), pattern = "^asv_")
  data_path <- two_datasets()
  before <- leftovers()
  convert(data_path)
  expect_equal(leftovers(), before)

  broken <- withr::local_tempdir()
  make_dataset(broken, "DS-X", drop = "asv")
  expect_error(asv_to_parquet(broken, file.path(withr::local_tempdir(), "parquet")), "asv.tsv")
  expect_equal(leftovers(), before)
})

test_that("a dataset with a missing file stops the run before anything is written for any dataset", {
  data_path <- withr::local_tempdir()
  make_dataset(data_path, "DS-A")
  make_dataset(data_path, "DS-X", drop = "asv")
  out <- file.path(withr::local_tempdir(), "parquet")

  expect_error(asv_to_parquet(data_path, out), "asv.tsv")
  expect_false(dir.exists(file.path(out, "occurrence")))
})

test_that("bad input is an error that names the problem", {
  out <- file.path(withr::local_tempdir(), "parquet")
  empty <- withr::local_tempdir()
  bad_id <- withr::local_tempdir()
  make_dataset(bad_id, "DS A")
  renamed <- withr::local_tempdir()
  make_dataset(renamed, "DS-A")
  file.rename(file.path(renamed, "DS-A.zip"), file.path(renamed, "DS-Z.zip"))

  expect_error(asv_to_parquet(file.path(empty, "nope"), out), "was not found")
  expect_error(asv_to_parquet(empty, out), "No ZIP files or dataset directories")
  expect_error(asv_to_parquet(bad_id, out), "Invalid dataset id")
  expect_error(asv_to_parquet(renamed, out), "DS-Z.*DS-A|DS-A.*DS-Z")
})

ecolongr_names <- list(
  occurrence = c("sample", "taxon", "quant"),
  asv = c(
    "taxon", "asv_id_alias", "dna_sequence", "scientific_name", "taxon_rank", "kingdom", "phylum",
    "class", "order", "family", "genus", "specific_epithet", "infraspecific_epithet", "otu",
    "date_identified", "identification_references", "associated_sequences",
    "previous_identifications", "identification_remarks"
  ),
  event = c(
    "sample", "event_date", "location_id", "verbatim_locality", "decimal_latitude",
    "coordinate_uncertainty_in_meters", "material_sample_id", "references",
    "pcr_primer_name_forward", "sop", "remarks"
  ),
  emof = c(
    "sample", "measurement_type", "measurement_type_id", "measurement_unit",
    "measurement_unit_id", "measurement_value"
  ),
  dataset = "dataset_name"
)

convert_as <- function(data_path, name_convert, env = parent.frame()) {
  out <- file.path(withr::local_tempdir(.local_envir = env), "parquet")
  asv_to_parquet(data_path, out, name_convert = name_convert)
  out
}

test_that("name_convert = 'ecolongr' gives exactly these columns per table, with dataset_id for the partition", {
  out <- convert_as(two_datasets(), "ecolongr")

  for (tbl in all_tables) {
    expect_setequal(names(read_pq(out, tbl)), c(ecolongr_names[[tbl]], "dataset_id"))
  }
})

test_that("the default and name_convert = 'none' keep the portal's names and datasetID", {
  data_path <- two_datasets()
  default <- convert(data_path)
  none <- convert_as(data_path, "none")

  expect_true(all(c("eventID", "taxonID", "organismQuantity", "datasetID") %in% names(read_pq(default, "occurrence"))))
  for (tbl in all_tables) expect_equal(names(read_pq(none, tbl)), names(read_pq(default, tbl)))
  expect_error(asv_to_parquet(data_path, tempfile(), name_convert = "camel"), "should be one of")
})

test_that("with 'ecolongr' the partition directories are dataset_id=<id> and no table has a datasetID one", {
  out <- convert_as(two_datasets(), "ecolongr")

  for (tbl in all_tables) expect_setequal(list.files(file.path(out, tbl)), c("dataset_id=DS-A", "dataset_id=DS-B"))
})

test_that("renaming changes no values: the ecolongr tables equal the default ones under the new names", {
  data_path <- two_datasets()
  plain <- convert(data_path)
  named <- convert_as(data_path, "ecolongr")
  renames <- c(eventID = "sample", taxonID = "taxon", organismQuantity = "quant", datasetID = "dataset_id")

  occ <- read_pq(plain, "occurrence")
  names(occ) <- renames[names(occ)]
  expect_equal(sorted(read_pq(named, "occurrence"), c("sample", "taxon")), sorted(occ, c("sample", "taxon")))

  ev <- read_pq(plain, "event")
  expect_equal(sorted(read_pq(named, "event"), "sample")$decimal_latitude, sorted(ev, "eventID")$decimalLatitude)
  expect_equal(sorted(read_pq(named, "event"), "sample")$sample, sorted(ev, "eventID")$eventID)
})

test_that("the key columns are named alike in every table, so the tables join without renaming", {
  out <- convert_as(two_datasets(), "ecolongr")
  occ <- read_pq(out, "occurrence")
  event <- read_pq(out, "event")
  asv <- read_pq(out, "asv")
  emof <- read_pq(out, "emof")

  expect_equal(nrow(merge(occ, event, by = c("dataset_id", "sample"))), nrow(occ))
  expect_equal(nrow(merge(occ, asv, by = c("dataset_id", "taxon"))), nrow(occ))
  expect_true(all(emof$sample %in% event$sample))
})

test_that("snake_case splits at case changes, keeps an acronym together and does not split before a digit", {
  data_path <- withr::local_tempdir()
  make_dataset(data_path, "DS-A", extra_event = list(DNAExtractionID = "x", sampleSize2Value = "1", ITSPrimer = "y"))
  event <- read_pq(convert_as(data_path, "ecolongr"), "event")

  expect_true(all(c("dna_extraction_id", "sample_size2_value", "its_primer") %in% names(event)))
})

test_that("two columns that get the same snake_case name are an error that names both, and nothing is written", {
  data_path <- withr::local_tempdir()
  make_dataset(data_path, "DS-A", extra_event = list(event_date = "dup"))
  out <- file.path(withr::local_tempdir(), "parquet")

  expect_error(asv_to_parquet(data_path, out, name_convert = "ecolongr"), "eventDate.*event_date|event_date.*eventDate")
  expect_false(dir.exists(file.path(out, "event")))
  expect_no_error(asv_to_parquet(data_path, out, name_convert = "none"))
})

test_that("a run with the other name_convert than the existing tables is an error that leaves them as they are", {
  data_path <- two_datasets()
  for (pair in list(c("none", "ecolongr"), c("ecolongr", "none"))) {
    out <- convert_as(data_path, pair[1])
    before <- read_pq(out, "occurrence")

    expect_error(asv_to_parquet(data_path, out, name_convert = pair[2]), "name_convert")
    expect_equal(read_pq(out, "occurrence"), before)
    expect_false(file.exists(file.path(out, ".staging")))
  }
})

test_that("a second run with the same name_convert is fine for both settings", {
  data_path <- two_datasets()
  for (nc in c("none", "ecolongr")) {
    out <- convert_as(data_path, nc)
    expect_no_error(asv_to_parquet(data_path, out, name_convert = nc))
  }
})

test_that("el_rda runs on the ecolongr tables as they are, without renaming", {
  skip_if_not_installed("ecolongr")
  out <- convert_as(two_datasets(), "ecolongr")
  occ <- read_pq(out, "occurrence")

  fit <- suppressMessages(ecolongr::el_rda(occ, k = 2, transform = "clr", pseudocount = 1))
  expect_equal(nrow(fit$CA$u), 7)
})

test_that("the real portal download gives the loaded counts back", {
  zip_file <- Sys.getenv("ASVOCCUR_TEST_ZIP")
  skip_if(!nzchar(zip_file) || !file.exists(zip_file), "set ASVOCCUR_TEST_ZIP to a portal zip")
  data_path <- withr::local_tempdir()
  file.copy(zip_file, data_path)
  id <- sub("\\.zip$", "", basename(zip_file))
  out <- convert(data_path)
  loaded <- load_data(data_path)
  occ <- read_pq(out, "occurrence")

  expect_equal(sorted(occ[c("eventID", "taxonID", "organismQuantity")], c("eventID", "taxonID")),
               triplets(loaded$counts[[id]]), ignore_attr = TRUE)
  expect_equal(nrow(read_pq(out, "asv")), nrow(loaded$asvs[[id]]))
  expect_equal(nrow(read_pq(out, "event")), nrow(loaded$events[[id]]))
})

test_that("the real portal download gets only snake_case names with 'ecolongr', and the known awkward ones right", {
  zip_file <- Sys.getenv("ASVOCCUR_TEST_ZIP")
  skip_if(!nzchar(zip_file) || !file.exists(zip_file), "set ASVOCCUR_TEST_ZIP to a portal zip")
  data_path <- withr::local_tempdir()
  file.copy(zip_file, data_path)
  out <- convert_as(data_path, "ecolongr")

  for (tbl in all_tables) {
    expect_match(names(read_pq(out, tbl)), "^[a-z][a-z0-9]*(_[a-z0-9]+)*$")
  }
  expect_true(all(c(
    "location_id", "institution_id", "material_sample_id", "coordinate_uncertainty_in_meters",
    "minimum_depth_in_meters", "pcr_primer_name_forward", "references", "sample"
  ) %in% names(read_pq(out, "event"))))
  expect_true(all(c("dna_sequence", "taxon", "specific_epithet") %in% names(read_pq(out, "asv"))))
})
