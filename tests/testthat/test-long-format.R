skip_if_not_installed("duckdb")
skip_if_not_installed("duckplyr")
skip_if_not_installed("DBI")

two_datasets <- function(env = parent.frame()) {
  data_path <- withr::local_tempdir(.local_envir = env)
  make_dataset(data_path, "DS-A", n_events = 3, asv_from = 1, n_asvs = 5)
  make_dataset(data_path, "DS-B", n_events = 4, asv_from = 3, n_asvs = 6)
  data_path
}

open_long <- function(data_path, env = parent.frame()) {
  store <- file.path(withr::local_tempdir(.local_envir = env), "store")
  load_long(data_path, store = store)
}

test_that("load_long() converts the downloads once and returns lazy tables", {
  data_path <- two_datasets()
  store <- file.path(withr::local_tempdir(), "store")

  loaded <- load_long(data_path, store = store)
  expect_named(loaded, c("counts", "asvs", "events", "emof", "datasets", "store"))
  for (tbl in c("counts", "asvs", "events", "emof")) expect_s3_class(loaded[[tbl]], "duckplyr_df")
  expect_equal(loaded$datasets$datasetID, c("DS-A", "DS-B"))

  stamp <- file.mtime(Sys.glob(file.path(store, "occurrence", "*", "*.parquet")))
  expect_message(again <- load_long(data_path, store = store), "stored copies.*DS-A, DS-B")
  expect_equal(file.mtime(Sys.glob(file.path(store, "occurrence", "*", "*.parquet"))), stamp)
  expect_equal(again$datasets, loaded$datasets)
})

test_that("load_long() converts all downloads again if one is missing from the store", {
  data_path <- two_datasets()
  store <- file.path(withr::local_tempdir(), "store")
  load_long(data_path, store = store)
  make_dataset(data_path, "DS-C", n_events = 2)

  expect_no_message(loaded <- load_long(data_path, store = store))
  expect_equal(loaded$datasets$datasetID, c("DS-A", "DS-B", "DS-C"))
})

test_that("load_long() opens an existing store without the downloads", {
  data_path <- two_datasets()
  store <- file.path(withr::local_tempdir(), "store")
  asv_to_parquet(data_path, store)

  loaded <- load_long(store = store)
  expect_equal(loaded$datasets$datasetID, c("DS-A", "DS-B"))
})

test_that("load_long() refuses a store written with ecolongr names", {
  store <- file.path(withr::local_tempdir(), "store")
  asv_to_parquet(two_datasets(), store, name_convert = "ecolongr")

  expect_error(load_long(store = store), "ecolongr")
})

test_that("merge_long() keeps the same counts as merge_data()", {
  data_path <- two_datasets()
  merged <- merge_long(open_long(data_path))
  expected <- merge_data(load_data(data_path))

  observed <- sorted(as.data.frame(merged$counts), c("eventID", "taxonID"))
  expect_equal(observed, triplets(expected$counts), ignore_attr = TRUE)
})

test_that("merge_long() has one row per ASV and in-memory events", {
  merged <- merge_long(open_long(two_datasets()))
  asvs <- as.data.frame(merged$asvs)

  expect_equal(sort(asvs$taxonID), asv_ids(8))
  expect_false("datasetID" %in% names(asvs))
  expect_true(all(is.na(asvs$species)))
  expect_false(inherits(merged$events, "duckplyr_df"))
  expect_equal(nrow(merged$events), 3 + 4)
  expect_null(merged$iasvs)
})

test_that("merge_long() merges only the datasets in ds", {
  merged <- merge_long(open_long(two_datasets()), ds = "DS-B")

  expect_equal(unique(merged$events$datasetID), "DS-B")
  expect_equal(merged$datasets$datasetID, "DS-B")
  expect_true(all(startsWith(as.data.frame(merged$counts)$eventID, "DS-B:")))
  expect_error(merge_long(open_long(two_datasets()), ds = "DS-X"), "DS-X")
})

test_that("merge_long() warns about ASVs with inconsistent taxonomy", {
  data_path <- withr::local_tempdir()
  make_dataset(data_path, "DS-A", n_events = 2, asv_from = 1, n_asvs = 3)
  make_dataset(data_path, "DS-B", as = "dir", n_events = 2, asv_from = 2, n_asvs = 3)
  asv_path <- file.path(data_path, "DS-B", "asv.tsv")
  asv <- data.table::fread(asv_path, sep = "\t", colClasses = "character")
  asv$scientificName[asv$taxonID == asv_ids(2)[2]] <- "Other name"
  write_tsv(asv, asv_path)

  loaded <- suppressWarnings(open_long(data_path))
  warnings <- capture_warnings(merged <- merge_long(loaded))
  expect_true(any(grepl("Inconsistent ASV taxonomy", warnings)))
  expect_equal(unique(merged$iasvs$taxonID), asv_ids(2)[2])
  expect_error(sum_clades_long(merged$counts, merged$asvs), "more than one row")
})

expect_same_clade_sums <- function(data_path) {
  merged <- merge_long(open_long(data_path))
  wide <- merge_data(load_data(data_path))
  # sum_by_clade() needs exactly the ASVs of the count matrix, in the same order.
  asvs <- wide$asvs[match(rownames(wide$counts), wide$asvs$taxonID)]
  expected <- sum_by_clade(wide$counts, asvs)

  observed <- as.data.frame(sum_clades_long(merged$counts, merged$asvs))
  expect_named(observed, c("rank", "clade", "eventID", "raw", "norm"))
  for (rank in names(expected$raw)) {
    obs <- observed[observed$rank == rank, ]
    raw <- expected$raw[[rank]]
    norm <- expected$norm[[rank]]
    for (i in seq_len(nrow(obs))) {
      expect_equal(as.numeric(obs$raw[i]), raw[obs$clade[i], obs$eventID[i]])
      expect_equal(obs$norm[i], norm[obs$clade[i], obs$eventID[i]])
    }
    expect_equal(sum(as.numeric(obs$raw)), sum(raw))
    expect_equal(nrow(obs), sum(raw != 0))
  }
}

test_that("sum_clades_long() matches sum_by_clade() at every rank", {
  expect_same_clade_sums(two_datasets())
})

test_that("species is genus and epithet, and NA without either, as in sum_by_clade()", {
  data_path <- withr::local_tempdir()
  make_dataset(data_path, "DS-A", as = "dir", n_events = 3, n_asvs = 5)
  asv_path <- file.path(data_path, "DS-A", "asv.tsv")
  asv <- data.table::fread(asv_path, sep = "\t", colClasses = "character")
  asv$genus <- c("Gen", "Gen", "Other", "", "Gen")
  asv$specificEpithet <- c("alpha", "beta", "alpha", "alpha", "")
  write_tsv(asv, asv_path)

  loaded <- suppressWarnings(open_long(data_path))
  expect_equal(
    sorted(as.data.frame(merge_long(loaded)$asvs), "taxonID")$species,
    c("Gen alpha", "Gen beta", "Other alpha", NA, NA)
  )
  suppressWarnings(expect_same_clade_sums(data_path))
})

test_that("sum_clades_long() respects a filtered table and stays lazy", {
  merged <- merge_long(open_long(two_datasets()))
  subset <- dplyr::filter(merged$counts, eventID == "DS-A:S1")

  summed <- sum_clades_long(subset, merged$asvs, ranks = "genus")
  expect_s3_class(summed, "duckplyr_df")
  result <- as.data.frame(summed)
  expect_equal(unique(result$eventID), "DS-A:S1")
  expect_equal(sum(result$norm), 1)
  expect_error(sum_clades_long(merged$counts, merged$asvs, ranks = "tribe"), "tribe")
})

test_that("values with quotes, tabs and line breaks are read as the portal wrote them", {
  data_path <- withr::local_tempdir()
  tricky <- c("\"Niskin, 5 m", "Station \"A\"\tdeep\nsecond line", "plain")
  make_dataset(data_path, "DS-Q", n_events = 3, extra_event = list(note = tricky))

  merged <- merge_long(open_long(data_path))
  expect_equal(nrow(merged$events), 3)
  expect_equal(merged$events$note, tricky)
  expect_true(all(is.na(merged$events$sop)))

  loaded <- load_data(data_path)
  expect_equal(loaded$events[["DS-Q"]]$note, tricky)
})

test_that("merge_long() ignores dataset-specific ASV attributes, as merge_data() does", {
  data_path <- withr::local_tempdir()
  make_dataset(data_path, "DS-A", n_events = 2, asv_from = 1, n_asvs = 3)
  make_dataset(data_path, "DS-B", as = "dir", n_events = 2, asv_from = 2, n_asvs = 3)
  asv_path <- file.path(data_path, "DS-B", "asv.tsv")
  asv <- data.table::fread(asv_path, sep = "\t", colClasses = "character")
  asv$asv_id_alias <- paste0("provider_", seq_len(nrow(asv)))
  asv$previousIdentifications <- "Bacteria|Other"
  write_tsv(asv, asv_path)

  merged <- expect_no_warning(merge_long(suppressWarnings(open_long(data_path))))
  asvs <- as.data.frame(merged$asvs)
  expect_null(merged$iasvs)
  expect_equal(sort(asvs$taxonID), asv_ids(4))
  expect_false(any(c("asv_id_alias", "previousIdentifications") %in% names(asvs)))
})
