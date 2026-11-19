test_that("counts are the union of ASVs and events, zero where an ASV is absent, ASVs sorted", {
  merged <- merge_data(two_loaded())
  expected <- matrix(
    c(1, 3, 0, 0,
      2, 4, 5, 0,
      0, 0, 0, 6),
    nrow = 3, byrow = TRUE, dimnames = list(c("t1", "t2", "t3"), c("a1", "a2", "b1", "b2"))
  )

  expect_equal(as.matrix(merged$counts), expected)
})

test_that("asvs has one row per ASV, keyed on taxonID", {
  asvs <- merge_data(two_loaded())$asvs

  expect_equal(asvs$taxonID, c("t1", "t2", "t3"))
  expect_equal(data.table::key(asvs), "taxonID")
})

test_that("events are stacked, keyed on eventID, and columns that look numeric are numeric again", {
  events <- merge_data(two_loaded())$events

  expect_equal(events$eventID, c("a1", "a2", "b1", "b2"))
  expect_equal(data.table::key(events), "eventID")
  expect_equal(events$decimalLatitude, c(55.5, 56, 60, 61.25))
  expect_equal(events$locationID, c(1, 2, 3, 4))
  expect_equal(events$note, c("x", "y", "z", NA))
})

test_that("emof from datasets with different measurements gets NA where one lacks a column", {
  emof <- merge_data(two_loaded())$emof

  expect_equal(emof[["salinity (psu)"]], c(1.5, 2.5, NA, NA))
  expect_equal(emof[["temperature (C)"]], c(NA, NA, 7, 8))
})

test_that("datasets are stacked and keyed on datasetID", {
  datasets <- merge_data(two_loaded())$datasets

  expect_equal(datasets$datasetID, c("DS-A", "DS-B"))
  expect_equal(data.table::key(datasets), "datasetID")
})

test_that("consistent datasets give no iasvs and no warning", {
  expect_no_warning(merged <- merge_data(two_loaded()))
  expect_null(merged$iasvs)
})

test_that("an ASV with a different name in two datasets warns and is listed in iasvs", {
  expect_warning(merged <- merge_data(two_loaded(t2_in_b = "Renamed")), "Inconsistent ASV taxonomy")

  expect_equal(merged$asvs$taxonID, c("t1", "t2", "t2", "t3"))
  expect_equal(merged$iasvs$datasetID, c("DS-A", "DS-B"))
  expect_equal(merged$iasvs$taxonID, c("t2", "t2"))
  expect_equal(merged$iasvs$scientificName, c("Taxon t2", "Renamed"))
})

test_that("ds selects the datasets to merge", {
  merged <- merge_data(two_loaded(), ds = "DS-B")

  expect_equal(colnames(merged$counts), c("b1", "b2"))
  expect_equal(merged$events$eventID, c("b1", "b2"))
  expect_equal(merged$asvs$taxonID, c("t2", "t3"))
})

test_that("a ds that was not loaded is an error that names it", {
  expect_error(merge_data(two_loaded(), ds = c("DS-A", "DS-Z")), "Invalid dataset ID\\(s\\) specified in ds: DS-Z")
})

test_that("an eventID in two datasets is an error that names it", {
  loaded <- two_loaded()
  colnames(loaded$counts[["DS-B"]]) <- c("a1", "b2")
  loaded$events[["DS-B"]] <- event_table(c("a1", "b2"), decimalLatitude = c(60, 61.25), locationID = c("3", "4"), note = c("z", NA))

  expect_error(merge_data(loaded), "Duplicated eventID\\(s\\) found.*a1")
})

test_that("input that is not a list of loaded lists is an error", {
  expect_error(merge_data(two_loaded()$counts), "hierarchical list")
  expect_error(merge_data("a"), "hierarchical list")
})

test_that("datasets read from files merge into the union of events and ASVs", {
  data_path <- withr::local_tempdir()
  make_dataset(data_path, "DS-A", n_events = 3, asv_from = 1, n_asvs = 5)
  make_dataset(data_path, "DS-B", n_events = 4, asv_from = 3, n_asvs = 6)
  merged <- merge_data(load_data(data_path))

  expect_equal(merged$events$eventID, c(paste0("DS-A:S", 1:3), paste0("DS-B:S", 1:4)))
  expect_equal(merged$asvs$taxonID, asv_ids(8))
  expect_null(merged$iasvs)
  expect_true(all(rownames(merged$counts) %in% asv_ids(8)))
  expect_equal(ncol(merged$counts), 7)
})
