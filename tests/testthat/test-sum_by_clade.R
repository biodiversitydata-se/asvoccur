# Four ASVs in two samples. t1, t2 share genus G1; t3 has genus G2; t4 has no genus.
# Only t1 and t3 have a specificEpithet.
clade_input <- function() {
  asvs <- data.table::data.table(
    taxonID = c("t1", "t2", "t3", "t4"),
    kingdom = "Bacteria",
    phylum = c("P1", "P1", "P2", NA),
    class = c("C1", "C1", "C2", NA),
    order = c("O1", "O1", "O2", NA),
    family = c("F1", "F1", "F2", NA),
    genus = c("G1", "G1", "G2", NA),
    specificEpithet = c("sp1", NA, "sp3", NA),
    otu = c("Otu1", NA, NA, NA)
  )
  counts <- sparse(asvs$taxonID, c("s1", "s2"), c(1, 2, 3, 4, 5, 0, 7, 8))
  list(counts = counts, asvs = asvs)
}

# s1 totals 10 and s2 totals 20.
genus_dimnames <- list(c("G1", "G2", "Unclassified"), c("s1", "s2"))

test_that("one result per rank, in the order kingdom to otu with species before otu", {
  x <- clade_input()
  summed <- sum_by_clade(x$counts, x$asvs)

  expect_named(summed, c("raw", "norm"))
  expect_equal(names(summed$raw), c("kingdom", "phylum", "class", "order", "family", "genus", "species", "otu"))
  expect_equal(names(summed$norm), names(summed$raw))
})

test_that("raw sums the counts of the ASVs in a clade, with NA as 'Unclassified', clades sorted", {
  x <- clade_input()
  raw <- sum_by_clade(x$counts, x$asvs)$raw

  expect_equal(as.matrix(raw$genus), matrix(c(3, 5, 3, 7, 4, 8), nrow = 3, byrow = TRUE, dimnames = genus_dimnames))
  expect_equal(as.matrix(raw$kingdom), matrix(c(10, 20), nrow = 1, dimnames = list("Bacteria", c("s1", "s2"))))
  expect_equal(as.matrix(raw$otu), matrix(c(1, 5, 9, 15), nrow = 2, byrow = TRUE, dimnames = list(c("Otu1", "Unclassified"), c("s1", "s2"))))
})

test_that("species is 'genus specificEpithet', and an ASV without an epithet is unclassified", {
  x <- clade_input()
  species <- sum_by_clade(x$counts, x$asvs)$raw$species

  expect_equal(as.matrix(species), matrix(c(1, 5, 3, 7, 6, 8), nrow = 3, byrow = TRUE,
                                          dimnames = list(c("G1 sp1", "G2 sp3", "Unclassified"), c("s1", "s2"))))
})

test_that("norm divides each sample's clade sums by the sample total, so every column sums to 1", {
  x <- clade_input()
  norm <- sum_by_clade(x$counts, x$asvs)$norm

  expect_equal(as.matrix(norm$genus), matrix(c(0.3, 0.25, 0.3, 0.35, 0.4, 0.4), nrow = 3, byrow = TRUE, dimnames = genus_dimnames))
  for (rank in names(norm)) expect_equal(unname(Matrix::colSums(norm[[rank]])), c(1, 1))
})

test_that("a sample with no counts normalises to 0, not NaN", {
  x <- clade_input()
  counts <- cbind(x$counts, s3 = c(0, 0, 0, 0))
  norm <- sum_by_clade(counts, x$asvs)$norm$genus

  expect_equal(unname(as.matrix(norm)[, "s3"]), c(0, 0, 0))
  expect_false(anyNA(as.matrix(norm)))
})

test_that("convert_to_dt gives data tables with the clade first, and warns once about memory", {
  x <- clade_input()

  expect_warning(summed <- sum_by_clade(x$counts, x$asvs, convert_to_dt = TRUE), "Converting sparse matrices to dense")
  genus <- summed$raw$genus
  expect_s3_class(genus, "data.table")
  expect_equal(names(genus), c("clade", "s1", "s2"))
  expect_equal(genus$clade, c("G1", "G2", "Unclassified"))
  expect_equal(genus$s2, c(5, 7, 8))
  expect_equal(summed$norm$genus$s1, c(0.3, 0.3, 0.4))
})

test_that("matrices above max_cells stay sparse and the warning names them", {
  x <- clade_input()
  warnings <- testthat::capture_warnings(summed <- sum_by_clade(x$counts, x$asvs, convert_to_dt = TRUE, max_cells = 5))

  expect_length(warnings, 2)
  expect_match(warnings[2], "Skipped converting large sparse matrices at: .*raw\\$genus")
  expect_no_match(warnings[2], "otu|kingdom")
  expect_s4_class(summed$raw$genus, "sparseMatrix")
  expect_s3_class(summed$raw$kingdom, "data.table")
  expect_s3_class(summed$raw$otu, "data.table")
})

test_that("counts whose rows differ from asvs$taxonID in content or order are an error", {
  x <- clade_input()

  expect_error(sum_by_clade(x$counts[4:1, ], x$asvs), "Mismatch detected")
  expect_error(sum_by_clade(x$counts[1:3, ], x$asvs), "Mismatch detected")
})

test_that("counts that are not a sparse matrix, or asvs that are not a data table, are an error", {
  x <- clade_input()

  expect_error(sum_by_clade(as.matrix(x$counts), x$asvs), "Input must be a sparse matrix")
  expect_error(sum_by_clade(x$counts, as.data.frame(x$asvs)), "Input must be a data table")
})

test_that("merged data from files sums to the same total at every rank", {
  data_path <- withr::local_tempdir()
  make_dataset(data_path, "DS-A", n_events = 3, asv_from = 1, n_asvs = 2)
  merged <- merge_data(load_data(data_path))
  summed <- sum_by_clade(merged$counts, merged$asvs)

  total <- sum(merged$counts)
  for (rank in names(summed$raw)) expect_equal(sum(summed$raw[[rank]]), total)
})
