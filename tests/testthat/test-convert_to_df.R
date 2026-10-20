test_that("a data table becomes a data.frame with the first column as row names", {
  dt <- data.table::data.table(id = c("a", "b"), x = c(1, 2), y = c("u", "v"))
  df <- convert_to_df(dt)

  expect_equal(class(df), "data.frame")
  expect_equal(rownames(df), c("a", "b"))
  expect_equal(names(df), c("x", "y"))
  expect_equal(df$x, c(1, 2))
})

test_that("a missing id becomes 'Unnamed-clades'", {
  dt <- data.table::data.table(id = c("a", NA), x = c(1, 2))

  expect_equal(rownames(convert_to_df(dt)), c("a", "Unnamed-clades"))
})

test_that("a list is converted element by element, keeping its names and nesting", {
  dt <- data.table::data.table(id = c("a", "b"), x = c(1, 2))
  out <- convert_to_df(list(first = dt, second = list(inner = dt)))

  expect_named(out, c("first", "second"))
  expect_equal(class(out$first), "data.frame")
  expect_equal(class(out$second$inner), "data.frame")
  expect_equal(rownames(out$second$inner), c("a", "b"))
})

test_that("sparse matrices are left as they are and a warning names where they are", {
  counts <- sparse(c("t1", "t2"), c("a", "b"), c(1, 2, 3, 4))
  merged <- list(counts = counts, events = event_table(c("a", "b"), x = 1:2))

  expect_warning(out <- convert_to_df(merged), "Skipped converting sparse counts at: merged\\$counts")
  expect_s4_class(out$counts, "sparseMatrix")
  expect_equal(class(out$events), "data.frame")
})

test_that("convert_counts = TRUE makes a data.frame with the ASVs as row names and warns about memory", {
  counts <- sparse(c("t1", "t2"), c("a", "b"), c(1, 2, 3, 4))

  expect_warning(df <- convert_to_df(counts, convert_counts = TRUE), "Converting sparse matrices to dense data.frame")
  expect_equal(class(df), "data.frame")
  expect_equal(rownames(df), c("t1", "t2"))
  expect_equal(df$b, c(3, 4))
})

test_that("a matrix above max_cells stays sparse with a warning that names it", {
  counts <- sparse(c("t1", "t2"), c("a", "b"), c(1, 2, 3, 4))

  expect_warning(out <- convert_to_df(counts, convert_counts = TRUE, max_cells = 3), "Skipped converting sparse counts")
  expect_s4_class(out, "sparseMatrix")
})

test_that("input that is already data.frames is returned as it is, with a message", {
  df <- data.frame(x = 1:2)

  expect_message(out <- convert_to_df(list(a = df)), "already converted")
  expect_identical(out, list(a = df))
})

test_that("an element that is neither table, matrix nor list is an error that names where", {
  dt <- data.table::data.table(id = "a", x = 1)

  expect_error(convert_to_df(list(ok = dt, bad = "text")), "Unsupported input type at .*\\$bad \\(class: character\\)")
})

test_that("the merged output of datasets from files converts, with the counts untouched by default", {
  data_path <- withr::local_tempdir()
  make_dataset(data_path, "DS-A", n_events = 3, n_asvs = 3)
  merged <- merge_data(load_data(data_path))

  expect_warning(out <- convert_to_df(merged), "Skipped converting sparse counts")
  expect_equal(rownames(out$asvs), asv_ids(3))
  expect_equal(rownames(out$events), paste0("DS-A:S", 1:3))
})
