# Internal helpers that find ASV portal datasets (ZIPs or folders) in a directory.

# Robust root detection: find the folder that contains occurrence.tsv + required TSVs.
# Works for:
# - flat datasets
# - datasets inside one or more nested folders
# - ZIPs re-packed on macOS (ignores __MACOSX and ._ files)
find_dataset_root <- function(dir_path) {

  occ <- list.files(
    dir_path,
    pattern = "^occurrence\\.tsv$",
    recursive = TRUE,
    full.names = TRUE
  )

  # Ignore macOS metadata artifacts
  occ <- occ[!grepl("(^|/)__MACOSX(/|$)", occ)]
  occ <- occ[!grepl("(^|/)\\._", occ)]

  if (length(occ) == 0) {
    stop("No occurrence.tsv found under: ", dir_path)
  }
  if (length(occ) > 1) {
    stop(
      "Multiple occurrence.tsv files found under: ", dir_path,
      "\nCannot determine unique dataset root."
    )
  }

  root <- dirname(occ)

  required <- c("asv.tsv", "event.tsv", "emof.tsv")
  missing <- required[!file.exists(file.path(root, required))]
  if (length(missing) > 0) {
    stop(
      "Dataset root found at ", root,
      " but missing required files: ",
      paste(missing, collapse = ", ")
    )
  }

  normalizePath(root, winslash = "/", mustWork = TRUE)
}

# Resolve source to a root directory containing TSVs.
# A ZIP is unzipped into a unique temp dir, returned as attribute "cleanup_dir" for the caller to remove;
# the temp dir is removed here if no dataset root is found in it.
resolve_root <- function(source_path) {
  if (grepl("\\.zip$", source_path, ignore.case = TRUE)) {
    td <- tempfile("asv_zip_")
    dir.create(td, showWarnings = FALSE, recursive = TRUE)
    found <- FALSE
    on.exit(if (!found) unlink(td, recursive = TRUE, force = TRUE))

    utils::unzip(source_path, exdir = td)

    root <- find_dataset_root(td)
    found <- TRUE
    attr(root, "cleanup_dir") <- td
    return(root)
  }

  if (dir.exists(source_path)) {
    root <- find_dataset_root(source_path)
    attr(root, "cleanup_dir") <- NULL
    return(root)
  }

  stop("Source is neither a ZIP nor a directory: ", source_path)
}

# The dataset id comes from the source path as supplied by the user: the dataset
# root may be a nested folder inside it, which would give the wrong id.
make_id <- function(source_path) {
  if (grepl("\\.zip$", source_path, ignore.case = TRUE)) {
    return(gsub("\\.zip$", "", basename(source_path), ignore.case = TRUE))
  }
  basename(source_path)
}

validate_ids <- function(ids) {
  bad <- ids[!grepl("^[A-Za-z0-9_\\-]+$", ids)]
  if (length(bad) > 0) {
    stop(
      "Invalid dataset id(s) detected: ",
      paste0("'", bad, "'", collapse = ", "),
      "\nPlease rename ZIP(s)/folder(s) to only use [A-Za-z0-9_-]."
    )
  }
}

list_sources <- function(data_path) {
  if (!dir.exists(data_path)) {
    stop(paste("Path of dataset directory:", data_path, "was not found."))
  }

  zip_sources <- list.files(data_path, pattern = "\\.zip$", full.names = TRUE)

  dir_sources <- list.dirs(data_path, recursive = FALSE, full.names = TRUE)
  dir_sources <- dir_sources[dir_sources != data_path]

  # Keep only directories that contain occurrence.tsv somewhere underneath (quick filter)
  dir_sources <- dir_sources[
    vapply(dir_sources, function(d) {
      length(list.files(d, pattern = "^occurrence\\.tsv$", recursive = TRUE)) > 0
    }, logical(1))
  ]

  # ---- De-duplicate: prefer ZIP archives over unzipped folders with same datasetID ----
  zip_ids <- gsub("\\.zip$", "", basename(zip_sources), ignore.case = TRUE)
  dir_ids <- basename(dir_sources)
  duplicate_ids <- intersect(zip_ids, dir_ids)

  if (length(duplicate_ids) > 0) {
    message(
      "Both ZIP archive(s) and unzipped folder(s) found for dataset(s): ",
      paste(duplicate_ids, collapse = ", "),
      ". Using ZIP archive(s) and ignoring folder(s)."
    )
    dir_sources <- dir_sources[!(dir_ids %in% duplicate_ids)]
    dir_ids <- basename(dir_sources) # refresh
  }

  # ---- Warn if reading any unzipped dataset folders ----
  if (length(dir_sources) > 0) {
    warning(
      "Reading unzipped dataset folder(s). Recommended workflow is to use ",
      "the original ZIP archives downloaded from the ASV portal.",
      call. = FALSE
    )
  }

  sources <- c(zip_sources, dir_sources)
  if (length(sources) == 0) {
    stop("No ZIP files or dataset directories found in ", data_path, ".")
  }
  sources
}
