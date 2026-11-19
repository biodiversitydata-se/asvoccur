write_tsv <- function(df, path) {
  utils::write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

asv_ids <- function(n) sprintf("ASV:%032x", seq_len(n))

# Writes the portal's files for one dataset into `root` and returns the dataset ID.
# ASVs with the same number have the same taxonID in every dataset, as in the portal.
write_dataset_files <- function(root, id, n_events = 3, asv_from = 1, n_asvs = 5,
                                empty_emof = FALSE, drop = character(), extra_event = NULL) {
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  events <- sprintf("%s:S%d", id, seq_len(n_events))
  nums <- seq(asv_from, length.out = n_asvs)
  taxa <- asv_ids(max(nums))[nums]

  event <- data.frame(
    dataset_pid = 1, ipt_resource_id = "sbdi-asv-1", eventID = events,
    datasetName = paste("Name of", id), eventDate = "2013-07-13T07:08",
    verbatimLocality = c("Arkona Basin", "Bothnian Bay", "Gulf of Finland")[(seq_len(n_events) - 1) %% 3 + 1],
    locationID = "L1", decimalLatitude = 55.185, coordinateUncertaintyInMeters = NA,
    materialSampleID = "https://example.org/sample", references = "ref",
    pcr_primer_name_forward = "515F",
    sop = "\"\"", remarks = "\"Niskin, 5 m", stringsAsFactors = FALSE
  )
  for (nm in names(extra_event)) event[[nm]] <- extra_event[[nm]]
  occurrence <- expand.grid(taxonID = taxa, eventID = events, stringsAsFactors = FALSE)
  occurrence <- occurrence[(seq_len(nrow(occurrence)) %% 3) != 0, ]
  occurrence$organismQuantity <- seq_len(nrow(occurrence)) * 7L
  occurrence <- data.frame(dataset_pid = 1, occurrence[, c("eventID", "taxonID", "organismQuantity")])
  asv <- data.frame(
    dataset_pid = 1, taxonID = taxa, asv_id_alias = sub("ASV:", "", taxa),
    DNA_sequence = paste0("ACGT", nums), scientificName = paste("Taxon", nums),
    taxonRank = "genus", kingdom = "Bacteria", phylum = "", class = "Gamma", order = "",
    family = "Fam", genus = "Gen", specificEpithet = "", infraspecificEpithet = "", otu = "",
    dateIdentified = "2026-01-01", identificationReferences = "ref", associatedSequences = "",
    previousIdentifications = "", identificationRemarks = "", stringsAsFactors = FALSE
  )
  emof <- data.frame(
    dataset_pid = 1, eventID = rep(events, 2),
    measurementType = rep(c("salinity", "temperature"), each = n_events),
    measurementTypeID = "", measurementUnit = rep(c("psu", "C"), each = n_events),
    measurementUnitID = "", measurementValue = seq_len(2 * n_events) + 0.5,
    stringsAsFactors = FALSE
  )
  if (empty_emof) emof <- emof[0, ]

  tables <- list(event = event, occurrence = occurrence, asv = asv, emof = emof)
  for (nm in setdiff(names(tables), drop)) write_tsv(tables[[nm]], file.path(root, paste0(nm, ".tsv")))
  writeLines("<eml/>", file.path(root, "eml.xml"))
  writeLines("README", file.path(root, "README.txt"))
  invisible(id)
}

# `as = "zip"` writes <id>.zip in `data_path`, `as = "dir"` an unzipped folder <id>.
# `wrapper` puts the files in a top-level folder and adds macOS metadata, as a re-packed zip has.
make_dataset <- function(data_path, id, as = c("zip", "dir"), wrapper = FALSE, ...) {
  as <- match.arg(as)
  dir.create(data_path, showWarnings = FALSE, recursive = TRUE)
  stage <- withr::local_tempdir(.local_envir = parent.frame())
  inner <- if (wrapper) file.path(stage, id, "inner") else file.path(stage, id)
  write_dataset_files(inner, id, ...)
  if (wrapper) {
    dir.create(file.path(stage, "__MACOSX", id, "inner"), recursive = TRUE)
    writeLines("junk", file.path(stage, "__MACOSX", id, "inner", "._occurrence.tsv"))
  }
  if (as == "zip") {
    top <- intersect(c(id, "__MACOSX"), list.files(stage, all.files = TRUE, no.. = TRUE))
    zip::zipr(file.path(data_path, paste0(id, ".zip")), files = top, root = stage)
  } else {
    file.copy(file.path(stage, id), data_path, recursive = TRUE)
  }
  invisible(id)
}

read_pq <- function(out_dir, table) {
  con <- DBI::dbConnect(duckdb::duckdb(shared_home = FALSE))
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE))
  DBI::dbGetQuery(con, sprintf(
    "SELECT * FROM read_parquet('%s/%s/*/*.parquet', hive_partitioning = true)", out_dir, table
  ))
}

triplets <- function(m) {
  s <- Matrix::summary(m)
  out <- data.frame(
    eventID = colnames(m)[s$j], taxonID = rownames(m)[s$i],
    organismQuantity = as.numeric(s$x), stringsAsFactors = FALSE
  )
  out[order(out$eventID, out$taxonID), ]
}

sorted <- function(df, by) {
  df <- df[do.call(order, unname(as.list(df[by]))), , drop = FALSE]
  rownames(df) <- NULL
  df
}

all_tables <- c("occurrence", "event", "asv", "emof", "dataset")
