#' Convert downloaded ASV occurrence data to Parquet tables
#'
#' Convert ZIP archives (or unzipped folders) downloaded from the Swedish ASV
#' portal, \url{https://asv-portal.biodiversitydata.se/}, to Parquet tables.
#' The conversion is done by DuckDB straight from the TSV files, so no dataset
#' is read into R and memory use does not grow with the size of a dataset.
#'
#' @param data_path Path to a directory containing one or more datasets, as for
#'   \code{\link{load_data}()}.
#' @param out_dir Directory to write the tables to. It is created if needed.
#' @param name_convert \code{"none"} (default) keeps the portal's column names.
#'   \code{"ecolongr"} names the keys as the \pkg{ecolongr} package expects, i.e.
#'   \code{eventID} to \code{sample}, \code{taxonID} to \code{taxon} and
#'   \code{organismQuantity} to \code{quant}, and writes all other columns, and the
#'   partition column \code{dataset_id}, in snake_case.
#'
#' @return \code{out_dir}, invisibly.
#'
#' @details
#' One directory per table is written to \code{out_dir}, each Hive-partitioned
#' by \code{datasetID}, e.g.
#' \code{out_dir/occurrence/datasetID=KTH-2013-Baltic-16S/data.parquet}:
#'
#' \itemize{
#'   \item \strong{occurrence}: \code{eventID}, \code{taxonID},
#'     \code{organismQuantity}; one row per non-zero count.
#'   \item \strong{asv}: ASV sequence and taxonomic assignment.
#'     Blank taxonomy cells are \code{NULL}.
#'   \item \strong{event}: Event/sample metadata. Coordinates, depths,
#'     elevations, coordinate uncertainty and sample size are numeric, all other
#'     columns text.
#'   \item \strong{emof}: Contextual measurements, in long format, with
#'     \code{measurementValue} as text.
#'   \item \strong{dataset}: \code{datasetName}.
#' }
#'
#' With \code{name_convert = "ecolongr"} the partition directories are
#' \code{dataset_id=...} and the occurrence table can be passed to the
#' \pkg{ecolongr} functions as it is. The names are fixed when the tables are
#' written, so a directory cannot hold tables of both kinds: the function stops
#' if \code{out_dir} was written with the other setting.
#'
#' The column \code{dataset_pid} is not kept, since it is not unique across
#' datasets. The \code{datasetID} is taken from the ZIP or folder name and must
#' equal the prefix of the \code{eventID}s in the dataset.
#'
#' All datasets in \code{out_dir} are read as one table with, for example,
#' \code{read_parquet('out_dir/occurrence/*/*.parquet', hive_partitioning = true)}
#' in DuckDB, or with \code{duckplyr::read_parquet_duckdb()}.
#' Running the function again replaces the datasets found in \code{data_path}
#' and leaves other datasets in \code{out_dir} as they are. Nothing is written to
#' \code{out_dir} if any dataset fails.
#' @export
asv_to_parquet <- function(data_path = "./datasets", out_dir, name_convert = c("none", "ecolongr")) {
  name_convert <- match.arg(name_convert)
  key <- if (name_convert == "ecolongr") "dataset_id" else "datasetID"
  sources <- list_sources(data_path)
  ds_ids <- vapply(sources, make_id, character(1), USE.NAMES = FALSE)
  validate_ids(ds_ids)

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  check_same_convention(out_dir, key, name_convert)
  staging <- file.path(normalizePath(out_dir, winslash = "/"), ".staging")
  unlink(staging, recursive = TRUE)
  on.exit(unlink(staging, recursive = TRUE), add = TRUE)

  con <- DBI::dbConnect(duckdb::duckdb(shared_home = FALSE))
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)

  for (i in seq_along(sources)) convert_source(con, sources[i], ds_ids[i], staging, key, name_convert)

  for (id in ds_ids) {
    for (tbl in parquet_tables) {
      target <- file.path(out_dir, tbl, paste0(key, "=", id))
      unlink(target, recursive = TRUE)
      dir.create(dirname(target), showWarnings = FALSE, recursive = TRUE)
      if (!file.rename(file.path(staging, tbl, basename(target)), target)) {
        stop("Could not move ", tbl, " of dataset ", id, " to ", target)
      }
    }
  }
  invisible(out_dir)
}

parquet_tables <- c("occurrence", "asv", "event", "emof", "dataset")

event_numeric_cols <- c(
  "decimalLatitude", "decimalLongitude", "coordinateUncertaintyInMeters",
  "minimumDepthInMeters", "maximumDepthInMeters",
  "minimumElevationInMeters", "maximumElevationInMeters", "sampleSizeValue"
)

convert_source <- function(con, source, id, staging, key, name_convert) {
  root <- resolve_root(source)
  on.exit(unlink(attr(root, "cleanup_dir"), recursive = TRUE, force = TRUE))

  event_file <- file.path(root, "event.tsv")
  prefixes <- DBI::dbGetQuery(con, sprintf(
    "SELECT DISTINCT split_part(eventID, ':', 1) AS p FROM %s", read_tsv_sql(con, event_file)
  ))$p
  if (!identical(prefixes, id)) {
    stop(
      "The eventIDs in '", basename(source), "' start with ",
      paste0("'", prefixes, "'", collapse = ", "), ", not with the dataset id '", id,
      "'. Do not rename the ZIP or folder."
    )
  }

  rename <- if (name_convert == "ecolongr") ecolongr_name else identity
  queries <- list(
    occurrence = tsv_select(con, file.path(root, "occurrence.tsv"), rename, drop = "dataset_pid", numeric = "organismQuantity"),
    asv = tsv_select(con, file.path(root, "asv.tsv"), rename, drop = "dataset_pid"),
    event = tsv_select(con, event_file, rename, drop = c("dataset_pid", "datasetName", "ipt_resource_id"), numeric = event_numeric_cols),
    emof = tsv_select(con, file.path(root, "emof.tsv"), rename, drop = "dataset_pid"),
    dataset = sprintf(
      "SELECT DISTINCT datasetName AS %s FROM %s",
      DBI::dbQuoteIdentifier(con, rename("datasetName")), read_tsv_sql(con, event_file)
    )
  )
  for (tbl in names(queries)) {
    dir <- file.path(staging, tbl, paste0(key, "=", id))
    dir.create(dir, recursive = TRUE)
    DBI::dbExecute(con, sprintf(
      "COPY (%s) TO %s (FORMAT parquet, COMPRESSION zstd)",
      queries[[tbl]], DBI::dbQuoteString(con, file.path(dir, "data.parquet"))
    ))
  }
}

# All columns are read as text so a column has the same type in every dataset;
# quoting is off because the files hold stray quote characters.
read_tsv_sql <- function(con, path) {
  sprintf(
    "read_csv(%s, delim = '\\t', quote = '', header = true, all_varchar = true)",
    DBI::dbQuoteString(con, path)
  )
}

tsv_select <- function(con, path, rename, drop = character(), numeric = character()) {
  from <- read_tsv_sql(con, path)
  cols <- DBI::dbGetQuery(con, sprintf("DESCRIBE SELECT * FROM %s", from))$column_name
  cols <- setdiff(cols, drop)
  new <- rename(cols)
  clash <- unique(new[duplicated(new)])
  if (length(clash) > 0) {
    stop(
      "In ", basename(path), ", ", paste0("'", cols[new %in% clash], "'", collapse = ", "),
      " would all be named ", paste0("'", clash, "'", collapse = ", "),
      ". Use name_convert = \"none\"."
    )
  }
  q <- vapply(cols, function(col) as.character(DBI::dbQuoteIdentifier(con, col)), character(1))
  qn <- vapply(new, function(col) as.character(DBI::dbQuoteIdentifier(con, col)), character(1))
  exprs <- ifelse(cols %in% numeric, sprintf("CAST(%s AS DOUBLE)", q), q)
  sprintf("SELECT %s FROM %s", paste(exprs, "AS", qn, collapse = ", "), from)
}

ecolongr_keys <- c(eventID = "sample", taxonID = "taxon", organismQuantity = "quant")

ecolongr_name <- function(x) {
  snake <- gsub("([A-Z]+)([A-Z][a-z])", "\\1_\\2", gsub("([a-z0-9])([A-Z])", "\\1_\\2", x))
  ifelse(x %in% names(ecolongr_keys), ecolongr_keys[x], tolower(snake))
}

check_same_convention <- function(out_dir, key, name_convert) {
  other <- setdiff(c("datasetID", "dataset_id"), key)
  found <- Sys.glob(file.path(out_dir, "*", paste0(other, "=*")))
  if (length(found) > 0) {
    stop(
      "'", out_dir, "' holds tables written with the other name_convert (",
      basename(found[1]), "). Use a new out_dir or remove the old tables."
    )
  }
}
