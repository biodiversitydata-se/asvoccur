#' Load downloaded ASV occurrence data as lazy long tables
#'
#' Open datasets downloaded from the Swedish ASV portal,
#' \url{https://asv-portal.biodiversitydata.se/}, as lazy long tables, without
#' building a taxa-by-samples matrix. The data are kept in Parquet files written by
#' \code{\link{asv_to_parquet}()} and queried with DuckDB, so they are not read
#' into memory. Combine the datasets with \code{\link{merge_long}()}.
#'
#' @param data_path Path to a directory with datasets downloaded from the portal,
#'   as for \code{\link{load_data}()}, or \code{NULL} to only open \code{store}.
#' @param store Directory with the Parquet tables. If any dataset in
#'   \code{data_path} is missing from \code{store}, all datasets in
#'   \code{data_path} are converted with \code{\link{asv_to_parquet}()} first.
#'   If all of them are in \code{store}, the stored copies are used and a message
#'   says so. A store can be prepared once and shared between users.
#' @param overwrite Logical. If \code{TRUE}, convert the datasets in
#'   \code{data_path} again even if they are in \code{store}. Use this after
#'   downloading a new version of a stored dataset, e.g. after a re-annotation.
#'
#' @return A list of lazy tables with one row per value and a \code{datasetID}
#'   column: \code{counts} (\code{eventID}, \code{taxonID},
#'   \code{organismQuantity}), \code{asvs}, \code{events} and \code{emof}, plus
#'   \code{datasets}, a data frame of dataset IDs and names, and \code{store}.
#'   \code{asvs} has a \code{species} column, which is genus and epithet
#'   combined, or \code{NA} if either is missing.
#'
#' @details
#' Requires the suggested packages \pkg{duckdb}, \pkg{duckplyr}, \pkg{DBI},
#' \pkg{dplyr} and \pkg{rlang}.
#'
#' The lazy tables work with \pkg{duckplyr} verbs, which run in DuckDB. They stop
#' with an error rather than read a whole table into R when a verb cannot be
#' translated. Use \code{as.data.frame()} to bring a small result into memory.
#'
#' @examples
#' \dontrun{
#' loaded <- load_long("./datasets", store = "./store")
#' merged <- merge_long(loaded)
#' }
#' @export
load_long <- function(data_path = NULL, store, overwrite = FALSE) {
  check_packages(long_packages)
  if (!is.null(data_path)) {
    sources <- list_sources(data_path)
    ids <- vapply(sources, make_id, character(1), USE.NAMES = FALSE)
    if (overwrite || !all(ids %in% stored_ids(store))) {
      asv_to_parquet(data_path, store)
    } else {
      message(
        "Using the stored copies in '", store, "' of: ", paste(ids, collapse = ", "),
        ". Use overwrite = TRUE to convert a newer download of a dataset."
      )
    }
  }
  if (length(Sys.glob(file.path(store, "*", "dataset_id=*"))) > 0) {
    stop(
      "'", store, "' was written with name_convert = \"ecolongr\". ",
      "load_long() needs the portal's column names: use another store."
    )
  }
  if (length(stored_ids(store)) == 0) {
    stop("No datasets found in '", store, "'.")
  }

  datasets <- as.data.frame(read_store(store, "dataset"))
  datasets <- datasets[order(datasets$datasetID), c("datasetID", "datasetName")]
  rownames(datasets) <- NULL

  list(
    counts = read_store(store, "occurrence"),
    asvs = read_asvs(store),
    events = read_store(store, "event"),
    emof = read_store(store, "emof"),
    datasets = datasets,
    store = store
  )
}

#' Merge long-format ASV occurrence datasets
#'
#' The long-table counterpart of \code{\link{merge_data}()}: combines datasets
#' opened with \code{\link{load_long}()} and checks that they can be combined.
#' Read counts, ASVs and contextual measurements stay lazy tables.
#'
#' @param loaded A list from \code{\link{load_long}()}.
#' @param ds An optional character vector of the datasets to merge. If excluded,
#'   all datasets are merged.
#'
#' @return A list with:
#' \itemize{
#'   \item \strong{counts}: Lazy table of read counts (\code{eventID},
#'     \code{taxonID}, \code{organismQuantity}), without zeros.
#'   \item \strong{asvs}: Lazy table with one row per ASV: sequence, taxonomy and
#'     \code{species}, as in \code{loaded$asvs}.
#'   \item \strong{emof}: Lazy table of contextual measurements, one row per
#'     measurement, with values as text.
#'   \item \strong{events}: Data frame of sample metadata, in memory.
#'   \item \strong{datasets}: Data frame of dataset IDs and names.
#'   \item \strong{iasvs}: Data frame of ASVs whose taxonomy differs between
#'     datasets, or \code{NULL} if there are none.
#' }
#'
#' @details
#' As in \code{\link{merge_data}()}, events found in more than one dataset give an
#' error, and ASVs whose taxonomy differs between datasets, typically because they
#' were downloaded before and after a re-annotation, give a warning. Such ASVs have
#' one row per version in \code{asvs}, which must be resolved before summing by
#' clade.
#'
#' Also as in \code{\link{merge_data}()}, \code{asvs} only has the attributes
#' that describe the ASV itself. Dataset-specific ones, such as
#' \code{asv_id_alias}, \code{associatedSequences} and
#' \code{previousIdentifications}, are in \code{loaded$asvs}, with one row per
#' dataset.
#'
#' @examples
#' \dontrun{
#' merged <- merge_long(load_long("./datasets", store = "./store"))
#' as.data.frame(head(merged$counts))
#' }
#' @export
merge_long <- function(loaded, ds = NULL) {
  check_packages(long_packages)
  col <- rlang::sym

  all_ids <- loaded$datasets$datasetID
  if (is.null(ds)) {
    ds <- all_ids
  } else {
    invalid <- setdiff(ds, all_ids)
    if (length(invalid) > 0) {
      stop(
        "Invalid dataset ID(s) specified in ds: ", paste(invalid, collapse = ", "),
        ". Don't include '.zip' in ID:s."
      )
    }
  }
  keep <- duckplyr::duckdb_tibble(datasetID = ds, .prudence = "stingy")
  in_ds <- function(x) dplyr::semi_join(x, keep, by = "datasetID")

  events <- as.data.frame(in_ds(loaded$events))
  duplicates <- unique(events$eventID[duplicated(events$eventID)])
  if (length(duplicates) > 0) {
    stop(
      "Duplicated eventID(s) found in datasets to be merged:\n",
      paste(duplicates, collapse = ", "), "\n",
      "Did you accidentally try to merge multiple copies of the same dataset? ",
      "Please resolve before proceeding."
    )
  }
  events <- events[order(events$eventID), , drop = FALSE]
  rownames(events) <- NULL

  asv_rows <- in_ds(loaded$asvs)
  # As in merge_data(), only attributes that describe the ASV itself are combined;
  # dataset-specific ones such as asv_id_alias stay in loaded$asvs.
  shared_cols <- intersect(c(asv_shared_cols, "species"), colnames(loaded$asvs))
  asvs <- dplyr::distinct(dplyr::select(asv_rows, dplyr::all_of(shared_cols)))

  conflicts <- dplyr::filter(
    dplyr::summarise(asvs, n = dplyr::n(), .by = "taxonID"),
    !!col("n") > 1
  )
  iasvs <- as.data.frame(dplyr::semi_join(
    dplyr::select(asv_rows, "datasetID", "taxonID", "scientificName"),
    conflicts,
    by = "taxonID"
  ))
  if (nrow(iasvs) > 0) {
    iasvs <- iasvs[order(iasvs$taxonID, iasvs$datasetID), , drop = FALSE]
    rownames(iasvs) <- NULL
    warning(
      "Inconsistent ASV taxonomy detected. This can occur when merging datasets ",
      "downloaded at different times, i.e. pre- and post-reannotation. Check ",
      "'merged$iasvs' for details, and resolve before proceeding with analysis.",
      call. = FALSE
    )
  } else {
    iasvs <- NULL
  }

  list(
    counts = dplyr::select(in_ds(loaded$counts), "eventID", "taxonID", "organismQuantity"),
    asvs = asvs,
    emof = dplyr::select(in_ds(loaded$emof), -"datasetID"),
    events = events,
    datasets = loaded$datasets[loaded$datasets$datasetID %in% ds, , drop = FALSE],
    iasvs = iasvs
  )
}

#' Sum long-format counts by clade within taxonomic ranks
#'
#' The long-table counterpart of \code{\link{sum_by_clade}()}: sums read counts
#' across ASVs within each clade, for each sample, and gives the share of each
#' clade in the sample. The result is a lazy table, so no matrix is built.
#'
#' @param counts Lazy table of read counts from \code{\link{merge_long}()},
#'   possibly filtered to selected samples.
#' @param asvs Lazy table of ASVs from \code{\link{merge_long}()}.
#' @param ranks Character vector of ranks to sum over. Defaults to all of them.
#'
#' @return A lazy table with columns \code{rank}, \code{clade}, \code{eventID},
#'   \code{raw} (summed reads) and \code{norm} (share of the reads in the sample).
#'
#' @details
#' ASVs with no value at a rank are summed as \code{"Unclassified"}, so at species
#' level all ASVs without an epithet form one group, whatever their genus. Labels
#' from the annotation, such as \code{unclassified.Acrididae}, are ordinary names and
#' are kept apart, so those reads can still be separated at genus level.
#'
#' Shares are computed from the table given, so filtering ASVs first changes the
#' totals. To keep the shares of the whole sample, filter the result instead.
#'
#' @examples
#' \dontrun{
#' merged <- merge_long(load_long("./datasets", store = "./store"))
#' orders <- sum_clades_long(merged$counts, merged$asvs, ranks = "order")
#' as.data.frame(orders)
#' }
#' @export
sum_clades_long <- function(counts, asvs, ranks = clade_ranks()) {
  check_packages(long_packages)
  unknown <- setdiff(ranks, clade_ranks())
  if (length(unknown) > 0) {
    stop(
      "Unknown rank(s): ", paste(unknown, collapse = ", "),
      ". Use one or more of: ", paste(clade_ranks(), collapse = ", "), "."
    )
  }

  # duckplyr translates columns given as symbols, but not the .data pronoun.
  col <- rlang::sym

  ids <- as.data.frame(dplyr::summarise(
    asvs,
    n_rows = dplyr::n(),
    n_ids = dplyr::n_distinct(!!col("taxonID"))
  ))
  if (ids$n_rows != ids$n_ids) {
    stop(
      "Some ASVs have more than one row in 'asvs', which would count their reads ",
      "more than once. See the iasvs element from merge_long(), and resolve before ",
      "summing by clade."
    )
  }

  one_rank <- function(rank) {
    taxa <- dplyr::transmute(
      asvs,
      taxonID = !!col("taxonID"),
      clade = dplyr::coalesce(!!col(rank), "Unclassified")
    )
    joined <- dplyr::inner_join(counts, taxa, by = "taxonID")
    summed <- dplyr::summarise(
      joined,
      raw = sum(!!col("organismQuantity"), na.rm = TRUE),
      .by = c("clade", "eventID")
    )
    totals <- dplyr::mutate(summed, total = sum(!!col("raw"), na.rm = TRUE), .by = "eventID")
    dplyr::transmute(
      totals,
      rank = !!rank,
      clade = !!col("clade"),
      eventID = !!col("eventID"),
      raw = !!col("raw"),
      norm = !!col("raw") / !!col("total")
    )
  }

  Reduce(dplyr::union_all, lapply(ranks, one_rank))
}

# Internal helpers ------------------------------------------------------------

long_packages <- c("duckdb", "duckplyr", "DBI", "dplyr", "rlang")

clade_ranks <- function() {
  c("kingdom", "phylum", "class", "order", "family", "genus", "species", "otu")
}

check_packages <- function(needed) {
  missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      "This function needs these packages: ", paste(missing, collapse = ", "),
      ". Install them with install.packages(c(",
      paste0('"', missing, '"', collapse = ", "), ")).",
      call. = FALSE
    )
  }
}

sql_string <- function(x) {
  paste0("'", gsub("'", "''", x, fixed = TRUE), "'")
}

# The ASV attributes that merge_data() combines across datasets.
asv_shared_cols <- c(
  "taxonID", "DNA_sequence", "scientificName", "taxonRank",
  "kingdom", "phylum", "order", "class", "family", "genus",
  "specificEpithet", "infraspecificEpithet", "otu",
  "identificationReferences", "identificationRemarks"
)

stored_ids <- function(store) {
  dirs <- Sys.glob(file.path(store, "dataset", "datasetID=*"))
  sub("^datasetID=", "", basename(dirs))
}

# duckplyr cannot translate paste(), so species is built in SQL. As in
# sum_by_clade(), an epithet without a genus gives no species.
read_asvs <- function(store) {
  duckplyr::read_sql_duckdb(
    sprintf(
      "SELECT *,
         CASE WHEN genus IS NOT NULL AND specificEpithet IS NOT NULL
           THEN genus || ' ' || specificEpithet END AS species
       FROM read_parquet(%s, hive_partitioning = true, hive_types_autocast = false)",
      sql_string(file.path(store, "asv", "*", "*.parquet"))
    ),
    prudence = "stingy"
  )
}

read_store <- function(store, table) {
  duckplyr::read_parquet_duckdb(
    file.path(store, table, "*", "*.parquet"),
    options = list(hive_partitioning = TRUE, hive_types_autocast = FALSE),
    prudence = "stingy"
  )
}
