#' Load downloaded ASV occurrence data
#'
#' Load Amplicon Sequence Variant (ASV) occurrence data from 'Darwin
#' Core (DwC)-like' archives downloaded from the Swedish ASP portal,
#' \url{https://asv-portal.biodiversitydata.se/}.
#' @param data_path Path of directory containing dataset (*.zip) files
#' @return A list of four sublists (\code{counts}, \code{asvs}, \code{events},
#'   \code{emof}) containing data table elements from each dataset, indexed by
#'   \code{datasetID}.
#' @usage load_data(data_path = './datasets');
#' @details Reads data from one or more compressed archives. Returns a list of
#'   sub lists, each of which contains data table objects (dt:s) from each
#'   included dataset:
#'
#' \itemize{
#'   \item \strong{counts}: List of dt:s representing read counts
#'     (\code{taxonID} [row] x \code{eventID} [col] matrix) from each dataset.
#'     \itemize{
#'       \item \code{`first-datasetID`} (dt)
#'       \item \code{`second-datasetID`} (dt)
#'     }
#'
#'   \item \strong{asvs}: List of dt:s representing `DNA_sequence` and taxonomy
#'     per \code{taxonID} from each dataset.
#'     \itemize{
#'       \item ...
#'     }
#'
#'  \item \strong{event}: List of dt:s representing basic event/sample metadata
#'     in a \code{eventID} [col] x parameter [row] matrix from each dataset.
#'     \itemize{
#'       \item ...
#'     }
#'
#'   \item \strong{emof}: List of dt:s representing additional contextual
#'     parameter values (\code{measurementValue}) in a \code{eventID} [row]
#'     x \code{measurementType} [col] matrix from each dataset.
#'     \itemize{
#'       \item ...
#'     }
#' }
#'
#' To access individual dataset tables:
#' \describe{
#'   \item{\code{loaded <- load_data(data_path = './datasets')}}{}
#'   \item{\code{View(loaded$counts$`first-datasetID`)}}{}
#'   \item{\code{# OR:}}{}
#'   \item{\code{View(loaded$counts[[1]])}}{}
#' }
#' @export
load_data <- function(data_path = './datasets') {
  
  # Locate datasets
  if (!dir.exists(data_path)) {
    stop(paste("Path of dataset directory:", data_path, "was not found."))
  }
  zip_files <- list.files(data_path, pattern = "\\.zip$", full.names = TRUE)
  if (length(zip_files) == 0) {
    stop(paste("No ZIP files found in", data_path, "."))
  }
  ds_ids <- gsub(".zip", "", basename(zip_files))
  # Detect e.g. '<datasetID> copy.zip' or '<datasetID> (1).zip'
  for (id in ds_ids)
    if (!grepl("^[A-Za-z0-9_\\-]+$", id))
      stop(paste("Invalid filename detected:", paste0("'",id,".zip'\n"),
                 "Please resolve before proceeding."))
  
  # Reads & reshapes occurrence.tsv into wide (taxonID x eventID) format
  get_counts <- function(zip) {
    occurrences <- fread(cmd = paste('unzip -p', zip, 'occurrence.tsv'))
    counts <- dcast(occurrences, taxonID ~ eventID,
                       value.var = "organismQuantity", fill = 0)
    setkey(counts, taxonID)
    # Set counts to integer
    counts[, names(counts)[-1] := lapply(.SD, as.integer), .SDcols = -1]
    return(counts)
  }
  
  # Reads ASV sequence and taxonomy from asv.tsv
  get_asvs <- function(zip) {
    asvs <- fread(cmd = paste('unzip -p', zip, 'asv.tsv'))
    asvs[, dataset_pid := NULL] # Col for admin use only
    # Replace "" with NA in taxonomy
    tax_cols <- c("kingdom", "phylum", "class", "order", "family", "genus",
                  "specificEpithet", "infraspecificEpithet", "otu")       
    asvs[, (tax_cols) := lapply(.SD, function(x) ifelse(x == "", NA, x)),
         .SDcols = tax_cols]
    setkey(asvs, taxonID)
    return(asvs)
  }
  
  # Reads and reshapes events.tsv
  get_events <- function(zip) {
    events <- fread(cmd = paste('unzip -p', zip, 'event.tsv'))
    events[, c("dataset_pid", "datasetName", "ipt_resource_id") := NULL]
    setkey(events, eventID)
    return(events)
  }
  
  get_datasets <- function(zip) {
    ds_cols <- c("eventID", "datasetName")
    datasets <- fread(cmd = paste('unzip -p', zip, 'event.tsv'), 
                      select = ds_cols)
    datasets[, datasetID := strsplit(eventID, ":")[[1]][1]]
    datasets[, eventID := NULL]
    setcolorder(datasets, c("datasetID", "datasetName"))
    datasets <- unique(datasets)
    return(datasets)
  }
  
  # Reads & reshapes emof.tsv
  # [eventID x measurementType (measurementUnit)]
  # and drops remaining fields, e.g.measurementMethod & measurementRemarks!
  get_emof <- function(zip) {
    emof <- fread(cmd = paste('unzip -p', zip, 'emof.tsv'))
    event_ids <- fread(cmd = paste('unzip -p', zip, 'event.tsv'), select = "eventID")
    # Handle datasets that have no contextual data
    if (nrow(emof) == 0) {
      message("Adding empty emof table for ", gsub(".zip", "", zip))
      emof <- data.table(eventID = event_ids$eventID)
    } else {
      emof <- dcast(emof, 
                    eventID ~ paste0(measurementType, " (", measurementUnit, ")"), 
                    value.var = "measurementValue")
      # Include events without data, if any
      emof <- merge(event_ids, emof, by = "eventID", all.x = TRUE)
    }
    setkey(emof, eventID)
    return(emof)
  }
  
  # Process data into data tables in (sub)lists, and return in parent list
  loaded <- list()
  loaded$counts <- setNames(lapply(zip_files, get_counts), ds_ids)
  loaded$asvs <- setNames(lapply(zip_files, get_asvs), ds_ids)
  loaded$events <- setNames(lapply(zip_files, get_events), ds_ids)
  loaded$datasets <- setNames(lapply(zip_files, get_datasets), ds_ids)
  loaded$emof <- setNames(lapply(zip_files, get_emof), ds_ids)
  return(loaded)
}


# Internal functions to check that input to downstream functions is 
# data table-based and matches the level of complexity accepted by the function.
get_input_category <- function(input) {
  is_dt <- function(input) {  # E.g. 'merged$asvs'
    return(inherits(input, "data.table"))}
  is_dt_lst <- function(input) {  # E.g. 'loaded$asvs', or 'merged'
    return(is.list(input) && all(sapply(input, is_dt)))}
  is_p_lst <- function(input) {  # 'Parent list' of sub lists E.g. 'loaded'
    return(is.list(input) && all(sapply(input, function(y) is_dt_lst(y))))}
  if (is_dt(input)) return('dt')
  if (is_dt_lst(input)) return('dt_lst')
  if (is_p_lst(input)) return('p_lst')
  return(FALSE)  # Non-dt-based input
}
check_input_category <- function(input, lowest_cat) {
  actual_cat <- get_input_category(input)
  errors <- list(
    dt = "Input must be a data table or a (possibly hierachical) list of data tables",
    dt_lst = "Input must be a (possibly hierachical) list of data tables",
    p_lst = "Input must be a hierarchical list of data tables."
  )
  if (lowest_cat == "dt") {
    if (actual_cat %in% c("dt", "dt_lst", "p_lst")) {
      return(TRUE) } else { stop(errors$dt) }}
  if (lowest_cat == "dt_lst") {
    if (actual_cat %in% c("dt_lst", "p_lst")) {
      return(TRUE) } else { stop(errors$dt_lst) }}
  if (lowest_cat == "p_lst") {
    if (actual_cat == "p_lst") {
      return(TRUE) } else { stop(errors$p_lst) }}
}

#' Merge data from different ASV occurrence datasets 
#' 
#' Merge data from different datasets previously loaded with
#' \code{\link[=load_data]{load_data()}} function.
#' @param loaded A multidimensional list of ASV occurrence data table
#' elements loaded with \code{\link[=load_data]{load_data()}}.
#' @param ds An optional character vector specifying the datasets
#'   to merge. If excluded, all datasets will be merged.
#' @return A list of data table elements (\code{counts}, \code{asvs},
#' \code{events}, \code{emof}) containing data merged from loaded datasets.
#' @usage merge_data(loaded, ds = NULL)
#' @details 
#' Takes the output from \code{\link[=load_data]{load_data()}} and merges data from
#' different ASV occurrence datasets into four data table objects:
#' 
#' \enumerate{
#'   \item \strong{counts}: Read counts in a \code{taxonID} [row] x 
#'      \code{eventID} [col] matrix.
#'   \item \strong{asvs}: \code{DNA_sequence}, and taxonomy columns per
#'     \code{taxonID}.
#'   \item \strong{events}: Basic event metadata in a \code{eventID} [row] x
#'     parameter [col] matrix.
#'   \item \strong{emof}: Additional contextual parameter values
#'     (\code{measurementValue}) in a \code{eventID} [row]
#'     x \code{measurementType} [col] matrix.
#' }
#' 
#' Merged tables include the UNION of unique rows and events (i.e. sequenced
#' samples) from loaded datasets. Table asvs includes non-dataset-specific data,
#' only, so that we get a single row per unique ASV. Resulting tables are
#' returned as elements of a list.
#' 
#' To access an individual table:
#' \describe{
#'   \item{\code{loaded <- load_data(data_path = './datasets')}}{}
#'   \item{\code{merged <- merge_data(loaded, ds = c(`first-datasetID`, `second-datasetID`))}}{}
#'   \item{\code{View(merged$counts)}}{}
#' }
#' @export
merge_data <- function(loaded, ds = NULL) {
  
  check_input_category(loaded, 'p_lst')
  
  # Check that specified data sets, if any, have actually been loaded
  if (is.null(ds)) {
    ds <- names(loaded$asvs)
  } else if (!is.vector(ds) || 
             length(ds) > length(intersect(ds, names(loaded$asvs)))) {
    invalid_names <- ds[!(ds %in% names(loaded$asvs))]
    if (length(invalid_names) > 0) {
      stop(paste("Invalid dataset ID(s) specified in ds:", 
                 paste(invalid_names, collapse = ", "),
                 " Don't include '.zip' in ID:s."))
    }
  }
  
  merged <- list()
  
  merged$counts <- Reduce(function(x, y)
    merge(x, y, by = "taxonID", all = TRUE), loaded$counts[ds])
  
  # Checks for overlaps in eventID between datasets to be merged
  detect_event_duplicates <- function(x, y) {
    duplicates <- intersect(x$eventID, y$eventID)
    if (length(duplicates) > 0) { 
      msg <- paste("Duplicated eventID(s) found in datasets to be merged:\n",
                   paste(duplicates, collapse = ", "), "\n",
                   "Did you accidentally try to merge multiple copies of",  
                   "the same dataset? Please resolve before proceeding.")
      stop(msg)
    }
  }
  # Reapplies numeric data type to cols (after merging as char)
  restore_numeric <- function(dt){
    dt[, names(dt) := lapply(.SD, function(col) {
      num_value <- suppressWarnings(as.numeric(col))
      ifelse(is.na(num_value), col, num_value)
    })]
  }
  
  merged$events <- restore_numeric(Reduce(function(x, y){
    detect_event_duplicates(x, y)
    # Convert all cols to char to ensure compatibility during merge
    x <- x[, lapply(.SD, as.character)]
    y <- y[, lapply(.SD, as.character)]
    rbindlist(list(x, y), use.names = TRUE, fill = TRUE)
    }, loaded$events[ds]))
  setkey(merged$events, eventID)
  
  merged$datasets <- restore_numeric(Reduce(function(x, y){
    rbindlist(list(x, y), use.names = TRUE, fill = TRUE)
  }, loaded$datasets[ds]))
  setkey(merged$datasets, datasetID)
  
  merged$emof <- restore_numeric(Reduce(function(x, y){
    detect_event_duplicates(x, y)
    # See events
    x <- x[, lapply(.SD, as.character)]
    y <- y[, lapply(.SD, as.character)]
    rbindlist(list(x, y), use.names = TRUE, fill = TRUE)
    }, loaded$emof[ds]))
  setkey(merged$emof, eventID)
  
  # We want 1 row/ASV, so only merge non-dataset-specific cols here
  merge_cols <- c("taxonID", "DNA_sequence", "scientificName", "taxonRank",
                  "kingdom" ,"phylum", "class", "order", "family", "genus",
                  "specificEpithet", "infraspecificEpithet", "otu",
                  "identificationReferences", "identificationRemarks")
  merged$asvs <- Reduce(function(x, y)
    merge(x[, ..merge_cols], y[, ..merge_cols], by = merge_cols, all = TRUE), 
    loaded$asvs[ds])
  setkey(merged$asvs, taxonID)
  
  # Identifies duplicated ASVs with inconsistent taxonomy across datasets,
  # likely due to merging datasets downloaded at different times
  # (pre- and post-reannotation).
  get_inconsistent_asvs <- function(merged_asvs, loaded_lst) {
    ids <- merged_asvs$taxonID[duplicated(merged_asvs$taxonID)]
    if (length(ids) == 0) return(NULL) # If no duplicates, stop here
    iasv_lst <- lapply(loaded_lst, function(dt) {
      dt[dt$taxonID %in% ids, c("taxonID", "scientificName"), drop = FALSE]
    })
    merged_iasvs <- rbindlist(iasv_lst, idcol = "datasetID")
    msg <- paste("Inconsistent ASV taxonomy detected. This can occur when",
                 "merging datasets downloaded at different times, i.e.",
                 "pre- and post-reannotation. Check 'View(merged$iasvs)' for",
                 "details, and resolve before proceeding with analysis.\n")
    warning(msg)
    return(merged_iasvs)
  }
  merged$iasvs <- get_inconsistent_asvs(merged$asvs, loaded$asvs)
  
  return(merged)
}

#' Sum counts by clade within taxonomic ranks
#'
#' Sum raw and normalized ASV counts by distinct clades at different taxonomic
#' ranks, for each sample in a \code{\link[=load_data]{loaded}} or
#' \code{\link[=merge_data]{merged}} ASV occurrence dataset.
#'
#' @param counts A data table with ASV read counts in a taxonID [row] x eventID
#'   [col] matrix, from a \code{\link[=load_data]{loaded}} or
#'   \code{\link[=merge_data]{merged}} ASV occurrence dataset.
#' @param asvs A data table containing the taxonomic classification of ASV:s
#'   included in the \code{counts} data table.
#' @return A list containing two sub-lists: `raw` and `norm`, each including
#'   data tables for summed ASV counts at each taxonomic rank.
#' @usage sum_by_clade(counts, asvs)
#' @details Sums raw and normalized read counts across ASVs within distinct
#'   clades, at specified taxonomic ranks, to provide higher-level views of the
#'   data. The function normalizes ASV read counts by total counts per sample,
#'   and returns a list of two sub list, each of which contains data tables
#'   (dt:s) of summed counts for each taxonomic rank:
#'
#' \itemize{
#'   \item \strong{raw}: List of dt:s showing raw read counts summed by 
#'   clade at different taxonomic ranks.
#'     \itemize{
#'       \item \code{kingdom} (dt)
#'       \item ...
#'       \item \code{species} (dt)
#'     }
#'   \item \strong{norm}: List of dt:s representing normalized read counts 
#'   summed by clade at different taxonomic ranks.
#'     \itemize{
#'       \item \code{kingdom} (dt)
#'       \item ...
#'       \item \code{species} (dt)
#'     }
#' }
#' #' To view an individual table:
#' \describe{
#'   \item{\code{loaded <- load_data(data_path = './datasets')}}{}
#'   \item{\code{merged <- merge_data(loaded)}}{}
#'   \item{\code{summed <- sum_by_clade(merged$counts, merged$asvs)}}{}
#'   \item{\code{View(summed$raw$family)}}{}
#' }
#' @export
sum_by_clade <- function(counts, asvs){
  
  check_input_category(counts, 'dt')
  check_input_category(asvs, 'dt')
  
  raw_counts <- counts
  count_cols <- names(counts)[-1]  # Drop taxonID
  norm_counts <- copy(raw_counts)
  norm_counts[, (count_cols) := lapply(.SD, function(x) x/sum(x, na.rm = TRUE)),
              .SDcols = count_cols]
  
  tax_cols <- c('taxonID', 'kingdom', 'phylum', 'class', 'order', 'family',
                'genus', 'specificEpithet', 'otu')
  taxa <- asvs[, ..tax_cols]
  taxa[, species := ifelse(is.na(specificEpithet), NA,
                           paste(genus, specificEpithet))]
  taxa[, specificEpithet := NULL]
  setcolorder(taxa, c(setdiff(names(taxa), 'otu'), 'otu'))  # Move otu last
  raw_counts <- merge(taxa, raw_counts, by = "taxonID")
  norm_counts <- merge(taxa, norm_counts, by = "taxonID")
  
  clade_sums <- list()
  ranks <- names(taxa)[-1]
  for (rank in ranks) {
    clade_sums$raw[[rank]] <-
      raw_counts[, lapply(.SD, sum, na.rm = TRUE),
                 by = setNames(list(get(rank)), rank), .SDcols = count_cols]
    clade_sums$norm[[rank]] <-
      norm_counts[, lapply(.SD, sum, na.rm = TRUE),
                  by = setNames(list(get(rank)), rank), .SDcols = count_cols]
  }
  return(clade_sums)
}

#' Convert data table(s) to data frame(s) with rownames
#'
#' Convert data table(s) into data frame(s), also transforming the first 
#' (assumed ID) field into row names. Can be used to e.g. convert a single 
#' \code{\link[=sum_by_clade]{summed}} data table, or all data table elements in
#' \code{\link[=load_data]{loaded}} or \code{\link[=merge_data]{merged}} lists.
#'
#' @param dt_obj A data table or (possibly hierarchical) list of data tables
#' where each dt has a unique ID in its first column.
#' @return A corresponding data frame or (possibly list of) list(s) of data
#' frames, in which the unique ID column has been transformed into row names.
#' @usage convert_to_df(dt_obj)
#' @details Converts data table(s) into data frame(s), also transforming the
#' first (assumed ID) field into row names. ASVs that lack an ID, i.e. a name
#' at some taxonomic rank, will be grouped into a row named 'Unnamed-clades' at
#' that level.
#'
#' Example usage:
#' \describe{
#'   \item{\code{loaded <- load_data(data_path = './datasets')}}{}
#'   \item{\code{merged <- merge_data(loaded)}}{}
#'   \item{\code{merged_df <- convert_to_df(merged)}}{}
#' }
#' @export
convert_to_df <- function(dt_obj) {
  
  check_input_category(dt_obj, 'dt')
  
  # Converts single dt to df, and first dt col to df rownames
  dt_to_df <- function(dt) {
    df <- as.data.frame(dt)
    df[[1]][is.na(df[[1]])] <- "Unnamed-clades"  # Make usable as rowname
    rownames(df) <- df[[1]]
    df[[1]] <- NULL
    return(df)
  }
  
  input_cat <- get_input_category(dt_obj)
  if (input_cat == 'dt') { return(dt_to_df(dt_obj)) }
  else if (input_cat == 'dt_lst') { return(lapply(dt_obj, dt_to_df)) }
  return(lapply(dt_obj, function(sub) { lapply(sub, dt_to_df) }))
}



