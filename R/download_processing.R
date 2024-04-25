#' Load downloaded ASV occurrence data
#'
#' Load Amplicon Sequence Variant (ASV) occurrence data from 'Darwin
#' Core (DwC)-like' archives downloaded from the Swedish ASP portal,
#' \url{https://asv-portal.biodiversitydata.se/}.
#' @param data_path Path of directory containing dataset (*.zip) files
#' @return A list of four sublists (\code{counts}, \code{asvs}, \code{events},
#'   \code{emof}) containing data.table elements from each dataset, indexed by
#'   \code{datasetID}.
#' @usage load_data(data_path = './datasets');
#' @details Reads data from one or more compressed archives. Returns a list of
#'   sub lists, each of which contains data.table objects (dt:s) from each
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
#'   \item \strong{asvs}: List of dt:s representing `asv_sequence` and taxonomy
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
  dirs <- gsub(".zip", "", basename(zip_files))
  
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
    setkey(asvs, taxonID)
    return(asvs)
  }
  
  # Reads and reshapes events.tsv
  get_events <- function(zip) {
    events <- fread(cmd = paste('unzip -p', zip, 'event.tsv')) # param x event
    events[, dataset_pid := NULL] # Col for admin use only
    setcolorder(events, c(setdiff(names(events), # Move last (& eventID first)
                                  "ipt_resource_id"), "ipt_resource_id"))
    return(events)
  }
  
  # Reads & reshapes emof.tsv
  # [eventID x measurementType (measurementUnit)]
  # and drops remaining fields, e.g.measurementMethod & measurementRemarks!
  get_emof <- function(zip) {
    emof <- fread(cmd = paste('unzip -p', zip, 'emof.tsv'))
    # Handle datasets that have no contextual data
    if (nrow(emof) == 0) {
      message("Adding empty emof table for ", gsub(".zip", "", zip))
      return(data.table("measurementType (measurementUnit)" = character()))
    }
    emof <- dcast(emof, 
                  eventID ~ paste0(measurementType, " (", measurementUnit, ")"), 
                  value.var = "measurementValue")
    return(emof)
  }
  
  # Process data into data.table:s in (sub)lists, and return in parent list
  loaded <- list()
  loaded$counts <- setNames(lapply(zip_files, get_counts), dirs)
  loaded$asvs <- setNames(lapply(zip_files, get_asvs), dirs)
  loaded$events <- setNames(lapply(zip_files, get_events), dirs)
  loaded$emof <- setNames(lapply(zip_files, get_emof), dirs)
  return(loaded)
}
#' Merge data from different ASV occurrence datasets 
#' 
#' Merge data from different datasets previously loaded with
#' \code{\link[=load_data]{load_data()}} function.
#' @param loaded A multidimensional list of ASV occurrence data.table
#' elements loaded with \code{\link[=load_data]{load_data()}}.
#' @param ds An optional character vector specifying the datasets
#'   to merge. If excluded, all datasets will be merged.
#' @return A list of data.table elements (\code{counts}, \code{asvs},
#' \code{events}, \code{emof}) containing data merged from loaded datasets.
#' @usage merge_data(loaded, ds = NULL)
#' @details 
#' Takes the output from \code{\link[=load_data]{load_data()}} and merges data from
#' different ASV occurrence datasets into four data.table objects:
#' 
#' \enumerate{
#'   \item \strong{counts}: Read counts in a \code{taxonID} [row] x 
#'      \code{eventID} [col] matrix.
#'   \item \strong{asvs}: \code{asv_sequence}, and taxonomy columns per
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
#'   \item{\code{merged <- merge_data(loaded)}}{}
#'   \item{\code{View(merged$counts)}}{}
#' }
#' @export
merge_data <- function(loaded, ds = NULL) {
  if (is.null(ds)) {
    ds <- names(loaded$asvs)
  } else if (!is.vector(ds) || 
             length(ds) > length(intersect(ds, names(loaded$asvs)))) {
    invalid_names <- ds[!(ds %in% names(loaded$asvs))]
    if (length(invalid_names) > 0) {
      stop(paste("Invalid dataset name(s) specified in ds:", 
                 paste(invalid_names, collapse = ", ")))
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
  
  merged$emof <- restore_numeric(Reduce(function(x, y){
    detect_event_duplicates(x, y)
    # See events
    x <- x[, lapply(.SD, as.character)]
    y <- y[, lapply(.SD, as.character)]
    rbindlist(list(x, y), use.names = TRUE, fill = TRUE)
    }, loaded$emof[ds]))
  
  # We want 1 row/ASV, so only merge non-dataset-specific cols here
  merge_cols <- c("taxonID", "asv_sequence", "scientificName", "taxonRank",
                  "kingdom" ,"phylum", "order", "class", "family", "genus",
                  "specificEpithet", "infraspecificEpithet", "otu",
                  "identificationReferences", "identificationRemarks")
  merged$asvs <- Reduce(function(x, y)
    merge(x[, ..merge_cols], y[, ..merge_cols], by = merge_cols, all = TRUE), 
    loaded$asvs[ds])
  
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