#' Load downloaded ASV occurrence data
#'
#' Load Amplicon Sequence Variant (ASV) occurrence data from 'Darwin
#' Core (DwC)-like' archives downloaded from the Swedish ASP portal,
#' \url{https://asv-portal.biodiversitydata.se/}.
#' @param data_path Path of directory containing dataset (*.zip) files
#' @return A list of three sublists (\code{counts}, \code{asvs}, \code{emof})
#'   containing data.table elements from each dataset, indexed by
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
#'   \item \strong{emof}: List of dt:s representing contextual parameter values
#'     (\code{measurementValue}) in a \code{measurementType}
#'     x \code{eventID} matrix from each dataset.
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
  build_counts <- function(zip) {
    occurrence <- fread(cmd = paste('unzip -p', zip, 'occurrence.tsv'))
    counts <- dcast(occurrence, taxonID ~ eventID,
                       value.var = "organismQuantity", fill = 0)
    setkey(counts, taxonID)
    # Set counts to integer
    counts[, names(counts)[-1] := lapply(.SD, as.integer), .SDcols = -1]
    return(counts)
  }
  
  # Reads ASV sequence and taxonomy from asv.tsv
  get_asvs <- function(zip) {
    asv <- fread(cmd = paste('unzip -p', zip, 'asv.tsv'))
    setkey(asv, taxonID)
    return(asv)
  }
  
  # Reads & reshapes emof.tsv into wide format
  # (measurementType [measurementUnit] x eventID)
  # and drops remaining fields, e.g.measurementMethod & measurementRemarks!
  build_emof <- function(zip) {
    emof <- fread(cmd = paste('unzip -p', zip, 'emof.tsv'))
    # Handle datasets that have no contextual data
    if (nrow(emof) == 0) {
      message("Adding empty emof table for ", gsub(".zip", "", zip))
      return(data.table("measurementType (measurementUnit)" = character()))
    }
    # Convert all cols to char, to not add unwanted decimals during dcast
    emof[, names(emof) := lapply(.SD, as.character)]

    emof <- dcast(emof, paste0(measurementType, " (", measurementUnit, ")")
                       ~ eventID, value.var = "measurementValue")
    setnames(emof, "measurementType", "measurementType (measurementUnit)")
    return(emof)
  }
  
  # Process data into data.table:s in (sub)lists, and return in parent list
  loaded <- list()
  loaded$counts <- setNames(lapply(zip_files, build_counts), dirs)
  loaded$asvs <- setNames(lapply(zip_files, get_asvs), dirs)
  loaded$emof <- setNames(lapply(zip_files, build_emof), dirs)
  return(loaded)
}
#' Merge data from different ASV occurrence datasets 
#' 
#' Merge data from different datasets previously loaded with
#' \code{load_data()} function.
#' @param loaded A list of three sublists (\code{counts}, \code{asvs},
#' \code{emof}) from \code{load_data()} function.
#' @return A list of three data.table elements (\code{counts}, \code{asvs},
#' \code{emof}) containing data merged from loaded datasets.
#' @usage 
#' merge_data(loaded);
#' @details 
#' Takes the output from \code{load_data()} and merges rows from
#' different ASV occurrence datasets into three data.table objects:
#' 
#' \enumerate{
#'   \item counts: Read counts in a \code{taxonID} [row] x \code{eventID} [col] matrix.
#'   \item asvs: \code{asv_sequence}, and taxonomy columns per \code{taxonID}.
#'   \item emof: Contextual parameter values (\code{measurementValue}) in a
#'     \code{measurementType} (\code{measurementUnit}) x \code{eventID} matrix.
#' }
#' 
#' Merged counts and emof tables include the UNION of unique rows
#' (i.e. ASVs or measurements) and events (i.e. sequenced samples)
#' from loaded datasets. Table asvs includes non-dataset-specific data, only, 
#' so that we get a single row per unique ASV. Resulting tables are
#' returned as elements of list.
#' 
#' To access an individual table:
#' \describe{
#'   \item{\code{loaded <- load_data(data_path = './datasets')}}{}
#'   \item{\code{merged <- merge_data(loaded)}}{}
#'   \item{\code{View(merged$counts)}}{}
#' }
#' @export
merge_data <- function(loaded) {
  merged <- list()
  merged$counts <- Reduce(function(x, y)
    merge(x, y, by = "taxonID", all = TRUE), loaded$counts)
  merged$emof <- Reduce(function(x, y)
    merge(x, y, by = "measurementType (measurementUnit)", all = TRUE), 
    loaded$emof)
  # We want 1 row/ASV, so only merge non-dataset-specific cols here
  merge_cols <- c("taxonID", "asv_sequence", "scientificName", "taxonRank",
                  "kingdom" ,"phylum", "order", "class", "family", "genus",
                  "specificEpithet", "infraspecificEpithet", "otu",
                  "identificationReferences", "identificationRemarks")
  merged$asvs <- Reduce(function(x, y) {
    merge(x[, ..merge_cols], y[, ..merge_cols], by = merge_cols, all = TRUE)
  }, loaded$asvs)
  return(merged)
}