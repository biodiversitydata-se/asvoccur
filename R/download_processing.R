#' Read and reshape downloaded ASV occurrence data
#'
#' Read and reshape Amplicon Sequence Variant (ASV) occurrence data from
#' 'Darwin Core (DwC)-like' archives downloaded from the Swedish ASP portal, 
#' \url{https://asv-portal.biodiversitydata.se/}.
#' @param data_path Path of directory containing dataset (*.zip) files
#' @return A list of three sublists (\code{asv_tables}, \code{asvs},
#' \code{emof_tables}) containing data.table elements from each dataset, 
#' indexed by \code{datasetID}.
#' @usage 
#' load_data(data_path = './datasets');
#' @details 
#' The function reads text files (\code{occurrence/asv/emof.tsv}) from compressed
#' archives (\code{[datasetID].zip}), and repackages data into three 
#' data.table objects per dataset:
#' 
#' \enumerate{
#'   \item ASV-table: Read counts in a \code{taxonID} [row] x \code{eventID} [col] matrix.
#'   \item asvs: \code{asv_sequence}, and taxonomy columns per \code{taxonID}.
#'   \item EMOF-table: Contextual parameter values (\code{measurementValue}) in a
#'     \code{measurementType} (\code{measurementUnit}) x \code{eventID} matrix.
#' }
#' 
#' Tables from different datasets are indexed by their respective \code{datasetID},
#' and organized into three lists (\code{asv_tables}, \code{asvs}, \code{emof_tables}).
#' These lists are returned as elements of a parent list.
#' 
#' To access individual dataset tables:
#' \describe{
#'   \item{\code{loaded <- load_data(data_path = './datasets')}}{}
#'   \item{\code{View(loaded$asv_tables$`some-datasetID`)}}{ # OR:}
#'   \item{\code{View(loaded$asvs[[1]])}}{}
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
  build_asv_table <- function(zip) {
    occurrence <- fread(cmd = paste('unzip -p', zip, 'occurrence.tsv'))
    asv_table <- dcast(occurrence, taxonID ~ eventID,
                       value.var = "organismQuantity", fill = 0)
    setkey(asv_table, taxonID)
    # Set counts to integer
    asv_table[, names(asv_table)[-1] := lapply(.SD, as.integer), .SDcols = -1]
    return(asv_table)
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
  build_emof_table <- function(zip) {
    emof <- fread(cmd = paste('unzip -p', zip, 'emof.tsv'))
    # Handle datasets that have no contextual data
    if (nrow(emof) == 0) {
      message("Adding empty emof table for ", gsub(".zip", "", zip))
      return(data.table("measurementType (measurementUnit)" = character()))
    }
    # Convert all cols to char, to not add unwanted decimals during dcast
    emof[, names(emof) := lapply(.SD, as.character)]

    emof_table <- dcast(emof, paste0(measurementType, " (", measurementUnit, ")")
                       ~ eventID, value.var = "measurementValue")
    setnames(emof_table, "measurementType", "measurementType (measurementUnit)")
    return(emof_table)
  }
  
  # Process data into data.table:s in (sub)lists, and return in parent list
  loaded <- list()
  loaded$asv_tables <- setNames(lapply(zip_files, build_asv_table), dirs)
  loaded$asvs <- setNames(lapply(zip_files, get_asvs), dirs)
  loaded$emof_tables <- setNames(lapply(zip_files, build_emof_table), dirs)
  return(loaded)
}
#' Merge data from different ASV occurrence datasets 
#' 
#' Merge ASV/EMOF tables and asvs from different datasets previously loaded with
#' \code{load_data()} function.
#' @param loaded A list of three sublists (\code{asv_tables}, \code{asvs},
#' \code{emof_tables}) from \code{load_data()} function.
#' @return A list of three data.table elements (\code{asv_table}, \code{asvs},
#' \code{emof_table}) containing data merged from loaded datasets.
#' @usage 
#' merge_data(loaded);
#' @details 
#' The function takes the output from \code{load_data()} and merges rows from
#' different ASV occurrence datasets into three data.table objects:
#' 
#' \enumerate{
#'   \item ASV-table: Read counts in a \code{taxonID} [row] x \code{eventID} [col] matrix.
#'   \item asvs: \code{asv_sequence}, and taxonomy columns per \code{taxonID}.
#'   \item EMOF-table: Contextual parameter values (\code{measurementValue}) in a
#'     \code{measurementType} (\code{measurementUnit}) x \code{eventID} matrix.
#' }
#' 
#' Merged ASV/EMOF tables include the UNION of unique rows
#' (i.e. ASVs or measurements) and events (i.e. sequenced samples)
#' from loaded datasets. Table asvs includes non-dataset-specific data, only, 
#' so that we get a single row per unique ASV. Resulting data.table:s are
#' returned as elements of list.
#' 
#' To access an individual table:
#' \describe{
#'   \item{\code{loaded <- load_data(data_path = './datasets')}}{}
#'   \item{\code{merged <- load_data(loaded)}}{}
#'   \item{\code{View(merged$asv_table)}}{}
#' }
#' @export
merge_data <- function(loaded) {
  merged <- list()
  merged$asv_table <- Reduce(function(x, y)
    merge(x, y, by = "taxonID", all = TRUE), loaded$asv_tables)
  merged$emof_table <- Reduce(function(x, y)
    merge(x, y, by = "measurementType (measurementUnit)", all = TRUE), 
    loaded$emof_table)
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