#' Load downloaded ASV occurrence data
#'
#' Load Amplicon Sequence Variant (ASV) occurrence data from 'Darwin
#' Core (DwC)-like' archives downloaded from the Swedish ASP portal,
#' \url{https://asv-portal.biodiversitydata.se/}.
#' @param data_path Path of directory containing dataset (*.zip) files
#' @return A list of four sublists (\code{counts}, \code{asvs}, \code{events},
#'   \code{emof}) containing sparse matrix or data table elements from each 
#'   dataset, indexed by \code{datasetID}.
#' @usage load_data(data_path = './datasets');
#' @details Reads data from one or more compressed archives. Returns a list of
#'   sub lists, each of which contains sparse matrix or data table 
#'   objects from each included dataset:
#'
#' \itemize{
#'   \item \strong{counts}: List of sparse matrices representing read counts
#'     (taxon [row] x event [col] sparse matrix) from each dataset.
#'     \itemize{
#'       \item \code{`first-datasetID`} (sparse matrix)
#'       \item \code{`second-datasetID`} (sparse matrix)
#'       \item ...
#'     }
#'
#'   \item \strong{asvs}: List of data tables containing the DNA sequence and
#'     taxonomic assignment of ASVs (taxon [row] x attribute [col]) from each dataset.
#'
#'     \itemize{
#'       \item \code{`first-datasetID`} (data table)
#'       \item \code{`second-datasetID`} (data table)
#'       \item ...
#'     }
#'
#'   \item \strong{events}: List of data tables representing basic event/sample
#'     metadata (event [row] x parameter [col]) from each dataset.
#'
#'     \itemize{
#'       \item \code{`first-datasetID`} (data table)
#'       \item \code{`second-datasetID`} (data table)
#'       \item ...
#'     }
#'
#'   \item \strong{emof}: List of data tables representing additional contextual
#'     parameter values (event [row] x measurementType [col]) from each dataset.
#'     \itemize{
#'       \item \code{`first-datasetID`} (data table)
#'       \item \code{`second-datasetID`} (data table)
#'       \item ...
#'     }
#' }
#'
#' To inspect individual matrices or data tables:
#' \describe{
#'   \item{\code{loaded <- load_data(data_path = './datasets')}}{}
#'   \item{\code{View(loaded$emof$`first-datasetID`)}}{}
#'   \item{\code{# OR (to show first 100 ASVs in first counts matrix):}}{}
#'   \item{\strong{Memory note:} `counts` is a sparse matrix; converting large
#'     matrices to dense format (e.g. `as.matrix()`) may exhaust RAM. Only convert
#'     small subsets for inspection.}{}
#'   \item{\code{View(as.matrix(loaded$counts[[1]][1:100,]))}}{}
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
  
  # Reads & reshapes occurrence.tsv into sparse (!) and wide (taxonID x eventID) format
  get_counts <- function(zip) {
    occurrences <- fread(utils::unzip(zip, files = "occurrence.tsv", 
                                      exdir = tempdir()))
    occurrences[, taxonID := as.factor(taxonID)]
    occurrences[, eventID := as.factor(eventID)]
    
    counts <- Matrix::sparseMatrix(
      i = as.integer(occurrences$taxonID),
      j = as.integer(occurrences$eventID),
      x = occurrences$organismQuantity,
      dims = c(nlevels(occurrences$taxonID), nlevels(occurrences$eventID)),
      dimnames = list(levels(occurrences$taxonID), levels(occurrences$eventID))
    )
    
    return(counts)
  }
  
  # Reads ASV sequence and taxonomy from asv.tsv
  get_asvs <- function(zip) {
    asvs <- fread(utils::unzip(zip, files = "asv.tsv", exdir = tempfile()))
    asvs[, dataset_pid := NULL] # Col for admin use only
    # Replace "" with NA in taxonomy
    tax_cols <- c("kingdom", "phylum", "order", "class", "family", "genus",
                  "specificEpithet", "infraspecificEpithet", "otu")       
    asvs[, (tax_cols) := lapply(.SD, function(x) ifelse(x == "", NA, x)),
         .SDcols = tax_cols]
    setkey(asvs, taxonID)
    return(asvs)
  }
  
  # Reads and reshapes events.tsv
  get_events <- function(zip) {
    events <- fread(utils::unzip(zip, files = "event.tsv", exdir = tempfile()))
    events[, c("dataset_pid", "datasetName", "ipt_resource_id") := NULL]
    setkey(events, eventID)
    return(events)
  }
  
  get_datasets <- function(zip) {
    ds_cols <- c("eventID", "datasetName")
    datasets <- fread(utils::unzip(zip, files = "event.tsv", 
                                   exdir = tempfile()), select = ds_cols)
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
    emof <- fread(utils::unzip(zip, files = "emof.tsv", exdir = tempfile()))
    event_ids <- fread(utils::unzip(zip, files = 'event.tsv', exdir = tempfile()), 
                       select = "eventID")

    # Handle datasets that have no contextual data
    if (nrow(emof) == 0) { warning("load_data(): Adding empty emof table for ", 
                                   gsub("\\.zip$", "", basename(zip)),
                                   call. = FALSE)
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
# data table- or matrix-based, and matches the level of complexity accepted by 
# the function.
get_input_category <- function(input) {
  
  is_dt <- function(input) {  # E.g. 'merged$asvs'
    return(inherits(input, "data.table"))}
  is_sp_mat <- function(input) {  # E.g. 'loaded$counts'
    return(inherits(input, "sparseMatrix"))}
  # Check if input is a list containing data.tables or sparse matrices
  is_dt_lst <- function(input) {  # E.g. 'loaded$asvs', or 'merged'
    return(is.list(input) && all(sapply(input, function(x) is_dt(x) || is_sp_mat(x))))}
  # Check if input is a parent list (list of lists) containing data.table or sparse matrix
  is_p_lst <- function(input) {  # 'Parent list' of sub lists E.g. 'loaded'
    return(is.list(input) && all(sapply(input, function(y) is_dt_lst(y) || is_sp_mat(y))))}
  
  if (is_dt(input)) return('dt')
  if (is_dt_lst(input)) return('dt_lst')
  if (is_p_lst(input)) return('p_lst')
  if (is_sp_mat(input)) return('sp_mat')
  return(FALSE)  # Non-dt or matrix-based input
}

check_input_category <- function(input, lowest_cat) {
  actual_cat <- get_input_category(input)
  errors <- list(
    dt = "Input must be a data table or a (possibly hierachical) list of data tables and matrices",
    sp_mat = "Input must be a sparse matrix or a (possibly hierachical) list of data tables and matrices",
    dt_lst = "Input must be a (possibly hierachical) list of data tables and matrices",
    p_lst = "Input must be a hierarchical list of data tables and matrices"
  )
  if (lowest_cat == "dt") {
    if (actual_cat %in% c("dt", "dt_lst", "p_lst")) {
      return(TRUE) } else { stop(errors$dt) }}
  if (lowest_cat == "sp_mat") {
    if (actual_cat %in% c("sp_mat", "dt_lst", "p_lst")) {
      return(TRUE) } else { stop(errors$sp_mat) }}
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
#' @param loaded A multidimensional list of ASV occurrence data elements loaded 
#' with \code{\link[=load_data]{load_data()}}.
#' @param ds An optional character vector specifying the datasets
#'   to merge. If excluded, all datasets will be merged.
#' @return A list of sparse matrix or data table elements: (\code{counts}, 
#' \code{asvs}, \code{events}, \code{emof}) containing data merged from loaded
#' datasets.
#' @usage merge_data(loaded, ds = NULL)
#' @details 
#' Takes the output from \code{\link[=load_data]{load_data()}} and merges data from
#' different ASV occurrence datasets into four sparse matrix or data table objects:
#' 
#' \enumerate{
#' 
#'   \item \strong{counts}: Read counts in a \code{taxonID} [row] x 
#'      \code{eventID} [col] sparse matrix.
#'      
#'   \item \strong{asvs}: DNA sequence and taxonomic assignment of ASVs in a  
#'      \code{taxonID} [row] x attribute [col] data table.
#'      
#'   \item \strong{events}: Basic event metadata in an \code{eventID} [row] x
#'     parameter [col] data table.
#'     
#'   \item \strong{emof}: Additional contextual parameter values
#'     (\code{measurementValue}) in an \code{eventID} [row]
#'     x \code{measurementType} [col] data table.
#' }
#' 
#' Merged matrices and tables include the UNION of unique rows and events (i.e. sequenced
#' samples) from loaded datasets. Table asvs includes non-dataset-specific data,
#' only, so that we get a single row per unique ASV. Resulting tables are
#' returned as elements of a list.
#' 
#' To inspect an individual matrix or data table:
#' \describe{
#'   \item{\code{loaded <- load_data(data_path = './datasets')}}{}
#'   \item{\code{merged <- merge_data(loaded, ds = c(`first-datasetID`, `second-datasetID`))}}{}
#'   \item{\code{View(merged$events)}}{}
#'   \item{\code{# OR (to show first 100 ASVs in merged counts matrix):}}{}
#'   \item{\code{View(as.matrix(merged$counts[1:100,]))}}{}
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
  
  # Checks for overlaps in eventID between datasets to be merged
  detect_event_duplicates <- function(x, y, is_sparse = FALSE) {
    event_x <- if (is_sparse) colnames(x) else x$eventID
    event_y <- if (is_sparse) colnames(y) else y$eventID
    duplicates <- intersect(event_x, event_y)
    if (length(duplicates) > 0) { 
      msg <- paste("Duplicated eventID(s) found in datasets to be merged:\n",
                   paste(duplicates, collapse = ", "), "\n",
                   "Did you accidentally try to merge multiple copies of",  
                   "the same dataset? Please resolve before proceeding.")
      stop(msg)
    }
  }
  
  merged <- list()
  
  # Iteratively merges sparse matrices with union of ASVs and events
  merge_sparse_counts <- function(matrices) {
    all_asvs <- unique(unlist(lapply(matrices, rownames)))
    # Add missing ASVs and sort identically across matrices
    padded_matrices <- lapply(matrices, function(mat) {
      missing_asvs <- setdiff(all_asvs, rownames(mat))
      if (length(missing_asvs) > 0) {  # Set missing counts to 0
        empty_mat <- Matrix::Matrix(0, nrow = length(missing_asvs), ncol = ncol(mat),
                            dimnames = list(missing_asvs, colnames(mat)))
        mat <- rbind(mat, empty_mat)
      }
      mat[order(rownames(mat)), , drop = FALSE]
    })
    Reduce(function(x, y) {
      detect_event_duplicates(x, y, is_sparse = TRUE)
      Matrix::cbind2(x, y)
    }, padded_matrices)
  }
  merged$counts <- merge_sparse_counts(loaded$counts[ds])
  
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
                  "kingdom" ,"phylum", "order", "class", "family", "genus",
                  "specificEpithet", "infraspecificEpithet", "otu",
                  "identificationReferences", "identificationRemarks")
  merged$asvs <- Reduce(function(x, y)
    merge(x[, ..merge_cols], y[, ..merge_cols], by = merge_cols, all = TRUE), 
    loaded$asvs[ds])
  setkey(merged$asvs, taxonID)
  
  # Identifies duplicated ASVs with inconsistent taxonomy across datasets,
  # likely due to merging datasets downloaded at different times
  # (pre- and post-annotation update).
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
#' @param counts ASV read counts in a taxon [row] x event [col] sparse matrix, 
#' from a \code{\link[=load_data]{loaded}} or
#'   \code{\link[=merge_data]{merged}} ASV occurrence dataset.
#' @param asvs A data table containing the DNA sequences and taxonomic assignment
#'   of ASV:s included in the \code{counts} matrix.
#' @return A list containing two sub-lists: `raw` and `norm`, each including
#'   data tables for summed ASV counts at each taxonomic rank.
#' @usage sum_by_clade(counts, asvs)
#' @details Sums raw and normalized read counts across ASVs within distinct
#'   clades, at specified taxonomic ranks, to provide higher-level views of the
#'   data. The function normalizes ASV read counts by total counts per sample,
#'   and returns a list of two sub list, each of which contains data tables
#'   of summed counts for each taxonomic rank:
#'
#' \itemize{
#'   \item \strong{raw}: List of data tables showing raw read counts summed by 
#'   clade at different taxonomic ranks.
#'     \itemize{
#'       \item \code{kingdom} (data table)
#'       \item ...
#'       \item \code{species} (data table)
#'     }
#'     
#'   \item \strong{norm}: List of data tables representing normalized read counts 
#'   summed by clade at different taxonomic ranks.
#'     \itemize{
#'       \item \code{kingdom} (data table)
#'       \item ...
#'       \item \code{species} (data table)
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
sum_by_clade <- function(counts, asvs) {
  
  check_input_category(counts, 'sp_mat')
  check_input_category(asvs, 'dt')
  
  if (!identical(rownames(counts), asvs$taxonID)) {
    stop("Mismatch detected: 
         Row names of 'counts' do not match 'taxonID' of 'asvs'. 
         Please ensure they are identical.")
  }
  
  tax_cols <- c('taxonID', 'kingdom', 'phylum', 'class', 'order', 'family', 
                'genus', 'specificEpithet', 'otu')
  taxa <- asvs[, ..tax_cols]
  taxa[, species := ifelse(is.na(specificEpithet), NA, 
                           paste(genus, specificEpithet))]
  taxa[, specificEpithet := NULL]
  setcolorder(taxa, c(setdiff(names(taxa), 'otu'), 'otu'))
  
  clade_sums_raw <- list()
  clade_sums_norm <- list()
  
  for (rank in names(taxa)[-1]) {
    clades <- ifelse(is.na(taxa[[rank]]), "Unclassified", taxa[[rank]])
    
    clade_levels <- unique(clades)
    clade_index <- match(clades, clade_levels)
    
    # Aggregate using matrix multiplication
    G <- Matrix::sparseMatrix(i = clade_index, j = seq_along(clades), x = 1,
                              dims = c(length(clade_levels), length(clades)))
    raw_matrix <- G %*% counts
    rownames(raw_matrix) <- clade_levels
    
    # Normalise
    norm_matrix <- raw_matrix %*% Matrix::Diagonal(x = 1 / Matrix::colSums(raw_matrix))
    colnames(norm_matrix) <- colnames(raw_matrix)
   
     # Sort clades
    clade_order <- order(rownames(raw_matrix))
    raw_matrix <- raw_matrix[clade_order, , drop = FALSE]
    norm_matrix <- norm_matrix[clade_order, , drop = FALSE]
    
    # Convert to dt:s
    clade_sums_raw[[rank]] <- data.table(clade = rownames(raw_matrix), as.matrix(raw_matrix))
    clade_sums_norm[[rank]] <- data.table(clade = rownames(norm_matrix), as.matrix(norm_matrix))
  }
  return(list(raw = clade_sums_raw, norm = clade_sums_norm))
}

#' Convert tabular ASV data to data.frame format
#'
#' Convert \code{data.table} objects to \code{data.frame} format for inspection
#' and compatibility with functions expecting base R data frames. Can be applied
#' to a single object or to the (possibly hierarchical) lists returned by
#' \code{\link[=load_data]{load_data()}} or
#' \code{\link[=merge_data]{merge_data()}}. Optionally attempts to convert sparse
#' count matrices.
#'
#' @param dt_obj A \code{data.table}, a \code{data.frame}, a sparse matrix (e.g.
#'   \code{dgCMatrix}), or a (possibly hierarchical) list containing these.
#' @param convert_counts Logical. If \code{TRUE}, attempts to convert sparse count
#'   matrices to \code{data.frame}. If \code{FALSE} (default), sparse matrices are
#'   left unchanged.
#' @param max_cells Maximum allowed number of cells (\code{nrow * ncol}) when
#'   converting a sparse matrix to a dense \code{data.frame}. If exceeded,
#'   conversion is skipped and a warning is issued.
#'
#' @return A \code{data.frame} or a (possibly hierarchical) list of
#'   \code{data.frame}s, with any sparse matrices that were not converted returned
#'   unchanged.
#'
#' @usage convert_to_df(dt_obj, convert_counts = FALSE, max_cells = 5e6)
#'
#' @details Converts \code{data.table} objects to \code{data.frame} and moves the
#' first column (assumed to be a unique ID) to row names. Missing IDs are replaced
#' with \code{"Unnamed-clades"}.
#'
#' If \code{convert_counts = TRUE}, sparse matrices are converted only when their
#' size does not exceed \code{max_cells}. Otherwise, they are left unchanged to
#' avoid excessive memory use.
#'
#' \strong{Memory note:} Converting sparse count matrices to dense
#' \code{data.frame}s can require large amounts of RAM. Only enable
#' \code{convert_counts} for small datasets, or when the resulting object size is
#' known to be manageable.
#'
#' Example usage:
#' \describe{
#'   \item{\code{loaded <- load_data(data_path = './datasets')}}{}
#'   \item{\code{merged <- merge_data(loaded)}}{}
#'   \item{\code{merged_df <- convert_to_df(merged)}}{}
#'   \item{\code{# Attempt counts conversion (may be skipped if too large):}}{}
#'   \item{\code{merged_df2 <- convert_to_df(merged, convert_counts = TRUE, max_cells = 1e8)}}{}
#' }
#'
#' @export
convert_to_df <- function(dt_obj, convert_counts = FALSE, max_cells = 5e6) {
  is_dt <- function(z) inherits(z, "data.table")
  is_df <- function(z) inherits(z, "data.frame") && !inherits(z, "data.table")
  is_spmat <- function(z) {
    inherits(z, "dgCMatrix") || inherits(z, "sparseMatrix") || inherits(z, "Matrix")
  }
  
  dt_to_df <- function(dt) {
    df <- as.data.frame(dt)
    id <- df[[1]]
    id[is.na(id)] <- "Unnamed-clades"
    rownames(df) <- id
    df[[1]] <- NULL
    df
  }
  
  skipped <- character(0)
  
  contains_convertible <- function(x) {
    if (is_dt(x) || is_spmat(x)) return(TRUE)
    if (is.list(x)) return(any(vapply(x, contains_convertible, logical(1))))
    FALSE
  }
  
  contains_only_df <- function(x) {
    if (is_df(x)) return(TRUE)
    if (is.list(x)) return(all(vapply(x, contains_only_df, logical(1))))
    FALSE
  }
  
  if (!contains_convertible(dt_obj) && contains_only_df(dt_obj)) {
    message("convert_to_df(): Input appears already converted (data.frame); nothing to do.")
    return(dt_obj)
  }
  
  walk <- function(obj, path = "") {
    if (is_dt(obj)) return(dt_to_df(obj))
    
    if (is_spmat(obj)) {
      if (!convert_counts) {
        skipped <<- c(skipped, path)
        return(obj)
      }
      n_cells <- as.numeric(nrow(obj)) * as.numeric(ncol(obj))
      if (n_cells > max_cells) {
        skipped <<- c(skipped, path)
        return(obj)
      }
      df <- as.data.frame(as.matrix(obj))
      rownames(df) <- rownames(obj)
      return(df)
    }
    
    if (is_df(obj)) return(obj)
    
    if (is.list(obj)) {
      nms <- names(obj)
      if (is.null(nms)) nms <- as.character(seq_along(obj))
      return(setNames(lapply(nms, function(nm) {
        walk(obj[[nm]], paste0(path, "$", nm))
      }), nms))
    }
    
    stop(
      sprintf(
        "convert_to_df(): Unsupported input type at %s (class: %s).",
        if (nzchar(path)) path else "<root>",
        paste(class(obj), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
  root <- tryCatch(deparse(substitute(dt_obj)), error = function(...) "dt_obj")
  result <- walk(dt_obj, root)
  
  if (length(skipped)) {
    warning(
      "convert_to_df(): Skipped converting sparse counts at: ",
      paste(sort(unique(skipped)), collapse = ", "),
      ". Counts left as sparse matrices. Use convert_counts = TRUE with larger max_cells or use disk-backed tools.",
      call. = FALSE
    )
  }
  
  result
}

