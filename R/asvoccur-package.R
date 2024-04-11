#' @title Tools for ASV occurrence data in SBDI
#' @name asvoccur
#' @description This package provides functions for handling ASV occurrence data in SBDI.
#' @importFrom data.table `:=`
#' @importFrom data.table `.SD`
#' @importFrom data.table data.table
#' @importFrom data.table dcast
#' @importFrom data.table fread
#' @importFrom data.table setkey
#' @importFrom data.table setnames
#' @importFrom stats setNames
NULL

# To avoid `No visible binding` Notes from R CMD check.
# Add any unquoted variables used by data.table functions 
# See 'non-standard evaluation (NSE)' & e.g.
# https://www.r-bloggers.com/no-visible-binding-for-global-variable/
.onLoad <- function(libname, pkgname){
  utils::globalVariables(c("..merge_cols", "taxonID"))  
}
