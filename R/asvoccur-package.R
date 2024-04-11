# Pull info from DESCRIPTION to asvoccur-package.Rd
#' @keywords internal
"_PACKAGE"

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
#' @importFrom data.table `:=`
#' @importFrom data.table `.SD`
#' @importFrom data.table data.table
#' @importFrom data.table dcast
#' @importFrom data.table fread
#' @importFrom data.table setkey
#' @importFrom data.table setnames
#' @importFrom stats setNames
## usethis namespace: end
NULL

# To avoid `No visible binding` Notes from R CMD check,
# add any unquoted variables used by data.table functions below
# See 'non-standard evaluation (NSE)' & e.g.
# https://www.r-bloggers.com/no-visible-binding-for-global-variable/
.onLoad <- function(libname, pkgname){
  utils::globalVariables(c("..merge_cols", "taxonID"))  
}
