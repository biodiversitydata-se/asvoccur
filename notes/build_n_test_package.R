# Some notes on how to build and test the package in RStudio
# To keep local version of this file, use:
# git update-index --assume-unchanged notes/build_n_test_package.R
# Undo with:
# git update-index --no-assume-unchanged notes/build_n_test_package.R

# install.packages("devtools")

# Load package into session for local testing
library(devtools)
load_all(".")  # Alt. use ⇧ + ⌘ + L OR Build tab | More | Load All

# Build help files with roxygen2:
# See https://cran.r-project.org/web/packages/roxygen2/vignettes/roxygen2.html
# Alt 1: set to automatically roxygenize when loading package
# Tools | Project Options | Build Tools | Generate doc. with roxygen
# Alt. 2 (if project setting has no effect) use:
# Build tab | More | Document

# View help files
?asvoccur
?asvoccur::load_data
?asvoccur::merge_data
?asvoccur::sum_by_clade

# Test functions
data_path <- '~/Downloads'
loaded <- load_data(data_path)
View(loaded$asvs[[1]])
merged <- merge_data(loaded)
View(merged$asvs)
clade_sums <- sum_by_clade(merged$counts, merged$asvs)
View(clade_sums$norm$class)
