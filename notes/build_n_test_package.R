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

View(loaded$events$`KTH-2013-Baltic-16S`)
View(loaded$asvs$`KTH-2013-Baltic-16S`)
View(loaded$emof$`KTH-2013-Baltic-16S`)
View(as.matrix(loaded$counts[[1]][1:100,]))

# View help files
?asvoccur
?asvoccur::load_data
?asvoccur::merge_data
?asvoccur::sum_by_clade
?asvoccur::convert_to_df

# Test functions
data_path <- '~/Downloads/16S/'
loaded <- load_data(data_path)
loaded_df <- convert_to_df(loaded)
loaded_df <- convert_to_df(loaded, convert_counts = TRUE, max_cells = 1e8)
rm(loaded_df)
gc()
View(loaded$asvs[[1]])
merged <- merge_data(loaded)
View(merged$asvs)
summed <- sum_by_clade(merged$counts, merged$asvs)
View(summed$norm$class)
summed_df <- convert_to_df(summed)
