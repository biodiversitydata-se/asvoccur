# Some notes on how to build and test the package in RStudio

# install.packages("devtools")
# install.packages("roxygen2")

# Load package into session for local testing
setwd(dir = '/some-path-to/asvoccur')
library(devtools)
load_all(".")  # Alt. use ⇧ + ⌘ + L OR Build tab | More | Load All

# Build help files:
# See https://cran.r-project.org/web/packages/roxygen2/vignettes/roxygen2.html
library(roxygen2)
roxygenise()
# Alt. set to automatically roxygenize when loading package
# Tools | Project Options | Build Tools | Generate doc. with roxygen

# Test help files
?asvoccur::load_data
help(asvoccur::load_data)

# Test functions
data_path <- '/some-path-to/datasets'
loaded <- load_data(data_path)
View(loaded$asv_tables[[1]])

merged <- merge_data(loaded)
View(merged$asv_table)







