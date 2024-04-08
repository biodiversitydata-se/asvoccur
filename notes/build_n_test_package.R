# install.packages("devtools")
# install.packages("roxygen2")

# Load package into session for local testing
setwd(dir = '/some-path-to/asvoccur')
library(devtools)
load_all(".")

# Build help files
# https://cran.r-project.org/web/packages/roxygen2/vignettes/roxygen2.html
library(roxygen2)
roxygenise()

# Test help files
?asvoccur::load_data
help(asvoccur::load_data)

# Test functions
data_path <- '/some-path-to/datasets'
loaded <- load_data(data_path)
View(loaded$emof_tables[[1]])



