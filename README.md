# asvoccur
Tools for ASV occurrence data in [SBDI](https://biodiversitydata.se/).

## Overview
The **asvoccur** R package, currently under development, processes ASV occurrence data and metadata downloaded from [the Swedish ASV portal](http://asv-portal.biodiversitydata.se/). It converts condensed DwC archives into ASV table format for easier downstream analysis in R, providing functions to load, merge, and aggregate ASV counts across taxonomic ranks.

## Install
```R
install.packages('devtools')
library(devtools)
install_github("biodiversitydata-se/asvoccur")
# or:
# install_github("biodiversitydata-se/asvoccur@develop")
library(asvoccur)

```
## Run
```R
data_path <- '/path-to/dataset-folder'
loaded <- load_data(data_path)
merged <- merge_data(loaded)
summed <- sum_by_clade(merged$counts, merged$asvs)
summed_df <- convert_to_df(summed)
```

## Get help
```R
?asvoccur::load_data
?asvoccur::merge_data
?asvoccur::sum_by_clade
?asvoccur::convert_to_df
```
