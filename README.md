# asvoccur
R tools for ASV occurrence data in [SBDI](https://biodiversitydata.se/).

## News
- **March 2025**: Switched counts table to sparse matrices to handle large datasets from the Insect Biome Atlas project, as previous implementation using `data.table` caused RAM exhaustion. Sparse matrices are currently applied only to counts, but further use for other objects (asvs, events, emof) may follow. This change is available in the **develop** branch.

## Overview
The **asvoccur** R package, currently under development, provides tools for unpacking and processing ASV occurrence data and metadata downloaded from [the Swedish ASV portal](http://asv-portal.biodiversitydata.se/). It enables users to convert condensed DwC archives into ASV table format for easier downstream analysis in R, by using functions that load, merge, and aggregate ASV counts across taxonomic ranks.

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
