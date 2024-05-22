# asvoccur
Tools for ASV occurrence data in SBDI

## Install
```R
install.packages('devtools')
library(devtools)
install_github("pragermh/asvoccur")
library(asvoccur)
```
## Run
```R
data_path <- '/path-to/dataset-folder'
loaded <- load_data(data_path)
merged <- merge_data(loaded)
clade_sums <- sum_by_clade(merged$counts, merged$asvs)
```

## Get help
```R
?asvoccur::load_data
?asvoccur::merge_data
?asvoccur::sum_by_clade
```
