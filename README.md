# asvoccur
R tools for ASV occurrence data in [SBDI](https://biodiversitydata.se/).

See [Releases](https://github.com/biodiversitydata-se/asvoccur/releases) for the changelog.

## Overview
The **asvoccur** R package, currently under development, provides tools for unpacking and processing ASV occurrence data and metadata downloaded from [the Swedish ASV portal](http://asv-portal.biodiversitydata.se/). It offers two ways to work with the data:

- **Long format** (`load_long()`, `merge_long()`, `sum_clades_long()`): keeps the data in the long format it is delivered in, with one row per ASV and sample, stores it as Parquet and queries it with DuckDB without reading it into memory. Suited for large datasets, and for analyses with [ecolongr](https://github.com/biodiversitydata-se/ecolongr).
- **Wide format** (`load_data()`, `merge_data()`, `sum_by_clade()`): builds ASV tables with one column per sample, as sparse matrices.

## Tutorial
See our [tutorial](https://biodiversitydata-se.github.io/sbdi-eDNA-smn/tutorial.html) for example workflows with ASV portal data in R: one track uses the wide format, and one uses the long format with ecolongr.

## Install
```R
install.packages('remotes')
remotes::install_github("biodiversitydata-se/asvoccur")
# or:
# remotes::install_github("biodiversitydata-se/asvoccur@develop")
library(asvoccur)
```
The long-format functions also need packages that are not installed with asvoccur:
```R
install.packages(c("duckplyr", "duckdb", "DBI"))
```
GitHub installations (local source builds) may occasionally fail. If so:
```R
remove.packages("asvoccur")
# restart R / RStudio
remotes::install_github("biodiversitydata-se/asvoccur", force = TRUE)
```

## Run
Download the datasets' ZIP archives from the ASV portal and put them in one folder.

### Long format
```R
library(asvoccur)
library(duckplyr)

loaded <- load_long('/path-to/dataset-folder', store = '/path-to/store')
merged <- merge_long(loaded)
orders <- sum_clades_long(merged$counts, merged$asvs, ranks = 'order')
as.data.frame(head(orders))
```
- `load_long()` converts the downloads to Parquet files in `store` the first time, and later only opens them. A store can be prepared once and shared between users.
- `merge_long()` combines the datasets with the same checks as `merge_data()`.
- `counts`, `asvs` and `emof` are lazy tables: use `dplyr` verbs on them, and `as.data.frame()` to bring a small result into R. `events` is an ordinary data frame.

The lazy tables can go straight into [ecolongr](https://github.com/biodiversitydata-se/ecolongr), so nothing is read into R until the distances are computed:
```R
library(ecolongr)

d <- merged$counts |>
  el_decostand('total', sample = 'eventID', taxon = 'taxonID', quant = 'organismQuantity') |>
  el_vegdist('bray', samples = merged$events$eventID,
             sample = 'eventID', taxon = 'taxonID', quant = 'organismQuantity')

# d is a dist object: use vegan on it.
ordination <- vegan::metaMDS(d)
```

### Wide format
```R
data_path <- '/path-to/dataset-folder'
loaded <- load_data(data_path)
merged <- merge_data(loaded)
summed <- sum_by_clade(merged$counts, merged$asvs)
summed_df <- convert_to_df(summed)
```

### Parquet tables
`asv_to_parquet()` writes the Parquet tables that `load_long()` uses: one directory per table (`occurrence`, `asv`, `event`, `emof`, `dataset`), partitioned by dataset. It can also be used on its own, for example from DuckDB with `read_parquet('store/occurrence/*/*.parquet', hive_partitioning = true)`.

`name_convert = 'ecolongr'` names the keys `sample`, `taxon` and `quant` and writes all other columns in snake_case, so the tables can be passed to ecolongr without renaming. `load_long()` needs the portal's names, so use a separate directory for such tables.

## Get help
```R
?asvoccur::load_long
?asvoccur::merge_long
?asvoccur::sum_clades_long
?asvoccur::asv_to_parquet
?asvoccur::load_data
?asvoccur::merge_data
?asvoccur::sum_by_clade
?asvoccur::convert_to_df
```
