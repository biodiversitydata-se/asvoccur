# asvoccur
R tools for ASV occurrence data in [SBDI](https://biodiversitydata.se/).

See [Releases](https://github.com/biodiversitydata-se/asvoccur/releases) for the changelog.

## Overview
The **asvoccur** R package, currently under development, provides tools for unpacking and processing ASV occurrence data and metadata downloaded from [the Swedish ASV portal](http://asv-portal.biodiversitydata.se/). It enables users to convert condensed DwC archives into ASV table format for easier downstream analysis in R, by using functions that load, merge, and aggregate ASV counts across taxonomic ranks.

## Tutorial
See our [tutorial](https://biodiversitydata-se.github.io/sbdi-days-eDNA/tutorial.html) for an example workflow using **asvoccur** with ASV portal data in R.

## Install
```R
install.packages('remotes')
remotes::install_github("biodiversitydata-se/asvoccur")
# or:
# remotes::install_github("biodiversitydata-se/asvoccur@develop")
library(asvoccur)
```
GitHub installations (local source builds) may occasionally fail. If so:
```R
remove.packages("asvoccur")
# restart R / RStudio
remotes::install_github("biodiversitydata-se/asvoccur", force = TRUE)
```

## Run
```R
data_path <- '/path-to/dataset-folder'
loaded <- load_data(data_path)
merged <- merge_data(loaded)
summed <- sum_by_clade(merged$counts, merged$asvs)
summed_df <- convert_to_df(summed)
```

### Large datasets: Parquet tables and ecolongr
`asv_to_parquet()` converts the downloads to Parquet tables with DuckDB, without reading the data into R.
The tables can be read lazily and passed to [ecolongr](https://github.com/biodiversitydata-se/ecolongr), which computes on long tables without making them dense.
`name_convert = 'ecolongr'` names the keys `sample`, `taxon` and `quant`, so no renaming is needed.
All other columns, and the partition column `dataset_id`, are written in snake_case.

```R
library(asvoccur)
library(ecolongr)
library(dplyr)

# One directory per table in parquet/, partitioned by dataset.
asv_to_parquet('/path-to/dataset-folder', out_dir = 'parquet', name_convert = 'ecolongr')

# A lazy table over all datasets: nothing is read into R until el_vegdist().
read_table <- function(table) {
  duckplyr::read_parquet_duckdb(
    file.path('parquet', table, '*', '*.parquet'),
    options = list(hive_partitioning = TRUE), prudence = 'stingy'
  )
}
counts <- read_table('occurrence')
samples <- read_table('event') |> pull(sample)

# Bray-Curtis distances on relative abundances.
# `samples` includes samples without counts and sets the order of d.
d <- counts |>
  el_decostand('total') |>
  el_vegdist('bray', samples = samples)

# d is a dist object: use vegan on it.
ordination <- vegan::metaMDS(d)
plot(ordination, type = 't')
```
Without `name_convert` the tables keep the portal's names, and `el_vegdist()` needs `sample = 'eventID', taxon = 'taxonID', quant = 'organismQuantity'`.
A directory written with one setting cannot be extended with the other: use a new `out_dir`.
In DuckDB the same tables are read with `read_parquet('parquet/occurrence/*/*.parquet', hive_partitioning = true)`.

## Get help
```R
?asvoccur::load_data
?asvoccur::asv_to_parquet
?asvoccur::merge_data
?asvoccur::sum_by_clade
?asvoccur::convert_to_df
```
