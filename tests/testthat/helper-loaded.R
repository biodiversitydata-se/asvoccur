# Hand-built objects shaped like the output of load_data(), with small numbers
# that can be checked by hand.

sparse <- function(rows, cols, values) {
  m <- matrix(values, nrow = length(rows), dimnames = list(rows, cols))
  Matrix::Matrix(m, sparse = TRUE)
}

asv_table <- function(taxon, scientific = paste("Taxon", taxon), genus = "Gen",
                      specific = NA_character_, phylum = "Phy") {
  dt <- data.table::data.table(
    taxonID = taxon, DNA_sequence = paste0("ACGT", taxon), scientificName = scientific,
    taxonRank = "genus", kingdom = "Bacteria", phylum = phylum, order = "Ord", class = "Cla",
    family = "Fam", genus = genus, specificEpithet = specific, infraspecificEpithet = NA_character_,
    otu = NA_character_, identificationReferences = "ref", identificationRemarks = NA_character_
  )
  data.table::setkey(dt, taxonID)
  dt
}

event_table <- function(events, ...) {
  dt <- data.table::data.table(eventID = events, ...)
  data.table::setkey(dt, eventID)
  dt
}

dataset_table <- function(id) data.table::data.table(datasetID = id, datasetName = paste("Name of", id))

# DS-A: ASVs t1, t2 in events a1, a2. DS-B: t2, t3 in b1, b2. t2 is in both.
two_loaded <- function(t2_in_b = "Taxon t2") {
  list(
    counts = list(
      "DS-A" = sparse(c("t1", "t2"), c("a1", "a2"), c(1, 2, 3, 4)),
      "DS-B" = sparse(c("t2", "t3"), c("b1", "b2"), c(5, 0, 0, 6))
    ),
    asvs = list(
      "DS-A" = rbind(asv_table("t1"), asv_table("t2")),
      "DS-B" = rbind(asv_table("t2", scientific = t2_in_b), asv_table("t3"))
    ),
    events = list(
      "DS-A" = event_table(c("a1", "a2"), decimalLatitude = c(55.5, 56), locationID = c("1", "2"), note = c("x", "y")),
      "DS-B" = event_table(c("b1", "b2"), decimalLatitude = c(60, 61.25), locationID = c("3", "4"), note = c("z", NA))
    ),
    datasets = list("DS-A" = dataset_table("DS-A"), "DS-B" = dataset_table("DS-B")),
    emof = list(
      "DS-A" = event_table(c("a1", "a2"), "salinity (psu)" = c(1.5, 2.5)),
      "DS-B" = event_table(c("b1", "b2"), "temperature (C)" = c(7, 8))
    )
  )
}
