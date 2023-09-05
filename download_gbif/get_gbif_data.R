library(tidyverse)
library(rgbif)
library(CoordinateCleaner)

setwd(here::here())

# List of species from GenBank download
seq_species <- read.csv("../download/metadata.csv") %>% # file with accession, gb_species, species (specific epithet, from Landis et al. 2021 Syst Bio supplemental)
  mutate(species = paste("Viburnum", species)) %>% # add genus name 
  dplyr::select(species)

# Find taxon keys for each species in the list
raw_taxon_keys <- map_dfr(seq_species$species, ~ name_backbone(name = ., verbose = TRUE)) %>% # name_backbone on each species and bind the rows
  filter(!is.na(species)) %>%
  filter(confidence >= 90) %>%
  distinct()
if (dir.exists("supplemental") == FALSE) {
  dir.create("supplemental")
}

# Split into exact and fuzzy matches, and manually check the fuzzy matches
exact_taxon_keys <- raw_taxon_keys %>%
  filter(matchType == "EXACT")
FUZZY_SPECIES_KEEP <- c( # Manually check names for reasonable matches, e.g., off by one or two letters
  "Viburnum rafinesqueanum Schult.",
  "Viburnum cotinifolium D.Don",
  "Viburnum betulaefolium Ward, 1885",
  "Viburnum colebrookianum Wall.",
  "Viburnum awafuki hort.",
  "Viburnum awafuki hort. ex Dippel",
  "Viburnum corylifolium Hook.fil. & Thomson",
  "Viburnum burjaeticum Regel & Herder",
  "Viburnum clemensae Kern",
  "Viburnum annamensis Fukuoka"
)
fuzzy_taxon_keys <- raw_taxon_keys %>%
  filter(matchType == "FUZZY") %>%
  filter(scientificName %in% FUZZY_SPECIES_KEEP) %>%
  mutate(verbatim_name = ifelse(verbatim_name == "Viburnum corylifolium", "Viburnum cotinifolium", 
                                ifelse(verbatim_name == "Viburnum cotinifolium", "Viburnum corylifolium", verbatim_name)))

# Save taxon keys and species keys (subset of taxon keys)
taxon_keys <- bind_rows(exact_taxon_keys, fuzzy_taxon_keys) # merge data frames with exact and fuzzy matches
write_csv(taxon_keys, "supplemental/taxon_keys.csv")
species_keys <- taxon_keys %>%
  pull(speciesKey)
write_csv(tibble(species_keys), "supplemental/species_keys.csv")

# Download occurrences
gbif_download <- occ_download(
  type = "and", 
  pred_in("taxonKey", species_keys),
  pred("hasCoordinate", TRUE),
  pred("hasGeospatialIssue", FALSE),
  pred_not(pred_in("issue",c("CONTINENT_INVALID",
                             "COUNTRY_COORDINATE_MISMATCH",
                             "COUNTRY_INVALID",
                             "COUNTRY_MISMATCH",
                             "PRESUMED_NEGATED_LATITUDE",
                             "PRESUMED_NEGATED_LONGITUDE",
                             "PRESUMED_SWAPPED_COORDINATE",
                             "TAXON_MATCH_FUZZY",
                             "COORDINATE_PRECISION_INVALID"))),
  pred_gt("decimalLatitude", -60),
  pred_lte("decimalLatitude",90),
  pred_gte("decimalLongitude",-180),
  pred_lte("decimalLongitude",180),
  pred("occurrenceStatus", "PRESENT"),
  format = "SIMPLE_CSV")
occ_download_wait(gbif_download) # wait for download to finish
gbif_download <- occ_download_get(gbif_download, path = getwd(), overwrite = TRUE) # save download
gbif_citation <- gbif_citation(gbif_download)
writeLines(gbif_citation$download, "supplemental/GBIF_citation.txt")

gbif_data <- occ_download_import(gbif_download) # load occurrence data

# Write GBIF data into csv
if (dir.exists("gbif_data") == FALSE) {
  dir.create("gbif_data")
}
write_csv(gbif_data, "gbif_data/raw.csv")

# Run tests to clean the GBIF data
clean_gbif_data <- gbif_data %>% 
  cc_dupl(., lon = "decimalLongitude", lat = "decimalLatitude") %>% # excludes duplicates
  cc_zero(., lon = "decimalLongitude", lat = "decimalLatitude") %>% # excludes zero coordinates
  cc_equ(., lon = "decimalLongitude", lat = "decimalLatitude") %>% # excludes points with identical lat/lon
  cc_sea(., lon = "decimalLongitude", lat = "decimalLatitude", ref = buffland, scale = 10, value = "clean", verbose = TRUE) %>% #excludes oceanic points
  cc_inst(., lon = "decimalLongitude", lat = "decimalLatitude", buffer = 500) %>% # excludes points near biodiversity institutions
  cc_cap(., lon = "decimalLongitude", lat = "decimalLatitude", buffer = 2500) %>% # excludes points near country capitals
  cc_cen(., lon = "decimalLongitude", lat = "decimalLatitude", buffer = 500) %>% # excludes points near country centroids
#  cc_outl(., lon = "decimalLongitude", lat = "decimalLatitude", method = "quantile", mltpl = 5) # excludes geographic outliers
  right_join(x = ., y = dplyr::select(taxon_keys, speciesKey, verbatim_name), by = "speciesKey") # add verbatim_name (the name we searched by) back in
write_csv(clean_gbif_data, "gbif_data/clean.csv")

# Review # of records per species
record_count <- clean_gbif_data %>%
  group_by(verbatim_name) %>%
  distinct(decimalLongitude, decimalLatitude, .keep_all = TRUE) %>%
  summarize(RecordCount = n())
write_csv(record_count, "supplemental/species_sample_counts.csv")
