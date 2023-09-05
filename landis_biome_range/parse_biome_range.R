library(dplyr)

here::i_am("biome_range/parse_biome_range.R")

read_nexus_data <- function(lines) {
  data <- lines %>%
    as_tibble() %>%
    filter(stringr::str_detect(value, "V_")) %>%
    mutate(value = stringr::str_replace(value, stringr::regex("(?<=[:punct:][:alnum:]) "), ",")) %>%
    mutate(value = stringr::str_squish(value)) %>% # removes whitespace at the start and end and replace all internal whitespace with a single space 
    mutate(Species = stringr::word(value, 1),
           Code = stringr::word(value, 2)) %>%
    dplyr::select(!value) %>%
    mutate(Species = stringr::str_remove(Species, "V_"),
           Code = stringr::str_remove_all(Code, stringr::regex("[\\(\\)]"))) %>%
    mutate(Code = stringr::str_split(Code, ",")) %>% # some species have multiple codes
    tidyr::unnest(Code)
  return(data)
}


biome_lines <- readLines(here::here("biome_range", "viburnum.biome.n4.nex"))
biome <- read_nexus_data(biome_lines) %>%
  mutate(Biome = case_when(
    Code == 0 ~ "Tropical",
    Code == 1 ~ "Lucidophyllous",
    Code == 2 ~ "Cloud",
    Code == 3 ~ "Temperate"
  )) %>%
  dplyr::select(!Code)
readr::write_csv(biome, here::here("biome_range", "landis_biomes.csv"))

range_lines <- readLines(here::here("biome_range", "viburnum.range.n6.nex"))
range_coding <- range_lines %>%
  as_tibble() %>%
  filter(stringr::str_detect(value, " :") & stringr::str_detect(value, "[\\[\\]]")) %>%
  mutate(value = stringr::str_remove_all(value, "[\\[\\]]")) %>%
  mutate(value = stringr::str_squish(value)) %>%
  mutate(Code = stringr::word(value, 1, 1, sep = " : "),
         Range = stringr::word(value, 2, 2, sep = " : ")) %>%
  mutate(Range = stringr::str_split(Range, stringr::fixed(" + "))) %>%
  dplyr::select(!value)
range <- read_nexus_data(range_lines) %>%
  left_join(., range_coding) %>%
  tidyr::unnest(Range) %>%
  dplyr::select(!Code)
readr::write_csv(range, here::here("biome_range", "landis_ranges.csv"))

