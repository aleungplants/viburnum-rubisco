library(dplyr)

# get taxa, then manually use figtree to label

here::i_am("csubst/analysis/get_convergent_clades.R")

tree <- ape::read.tree(here::here("csubst", "radseq_pruned_csubst.tre")) %>%
  TreeTools::Preorder()
sites <- readr::read_csv(here::here("csubst", "analysis", "csubst_site_merged.csv")) %>%
  filter(OCNany2spe > 0.9) %>%
  mutate(codon_site_alignment = as.numeric(codon_site_alignment) + 31)

branches <- sites %>%
  dplyr::select(branch_id_1, branch_id_2)
# 
# if (dir.exists(here::here("analysis/clade_trees")) == FALSE) {
#   dir.create("analysis/clade_trees")
# }
# 
# purrr::pwalk(branches, ~ print(c(.x, .y)))
# 
# purrr::pwalk(branches, ~ file.copy(from = here::here("radseq_pruned_csubst.tre"), 
#                                   to = here::here("analysis", "clade_trees", paste("radseq_pruned_", .x, "_", .y, ".tre", sep = ""))))

nodes <- sites %>%
  dplyr::select(branch_id_1, branch_id_2, node1, node2, codon_site_alignment, clade1, clade2, aa_1, aa_2, aa_1_anc, aa_2_anc) %>%
  rowwise() %>%
  mutate(taxa1 = ifelse(stringr::str_detect(node1, stringr::regex("[:digit:]")), list(TreeTools::Subtree(tree, as.numeric(node1))$tip.label), list(node1)),
         taxa2 = ifelse(stringr::str_detect(node2, stringr::regex("[:digit:]")), list(TreeTools::Subtree(tree, as.numeric(node2))$tip.label), list(node2))) %>%
  mutate(taxa1 = paste(taxa1, collapse = " "),
         taxa2 = paste(taxa2, collapse = " "))  %>%
  arrange(branch_id_1, branch_id_2)
readr::write_csv(nodes, here::here("csubst", "analysis", "convergent_clades.csv"))

nodes_subset <- nodes %>%
  filter(stringr::str_detect(node1, stringr::regex("[:digit:]"))) %>%
  filter(stringr::str_detect(node2, stringr::regex("[:digit:]")))
readr::write_csv(nodes_subset, here::here("csubst", "analysis", "convergent_clades_subset.csv"))

# make xlsx table

nodes_clean <- nodes %>%
  select(!contains(c("branch"))) %>%
  rename(site = codon_site_alignment) %>%
  mutate(substitution1 = stringr::str_c(aa_1_anc, site, aa_1),
         substitution2 = stringr::str_c(aa_2_anc, site, aa_2),
         .after = "site") %>%
  select(!contains("aa"))
nodes_clean1 <- nodes_clean %>%
  select(site, contains("1")) %>%
  rename_with(~ stringr::str_remove(., "1"))
nodes_clean2 <- nodes_clean %>%
  select(site, contains("2")) %>%
  rename_with(~ stringr::str_remove(., "2"))

full_specific_epithet <- readr::read_csv(here::here("download", "metadata.csv")) %>%
  select(species) %>%
  mutate(taxa = stringr::str_sub(species, 1, 10), .before = "species") 
landis_biomes <- readr::read_csv(here::here("landis_biome_range", "landis_biomes.csv")) %>%
  rename_with(~stringr::str_to_lower(.)) %>%
  mutate(species = stringr::str_replace(species, "_", "-"))

nodes_biome <- bind_rows(nodes_clean1, nodes_clean2) %>%
  distinct() %>%
  mutate(clade = stringr::str_to_title(clade),
         taxa = stringr::str_split(taxa, " "), 
         .keep = "unused") %>%
  tidyr::unnest(taxa) %>%
  left_join(full_specific_epithet) %>%
  mutate(node = ifelse(!stringr::str_detect(node, "[:digit:]"), species, node)) %>%
  select(!taxa) %>%
  left_join(landis_biomes, relationship = "many-to-many")

nodes_table <- nodes_biome %>%
  group_by(node, site) %>% 
  arrange(site, substitution, node, biome) %>%
  tidyr::nest(species_biome = c(species, biome)) %>%
  rowwise() %>%
  mutate(species = purrr::map(select(species_biome, species), ~stringr::str_c(., collapse = ", ")) %>% as.character(),
         biome = purrr::map(select(species_biome, biome), ~stringr::str_c(., collapse = ", ")) %>% as.character(),
         .keep = "unused") %>%
  relocate(clade, .after = "species") %>%
  arrange(site, substitution, node) %>%
  ungroup() %>%
  rename_with(~stringr::str_to_title(.))

nodes_table_unfill <- nodes_table %>%
  mutate(Site = case_when(Site == lag(Site) ~ NA,
                          .default = Site)) 

writexl::write_xlsx(nodes_table, here::here("analyze", "convergent_clades.xlsx"))
writexl::write_xlsx(nodes_table_unfill, here::here("analyze", "convergent_clades_unfill.xlsx"))



