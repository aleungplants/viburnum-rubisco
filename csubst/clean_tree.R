library(ape)
library(ggtree)
library(dplyr)

print("Reading tree.")
tree <- read.tree(here::here("radseq_pruned_csubst.tre"))
print("Writing clean tree.")
write.tree(tree, file = here::here("radseq_pruned_csubst.tre"))
