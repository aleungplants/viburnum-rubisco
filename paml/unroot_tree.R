setwd(here::here())
print("Reading tree.")
tree <- ape::read.tree("radseq_pruned.tre")
print("Writing unrooted tree.")
tree$edge.length <- NULL
print("Removing branch lengths for PAML.")
ape::write.tree(ape::unroot(tree), file = "radseq_pruned_unrooted.tre")
