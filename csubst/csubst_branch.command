#!/bin/bash

# cd the folder containing this file
cd "$(dirname "$0")"

# get rooted tree file and alignment file
cp ../paml/radseq_pruned.tre radseq_pruned.tre
cp ../align/clean_aligned_checked.fasta clean_aligned_checked.fasta
python clean_fasta_csubst.py # fix length of taxa names (10, because the tree tip labels also has this length)

# Run R script for fixing newick file format
Rscript clean_tree.R

# analyze Landis et al. radseq rooted tree with aligned rbcL sequences

mkdir -p branch
cd branch

csubst analyze \
--alignment_file ../clean_aligned_csubst.fasta \
--rooted_tree_file ../radseq_pruned_csubst.tre \
--genetic_code 11 \
--iqtree_redo yes \
--threads 3

# write branches of interest, where OCN > 1 (i.e., at least one convergent substitution)
Rscript branches_of_interest.R

# run site analysis on those branches
csubst_site.command


