#!/bin/bash

# cd the folder containing this file
cd "$(dirname "$0")"

mkdir -p site
cd site

# loop through branches of interest file generated from the csubst branch analayses
while IFS=, read -r line; do
  csubst site \
  --alignment_file ../clean_aligned_csubst.fasta \
  --rooted_tree_file ../radseq_pruned_csubst.tre \
  --cb2 ../branch/csubst_cb_2.tsv \
  --branch_id $line \
  --genetic_code 11 \
  --threads 3
done < ../branch/branches_of_interest.txt

# csubst site --alignment_file ../clean_aligned_csubst.fasta --rooted_tree_file ../radseq_pruned_csubst.tre --cb2 ../branch/csubst_cb_2.tsv --branch_id 282,293 --genetic_code 11 --threads 3