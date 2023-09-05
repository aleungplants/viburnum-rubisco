import os
from Bio import AlignIO
from Bio import Phylo
import shutil

names = []
with open("clean_aligned.phylip", "r") as infile:
    data = AlignIO.read(infile, "phylip")
    for record in data:
        if len(record.name) > 10:
            record.name = record.name[0:10]
        names.append(record.name)

os.chdir("..")
tree = Phylo.read("prune_tree/radseq_cpdna_stage2.tre", "nexus") # had to rename rigidum to rugosum
print(tree.get_terminals())
for tree_tip in tree.get_terminals():
    FOUND = False
    for name in names:
        if name in tree_tip.name:
            tree_tip.name = name
            FOUND = True
    if FOUND == False:
        print("Pruning", tree_tip.name)
        tree.prune(target = tree_tip.name)

Phylo.draw_ascii(tree)
Phylo.write(tree, "prune_tree/radseq_pruned.tre", "newick")
shutil.copyfile("prune_tree/radseq_pruned.tre", "paml/radseq_pruned.tre")
