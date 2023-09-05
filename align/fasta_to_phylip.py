import os
import pathlib
from Bio import AlignIO

os.chdir(pathlib.Path(__file__).parent.resolve())

with open("clean_aligned_checked.fasta", "r") as infile: # removed first 31 codons due to missing coding sequencing, also manually aligned codon 475
    data = AlignIO.read(infile, "fasta")
    print("Read file with", len(data), "records")

# Save files for (1) tree making, (2) tree pruning, and (3) PAML
os.chdir("..") # go up one dir
outpaths = ["iqtree/clean_aligned.phylip", "prune_tree/clean_aligned.phylip", "paml/clean_aligned.phylip"]
for path in outpaths:
    with open(path, "w") as outfile:
        AlignIO.write(data, outfile, "phylip")
        print("Saved", path, "with", len(data), "records")

# Fix formatting for PAML: add "I" to line 1 and add spaces after taxon names
for path in outpaths:
    outpath = path.split(".")
    outpath = "_paml.".join(outpath)
    with open(path, "r") as infile, open(outpath, "w") as outfile:
        for line_number, line in enumerate(infile):
            if line_number == 0:
                line = line.replace("\n", ' I\n')
            else:
                line = line.split(" ", maxsplit = 1) # split removes the delimiter
                line = "  ".join(line) # add back delimited (space) with an extra space - PAML wants two spaces between taxon name and sequence
            outfile.write(line)
    print("Saved", outpath, "in PAML-friendly format")
