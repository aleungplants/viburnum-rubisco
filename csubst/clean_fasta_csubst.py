import os
import pathlib
from Bio import AlignIO

os.chdir(pathlib.Path(__file__).parent.resolve())

with open("clean_aligned_checked.fasta", "r") as infile:
    data = AlignIO.read(infile, "fasta")
    print("Read file with", len(data), "records")

for record in data:
    record.id = record.id[0:10] # shorten name to match tree
    record.description = "" # remove description, AlignIO.write gives id and description to the taxon
    # record.seq = record.seq[0:1425]  # remove stop codon- not needed anymore because of manual check of the alignment

OUTPATH = "clean_aligned_csubst.fasta"
with open(OUTPATH, "w") as outfile:
    AlignIO.write(data, outfile, "fasta")
    print("Saved", OUTPATH, "with", len(data), "records")
