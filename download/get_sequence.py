import csv
from Bio import Entrez

# File containing accessions from Landis et al. 2021 Syst Biol and species names. Names are needed because GenBank names are not always accurate.
reader = csv.DictReader(open("accessions.csv"))
accessions = {}
for row in reader:
    accessions[row["species"]] = row["accession"]

# Download from the list of accessions
GB_FILENAME = "sequence.gb"
Entrez.email = "arthur.aqualung@gmail.com"

open(GB_FILENAME, "w")
for species, accession in accessions.items():
    record = Entrez.efetch(db = "nucleotide", id = accession, rettype = "gb", retmode = "text")
    print("Writing", accession, "to", GB_FILENAME)
    with open(GB_FILENAME, "a") as f:
        f.write(record.read())
print("Wrote .gb file")

