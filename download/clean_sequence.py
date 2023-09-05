import csv
from Bio import SeqIO

# File containing accessions from Landis et al. 2021 Syst Biol and species names. Names are needed because GenBank names are not always accurate.

reader = csv.DictReader(open("accessions.csv"))
accessions = {}
for row in reader:
    accessions[row["accession"]] = row["species"]

# Read GenBank records as a SeqIO object.

records = list(SeqIO.parse("sequence.gb", "genbank"))

clean_records = []
record_metadata = [] # to keep track of the original .gb data

print("Opened .gb file with", len(records), "records")
for i in reversed(range(len(records))):
    for j in range(len(records[i].features)):
        if records[i].features[j].type == "CDS" and "gene" in records[i].features[j].qualifiers:
            if records[i].features[j].qualifiers["gene"] == ["rbcL"]: # check all features for rbcL gene
                cds = records[i].features[j].location # get the location data for the rbcL CDS
                if cds.strand == -1:  # 1 is positive, -1 is negative
                    cds_seq = cds.extract(records[i].seq).complement()
                    records[i].seq = cds_seq
                if cds.strand == 1:
                    cds_seq = cds.extract(records[i].seq)
                    records[i].seq = cds_seq
    records[i].id = accessions[records[i].name]
    # if len(records[i].seq) < 1150:
    #     print(records[i].id, "sequence is too short, with length of", len(records[i].seq))
    #     continue
    records[i].description = ""
    clean_records.append(records[i])
    record_metadata.append(",".join([records[i].name, records[i].annotations["organism"], records[i].id]))

# Save files.
with open("clean.fasta", "w") as outfile:
    SeqIO.write(clean_records, outfile, "fasta")
print("Saved .fasta file with", len(clean_records), "sequences")
with open("metadata.csv", "w") as outfile:
    outfile.write("accession,gb_species,species\n")
    for k in range(len(record_metadata)):
        outfile.write(record_metadata[k] + "\n")
print("Saved metadata file with", len(record_metadata), "accessions")
