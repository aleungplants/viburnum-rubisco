import pandas as pd
from Bio import AlignIO


with open("../paml/clean_aligned.phylip", "r") as infile:
    data = AlignIO.read(infile, "phylip")

omega = pd.read_csv("paml_site_omega.csv")
sites = omega.loc[(omega["PosteriorProb"] > 0.9) & (omega["NSSitesModel"] == 2)]["Site"] # AIC for model 2 < model 8
sites = list(sites)
print(sites)

aa = []
for i in range(len(data)):
    species_aa = [data[i].name]
    for j in range(len(sites)):
        start = 3 * (sites[j] - 1 - 31) # missing 31 sites. we printed the other data files based on the actual position, but this reads the alignment file so we need to shift it back
        end = 3 * (sites[j] - 1 - 31) + 3
        try:
            species_aa.append(data[i].seq[start:end].translate(gap = "-", table = 11))
        except:
            print("Error with translating at taxon", i,"and AA position", sites[j],
                  ", in", data[i].name, ".",
                  "The codon was", data[i].seq[start:end])
            species_aa.append("-")
            continue
    aa.append(species_aa)

headers = sites
headers.insert(0, "ShortName")
aa = pd.DataFrame(aa, columns = headers)
aa.to_csv("aa_at_sites.csv", index = False)
