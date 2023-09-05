Data and code for analyses in Leung et al. unpublished manuscript
================

Here is an explanation for the files provided in each folder. As much as
possible, I wrote the portable code, i.e., no need to enter any working
directories. However, to do so R scripts will require installing the
`here` package. Filepaths in Python scripts should work relative to the
location of the script. Data files \>100 MB in size that can be
generated with the scripts are not included.

/download:

- `accessions.csv` was manually taken from the supplemental information
  of Landis et al. (2021) Systematic Biology. Note that if you’re going
  to look there for accession numbers, some are missing (`missing.csv`),
  so I had to manually search on GenBank for missing accession numbers
  for species in their tree (I manually checked the citation in the
  GenBank entry that it was from that paper).
- `get_sequence.py` uses `accession.csv` to generate a `sequence.gb`
  file.
- `clean_sequence.py` uses `sequence.gb` to generate a .fasta file and
  makes sure the sequences are only rbcL sequences (and not other
  genes).
- We also use species names from Landis et al. and not the species from
  the GenBank description; this is documented in `metadata.csv`.

`/align`:

- `clean_aligned.fasta` was generated from `/download/clean.fasta` using
  `muscle5.1.win64.exe` (not provided). Batch (.bat) files were used for
  reproducibility; must be run on Windows.
- I manually checked the alignment in Mega 11; this is
  `clean_aligned_checked.fasta`.
- I also generated an `ensemble.efa` file, which tells you how much
  dispersion (D) between replicates. D was close to zero (1.757e-5).
- `fasta_to_phylip.py` converts the aligned .fasta file to .phylip,
  which is needed for PAML. Also copies to the /paml, /prune_tree, and
  /iqtree folders. Also formats the file for PAML-specific phylip
  format.

`/prune_tree`:

- `prune_tree.py` prunes the Landis et al. `radseq_cpdna_stage2.tre`
  tree to the tips for which we have rbcL sequences (not all do).

*Because all the rbcL sequences start at nucleotide position 32, the
position numbers are shifted.*

`/paml`:

- `unroot_tree.R` is used to unroot the Landis et al. tree, which is
  needed for PAML
- `clean_aligned_paml.phylip` is used with the unrooted RADseq tree in
  `codeml.exe` (not provided) to conduct analyses of selection.
- See `codeml.ctl` for settings.

`/paml_site_reconstruction`:

- The best fitting codon model from `/paml` analyses was used to
  reconstruct ancestral amino acids, see `codeml.ctl` for settings.

`/csubst`:

- `csubst_branch.command` is a macOS command file for running csubst. No
  executable is needed because CSUBST should install to your PATH (i.e.,
  just need to type csubst function in Terminal). It will also run some
  other scripts as needed. See comments for explanation of most of the
  command file.
  - `clean_tree.R` is needed because CSUBST doesn’t play well with the
    tree annotated with other information.
  - output is stored in `/csubst/branch`
- `csubst_site.command` is run after the branch analyses (which
  identifies branches with at least one convergent amino acid
  substitutions), to identify which codon site the substituion(s) is/are
  at.
  - output is stored in `/csubst/site`
- `analysis/merge_site_data.R` merges the site data across all branches
  of interest.
- `analysis/get_convergent_clades.R` will use the branch and node
  numbers to identify which clades/taxa are descendents of the branches
  where convergent evolution occurred. Also if those branches correspond
  to named clades (e.g., Succotinus), then it will note that too.

`/download_gbif`:

- `/wc2.1_30s_bio` should contain the WorldClim Bioclimatic variable
  raster files. Left empty on purpose because you can get these files at
  [https://worldclim.org/data/worldclim21.html](). Used in the main
  `/analysis/analyze.rmd` file.
- `get_gbif_data.r` will grab the GBIF occurrences. I provide the GBIF
  download in a zip file. The script saves the citation and other useful
  info in the `supplemental` folder.

`/edwards_climate` just contains the highly curated climate data for 120
*Viburnum* species from Edwards et al. 2017 American Naturalist
Tropical-Temperate transition paper.

`/landis_biome_range`:

- `viburnum.biome.n4.nex` and `viburnum.range.n6.nex` are from Landis et
  al.
- `parse_biome_range.R` turns those data into .csv files; I only ended
  up using the biome data.

`/analyze`:

- `analyze.rmd` is the main analysis file I used to process data and
  generate figures. This folder contains .csv files of processed data
  used to make the figures.
- `aa_at_sites.py` is needed because although I parsed the PAML outputs
  in `analyze.rmd` to get the positively selected sites (where in the
  protein), I needed to know what the amino acid at that site is; thus,
  needed to use the alignment file to get it.

`/structure`:

- pymol.txt has the PyMol functions I used to generate the figures.
  Outputed images are combined in `analyze.rmd`.
- I used 8RUC; easy to get from PDB so I don’t include it here.

`/iqtree` has the IQTree .bat file for generating the gene tree.
