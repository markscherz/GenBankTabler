# GenBankTabler
Tool to generate tables of sequences available in GenBank for a given taxon or species list.

Note that this code is generated entirely in ChatGPT. I hate ChatGPT, so I have reservations, but it is hard to deny that it is effective, considering that my python skills are nearly zero.

# Requirements
`bioconda`

`pandas`

# Usage
`python genbank_taxon_tool.py --csv_file species_list.csv --email yourncbiemail@email.com --min_marker_count 5 --merge_on both`

where `species_list.csv` is a csv file containing one species name per line

or

`python genbank_taxon_tool.py --taxon "taxon" --email yourncbiemail@email.com --min_marker_count 5 --merge_on both`

where "taxon" is your taxon of interest





