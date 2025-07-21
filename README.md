# GenBankTabler
Tool to generate tables of sequences available in GenBank for a given taxon or species list.

Note that this code is generated entirely in ChatGPT. I hate ChatGPT, so I have reservations, but it is hard to deny that it is effective, considering that my python skills are nearly zero.

# Requirements
`bioconda`

`pandas`

so depending how you are running python: 

`conda install -c conda-forge biopython`

or 

`pip install pandas biopython`

# Usage
`python GenBankTabler.py --csv_file species_list.csv --email yourNCBIemail@email.com --min_marker_count 5 --merge_on both`

where `species_list.csv` is a csv file containing one species name per line

or

`python GenBankTabler.py --taxon "taxon" --email yourNCBIemail@email.com --min_marker_count 5 --merge_on both`

where "taxon" is your taxon of interest





