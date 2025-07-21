# GenBankTabler
Tool to generate tables of sequences available in GenBank for a given taxon or species list.

Note that this code is generated entirely in ChatGPT. I hate ChatGPT, so I have reservations, but it is hard to deny that it is effective, considering that my python skills are nearly zero.

The general idea:
1. User supplies their NCBI user email and either a --csv_file containing a list of species with one pers line, or --taxon "taxon" 
2. For each species in the list (-csv_file version), or for an entire taxon (--taxon), all published sequences from GenBank are downloaded that match the search criteria
3. These are split into tables by marker, populating various columns with details from the GenBank metadata
4. Tables are merged into a single table by Isolate and specimen_voucher (options to use just one or the other are included), retaining information if sequence comes from type
5. The resulting table is output. Currently all marker tables are also generated along the way.

**Warning:** this script does not check how many files it is about to download, and has not been tested on cases where there are genome-level datasets available. It is not recommend that you run it with groups containing thousands of species or tens of thousands of sequences in GenBank. I have no idea what would happen, but it would probably not be pretty.

**Known deficiencies:** 
+ Resulting tables have to be checked manually. Markers that deviate by small name differences (e.g. Rag-1 and Rag1) are not merged into a single column.
+ If the script cannot identify a marker, it will put it into an "Unknown" column. If there are many unknown markers for the same individual, they will appear as duplicated rows of that individual. 


# Requirements
`bioconda`

`pandas`

so depending how you are running python: 

`conda install -c conda-forge biopython`

or 

`pip install pandas biopython`

# Usage
`python GenBankTabler.py --csv_file species_list.csv --email yourNCBIemail@email.com --min_marker_count 5 --merge_on both`

where `species_list.csv` is a csv file containing one species name per line. --min_marker_count is the minimum number of sequences per marker to retain that marker. 

or

`python GenBankTabler.py --taxon "taxon" --email yourNCBIemail@email.com --min_marker_count 5 --merge_on both`

where "taxon" is your taxon of interest





