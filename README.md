# GenBankTabler
Tool to generate tables of sequences available in GenBank for a given taxon or species list.

Note that this code is generated entirely in ChatGPT. I hate ChatGPT, so I have reservations, but it is hard to deny that it is effective, considering that my python skills are nearly zero.

The general idea:
1. User supplies their NCBI user email and either a --csv_file containing a list of species with one pers line, or --taxon "taxon" 
2. For each species in the list (-csv_file version), or for an entire taxon (--taxon), all published sequences from GenBank are downloaded that match the search criteria (or sequences up to a given `--max_seqs`, default 10,000). If desired, the user can specify to only keep markers with a certain number of sequences (`--min_marker_count`), which is useful for keeping only markers that will be phylogenetically informative. 
3. These are split into tables by marker, populating various columns with details from the GenBank metadata
4. Tables are merged into a single table by Isolate and specimen_voucher (options to use just one or the other are included), retaining information if sequence comes from type specimen.
5. The resulting table is output. Currently all marker tables are also generated along the way.
6. If desired, specimens that do not have voucher or isolate numbers are output to separate CSV files of unvouchered sequences, which can be checked by hand.

**Known deficiencies:** 
+ Resulting tables have to be checked manually. Markers that deviate by small name differences (e.g. Rag-1 and Rag1) are not merged into a single column.
+ If the script cannot identify a marker, it will put it into an "Unknown" column. If there are many unknown markers for the same individual, they will appear as duplicated rows of that individual.
+ does not yet merge cases where voucher and isolate are identical. 


# Requirements
`bioconda`

`pandas`

so depending how you are running python: 

`conda install -c conda-forge biopython`

or 

`pip install pandas biopython`

# Usage
`python GenBankTabler.py --csv_file species_list.csv --email yourNCBIemail@email.com --min_marker_count 5 --merge_on both --include_unlinked --max_seqs 10000`

where `species_list.csv` is a csv file containing one species name per line. --min_marker_count is the minimum number of sequences per marker to retain that marker. the `--include_unlinke` flag tells it to output tables of voucherless specimens, which can be used at your discretion. If you do not specify this flag, those specimens are dropped altogether. 

or

`python GenBankTabler.py --taxon "taxon" --email yourNCBIemail@email.com --min_marker_count 5 --merge_on both --include_unlinked --max_seqs 10000`

where "taxon" is your taxon of interest

# Args
+ `--csv_file` the species list you want to get sequences for
+ `--taxon` the taxon you want to get all sequences for
+ `--email` your Entrez email
+ `--min_marker_count` The minimum number of sequences that should be present for a marker to be retained.
+ `--merge_on` What to merge specimens based on. Either Isolate, Voucher, or Both. Default: Both
+ `--include_unlinked` optional call to keep sequences that do not have vouchers. These are outpuat as separate CSV files for the user to inspect. 
+ `--max_seqs` The maximum number of sequences that should be downloaded. Default: 10,000 
+ `--auto-confirm` optional argument to bypass sanity check for the number of sequences you are going to download. Use with caution.
