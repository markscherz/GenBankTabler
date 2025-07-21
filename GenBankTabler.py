import os
import pandas as pd
from Bio import Entrez, SeqIO
from collections import defaultdict
import argparse

# -------------------------
# Function: Fetch sequences from GenBank for a given taxonomic query
# -------------------------
def fetch_genbank_records(query, email, max_seqs, auto_confirm):
    Entrez.email = email
    print(f"Searching GenBank for: {query}")

    handle = Entrez.esearch(db="nucleotide", term=f"{query}[Organism]", retmax=1000000)
    record = Entrez.read(handle)
    handle.close()
    id_list = record["IdList"]

    total_found = len(id_list)
    if max_seqs is not None:
        id_list = id_list[:max_seqs]

    print(f"Found {total_found} sequences for {query}. Will download {len(id_list)}.")

    if not auto_confirm:
        proceed = input("Proceed with download? (y/n): ")
        if proceed.lower() != 'y':
            print("Aborting.")
            exit()

    all_seqs = []
    if not id_list:
        return all_seqs

    for start in range(0, len(id_list), 500):
        end = min(start + 500, len(id_list))
        ids = id_list[start:end]
        handle = Entrez.efetch(db="nucleotide", id=ids, rettype="gb", retmode="text")
        records = list(SeqIO.parse(handle, "genbank"))
        all_seqs.extend(records)
        handle.close()

    return all_seqs

# -------------------------
# Function: Extract marker name directly from gene/product qualifiers, prefer gene
# -------------------------
def extract_marker(record):
    is_mito = False
    for feature in record.features:
        if feature.type == "source" and "organelle" in feature.qualifiers:
            if "mitochondrion" in feature.qualifiers["organelle"][0].lower():
                is_mito = True
                break

    for feature in record.features:
        if feature.type in ["gene", "CDS", "misc_feature", "rRNA", "tRNA", "ncRNA"]:
            if "gene" in feature.qualifiers:
                marker_name = feature.qualifiers["gene"][0].strip().upper()
                return ("_MITO_" + marker_name) if is_mito else marker_name

    for feature in record.features:
        if feature.type in ["gene", "CDS", "misc_feature", "rRNA", "tRNA", "ncRNA"]:
            if "product" in feature.qualifiers:
                marker_name = feature.qualifiers["product"][0].strip().upper()
                return ("_MITO_" + marker_name) if is_mito else marker_name

    return "Unknown"

# -------------------------
# Function: Extract metadata from GenBank record
# -------------------------
def extract_metadata(record, fallback_species):
    metadata = {
        "Isolate": record.annotations.get("isolate", ""),
        "Specimen_Voucher": record.annotations.get("specimen_voucher", ""),
        "Species": record.annotations.get("organism", fallback_species),
        "Geo_Loc_Name": "",
        "Publication_Title": "",
        "GenBank_Accession": record.id,
        "Sequence": str(record.seq),
        "Type_Material": ""
    }

    for feature in record.features:
        if feature.type == "source":
            metadata["Geo_Loc_Name"] = feature.qualifiers.get("country", [""])[0]
            metadata["Isolate"] = feature.qualifiers.get("isolate", [metadata["Isolate"]])[0]
            metadata["Specimen_Voucher"] = feature.qualifiers.get("specimen_voucher", [metadata["Specimen_Voucher"]])[0]
            if "type_material" in feature.qualifiers:
                metadata["Type_Material"] = feature.qualifiers["type_material"][0]

    if "references" in record.annotations and record.annotations["references"]:
        metadata["Publication_Title"] = record.annotations["references"][0].title

    return metadata

# -------------------------
# Function: Process species list file
# -------------------------
def process_species_list(csv_path, email, min_marker_count=0, include_unlinked=False, max_seqs=10000, auto_confirm=False):
    species_df = pd.read_csv(csv_path, header=None, names=["Species"])
    marker_tables = defaultdict(list)

    for species in species_df["Species"]:
        records = fetch_genbank_records(species, email, max_seqs, auto_confirm)
        for record in records:
            marker = extract_marker(record)
            metadata = extract_metadata(record, species)
            marker_tables[marker].append(metadata)

    return finalize_marker_tables(marker_tables, min_marker_count, include_unlinked)

# -------------------------
# Function: Process a single taxon (any level)
# -------------------------
def process_taxon(taxon_name, email, min_marker_count=0, include_unlinked=False, max_seqs=10000, auto_confirm=False):
    records = fetch_genbank_records(taxon_name, email, max_seqs, auto_confirm)
    marker_tables = defaultdict(list)

    for record in records:
        species_name = record.annotations.get("organism", taxon_name)
        marker = extract_marker(record)
        metadata = extract_metadata(record, species_name)
        marker_tables[marker].append(metadata)

    return finalize_marker_tables(marker_tables, min_marker_count, include_unlinked)

# -------------------------
# Finalize and sort marker tables
# -------------------------
def finalize_marker_tables(marker_tables, min_marker_count, include_unlinked):
    sorted_markers = sorted(marker_tables.keys(), key=lambda x: (not x.startswith("_MITO_"), x))
    marker_dfs = {}

    for marker in sorted_markers:
        records = marker_tables[marker]
        df = pd.DataFrame(records)

        has_keys = df[["Isolate", "Specimen_Voucher"]].notna().any(axis=1) & ~(df[["Isolate", "Specimen_Voucher"]] == "").all(axis=1)
        df_linked = df[has_keys].copy()
        df_unlinked = df[~has_keys].copy()

        if min_marker_count > 0 and len(df_linked) < min_marker_count:
            continue

        if not df_linked.empty:
            df_linked.to_csv(f"marker_{marker.replace('/', '-')}.csv", index=False)
            marker_dfs[marker] = df_linked

        if include_unlinked and not df_unlinked.empty:
            df_unlinked_out = df_unlinked[["Species", "GenBank_Accession", "Sequence"]]
            df_unlinked_out.to_csv(f"voucherless_{marker.replace('/', '-')}.csv", index=False)

    return marker_dfs

# -------------------------
# Merge marker tables
# -------------------------
def merge_tables(marker_dfs, merge_on, include_unlinked=False):
    merge_keys = merge_on if isinstance(merge_on, list) else [merge_on]
    merged_df = None

    for marker, df in marker_dfs.items():
        df = df.copy()
        has_merge_keys = df[merge_keys].notnull().any(axis=1)
        df_with_keys = df[has_merge_keys]

        cols = df_with_keys.columns.tolist()
        fixed_order = [col for col in ["Isolate", "Specimen_Voucher", "Species", "Type_Material"] if col in cols]
        remaining = [col for col in cols if col not in fixed_order]
        df_with_keys = df_with_keys[fixed_order + remaining]
        df_with_keys.columns = [col if col in merge_keys else f"{col}_{marker}" for col in df_with_keys.columns]

        if merged_df is None:
            merged_df = df_with_keys
        else:
            merged_df = pd.merge(merged_df, df_with_keys, on=merge_keys, how="outer")

    if merged_df is None:
        print("No sequences with mergeable voucher/isolate information found.")
        return pd.DataFrame()

    species_cols = [col for col in merged_df.columns if col.startswith("Species")]
    if species_cols:
        merged_df["Species"] = merged_df[species_cols].bfill(axis=1).iloc[:, 0]
        merged_df.drop(columns=species_cols, inplace=True)

    type_material_cols = [col for col in merged_df.columns if col.startswith("Type_Material")]
    if type_material_cols:
        merged_df["Type_Material"] = merged_df[type_material_cols].bfill(axis=1).iloc[:, 0]
        merged_df.drop(columns=type_material_cols, inplace=True)
    else:
        merged_df["Type_Material"] = ""

    cols = merged_df.columns.tolist()
    reordered = ["Isolate", "Specimen_Voucher", "Species", "Type_Material"] + [col for col in cols if col not in ["Isolate", "Specimen_Voucher", "Species", "Type_Material"]]
    return merged_df[reordered]

# -------------------------
# Main
# -------------------------
def main():
    parser = argparse.ArgumentParser(description="Download and process GenBank sequences by species or taxon")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--csv_file", help="Path to CSV file with species names (one per line)")
    group.add_argument("--taxon", help="Name of taxonomic group to search (e.g., genus, family)")

    parser.add_argument("--email", required=True, help="Email for NCBI Entrez")
    parser.add_argument("--min_marker_count", type=int, default=0, help="Remove marker tables with fewer than this number of entries")
    parser.add_argument("--merge_on", choices=["Isolate", "Specimen_Voucher", "both"], default="both",
                        help="Column(s) to merge marker tables on")
    parser.add_argument("--output", default="merged_output.csv", help="Output file for merged data")
    parser.add_argument("--include_unlinked", action="store_true", help="Include sequences without isolate or voucher info as separate tables")
    parser.add_argument("--max_seqs", type=int, default=10000, help="Maximum number of sequences to download per search. Use -1 for unlimited.")
    parser.add_argument("--auto_confirm", action="store_true", help="Skip confirmation prompt before downloading sequences")

    args = parser.parse_args()
    max_seqs = None if args.max_seqs == -1 else args.max_seqs

    if args.csv_file:
        marker_dfs = process_species_list(args.csv_file, args.email, args.min_marker_count, args.include_unlinked, max_seqs, args.auto_confirm)
    else:
        marker_dfs = process_taxon(args.taxon, args.email, args.min_marker_count, args.include_unlinked, max_seqs, args.auto_confirm)

    merge_columns = ["Isolate", "Specimen_Voucher"] if args.merge_on == "both" else [args.merge_on]

    merged_df = merge_tables(marker_dfs, merge_columns, include_unlinked=args.include_unlinked)
    if not merged_df.empty:
        merged_df.to_csv(args.output, index=False)
        print(f"Merged output saved to {args.output}")

if __name__ == "__main__":
    main()
