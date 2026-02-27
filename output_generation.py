from Bio import SeqIO
import pandas as pd
import os
from utils import build_gene_entry

def setup_output_columns(df, method_info, type1_validation_status=None, family_map=None):
    if df.empty:
        return df

    # Column mapping for cleaner code
    score_columns = {
        'mmseqs': ('MMseqs_Evalue', 'MMseqs_Bitscore'),
        'mmseqs2_rna': ('MMseqs_Evalue', 'MMseqs_Bitscore'),  # RNA uses same columns
        'hmm': ('HMM_Evalue', 'HMM_Bitscore'),
        'foldseek': ('Foldseek_Evalue', 'Foldseek_Bitscore')
    }

    # Methods to check for target, in priority order
    target_methods = ['mmseqs2', 'mmseqs2_dna', 'mmseqs2_rna', 'hmm', 'foldseek']

    # Add base columns (including Type_I_Status for validation)
    base_cols = ["TA_Family", "domain", "MMseqs_hit", "Method",
                 "MMseqs_Evalue", "MMseqs_Bitscore", "HMM_Evalue", "HMM_Bitscore", "Foldseek_Evalue", "Foldseek_Bitscore",
                 "Type_I_Status"]

    for col in base_cols:
        df[col] = '-'

    if family_map is None:
        family_map = {}

    # Fill method information
    for i in range(len(df)):
        protein_id = df.iloc[i]['Protein_ID']

        # Fill Type I validation status
        if type1_validation_status and protein_id in type1_validation_status:
            df.iloc[i, df.columns.get_loc('Type_I_Status')] = type1_validation_status[protein_id]

        if protein_id in method_info:
            info = method_info[protein_id]
            df.iloc[i, df.columns.get_loc('Method')] = '/'.join(info.get('methods', []))

            # Fill scores using mapping
            for method, (evalue_col, bitscore_col) in score_columns.items():
                evalue_key = f'{method}_evalue'
                if evalue_key in info:
                    df.iloc[i, df.columns.get_loc(evalue_col)] = info[evalue_key]

                bitscore_key = f'{method}_bitscore'
                if bitscore_key in info:
                    df.iloc[i, df.columns.get_loc(bitscore_col)] = info[bitscore_key]

            # Resolve TA family name from the best available target hit
            for method in target_methods:
                target_key = f'{method}_target'
                if target_key in info:
                    target_id = info[target_key]
                    family_name = _resolve_family_name(target_id, family_map)
                    if family_name:
                        df.iloc[i, df.columns.get_loc('TA_Family')] = family_name
                        break

    if 'Contig' in df.columns:
        cols = [c for c in df.columns if c != 'Contig']
        product_idx = cols.index('Product') if 'Product' in cols else len(cols)-1
        cols.insert(product_idx + 1, 'Contig')
        df = df[cols]

    return df


def _resolve_family_name(target_id, family_map):
    """Resolve a database target ID to a human-readable protein family name.

    Checks the lookup table first, then falls back to stripping the
    toxin_/antitoxin_ prefix from the target ID.  If the stripped result
    is purely numeric (e.g. HMM profile numbers), return None so the
    caller can try the next method in the priority list.
    """
    if target_id in family_map:
        return family_map[target_id]
    # Fall back: strip prefix
    stripped = None
    if target_id.startswith('toxin_'):
        stripped = target_id[6:]
    elif target_id.startswith('antitoxin_'):
        stripped = target_id[10:]
    else:
        stripped = target_id
    # Don't return bare numbers (HMM profile IDs like toxin_169)
    if stripped and not stripped.isdigit():
        return stripped
    return None

def create_all_outputs(search_results, neighborhoods, config, args, type1_validation_status=None):
    """Generate all output files."""
    toxin_genes, antitoxin_genes, method_info = search_results[0], search_results[1], search_results[2]

    # Import here to avoid circular dependency
    from result_processing import pair_ta_systems, filter_same_type_TA, validate_complete_pairs

    df_pairs, orphan_genes = pair_ta_systems(toxin_genes, antitoxin_genes, neighborhoods, args.max_pairing_distance)

    if not df_pairs.empty:
        # Add annotations BEFORE filtering
        family_map = config.get('family_map', {})
        df_pairs = setup_output_columns(df_pairs, method_info, type1_validation_status, family_map=family_map)

        # Apply validation filters
        df_pairs = filter_same_type_TA(df_pairs)
        df_pairs = validate_complete_pairs(df_pairs)

        if not df_pairs.empty:
            df_pairs.to_csv(args.outfile, index=False)
            print(f"TA pairs saved to {args.outfile}")
        else:
            print("No valid TA pairs after validation")

    # Orphan output
    if args.predict_orphan and orphan_genes:
        ptt = neighborhoods['ptt']
        has_contig = "Contig" in ptt.columns
        orphan_data = []
        for i, gene_id in enumerate(orphan_genes, 1):
            gene_row = ptt[ptt['PID'] == gene_id].iloc[0]
            gene_type = 'Toxin' if gene_id in toxin_genes else 'Antitoxin'
            orphan_entry = build_gene_entry(gene_row, gene_type, i, has_contig)
            # Rename TA_ID to Orphan_TA_ID for orphans
            orphan_entry['Orphan_TA_ID'] = orphan_entry.pop('TA_ID')
            orphan_data.append(orphan_entry)

        df_orphans = pd.DataFrame(orphan_data)
        family_map = config.get('family_map', {})
        df_orphans = setup_output_columns(df_orphans, method_info, type1_validation_status, family_map=family_map)
        orphan_output = args.outfile.replace('.csv', '_orphan.csv')
        df_orphans.to_csv(orphan_output, index=False)
        print(f"Orphans saved to {orphan_output}")

    # Neighborhood output
    if args.save_unidentified:
        save_neighborhoods_and_annotation(neighborhoods, toxin_genes.union(antitoxin_genes), orphan_genes,
                                         toxin_genes, antitoxin_genes, config, args)

def save_neighborhoods_and_annotation(neighborhoods, all_identified, orphan_genes, toxin_genes, antitoxin_genes, config, args):
    unidentified = set(neighborhoods['proteins_in_neighborhoods']) - all_identified
    target_neighborhoods = set()

    for gene in unidentified.union(orphan_genes):
        if gene in neighborhoods['gene_to_neighborhood']:
            target_neighborhoods.add(neighborhoods['gene_to_neighborhood'][gene])

    if not target_neighborhoods:
        print("No neighborhoods to save")
        return

    ptt = neighborhoods['ptt']
    genome_seqs = {r.id: str(r.seq) for r in SeqIO.parse(config['fna_file'], "fasta")}
    annotation_data = []

    with open(args.save_unidentified, 'w') as fasta_f:
        for nbh_id in sorted(target_neighborhoods):
            genes = neighborhoods['neighborhood_to_genes'][nbh_id]
            gene_coords = []

            for gene_id in genes:
                gene_row = ptt[ptt['PID'] == gene_id]
                if len(gene_row) > 0:
                    location = gene_row.iloc[0]['Location']
                    strand = gene_row.iloc[0]['Strand']
                    start, end = map(int, location.split('..'))
                    gene_coords.append((start, end, gene_id, strand))

            if gene_coords:
                gene_coords.sort()
                nbh_start = max(1, min(c[0] for c in gene_coords) - 100)
                nbh_end = max(c[1] for c in gene_coords) + 100

                contig_id = list(genome_seqs.keys())[0]
                if nbh_end <= len(genome_seqs[contig_id]):
                    seq = genome_seqs[contig_id][nbh_start-1:nbh_end]
                    fasta_f.write(f">{nbh_id} {nbh_start}..{nbh_end} {contig_id}\n{seq}\n")

                    for start, end, gene_id, strand in gene_coords:
                        rel_start = start - nbh_start + 1
                        rel_end = end - nbh_start + 1

                        annotated_id = gene_id
                        if gene_id in toxin_genes:
                            annotated_id += "_T"
                        elif gene_id in antitoxin_genes:
                            annotated_id += "_AT"

                        annotation_data.append({
                            'Neighborhood': nbh_id, 'Annotated_ID': annotated_id,
                            'Location': f"{rel_start}..{rel_end}", 'Strand': strand
                        })

    if annotation_data:
        ann_df = pd.DataFrame(annotation_data)
        ann_file = args.save_unidentified.replace('.fna', '.tsv')
        ann_df.to_csv(ann_file, sep='\t', index=False)
        print(f"Saved {len(target_neighborhoods)} neighborhoods to {args.save_unidentified}")
        print(f"Saved annotations to {ann_file}")
