from Bio import SeqIO
import pandas as pd
import os
from utils import build_gene_entry

def setup_output_columns(df, method_info, type1_validation_status=None):
    if df.empty:
        return df

    # Column mapping for cleaner code
    score_columns = {
        'mmseqs': ('MMseqs_Evalue', 'MMseqs_Bitscore'),
        'hmm': ('HMM_Evalue', 'HMM_Bitscore'),
        'foldseek': ('Foldseek_Evalue', 'Foldseek_Bitscore')
    }

    # Add base columns (including Type_I_Status for validation)
    base_cols = ["domain", "MMseqs_hit", "Method",
                 "MMseqs_Evalue", "MMseqs_Bitscore", "HMM_Evalue", "HMM_Bitscore", "Foldseek_Evalue", "Foldseek_Bitscore",
                 "Type_I_Status"]

    for col in base_cols:
        df[col] = '-'

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

    if 'Contig' in df.columns:
        cols = [c for c in df.columns if c != 'Contig']
        product_idx = cols.index('Product') if 'Product' in cols else len(cols)-1
        cols.insert(product_idx + 1, 'Contig')
        df = df[cols]

    return df

def create_all_outputs(search_results, neighborhoods, config, args, type1_validation_status=None):
    """Generate all output files."""
    toxin_genes, antitoxin_genes, method_info = search_results[0], search_results[1], search_results[2]
    self_ta_genes = search_results[3] if len(search_results) > 3 else set()

    # Import here to avoid circular dependency
    from result_processing import pair_ta_systems, filter_same_type_TA, validate_complete_pairs

    # Pass self_ta_genes to handle pairing logic
    df_pairs, orphan_genes, final_self_ta_genes = pair_ta_systems(toxin_genes, antitoxin_genes, neighborhoods, self_ta_genes)

    if not df_pairs.empty:
        # Add annotations BEFORE filtering
        df_pairs = setup_output_columns(df_pairs, method_info, type1_validation_status)

        # Apply validation filters
        df_pairs = filter_same_type_TA(df_pairs)
        df_pairs = validate_complete_pairs(df_pairs)

        if not df_pairs.empty:
            df_pairs.to_csv(args.output, index=False)
            print(f"TA pairs saved to {args.output}")
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
        df_orphans = setup_output_columns(df_orphans, method_info, type1_validation_status)
        orphan_output = args.output.replace('.csv', '_orphan.csv')
        df_orphans.to_csv(orphan_output, index=False)
        print(f"Orphans saved to {orphan_output}")

    # Self-TA output (genes that are their own antitoxin) - use final_self_ta_genes (unpaired only)
    if final_self_ta_genes:
        ptt = neighborhoods['ptt']
        has_contig = "Contig" in ptt.columns
        self_ta_data = []
        for i, gene_id in enumerate(final_self_ta_genes, 1):
            gene_rows = ptt[ptt['PID'] == gene_id]
            if len(gene_rows) == 0:
                continue
            gene_row = gene_rows.iloc[0]
            # Self-TAs are both toxin and antitoxin
            self_ta_entry = build_gene_entry(gene_row, 'Self-TA', i, has_contig)
            self_ta_entry['Self_TA_ID'] = self_ta_entry.pop('TA_ID')
            self_ta_data.append(self_ta_entry)

        if self_ta_data:
            df_self_ta = pd.DataFrame(self_ta_data)
            df_self_ta = setup_output_columns(df_self_ta, method_info, type1_validation_status)

            # Create filename with "self_" prefix as requested
            base_name = os.path.basename(args.output)
            dir_name = os.path.dirname(args.output)
            self_ta_output = os.path.join(dir_name, f"self_{base_name}")

            df_self_ta.to_csv(self_ta_output, index=False)
            print(f"Self-TA genes saved to {self_ta_output}")

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
