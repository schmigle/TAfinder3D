import pandas as pd
import os

# Consolidated helper functions
def get_gene_id(ptt, idx, use_locus_tag):
    return ptt.loc[idx, "Synonym"] if use_locus_tag else ptt.loc[idx, "PID"]

def is_ncrna(gene_id):
    return gene_id.startswith('Nucl')

def calculate_distance(location_up, location_down):
    location_up_right = float(location_up.split("..")[1])
    location_down_left = float(location_down.split("..")[0])
    return location_down_left - location_up_right

def calculate_overlap(coord_query, coord_subject):
    start_1, end_1 = coord_query.split("..")
    start_2, end_2 = coord_subject.split("..")
    start_1, start_2, end_1, end_2 = int(start_1), int(start_2), int(end_1), int(end_2)
    overlap = min(end_1, end_2) - max(start_1, start_2) + 1
    return overlap / (end_1 - start_1 + 1)

def determine_cutoff_params(args, search_type):
    if args.use_evalue:
        # Use e-value cutoffs
        return args.evalue, True, str(args.evalue)

    # Use bitscore cutoff
    bitscore_map = {
        'mmseqs': args.mmseqs_bitscore,
        'hmm': args.hmm_bitscore,
        'foldseek': args.foldseek_bitscore
    }
    return bitscore_map[search_type], False, "1000"

def find_overlapping_gene_in_ptt(ptt, query_start, query_end, has_contig, query_contig=None):
    for idx in range(len(ptt)):
        # If we have contig info, filter by it first
        if has_contig and query_contig:
            if ptt.loc[idx, "Contig"] != query_contig:
                continue

        loc = ptt.loc[idx, "Location"]
        start, end = map(int, loc.split(".."))
        # Check for overlap
        if max(start, query_start) <= min(end, query_end):
            gene_id = ptt.loc[idx, "PID"]
            strand = ptt.loc[idx, "Strand"]
            contig = ptt.loc[idx, "Contig"] if has_contig else None
            return gene_id, start, end, strand, contig
    return None, None, None, None, None

def update_method_info(results, gene_id, method_name, evalue, bitscore, use_evalue):
    """Update method_info dict with new results, avoiding duplicates."""
    if gene_id in results['method_info']:
        # Don't add duplicate method name
        if method_name not in results['method_info'][gene_id]['methods']:
            results['method_info'][gene_id]['methods'].append(method_name)
        # Update scores if better
        evalue_key = f'{method_name.lower()}_evalue'
        bitscore_key = f'{method_name.lower()}_bitscore'
        if use_evalue:
            if evalue < results['method_info'][gene_id].get(evalue_key, float('inf')):
                results['method_info'][gene_id][evalue_key] = evalue
                results['method_info'][gene_id][bitscore_key] = bitscore
        else:
            if bitscore > results['method_info'][gene_id].get(bitscore_key, 0):
                results['method_info'][gene_id][evalue_key] = evalue
                results['method_info'][gene_id][bitscore_key] = bitscore
    else:
        method_key = method_name.lower()
        results['method_info'][gene_id] = {
            'methods': [method_name],
            f'{method_key}_evalue': evalue,
            f'{method_key}_bitscore': bitscore
        }

def build_gene_entry(gene_row, gene_type, ta_id, has_contig):
    """Build a single gene entry dict for output."""
    entry = {
        'TA_ID': ta_id, 'Type': gene_type, 'Locus_tag': gene_row['Synonym'],
        'Coordinates': gene_row['Location'], 'Strand': gene_row['Strand'],
        'Length': gene_row['Length'], 'Protein_ID': gene_row['PID'],
        'Product': gene_row['Product']
    }
    if has_contig:
        entry['Contig'] = gene_row['Contig']
    return entry

# Consolidated parsing functions
def parse_search_output(output_file, search_type, cutoff_value, use_evalue=False, havalue=None, maximum_length=500):
    if not os.path.exists(output_file) or os.path.getsize(output_file) == 0:
        return pd.DataFrame()

    try:
        if search_type == 'mmseqs':
            # MMseqs2 easy-search output format:
            # query, target, pident, alnlen, mismatch, gapopen, qstart, qend, tstart, tend, evalue, bits
            df = pd.read_csv(output_file, header=None, sep="\t")
            # Columns: 0=query, 1=target, 2=pident, 3=alnlen, 10=evalue, 11=bitscore

            # Apply cutoff based on mode
            if use_evalue:
                df = df[df[10] < cutoff_value]  # E-value filter
            else:
                df = df[df[11] >= cutoff_value]  # Bitscore filter

            # Apply length and identity filters
            df = df[(df[3] >= 30) & (df[3] <= maximum_length) & (df[2] >= 30)]

            # Calculate Ha-value if needed
            if havalue:
                # Ha-value = (identity/100) * alignment_length / query_length
                # Columns: 2=pident, 11=bits (keeping for compatibility)
                df[14] = (df[2] / 100) * df[5] / df[3]  # Approximation
                df = df[df[14] >= havalue]

            # Return relevant columns: query, target, pident, alnlen, mismatch, gapopen, evalue, bitscore
            return df[[0, 1, 2, 3, 4, 5, 10, 11]] if df.shape[0] > 0 else pd.DataFrame()

        elif search_type == 'hmm':
            df = pd.read_csv(output_file, skiprows=2, header=None, sep=r"\s+", comment="#")
            if df.empty:
                return pd.DataFrame()
            # Columns: 0=target, 2=tlen, 12=ievalue, 13=bitscore
            df = df[[0, 2, 3, 4, 17, 18, 12, 13]]

            # Apply cutoff based on mode
            if use_evalue:
                df = df[(df[2] >= 30) & (df[2] <= maximum_length) & (df[12] < cutoff_value) & (df[0] != "")]
            else:
                df = df[(df[2] >= 30) & (df[2] <= maximum_length) & (df[13] >= cutoff_value) & (df[0] != "")]
            return df

        elif search_type == 'mmseqs_rna':
            # MMseqs2 RNA search - same format as protein search
            df = pd.read_csv(output_file, header=None, sep="\t")

            # Apply cutoff based on mode
            if use_evalue:
                df = df[df[10] < cutoff_value]  # E-value filter
            else:
                df = df[df[11] >= cutoff_value]  # Bitscore filter

            # Calculate Ha-value for compatibility
            df[14] = (df[2] / 100) * df[5] / df[4]
            # Include columns 7,8 (qstart, qend) for coordinate matching in DNA searches
            return df[[0, 1, 2, 3, 4, 5, 7, 8, 10, 11, 14]]

    except Exception:
        # Silently return empty DataFrame - likely means no hits found
        return pd.DataFrame()
