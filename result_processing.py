from Bio.Seq import Seq
import pandas as pd
from utils import is_ncrna, build_gene_entry

def has_antisense_complementarity(rna_seq, mrna_seq, rna_strand, prot_strand, window_size=20, max_mismatches=1):
    """Check if RNA and mRNA have complementarity in antisense orientation.

    For Type I systems, RNA on one strand should be complementary to mRNA on opposite strand.
    The RNA reverse complement should match regions of the mRNA with up to max_mismatches.
    """
    try:
        rna_rc = str(Seq(rna_seq).reverse_complement())
    except (ValueError, TypeError, AttributeError):
        return False

    for i in range(len(rna_rc) - window_size + 1):
        rna_window = rna_rc[i:i + window_size]
        for j in range(len(mrna_seq) - window_size + 1):
            mrna_window = mrna_seq[j:j + window_size]
            mismatches = sum(1 for a, b in zip(rna_window, mrna_window) if a != b)
            if mismatches <= max_mismatches:
                return True

    return False

def validate_type1_complementarity(neighborhoods, rna_results, genome_seq, toxin_genes, antitoxin_genes, min_complementarity=30):
    """Validate Type I TA systems: RNA and protein on opposite strands within 150bp with complementarity.

    Type I: RNA antitoxin must be within 150bp of protein toxin on opposite strand with antisense complementarity.
    Type III: Treated as regular TA pairs (no special validation needed).

    Args:
        neighborhoods: Neighborhood data structure
        rna_results: RNA search results
        genome_seq: Genome sequence string
        toxin_genes: Set of identified toxin genes (from merged search results)
        antitoxin_genes: Set of identified antitoxin genes (from merged search results)
        min_complementarity: Minimum complementarity window size
    """
    ptt = neighborhoods['ptt']
    validation_status = {}  

    total_type1_rnas = sum(1 for g in antitoxin_genes if g.startswith('Nucl_type1_at_'))
    total_type3_rnas = sum(1 for g in antitoxin_genes if g.startswith('Nucl_type3_at_'))

    print(f"Validating Type I/III TA systems (opposite strands + complementarity within 150bp)...")
    print(f"  Found {total_type1_rnas} Type I RNA antitoxins, {total_type3_rnas} Type III RNA antitoxins to validate")

    for nbh_id, genes in list(neighborhoods['neighborhood_to_genes'].items()):
        type1_rna_antitoxins = [g for g in genes if g.startswith('Nucl_type1_at_') and g in antitoxin_genes]
        protein_toxins = [g for g in genes if not is_ncrna(g) and g in toxin_genes]

        if not type1_rna_antitoxins or not protein_toxins:
            continue  

        for rna_gene in type1_rna_antitoxins:
            rna_row = ptt[ptt['PID'] == rna_gene].iloc[0]
            for prot_gene in protein_toxins:
                prot_row = ptt[ptt['PID'] == prot_gene].iloc[0]

                if rna_row['Strand'] == prot_row['Strand']:
                    continue

                try:
                    rna_start, rna_end = map(int, rna_row['Location'].split('..'))
                    prot_start, prot_end = map(int, prot_row['Location'].split('..'))

                    if rna_start < prot_start:
                        distance = prot_start - rna_end  # Gap between genes
                    else:
                        distance = rna_start - prot_end

                    # Allow negative distance (overlap) or small positive distance
                    if distance > 150:
                        continue  # Too far apart

                    # 3. Now check complementarity (only for nearby opposite-strand pairs)
                    rna_seq = genome_seq[rna_start-1:rna_end]
                    mrna_seq = genome_seq[prot_start-1:prot_end]

                    if has_antisense_complementarity(rna_seq, mrna_seq, rna_row['Strand'],
                                                    prot_row['Strand'], window_size=min_complementarity):
                        # Type I system confirmed! Mark both genes as validated
                        status = "validated" if rna_gene in rna_results.get('all_rna_genes', set()) else "putative"
                        validation_status[rna_gene] = status
                        validation_status[prot_gene] = status
                    # If no complementarity, just don't mark as validated (but allow pairing)

                except (ValueError, IndexError, KeyError):
                    continue

    validated_count = sum(1 for v in validation_status.values() if v == "validated")
    putative_count = sum(1 for v in validation_status.values() if v == "putative")

    # Count how many RNA antitoxins were validated
    validated_rna_ats = sum(1 for g, v in validation_status.items() if g.startswith('Nucl_type1_at_') or g.startswith('Nucl_type3_at_'))

    print(f"  Complementarity validation: {validated_count} genes validated, {putative_count} putative")
    print(f"  {validated_rna_ats} out of {total_type1_rnas + total_type3_rnas} RNA antitoxins passed complementarity check")

    return validation_status


def merge_search_results(*search_results, args=None):
    # Track which genes are toxins/antitoxins and their normalized scores
    gene_assignments = {}  # gene_id -> {'type': 'toxin'/'antitoxin', 'method': str, 'score': float, 'normalized_score': float}
    combined_method_info = {}

    def normalize_score(method, score):
        """Normalize bitscores for cross-method comparison.
        - MMseqs2 protein: multiply by 0.675
        - Foldseek: multiply by 0.2
        - MMseqs2 nucleotide (DNA/RNA): multiply by 0.675
        - HMM: use raw bitscore"""
        if 'foldseek' in method.lower():
            return score * 0.2
        elif 'mmseqs' in method.lower():
            # All MMseqs2 searches (protein, DNA, RNA) use same normalization
            return score * 0.675
        else:
            # HMM uses raw bitscore
            return score

    def get_best_method_and_score(info):
        """Get method with highest normalized bitscore from method_info."""
        best_method = None
        best_raw_score = 0
        best_normalized_score = 0

        for method in info.get('methods', []):
            # Get bitscore for this method
            if 'mmseqs2' in method.lower() and 'dna' not in method.lower() and 'rna' not in method.lower():
                score_key = 'mmseqs2_bitscore'
            elif 'dna' in method.lower():
                score_key = 'mmseqs2_dna_bitscore'
            elif 'hmm' in method.lower():
                score_key = 'hmm_bitscore'
            elif 'foldseek' in method.lower():
                score_key = 'foldseek_bitscore'
            else:
                score_key = f"{method.lower().replace('_', '')}_bitscore"

            raw_score = info.get(score_key, 0)
            normalized = normalize_score(method, raw_score)

            if normalized > best_normalized_score:
                best_method = method
                best_raw_score = raw_score
                best_normalized_score = normalized

        return best_method, best_raw_score, best_normalized_score

    # Process all results
    for result in search_results:
        # Merge method info first
        for gene_id, info in result['method_info'].items():
            if gene_id in combined_method_info:
                existing_methods = set(combined_method_info[gene_id]['methods'])
                new_methods = set(info['methods'])
                combined_method_info[gene_id]['methods'] = list(existing_methods.union(new_methods))
                combined_method_info[gene_id].update({k: v for k, v in info.items() if k != 'methods'})
            else:
                combined_method_info[gene_id] = info.copy()

        # Track toxin assignments
        for gene_id in result['toxins']:
            if gene_id not in combined_method_info:
                continue

            method, raw_score, normalized_score = get_best_method_and_score(combined_method_info[gene_id])

            if gene_id not in gene_assignments:
                gene_assignments[gene_id] = {'type': 'toxin', 'method': method, 'score': raw_score, 'normalized_score': normalized_score}
            else:
                # Already assigned - check if this has higher normalized score
                existing = gene_assignments[gene_id]
                if normalized_score > existing['normalized_score']:
                    gene_assignments[gene_id] = {'type': 'toxin', 'method': method, 'score': raw_score, 'normalized_score': normalized_score}

        # Track antitoxin assignments
        for gene_id in result['antitoxins']:
            if gene_id not in combined_method_info:
                continue

            method, raw_score, normalized_score = get_best_method_and_score(combined_method_info[gene_id])

            if gene_id not in gene_assignments:
                gene_assignments[gene_id] = {'type': 'antitoxin', 'method': method, 'score': raw_score, 'normalized_score': normalized_score}
            else:
                # Already assigned - check if this has higher normalized score
                existing = gene_assignments[gene_id]
                if normalized_score > existing['normalized_score']:
                    gene_assignments[gene_id] = {'type': 'antitoxin', 'method': method, 'score': raw_score, 'normalized_score': normalized_score}

    # Build final sets based on resolved assignments
    all_toxins = set()
    all_antitoxins = set()

    for gene_id, assignment in gene_assignments.items():
        if assignment['type'] == 'toxin':
            all_toxins.add(gene_id)
        else:
            all_antitoxins.add(gene_id)

    return all_toxins, all_antitoxins, combined_method_info

def filter_same_type_TA(df_result):
    """Remove pairs where both members have different TA types."""
    if df_result.empty or 'TA_type' not in df_result.columns:
        return df_result

    ta_id_to_remove = []
    for i in range(0, df_result.shape[0], 2):
        if i+1 >= df_result.shape[0]:
            break
        if df_result.loc[i, 'TA_type'] != df_result.loc[i+1, 'TA_type']:
            ta_id_to_remove.append(df_result.loc[i, 'TA_ID'])

    df_result = df_result[~df_result.TA_ID.isin(ta_id_to_remove)]
    df_result.index = pd.RangeIndex(start=0, stop=df_result.shape[0])

    if df_result.shape[0] > 0:
        # Handle odd number of rows (incomplete pairs)
        num_complete_pairs = df_result.shape[0] // 2
        ta_ids = []
        for i in range(num_complete_pairs):
            ta_ids.extend([i + 1, i + 1])
        if df_result.shape[0] % 2 == 1:  # Odd number
            ta_ids.append(num_complete_pairs + 1)
        df_result['TA_ID'] = ta_ids
    return df_result

def validate_complete_pairs(df_result):
    """Ensure every TA_ID has exactly one Toxin and one Antitoxin."""
    if len(df_result) == 0:
        return df_result

    valid_ta_ids = []
    for ta_id in df_result['TA_ID'].unique():
        pair = df_result[df_result['TA_ID'] == ta_id]
        if len(pair) == 2:
            types = set(pair['Type'].values)
            # Must have both 'Toxin' and 'Antitoxin'
            if types == {'Toxin', 'Antitoxin'}:
                valid_ta_ids.append(ta_id)

    initial_rows = len(df_result)
    df_result = df_result[df_result['TA_ID'].isin(valid_ta_ids)]
    removed_rows = initial_rows - len(df_result)
    print(f"Removed {removed_rows // 2} incomplete/invalid pairs ({removed_rows} rows) → kept {len(df_result) // 2} valid TA pairs ({len(df_result)} rows)")

    # Reset index and renumber TA_IDs sequentially
    df_result = df_result.sort_values('TA_ID').reset_index(drop=True)
    if len(df_result) > 0:
        ta_ids = []
        for i in range(len(df_result) // 2):
            ta_ids.extend([i + 1, i + 1])
        df_result['TA_ID'] = ta_ids
    return df_result

def pair_ta_systems(toxin_genes, antitoxin_genes, neighborhoods, max_pairing_distance=150):
    ptt = neighborhoods['ptt']
    has_contig = "Contig" in ptt.columns

    ta_pairs = []
    paired_genes = set()
    ta_id = 1

    for nbh_id, genes in neighborhoods['neighborhood_to_genes'].items():
        nbh_toxins = set(genes).intersection(toxin_genes)
        nbh_antitoxins = set(genes).intersection(antitoxin_genes)

        for tox_gene in nbh_toxins:
            for at_gene in nbh_antitoxins:
                # Check RNA status
                tox_is_rna = is_ncrna(tox_gene)
                at_is_rna = is_ncrna(at_gene)

                if tox_is_rna and not at_is_rna:
                    continue  

                tox_row = ptt[ptt['PID'] == tox_gene].iloc[0]
                at_row = ptt[ptt['PID'] == at_gene].iloc[0]

                tox_loc = tox_row['Location']
                at_loc = at_row['Location']

                tox_start, tox_end = map(int, tox_loc.split('..'))
                at_start, at_end = map(int, at_loc.split('..'))

                if tox_start < at_start:
                    intergenic_distance = at_start - tox_end
                else:
                    intergenic_distance = tox_start - at_end

                # Allow overlapping genes (negative distance is OK)
                if intergenic_distance > max_pairing_distance:
                    continue

                # Determine TA system type based on gene characteristics
                ta_type = None
                if not tox_is_rna and at_is_rna:
                    # Protein toxin + RNA antitoxin
                    if at_gene.startswith('Nucl_type1_at_'):
                        ta_type = 'Type_I'
                    elif at_gene.startswith('Nucl_type3_at_'):
                        ta_type = 'Type_III'
                elif tox_is_rna and at_is_rna:
                    # RNA toxin + RNA antitoxin (Type VIII)
                    if tox_gene.startswith('Nucl_type8_') or at_gene.startswith('Nucl_type8_'):
                        ta_type = 'Type_VIII'
                else:
                    # Protein toxin + protein antitoxin (Type II, IV, V, VI, VII, etc.)
                    ta_type = 'Type_II_VII'  # Generic for protein-protein systems

                # Build entries using helper function
                tox_entry = build_gene_entry(tox_row, 'Toxin', ta_id, has_contig)
                at_entry = build_gene_entry(at_row, 'Antitoxin', ta_id, has_contig)

                # Add TA_type to both entries
                if ta_type:
                    tox_entry['TA_type'] = ta_type
                    at_entry['TA_type'] = ta_type

                ta_pairs.append(tox_entry)
                ta_pairs.append(at_entry)

                paired_genes.update([tox_gene, at_gene])
                ta_id += 1

    orphan_genes = (toxin_genes.union(antitoxin_genes)) - paired_genes
    print(f"Found {len(ta_pairs)//2} TA pairs, {len(orphan_genes)} orphans")

    df_result = pd.DataFrame(ta_pairs)

    return df_result, orphan_genes
