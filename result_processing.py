from Bio.Seq import Seq
import pandas as pd
from utils import is_ncrna, build_gene_entry

# RNA-mRNA complementarity checking for Type I systems
def has_rna_mrna_complementarity(rna_seq, mrna_seq, window_size=30, max_mismatches=2):
    try:
        rna_rc = str(Seq(rna_seq).reverse_complement())
    except (ValueError, TypeError, AttributeError):
        return False

    # Slide window across rna_rc
    for i in range(len(rna_rc) - window_size + 1):
        window = rna_rc[i:i + window_size]

        # Check this window against all positions in mrna_seq
        for j in range(len(mrna_seq) - window_size + 1):
            mrna_window = mrna_seq[j:j + window_size]

            # Count mismatches
            mismatches = sum(1 for a, b in zip(window, mrna_window) if a != b)

            # Early return if match found
            if mismatches <= max_mismatches:
                return True

    return False

def validate_type1_complementarity(neighborhoods, rna_results, genome_seq, min_complementarity=30):
    ptt = neighborhoods['ptt']
    validation_status = {}  # gene_id -> "validated" or "putative"
    neighborhoods_to_remove = set()

    print(f"\nValidating Type I TA systems (RNA antitoxin + protein toxin complementarity)...")

    for nbh_id, genes in list(neighborhoods['neighborhood_to_genes'].items()):
        # Find RNA antitoxins and protein genes in this neighborhood
        rna_antitoxins = [g for g in genes if is_ncrna(g)]
        proteins = [g for g in genes if not is_ncrna(g)]

        if not rna_antitoxins or not proteins:
            continue  # No Type I system possible

        # Check each RNA-protein pair
        for rna_gene in rna_antitoxins:
            rna_row = ptt[ptt['PID'] == rna_gene].iloc[0]
            for prot_gene in proteins:
                prot_row = ptt[ptt['PID'] == prot_gene].iloc[0]

                # Only check opposite strand pairs
                if rna_row['Strand'] == prot_row['Strand']:
                    continue

                # Extract sequences
                try:
                    rna_start, rna_end = map(int, rna_row['Location'].split('..'))
                    prot_start, prot_end = map(int, prot_row['Location'].split('..'))
                    rna_seq = genome_seq[rna_start-1:rna_end]
                    mrna_seq = genome_seq[prot_start-1:prot_end]

                    # Check complementarity (with 2bp mismatch tolerance)
                    if has_rna_mrna_complementarity(rna_seq, mrna_seq, window_size=min_complementarity):
                        # Has complementarity! Determine status based on MMseqs2 detection
                        original_rna_id = rna_gene.replace('Nucl_', '') if rna_gene.startswith('Nucl_') else rna_gene
                        status = "validated" if original_rna_id in rna_results.get('validated_rna_genes', set()) else "putative"
                        validation_status[rna_gene] = status
                        validation_status[prot_gene] = status
                    else:
                        # No complementarity - mark neighborhood for removal
                        neighborhoods_to_remove.add(nbh_id)

                except (ValueError, IndexError, KeyError):
                    continue

    # Remove neighborhoods without complementarity
    for nbh_id in neighborhoods_to_remove:
        if nbh_id in neighborhoods['neighborhood_to_genes']:
            genes = neighborhoods['neighborhood_to_genes'][nbh_id]
            # Remove from all structures
            del neighborhoods['neighborhood_to_genes'][nbh_id]
            for gene in genes:
                if gene in neighborhoods['gene_to_neighborhood']:
                    del neighborhoods['gene_to_neighborhood'][gene]
                if gene in neighborhoods['proteins_in_neighborhoods']:
                    neighborhoods['proteins_in_neighborhoods'].remove(gene)

    validated_count = sum(1 for v in validation_status.values() if v == "validated")
    putative_count = sum(1 for v in validation_status.values() if v == "putative")
    print(f"Type I validation: {validated_count} validated, {putative_count} putative, {len(neighborhoods_to_remove)} neighborhoods removed")

    return validation_status

def detect_method_specific_self_ta(combined_method_info, args):
    self_ta_candidates = {}

    for gene_id, info in combined_method_info.items():
        methods = info.get('methods', [])

        # Check each method for identical T/AT scores
        for method in methods:
            method_key = method.lower().replace('_', '').replace('2', '')

            # Get the bitscore key for this method
            if 'foldseek' in method.lower():
                score_key = 'foldseek_bitscore'
                base_cutoff = args.foldseek_bitscore
            elif 'hmm' in method.lower():
                score_key = 'hmm_bitscore'
                base_cutoff = args.hmm_bitscore
            elif 'mmseqs' in method.lower():
                if 'dna' in method.lower():
                    score_key = 'mmseqs2_dna_bitscore'
                else:
                    score_key = 'mmseqs2_bitscore'
                base_cutoff = args.mmseqs_bitscore
            else:
                continue  # Skip unknown methods

            # Check if this gene has a score for this method
            if score_key in info:
                raw_score = info[score_key]
                higher_cutoff = base_cutoff * 3.0

                # Check if score meets 3x cutoff
                if raw_score >= higher_cutoff:
                    # This is a potential self-TA candidate
                    self_ta_candidates[gene_id] = {
                        'method': method,
                        'score': raw_score,
                        'cutoff': higher_cutoff
                    }
                    break  # One qualifying method is enough

    return self_ta_candidates

def merge_search_results(*search_results, args=None):
    # Track which genes are toxins/antitoxins and their normalized scores
    gene_assignments = {}  # gene_id -> {'type': 'toxin'/'antitoxin', 'method': str, 'score': float, 'normalized_score': float}
    combined_method_info = {}

    # Track method-specific scores for self-TA detection
    method_specific_scores = {}  # gene_id -> {method -> {'toxin_score': float, 'antitoxin_score': float}}

    def normalize_score(method, score):
        if 'foldseek' in method.lower():
            return score / 10.0
        elif 'mmseqs' in method.lower():
            return score / 2.0
        else:  # HMM
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

            # Track method-specific scores from this result's method_info
            # This captures the actual method that found this gene
            for method in info.get('methods', []):
                if gene_id not in method_specific_scores:
                    method_specific_scores[gene_id] = {}
                if method not in method_specific_scores[gene_id]:
                    method_specific_scores[gene_id][method] = {}

                # Get the score for this specific method
                if 'foldseek' in method.lower():
                    score_key = 'foldseek_bitscore'
                elif 'hmm' in method.lower():
                    score_key = 'hmm_bitscore'
                elif 'dna' in method.lower():
                    score_key = 'mmseqs2_dna_bitscore'
                else:
                    score_key = 'mmseqs2_bitscore'

                if score_key in info:
                    method_specific_scores[gene_id][method]['score'] = info[score_key]

        # Track toxin assignments
        for gene_id in result['toxins']:
            if gene_id not in combined_method_info:
                continue

            # Mark this gene as toxin for each method in this result
            if gene_id in result['method_info']:
                for method in result['method_info'][gene_id].get('methods', []):
                    if gene_id in method_specific_scores and method in method_specific_scores[gene_id]:
                        method_specific_scores[gene_id][method]['is_toxin'] = True

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

            # Mark this gene as antitoxin for each method in this result
            if gene_id in result['method_info']:
                for method in result['method_info'][gene_id].get('methods', []):
                    if gene_id in method_specific_scores and method in method_specific_scores[gene_id]:
                        method_specific_scores[gene_id][method]['is_antitoxin'] = True

            method, raw_score, normalized_score = get_best_method_and_score(combined_method_info[gene_id])

            if gene_id not in gene_assignments:
                gene_assignments[gene_id] = {'type': 'antitoxin', 'method': method, 'score': raw_score, 'normalized_score': normalized_score}
            else:
                # Already assigned - check if this has higher normalized score
                existing = gene_assignments[gene_id]
                if normalized_score > existing['normalized_score']:
                    gene_assignments[gene_id] = {'type': 'antitoxin', 'method': method, 'score': raw_score, 'normalized_score': normalized_score}

    # Detect method-specific self-TA candidates
    # A gene is self-TA if the SAME method found it as both toxin AND antitoxin
    self_ta_candidates = {}
    for gene_id, method_scores in method_specific_scores.items():
        for method, scores in method_scores.items():
            is_toxin = scores.get('is_toxin', False)
            is_antitoxin = scores.get('is_antitoxin', False)
            score = scores.get('score', 0)

            # Check if this method found the gene as both toxin and antitoxin
            if is_toxin and is_antitoxin:
                # Determine method cutoff
                if args:
                    if 'foldseek' in method.lower():
                        base_cutoff = args.foldseek_bitscore
                    elif 'hmm' in method.lower():
                        base_cutoff = args.hmm_bitscore
                    else:  # MMseqs variants
                        base_cutoff = args.mmseqs_bitscore

                    higher_cutoff = base_cutoff * 3.0

                    # Check if score meets 3x cutoff
                    if score >= higher_cutoff:
                        self_ta_candidates[gene_id] = {
                            'method': method,
                            'score': score,
                            'cutoff': higher_cutoff
                        }
                        # Mark in method info for output
                        combined_method_info[gene_id]['is_self_ta'] = True
                        break  # One qualifying method is enough

    # Build final sets based on resolved assignments
    all_toxins = set()
    all_antitoxins = set()

    for gene_id, assignment in gene_assignments.items():
        # Skip self-TA genes from regular assignment (they go in both sets)
        if gene_id in self_ta_candidates:
            continue
        if assignment['type'] == 'toxin':
            all_toxins.add(gene_id)
        else:
            all_antitoxins.add(gene_id)

    # Add self-TA genes to BOTH sets so they can pair with neighbors
    # The pairing function will assign them the opposite role of their partner
    for gene_id in self_ta_candidates:
        all_toxins.add(gene_id)
        all_antitoxins.add(gene_id)

    # Report results
    if self_ta_candidates:
        print(f"  Found {len(self_ta_candidates)} potential self-TA genes (identical T/AT scores ≥3x cutoff)")

    return all_toxins, all_antitoxins, combined_method_info, set(self_ta_candidates.keys())

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

    initial_pairs = len(df_result)
    df_result = df_result[df_result['TA_ID'].isin(valid_ta_ids)]
    print(f"Removed {initial_pairs - len(df_result)} incomplete/invalid pairs (kept {len(df_result)} valid rows = {len(df_result)//2} pairs)")

    # Reset index and renumber TA_IDs sequentially
    df_result = df_result.sort_values('TA_ID').reset_index(drop=True)
    if len(df_result) > 0:
        ta_ids = []
        for i in range(len(df_result) // 2):
            ta_ids.extend([i + 1, i + 1])
        df_result['TA_ID'] = ta_ids
    return df_result

def pair_ta_systems(toxin_genes, antitoxin_genes, neighborhoods, self_ta_genes=None):
    ptt = neighborhoods['ptt']
    has_contig = "Contig" in ptt.columns

    if self_ta_genes is None:
        self_ta_genes = set()

    ta_pairs = []
    paired_genes = set()
    paired_self_ta_genes = set()  # Track self-TA genes that got paired
    ta_id = 1

    for nbh_id, genes in neighborhoods['neighborhood_to_genes'].items():
        nbh_toxins = set(genes).intersection(toxin_genes)
        nbh_antitoxins = set(genes).intersection(antitoxin_genes)

        for tox_gene in nbh_toxins:
            for at_gene in nbh_antitoxins:
                # Check if toxin is RNA
                tox_is_rna = is_ncrna(tox_gene)
                at_is_rna = is_ncrna(at_gene)

                if tox_is_rna and not at_is_rna:
                    continue  # Skip this pairing

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

                if intergenic_distance > 200:
                    continue

                tox_is_self_ta = tox_gene in self_ta_genes
                at_is_self_ta = at_gene in self_ta_genes

                # Case 1: Both are self-TA → mark as "ambiguous"
                if tox_is_self_ta and at_is_self_ta:
                    tox_type = 'Ambiguous'
                    at_type = 'Ambiguous'
                    paired_self_ta_genes.add(tox_gene)
                    paired_self_ta_genes.add(at_gene)
                # Case 2: Only toxin gene is self-TA → partner is antitoxin, so self-TA takes toxin role
                elif tox_is_self_ta:
                    tox_type = 'Toxin'
                    at_type = 'Antitoxin'
                    paired_self_ta_genes.add(tox_gene)
                # Case 3: Only antitoxin gene is self-TA → partner is toxin, so self-TA takes antitoxin role
                elif at_is_self_ta:
                    tox_type = 'Toxin'
                    at_type = 'Antitoxin'
                    paired_self_ta_genes.add(at_gene)
                # Case 4: Neither is self-TA → normal pairing
                else:
                    tox_type = 'Toxin'
                    at_type = 'Antitoxin'

                # Build entries using helper function
                ta_pairs.append(build_gene_entry(tox_row, tox_type, ta_id, has_contig))
                ta_pairs.append(build_gene_entry(at_row, at_type, ta_id, has_contig))

                paired_genes.update([tox_gene, at_gene])
                ta_id += 1

    # Calculate final unpaired self-TA genes (those that didn't get paired)
    final_self_ta_genes = self_ta_genes - paired_self_ta_genes

    orphan_genes = (toxin_genes.union(antitoxin_genes)) - paired_genes
    print(f"Found {len(ta_pairs)//2} TA pairs, {len(orphan_genes)} orphans")
    if final_self_ta_genes:
        print(f"Found {len(final_self_ta_genes)} unpaired self-TA genes")

    df_result = pd.DataFrame(ta_pairs)

    return df_result, orphan_genes, final_self_ta_genes
