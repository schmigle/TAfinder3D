#!/usr/bin/env python3
"""
TAfinder - Modular Version
Identifies toxin-antitoxin systems in bacterial genomes using sequence and structure-based searches.

This is a refactored version that splits the original monolithic script into modular components.
The original tafinder_foldseek_simplified.py remains unchanged.
"""

from Bio import SeqIO
import os

# Import configuration and utilities
from config import parse_args, initialize_environment
from utils import is_ncrna, ensure_prostt5_weights

# Import step modules
from input_and_neighborhoods import (
    process_input_files,
    build_neighborhoods_with_rna_support
)
from combined_search import (
    run_rna_mmseqs_searches,
    run_mmseqs_searches,
    run_dna_mmseqs_searches,
    run_hmm_searches,
    run_foldseek_searches
)
from result_processing import (
    validate_type1_complementarity,
    merge_search_results
)
from output_generation import create_all_outputs


def main():
    args = parse_args()
    config = initialize_environment(args)

    # Ensure ProstT5 weights are available (only if Foldseek enabled and no pre-built query DB)
    if not args.no_foldseek and not args.foldseek_query_db:
        ensure_prostt5_weights(config['prostt5_weights'])

    print("Processing input files...")
    process_input_files(args, config)

    # STEP 1: Run RNA MMseqs2 searches FIRST to identify RNA components
    print("\n=== STEP 1: RNA Detection ===")
    rna_results = run_rna_mmseqs_searches(config['fna_file'], config['ptt_file'], config, args)

    # STEP 2: Build neighborhoods with RNA support (including MMseqs2-detected RNAs)
    print("\n=== STEP 2: Building Neighborhoods ===")
    neighborhoods = build_neighborhoods_with_rna_support(config['ptt_file'], args, rna_results)

    if args.search_all:
        search_file = config['faa_file']
        print(f"Searching all proteins (--search_all enabled)")
    else:
        search_file = f"{args.tmp_dir}/neighborhood_proteins.faa"
        with open(search_file, 'w') as f:
            for record in SeqIO.parse(config['faa_file'], "fasta"):
                if record.id in neighborhoods['proteins_in_neighborhoods']:
                    SeqIO.write(record, f, "fasta")

    print("\n=== STEP 3: DNA Searches ===")

    cds_nuc_file = f"{args.tmp_dir}/sequence.ffa"
    if not os.path.exists(cds_nuc_file):
        ptt = neighborhoods['ptt']
        genome_seqs = {r.id: r.seq for r in SeqIO.parse(config['fna_file'], "fasta")}
        with open(cds_nuc_file, 'w') as f:
            for idx in range(len(ptt)):
                gene_id = ptt.loc[idx, 'PID']
                if is_ncrna(gene_id):
                    continue
                loc = ptt.loc[idx, 'Location']
                strand = ptt.loc[idx, 'Strand']
                contig = ptt.loc[idx, 'Contig'] if 'Contig' in ptt.columns else list(genome_seqs.keys())[0]
                start, end = map(int, loc.split('..'))
                seq = genome_seqs[contig][start-1:end]
                if strand == '-':
                    seq = seq.reverse_complement()
                f.write(f">{gene_id}\n{str(seq)}\n")
    dna_results = run_dna_mmseqs_searches(cds_nuc_file, config, args)
    if dna_results['toxins'] or dna_results['antitoxins']:
        print(f"  DNA results to merge: {len(dna_results['toxins'])} toxins, {len(dna_results['antitoxins'])} antitoxins")
        print(f"  DNA toxin IDs: {list(dna_results['toxins'])[:5]}")
        print(f"  DNA antitoxin IDs: {list(dna_results['antitoxins'])[:5]}")

    # STEP 4: Run protein MMseqs2 searches
    print("\n=== STEP 4: Protein MMseqs2 Searches ===")
    mmseqs_results = run_mmseqs_searches(search_file, config, args)

    # STEP 5: Run HMM searches
    print("\n=== STEP 5: HMM Searches ===")
    hmm_results = run_hmm_searches(search_file, config, args)

    # Add RNA results to the tracking
    # RNA genes are now in rna_results with their new gene IDs
    rna_mmseqs_results = {
        'toxins': set(rna_results['rna_toxin_locations'].keys()),
        'antitoxins': set(rna_results['rna_antitoxin_locations'].keys()),
        'method_info': {}
    }

    # Add method info for RNA hits with actual scores
    for gene in rna_results['rna_toxin_locations'].keys():
        evalue, bitscore = rna_results['rna_scores'].get(gene, (0.0, 0.0))
        rna_mmseqs_results['method_info'][gene] = {
            'methods': ['MMseqs2_RNA'],
            'mmseqs2_rna_evalue': evalue,
            'mmseqs2_rna_bitscore': bitscore
        }
    for gene in rna_results['rna_antitoxin_locations'].keys():
        evalue, bitscore = rna_results['rna_scores'].get(gene, (0.0, 0.0))
        rna_mmseqs_results['method_info'][gene] = {
            'methods': ['MMseqs2_RNA'],
            'mmseqs2_rna_evalue': evalue,
            'mmseqs2_rna_bitscore': bitscore
        }

    # STEP 6: Run Foldseek on unidentified
    print("\n=== STEP 6: Foldseek Structural Searches ===")
    identified_genes = mmseqs_results['toxins'].union(mmseqs_results['antitoxins']).union(
        hmm_results['toxins']).union(hmm_results['antitoxins']).union(
        rna_mmseqs_results['toxins']).union(rna_mmseqs_results['antitoxins']).union(
        dna_results['toxins']).union(dna_results['antitoxins'])

    # Determine Foldseek candidates based on search mode
    if args.search_all:
        # Search all proteins (excluding ncRNAs) when --search_all is enabled
        all_proteins = set(neighborhoods['ptt']['PID'].tolist())
        all_proteins = {p for p in all_proteins if not is_ncrna(p)}
        candidate_proteins = all_proteins - identified_genes
    else:
        # Only search proteins in neighborhoods
        candidate_proteins = set(neighborhoods['proteins_in_neighborhoods']) - identified_genes

    foldseek_results = run_foldseek_searches(candidate_proteins, config, args)

    # Merge all results (including RNA and DNA results)
    final_results = merge_search_results(rna_mmseqs_results, mmseqs_results, hmm_results, foldseek_results, dna_results, args=args)

    print(f"\n=== Final Results ===")
    print(f"Total toxins found: {len(final_results[0])}")
    print(f"Total antitoxins found: {len(final_results[1])}")
    if len(final_results) > 3 and final_results[3]:
        print(f"Self-TA genes found: {len(final_results[3])}")

    # STEP 7: Validate Type I complementarity (AFTER all searches complete)
    print("\n=== STEP 7: Type I Validation ===")
    toxin_genes, antitoxin_genes = final_results[0], final_results[1]
    type1_validation_status = validate_type1_complementarity(
        neighborhoods, rna_results, neighborhoods['genome_seq'],
        toxin_genes, antitoxin_genes, min_complementarity=30
    )

    # Create outputs
    print("\n=== Creating Outputs ===")
    create_all_outputs(final_results, neighborhoods, config, args, type1_validation_status)

    print("Analysis complete!")


if __name__ == "__main__":
    main()
