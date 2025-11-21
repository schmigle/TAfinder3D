from Bio import SeqIO
import pandas as pd
import os
import subprocess
import shutil
from config import MAX_NCRNA_DISTANCE
from utils import (determine_cutoff_params, parse_search_output,
                   find_overlapping_gene_in_ptt, update_method_info, is_ncrna)

def run_single_rna_mmseqs_search(fna_file, db_path, output_file, args, ptt, has_contig, cutoff_value, use_evalue, evalue_param):
    locations = {}
    all_genes = set()

    cmd = ["mmseqs", "easy-search", fna_file, db_path, output_file, args.tmp_dir,
           "--format-output", "query,target,pident,qlen,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits",
           "--threads", str(args.threads), "--search-type", "3",  # 3 = nucleotide search
           "-e", evalue_param]

    try:
        subprocess.run(cmd, capture_output=True, text=True, check=True)
        df = parse_search_output(output_file, 'mmseqs_rna', cutoff_value, use_evalue)

        if not df.empty:
            for _, row in df.iterrows():
                query_name = str(row[0])  # Query is contig/sequence name
                query_start, query_end = int(row[7]), int(row[8])
                # Find overlapping genes in PTT
                gene_id, start, end, strand, contig = find_overlapping_gene_in_ptt(ptt, query_start, query_end, has_contig, query_name)
                if gene_id:
                    locations[gene_id] = (start, end, strand, contig)
                    all_genes.add(gene_id)
    except subprocess.CalledProcessError:
        pass

    return locations, all_genes

def run_rna_mmseqs_searches(fna_file, ptt_file, config, args):
    print("Running RNA MMseqs2 searches...")

    rna_results = {
        'rna_toxin_locations': {},      # gene_id -> (start, end, strand, contig)
        'rna_antitoxin_locations': {},  # gene_id -> (start, end, strand, contig)
        'validated_rna_genes': set(),   # Only RNA toxins with paired antitoxins
        'all_rna_genes': set(),         # All detected RNAs (for tracking)
        'type8_genes': set(),           # Type VIII RNA genes
        'type3_genes': set()            # Type III RNA genes (neighbor rule)
    }

    if args.no_mmseqs:
        return rna_results

    # Read PTT to get coordinates for matching
    ptt = pd.read_csv(ptt_file, sep="\t", skiprows=2)
    has_contig = "Contig" in ptt.columns

    # Determine cutoff value and mode
    cutoff_value, use_evalue, evalue_param = determine_cutoff_params(args, 'mmseqs')

    # STEP 1: Type VIII RNA systems (searched first, separately)
    print("  Searching Type VIII RNA databases...")
    type8_tox_locations = {}
    type8_at_locations = {}

    if 'mmseqs_type8_tox_db' in config:
        type8_tox_output = f"{args.tmp_dir}/mmseqs_type8_T.out"
        type8_tox_locations, type8_tox_genes = run_single_rna_mmseqs_search(
            fna_file, config['mmseqs_type8_tox_db'], type8_tox_output, args, ptt, has_contig, cutoff_value, use_evalue, evalue_param
        )
        rna_results['all_rna_genes'].update(type8_tox_genes)

    if 'mmseqs_type8_antitox_db' in config:
        type8_at_output = f"{args.tmp_dir}/mmseqs_type8_AT.out"
        type8_at_locations, type8_at_genes = run_single_rna_mmseqs_search(
            fna_file, config['mmseqs_type8_antitox_db'], type8_at_output, args, ptt, has_contig, cutoff_value, use_evalue, evalue_param
        )
        rna_results['all_rna_genes'].update(type8_at_genes)

    # Validate Type VIII pairs (neighbor rule: within 200bp)
    for tox_gene, (tox_start, tox_end, _, tox_contig) in type8_tox_locations.items():
        for at_gene, (at_start, at_end, _, at_contig) in type8_at_locations.items():
            if tox_contig != at_contig:
                continue
            if tox_start < at_start:
                intergenic_distance = at_start - tox_end
            else:
                intergenic_distance = tox_start - at_end
            if intergenic_distance <= MAX_NCRNA_DISTANCE:
                rna_results['validated_rna_genes'].add(tox_gene)
                rna_results['validated_rna_genes'].add(at_gene)
                rna_results['type8_genes'].add(tox_gene)
                rna_results['type8_genes'].add(at_gene)

    # Run MMseqs2 for RNA antitoxins only (Type I/III toxins are proteins)
    rna_at_output = f"{args.tmp_dir}/mmseqs_rna_AT.out"
    rna_results['rna_antitoxin_locations'], at_genes = run_single_rna_mmseqs_search(
        fna_file, config['mmseqs_rna_antitox_db'], rna_at_output, args, ptt, has_contig, cutoff_value, use_evalue, evalue_param
    )
    rna_results['all_rna_genes'].update(at_genes)

    # Merge Type VIII locations into main results
    rna_results['rna_toxin_locations'].update(type8_tox_locations)
    rna_results['rna_antitoxin_locations'].update(type8_at_locations)

    # Type I/III RNA antitoxins are validated later when paired with protein toxins
    # (validation happens in main after protein searches complete)

    return rna_results

def run_mmseqs_searches(search_file, config, args):
    results = {'toxins': set(), 'antitoxins': set(), 'method_info': {}}

    if args.no_mmseqs:
        return results

    mmseqs_configs = [
        ('T', config['mmseqs_tox_db']),
        ('AT', config['mmseqs_antitox_db'])
    ]

    # Determine cutoff value and mode
    cutoff_value, use_evalue, evalue_param = determine_cutoff_params(args, 'mmseqs')

    # # Check if search file exists and has content
    # if not os.path.exists(search_file):
    #     print(f"  ERROR: Search file does not exist: {search_file}")
    #     return results

    # file_size = os.path.getsize(search_file)
    # if file_size == 0:
    #     print(f"  ERROR: Search file is empty: {search_file}")
    #     return results

    # Count sequences in file
    seq_count = sum(1 for line in open(search_file) if line.startswith('>'))
    file_size = os.path.getsize(search_file)
    print(f"  Search file: {search_file} ({seq_count} sequences, {file_size} bytes)")

    for target, db_path in mmseqs_configs:
        output_file = f"{args.tmp_dir}/mmseqs_{target}.out"

        # MMseqs2 easy-search command
        # Format: query target output tmpDir
        # Output format: query, target, pident, alnlen, mismatch, gapopen, qstart, qend, tstart, tend, evalue, bits
        cmd = ["mmseqs", "easy-search", search_file, db_path, output_file, args.tmp_dir,
               "--format-output", "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits",
               "--threads", str(args.threads), "--gpu", "1",
               "-e", evalue_param]

        try:
            subprocess.run(cmd, capture_output=True, text=True, check=True)
            df = parse_search_output(output_file, 'mmseqs', cutoff_value, use_evalue, args.havalue, args.maximum_length)

            if not df.empty:
                gene_set = results['toxins'] if target == 'T' else results['antitoxins']
                for _, row in df.iterrows():
                    gene_id = row.iloc[0]
                    gene_set.add(gene_id)
                    update_method_info(results, gene_id, 'MMseqs2', row.iloc[6], row.iloc[7], use_evalue)
        except subprocess.CalledProcessError as e:
            print(f"  MMseqs2 {target} search failed: {e}")
            continue

    print(f"  Protein MMseqs2 found {len(results['toxins'])} toxins, {len(results['antitoxins'])} antitoxins")
    return results

def run_dna_mmseqs_searches(cds_nuc_file, config, args):
    results = {'toxins': set(), 'antitoxins': set(), 'method_info': {}}

    if args.no_mmseqs:
        return results

    if 'mmseqs_dna_tox_db' not in config or 'mmseqs_dna_antitox_db' not in config:
        return results

    print("Running MMseqs2 DNA searches...")

    # Determine cutoff value and mode
    cutoff_value, use_evalue, evalue_param = determine_cutoff_params(args, 'mmseqs')

    dna_configs = [
        ('T', config['mmseqs_dna_tox_db']),
        ('AT', config['mmseqs_dna_antitox_db'])
    ]

    for target, db_path in dna_configs:
        output_file = f"{args.tmp_dir}/mmseqs_dna_{target}.out"

        # Search CDS nucleotide sequences against DNA database
        cmd = ["mmseqs", "easy-search", cds_nuc_file, db_path, output_file, args.tmp_dir,
               "--format-output", "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits",
               "--threads", str(args.threads), "--search-type", "3",  # nucleotide search
               "-e", evalue_param]

        try:
            subprocess.run(cmd, capture_output=True, text=True, check=True)
            df = parse_search_output(output_file, 'mmseqs_rna', cutoff_value, use_evalue)

            if not df.empty:
                gene_set = results['toxins'] if target == 'T' else results['antitoxins']
                for _, row in df.iterrows():
                    # Query ID is now the ORF name directly (e.g., ORF123)
                    gene_id = str(row[0])
                    gene_set.add(gene_id)
                    update_method_info(results, gene_id, 'MMseqs2_DNA', row[10], row[11], use_evalue)

        except subprocess.CalledProcessError:
            continue

    print(f"  DNA search found {len(results['toxins'])} toxins, {len(results['antitoxins'])} antitoxins")
    return results

def run_hmm_searches(search_file, config, args):
    results = {'toxins': set(), 'antitoxins': set(), 'method_info': {}}

    if args.no_hmm:
        return results

    hmm_configs = [
        ('T', config['hmm_tox_file']),
        ('AT', config['hmm_antitox_file'])
    ]

    # Determine cutoff value and mode
    cutoff_value, use_evalue, evalue_param = determine_cutoff_params(args, 'hmm')

    for target, hmm_file in hmm_configs:
        output_file = f"{args.tmp_dir}/hmmsearch_{target}.domtbl"
        cmd = ["hmmsearch", "--cpu", str(args.threads), "--noali", "-E", evalue_param,
               "--domE", evalue_param, "-o", f"{args.tmp_dir}/hmmsearch.txt",
               "--domtblout", output_file, hmm_file, search_file]

        try:
            subprocess.run(cmd, capture_output=True, text=True, check=True)

            df = parse_search_output(output_file, 'hmm', cutoff_value, use_evalue, None, args.maximum_length)

            if not df.empty:
                gene_set = results['toxins'] if target == 'T' else results['antitoxins']
                for _, row in df.iterrows():
                    gene_id = row.iloc[0]
                    gene_set.add(gene_id)
                    update_method_info(results, gene_id, 'HMM', row.iloc[6], row.iloc[7], use_evalue)
        except subprocess.CalledProcessError as e:
            print(f"  HMM {target} search failed: {e}")
            continue

    return results

def run_foldseek_searches(candidate_proteins, config, args):
    """Run Foldseek searches on candidates using createdb -> search -> convertalis pipeline."""
    results = {'toxins': set(), 'antitoxins': set(), 'method_info': {}}

    if args.no_foldseek or not candidate_proteins:
        return results

    foldseek_path = shutil.which('foldseek')
    if not foldseek_path:
        print("Foldseek not found in PATH")
        return results

    # Setup Foldseek weights
    if not os.path.exists(config['prostt5_weights']):
        try:
            subprocess.run([foldseek_path, "databases", "ProstT5", config['prostt5_weights'], args.tmp_dir],
                         check=True, capture_output=True)
        except subprocess.CalledProcessError:
            return results

    # Create candidate FASTA
    candidate_faa = f"{args.tmp_dir}/candidates.faa"
    protein_seqs = {r.id: str(r.seq) for r in SeqIO.parse(config['faa_file'], "fasta") if r.id in candidate_proteins}

    if not protein_seqs:
        return results

    with open(candidate_faa, 'w') as f:
        for prot_id, seq in protein_seqs.items():
            f.write(f">{prot_id}\n{seq}\n")

    # Step 1: Create query database with ProstT5
    foldseek_query_db = f"{args.tmp_dir}/foldseek_query_db"
    try:
        subprocess.run([foldseek_path, "createdb", candidate_faa, foldseek_query_db,
                       "--prostt5-model", config['prostt5_weights'],
                       "--gpu", "1"],
                      capture_output=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Foldseek createdb failed: {e}")
        return results

    # Step 2: Make padded sequence database
    foldseek_query_db_padded = f"{args.tmp_dir}/foldseek_query_db_padded"
    try:
        subprocess.run([foldseek_path, "makepaddedseqdb", foldseek_query_db, foldseek_query_db_padded],
                      capture_output=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Foldseek makepaddedseqdb failed: {e}")
        return results

    # Determine cutoff value and mode
    cutoff_value, use_evalue, evalue_param = determine_cutoff_params(args, 'foldseek')

    # Step 3 & 4: Search both databases (already padded) and convert results
    for target, db_path in [('toxins', config['foldseek_tox_db']), ('antitoxins', config['foldseek_antitox_db'])]:
        foldseek_result_db = f"{args.tmp_dir}/foldseek_{target}_result"
        output_file = f"{args.tmp_dir}/foldseek_{target}.tsv"

        try:
            # Step 3: Run foldseek search on padded databases
            subprocess.run([foldseek_path, "search", foldseek_query_db_padded, db_path, foldseek_result_db, args.tmp_dir,
                           "--gpu", "1",
                           "-e", evalue_param,
                           "-s", str(args.foldseeksensitivity)],
                          capture_output=True, check=True)

            # Step 4: Convert alignment to TSV format
            subprocess.run([foldseek_path, "convertalis", foldseek_query_db_padded, db_path, foldseek_result_db, output_file,
                           "--format-output", "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits"],
                          capture_output=True, check=True)

            # Parse results
            if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
                df = pd.read_csv(output_file, sep="\t", header=None,
                               names=["query", "target", "pident", "alnlen", "mismatch", "gapopen",
                                     "qstart", "qend", "tstart", "tend", "evalue", "bitscore"])

                # Apply cutoff based on mode
                if use_evalue:
                    df = df[df['evalue'] < cutoff_value]
                else:
                    df = df[df['bitscore'] >= cutoff_value]

                # Keep best hit per query
                df = df.sort_values('evalue').groupby('query').first().reset_index()

                gene_set = results['toxins'] if target == 'toxins' else results['antitoxins']
                for _, row in df.iterrows():
                    gene_id = row['query']
                    gene_set.add(gene_id)
                    results['method_info'][gene_id] = {
                        'methods': ['Foldseek'],
                        'foldseek_evalue': row['evalue'],
                        'foldseek_bitscore': row['bitscore']
                    }
        except subprocess.CalledProcessError as e:
            print(f"Foldseek {target} search failed: {e}")
            continue

    # Handle dual identifications
    dual_identified = results['toxins'].intersection(results['antitoxins'])
    for gene_id in dual_identified:
        info = results['method_info'][gene_id]
        tox_score = info.get('foldseek_bitscore', 0)
        # For simplicity, remove from antitoxins if found in both
        results['antitoxins'].discard(gene_id)

    print(f"Foldseek found {len(results['toxins'])} toxins, {len(results['antitoxins'])} antitoxins")
    return results
