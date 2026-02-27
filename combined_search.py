from Bio import SeqIO
import pandas as pd
import os
import subprocess
import shutil
from config import MAX_NCRNA_DISTANCE
from utils import (determine_cutoff_params, parse_search_output,
                   find_overlapping_gene_in_ptt, update_method_info, is_ncrna)

def build_search_flags(args, evalue):
    """Build common search flags (without GPU)."""
    return [
        "--format-output", "query,target,evalue,bits",
        "--threads", str(args.threads),
        "-e", str(evalue),
        "--max-seqs", "10"
    ]

def run_single_rna_mmseqs_search(fna_file, db_path, output_file, args, ptt, has_contig, cutoff_value, use_evalue, evalue_param):
    locations = {}
    all_genes = set()

    cmd = ["mmseqs", "easy-search", fna_file, db_path, output_file, args.tmp_dir,
           "--format-output", "query,target,evalue,bits",
           "--threads", str(args.threads), "--search-type", "3",  # 3 = nucleotide search
           "-e", evalue_param]

    try:
        subprocess.run(cmd, capture_output=True, text=True, check=True)
        df = parse_search_output(output_file, 'mmseqs_rna', cutoff_value, use_evalue)

        if not df.empty:
            for _, row in df.iterrows():
                query_name = str(row[0])  # Query is contig/sequence name
                query_start, query_end = int(row[7]), int(row[8])
                gene_id, start, end, strand, contig = find_overlapping_gene_in_ptt(ptt, query_start, query_end, has_contig, query_name)
                if gene_id:
                    locations[gene_id] = (start, end, strand, contig)
                    all_genes.add(gene_id)

    except subprocess.CalledProcessError as e:
        print(f"  Single RNA MMseqs2 search failed")
        if e.stderr:
            print(f"  Error: {e.stderr}")
        if e.stdout:
            print(f"  Output: {e.stdout}")

    return locations, all_genes

def run_rna_mmseqs_searches(fna_file, ptt_file, config, args):
    print("Running nucleotide MMseqs2 searches...")

    rna_results = {
        'rna_toxin_locations': {},      # gene_id -> (start, end, strand, contig)
        'rna_antitoxin_locations': {},  # gene_id -> (start, end, strand, contig)
        'validated_rna_genes': set(),   # Only RNA toxins with paired antitoxins
        'all_rna_genes': set(),         # All detected RNAs (for tracking)
        'type8_genes': set(),           # Type VIII RNA genes
        'type3_genes': set(),           # Type III RNA genes (neighbor rule)
        'new_ptt_rows': [],             # New PTT entries to add for RNA genes
        'rna_scores': {}                # gene_id -> (evalue, bitscore)
    }

    # Read PTT to get contig info
    ptt = pd.read_csv(ptt_file, sep="\t", skiprows=2)
    has_contig = "Contig" in ptt.columns

    # Determine cutoff value and mode
    cutoff_value, use_evalue, evalue_param = determine_cutoff_params(args, 'mmseqs')

    print("  Using unified nucleotide database...")
    nucleotide_output = f"{args.tmp_dir}/mmseqs_nucleotide.out"

    # RNA search needs coordinates for PTT mapping
    cmd = ["mmseqs", "easy-search", fna_file, config['mmseqs_nucleotide_db'], nucleotide_output, args.tmp_dir,
           "--search-type", "3", "--format-output", "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,evalue,bits",
           "--threads", str(args.threads), "-e", evalue_param]

    try:
        result = subprocess.run(cmd, check=True)

        df = parse_search_output(nucleotide_output, 'mmseqs_rna', cutoff_value, use_evalue)

        if not df.empty:
            type8_tox_locations = {}
            type8_at_locations = {}
            type1_counter = 0
            type3_counter = 0
            type8_counter = 0

            for _, row in df.iterrows():
                query_name = str(row[0])  # Contig/sequence name
                target_id = str(row[1])   # Target has prefix
                query_start, query_end = int(row[6]), int(row[7])  # qstart, qend now at columns 6, 7
                evalue = float(row[8])  # E-value at column 8
                bitscore = float(row[9])  # Bitscore at column 9

                # Determine strand from coordinates
                if query_start < query_end:
                    strand = '+'
                    start, end = query_start, query_end
                else:
                    strand = '-'
                    start, end = query_end, query_start

                # Contig is the query name
                contig = query_name
                length = (end - start + 1) // 3  # Approximate length in AA for consistency

                # Parse prefix to determine type and create new PTT entry
                if target_id.startswith('toxin8_'):
                    # Type VIII RNA toxin
                    gene_id = f"Nucl_type8_tox_{type8_counter}"
                    type8_counter += 1
                    type8_tox_locations[gene_id] = (start, end, strand, contig)
                    rna_results['all_rna_genes'].add(gene_id)
                    rna_results['rna_scores'][gene_id] = (evalue, bitscore)

                    # Create PTT row
                    ptt_row = {
                        'Location': f"{start}..{end}",
                        'Strand': strand,
                        'Length': length,
                        'PID': gene_id,
                        'Gene': '-',
                        'Synonym': '-',
                        'Code': '-',
                        'COG': '-',
                        'Product': f"Type VIII RNA toxin ({target_id})",
                        'Contig': contig if has_contig else ptt.iloc[0]['Contig']
                    }
                    rna_results['new_ptt_rows'].append(ptt_row)

                elif target_id.startswith('antitoxin8_'):
                    # Type VIII RNA antitoxin
                    gene_id = f"Nucl_type8_at_{type8_counter}"
                    type8_counter += 1
                    type8_at_locations[gene_id] = (start, end, strand, contig)
                    rna_results['all_rna_genes'].add(gene_id)
                    rna_results['rna_scores'][gene_id] = (evalue, bitscore)

                    ptt_row = {
                        'Location': f"{start}..{end}",
                        'Strand': strand,
                        'Length': length,
                        'PID': gene_id,
                        'Gene': '-',
                        'Synonym': '-',
                        'Code': '-',
                        'COG': '-',
                        'Product': f"Type VIII RNA antitoxin ({target_id})",
                        'Contig': contig if has_contig else ptt.iloc[0]['Contig']
                    }
                    rna_results['new_ptt_rows'].append(ptt_row)

                elif target_id.startswith('antitoxin1_'):
                    # Type I RNA antitoxin
                    gene_id = f"Nucl_type1_at_{type1_counter}"
                    type1_counter += 1
                    rna_results['rna_antitoxin_locations'][gene_id] = (start, end, strand, contig)
                    rna_results['all_rna_genes'].add(gene_id)
                    rna_results['rna_scores'][gene_id] = (evalue, bitscore)

                    ptt_row = {
                        'Location': f"{start}..{end}",
                        'Strand': strand,
                        'Length': length,
                        'PID': gene_id,
                        'Gene': '-',
                        'Synonym': '-',
                        'Code': '-',
                        'COG': '-',
                        'Product': f"Type I RNA antitoxin ({target_id})",
                        'Contig': contig if has_contig else ptt.iloc[0]['Contig']
                    }
                    rna_results['new_ptt_rows'].append(ptt_row)

                elif target_id.startswith('antitoxin3_'):
                    # Type III RNA antitoxin
                    gene_id = f"Nucl_type3_at_{type3_counter}"
                    type3_counter += 1
                    rna_results['rna_antitoxin_locations'][gene_id] = (start, end, strand, contig)
                    rna_results['all_rna_genes'].add(gene_id)
                    rna_results['type3_genes'].add(gene_id)
                    rna_results['rna_scores'][gene_id] = (evalue, bitscore)

                    ptt_row = {
                        'Location': f"{start}..{end}",
                        'Strand': strand,
                        'Length': length,
                        'PID': gene_id,
                        'Gene': '-',
                        'Synonym': '-',
                        'Code': '-',
                        'COG': '-',
                        'Product': f"Type III RNA antitoxin ({target_id})",
                        'Contig': contig if has_contig else ptt.iloc[0]['Contig']
                    }
                    rna_results['new_ptt_rows'].append(ptt_row)

                # Skip other prefixes (toxin_, antitoxin_) - handled by DNA search

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

            # Merge Type VIII locations into main results
            rna_results['rna_toxin_locations'].update(type8_tox_locations)
            rna_results['rna_antitoxin_locations'].update(type8_at_locations)

            print(f"  Found {type1_counter} Type I, {type3_counter} Type III, {len(type8_tox_locations)} Type VIII toxin, {len(type8_at_locations)} Type VIII antitoxin RNA hits")

    except subprocess.CalledProcessError as e:
        print("  Warning: Nucleotide database search failed")
        if e.stderr:
            print(f"  Error: {e.stderr}")
        if e.stdout:
            print(f"  Output: {e.stdout}")
    return rna_results

def run_mmseqs_searches(search_file, config, args):
    results = {'toxins': set(), 'antitoxins': set(), 'method_info': {}}

    if args.no_mmseqs:
        return results

    seq_count = sum(1 for line in open(search_file) if line.startswith('>'))
    file_size = os.path.getsize(search_file)
    print(f"  Search file: {search_file} ({seq_count} sequences, {file_size} bytes)")

    output_file = f"{args.tmp_dir}/mmseqs_unified.out"

    common_flags = build_search_flags(args, args.mmseqs_evalue)
    gpu_flag = ["--gpu", "1"] if config['gpu_available'] else []
    sensitivity_flag = ["-s", str(args.mmseqs_sensitivity)] if args.mmseqs_sensitivity else []
    cmd = ["mmseqs", "easy-search", search_file, config['mmseqs_protein_db'], output_file, args.tmp_dir] + common_flags + gpu_flag + sensitivity_flag

    try:
        subprocess.run(cmd, check=True)

        df = pd.read_csv(output_file, sep="\t", header=None,
                        names=["query", "target", "evalue", "bitscore"])

        # Filter by bitscore threshold
        df = df[df['bitscore'] >= args.mmseqs_bitscore]

        if not df.empty:
            # Keep best bitscore per query
            df = df.sort_values('bitscore', ascending=False).groupby('query').first().reset_index()

            for _, row in df.iterrows():
                query_id = row['query']
                target_id = row['target']

                # Determine if hit is toxin or antitoxin based on prefix
                if target_id.startswith('toxin_'):
                    results['toxins'].add(query_id)
                    update_method_info(results, query_id, 'MMseqs2', row['evalue'], row['bitscore'], False, target=target_id)
                elif target_id.startswith('antitoxin_'):
                    results['antitoxins'].add(query_id)
                    update_method_info(results, query_id, 'MMseqs2', row['evalue'], row['bitscore'], False, target=target_id)

    except subprocess.CalledProcessError as e:
        print(f"  MMseqs2 search failed")
        if e.stderr:
            print(f"  Error: {e.stderr}")
        if e.stdout:
            print(f"  Output: {e.stdout}")

    print(f"  Protein MMseqs2 found {len(results['toxins'])} toxins, {len(results['antitoxins'])} antitoxins")
    return results

def run_dna_mmseqs_searches(cds_nuc_file, config, args):
    print("Running MMseqs2 DNA searches...")

    results = {'toxins': set(), 'antitoxins': set(), 'method_info': {}}

    output_file = f"{args.tmp_dir}/mmseqs_dna.out"

    common_flags = build_search_flags(args, args.mmseqs_evalue)
    sensitivity_flag = ["-s", str(args.mmseqs_sensitivity)] if args.mmseqs_sensitivity else []

    cmd = ["mmseqs", "easy-search", cds_nuc_file, config['mmseqs_nucleotide_db'], output_file, args.tmp_dir,
           "--search-type", "3"] + common_flags + sensitivity_flag

    try:
        subprocess.run(cmd, check=True)

        df = pd.read_csv(output_file, sep="\t", header=None,
                        names=["query", "target", "evalue", "bitscore"])

        # Filter by bitscore threshold
        df = df[df['bitscore'] >= args.mmseqs_bitscore]

        if not df.empty:
            # Keep best bitscore per query
            df = df.sort_values('bitscore', ascending=False).groupby('query').first().reset_index()

            for _, row in df.iterrows():
                query_id = row['query']
                target_id = row['target']

                # Determine if hit is toxin or antitoxin based on prefix
                # DNA toxins use toxin_ prefix, DNA antitoxins use antitoxin_ prefix
                if target_id.startswith('toxin_'):
                    results['toxins'].add(query_id)
                    update_method_info(results, query_id, 'MMseqs2_DNA', row['evalue'], row['bitscore'], False, target=target_id)
                elif target_id.startswith('antitoxin_'):
                    results['antitoxins'].add(query_id)
                    update_method_info(results, query_id, 'MMseqs2_DNA', row['evalue'], row['bitscore'], False, target=target_id)

    except subprocess.CalledProcessError as e:
        print(f"  DNA search failed")
        if e.stderr:
            print(f"  Error: {e.stderr}")
        if e.stdout:
            print(f"  Output: {e.stdout}")

    print(f"  DNA search found {len(results['toxins'])} toxins, {len(results['antitoxins'])} antitoxins")
    return results

def run_hmm_searches(search_file, config, args):
    results = {'toxins': set(), 'antitoxins': set(), 'method_info': {}}

    if args.no_hmm:
        print("  HMM searches disabled (--no_hmm)")
        return results

    print("  Running HMM searches...")

    # Use unified HMM database with e-value cutoff (bitscore filtering happens post-search)
    output_file = f"{args.tmp_dir}/hmmsearch_unified.domtbl"
    cmd = ["hmmsearch", "--cpu", str(args.threads), "--noali", "-E", str(args.hmm_evalue),
           "-o", f"{args.tmp_dir}/hmmsearch.txt",
           "--domtblout", output_file, config['hmm_combined_file'], search_file]

    try:
        subprocess.run(cmd, check=True)

        # Parse domtblout format
        # Columns: target_name, accession, tlen, query_name, accession, qlen, evalue, score, bias, #, of, c-Evalue, i-Evalue, score, bias, from, to, from, to, from, to, acc, description
        hits = []
        with open(output_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.split()
                if len(parts) < 23:
                    continue
                query_id = parts[0]  # Query protein
                target_id = parts[3]  # HMM profile name
                score = float(parts[7])  # Full sequence score
                evalue = float(parts[6])  # Full sequence e-value

                hits.append({
                    'query': query_id,
                    'target': target_id,
                    'score': score,
                    'evalue': evalue
                })

        for hit in hits:
            query_id = hit['query']
            target_id = hit['target']

            # Apply bitscore filter
            if hit['score'] < args.hmm_bitscore:
                continue

            # Determine if hit is toxin or antitoxin based on prefix
            if target_id.startswith('toxin_'):
                results['toxins'].add(query_id)
                update_method_info(results, query_id, 'HMM', hit['evalue'], hit['score'], False, target=target_id)
            elif target_id.startswith('antitoxin_'):
                results['antitoxins'].add(query_id)
                update_method_info(results, query_id, 'HMM', hit['evalue'], hit['score'], False, target=target_id)

    except subprocess.CalledProcessError as e:
        print(f"  HMM search failed")
        if e.stderr:
            print(f"  Error: {e.stderr}")

    print(f"  HMM found {len(results['toxins'])} toxins, {len(results['antitoxins'])} antitoxins")
    return results

def run_foldseek_searches(candidate_proteins, config, args):
    """Run Foldseek searches using easy-search."""
    import time
    start_time = time.time()

    results = {'toxins': set(), 'antitoxins': set(), 'method_info': {}}

    if args.no_foldseek or not candidate_proteins:
        return results

    foldseek_path = shutil.which('foldseek')
    if not foldseek_path:
        print("Foldseek not found in PATH")
        return results

    # Create candidate FASTA file (always needed for easy-search)
    candidate_faa = f"{args.tmp_dir}/candidates.faa"
    protein_seqs = {r.id: str(r.seq) for r in SeqIO.parse(config['faa_file'], "fasta") if r.id in candidate_proteins}

    if not protein_seqs:
        return results

    with open(candidate_faa, 'w') as f:
        for prot_id, seq in protein_seqs.items():
            f.write(f">{prot_id}\n{seq}\n")

    output_file = f"{args.tmp_dir}/foldseek_unified.tsv"

    common_flags = build_search_flags(args, args.foldseek_evalue)
    gpu_flag = ["--gpu", "1"] if config['gpu_available'] else []
    sensitivity_flag = ["-s", str(args.foldseek_sensitivity)] if args.foldseek_sensitivity else []

    # Check if pre-computed query database is provided
    if args.foldseek_query_db:
        print(f"  Using pre-built Foldseek query database: {args.foldseek_query_db}")
        search_cmd = [foldseek_path, "easy-search", args.foldseek_query_db, config['foldseek_db'], output_file, args.tmp_dir] + common_flags + gpu_flag + sensitivity_flag
    else:
        create_cmd = [foldseek_path, "createdb", candidate_faa, str(args.tmp_dir + "/foldseek_query_db"), "--prostt5-model", config['prostt5_weights']] + gpu_flag
        search_cmd = [foldseek_path, "easy-search", str(args.tmp_dir + "/foldseek_query_db"), config['foldseek_db'], output_file, args.tmp_dir,
               "--prostt5-model", config['prostt5_weights']] + common_flags + gpu_flag + sensitivity_flag

    try:
        search_start = time.time()
        if not args.foldseek_query_db:
            subprocess.run(create_cmd, check=True)
        subprocess.run(search_cmd, check=True)
        print(f"  Foldseek completed in {time.time() - search_start:.2f}s")

        # Parse results and track bitscores for dual hits
        if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
            df = pd.read_csv(output_file, sep="\t", header=None,
                           names=["query", "target", "evalue", "bitscore"])

            # Filter by bitscore threshold
            df = df[df['bitscore'] >= args.foldseek_bitscore]

            toxin_scores = {}  
            antitoxin_scores = {}  

            for _, row in df.iterrows():
                query_id = row['query']
                target_id = row['target']
                evalue = row['evalue']
                bitscore = row['bitscore']

                # Track best bitscore for toxins and antitoxins separately
                if target_id.startswith('toxin_'):
                    if query_id not in toxin_scores or bitscore > toxin_scores[query_id][2]:
                        toxin_scores[query_id] = (target_id, evalue, bitscore)
                elif target_id.startswith('antitoxin_'):
                    if query_id not in antitoxin_scores or bitscore > antitoxin_scores[query_id][2]:
                        antitoxin_scores[query_id] = (target_id, evalue, bitscore)

            # Handle dual identifications: keep the one with higher bitscore
            for gene_id in set(toxin_scores.keys()).intersection(antitoxin_scores.keys()):
                tox_score = toxin_scores[gene_id][2]
                at_score = antitoxin_scores[gene_id][2]

                if tox_score >= at_score:
                    del antitoxin_scores[gene_id]
                else:
                    del toxin_scores[gene_id]

            # Populate results
            for gene_id, (target_id, evalue, bitscore) in toxin_scores.items():
                results['toxins'].add(gene_id)
                results['method_info'][gene_id] = {
                    'methods': ['Foldseek'],
                    'foldseek_evalue': evalue,
                    'foldseek_bitscore': bitscore,
                    'foldseek_target': target_id
                }

            for gene_id, (target_id, evalue, bitscore) in antitoxin_scores.items():
                results['antitoxins'].add(gene_id)
                results['method_info'][gene_id] = {
                    'methods': ['Foldseek'],
                    'foldseek_evalue': evalue,
                    'foldseek_bitscore': bitscore,
                    'foldseek_target': target_id
                }

    except subprocess.CalledProcessError as e:
        print(f"Foldseek easy-search failed")
        if e.stderr:
            print(f"  Error: {e.stderr}")
        if e.stdout:
            print(f"  Output: {e.stdout}")

    print(f"Foldseek found {len(results['toxins'])} toxins, {len(results['antitoxins'])} antitoxins")
    return results
