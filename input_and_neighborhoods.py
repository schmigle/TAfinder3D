from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import numpy as np
import os
import re
import subprocess
from utils import get_gene_id, is_ncrna, calculate_distance

def extract_cds_sequences(fna_file, coords_file, output_file):
    records = {rec.id: rec.seq for rec in SeqIO.parse(fna_file, "fasta")}

    with open(coords_file) as c, open(output_file, "w") as o:
        for line in c:
            parts = line.strip().split("\t")
            if len(parts) == 4:
                contig = list(records.keys())[0]
                start, stop, strand, name = parts
            elif len(parts) >= 5:
                contig, start, stop, strand, name = parts[:5]
            else:
                continue

            s, e = int(start), int(stop)
            seq = records[contig][s-1:e]
            if strand == "-":
                seq = seq.reverse_complement()

            o.write(f">{name} {contig}:{s}-{e}({strand})\n")
            o.write("\n".join(str(seq)[i:i+60] for i in range(0, len(seq), 60)) + "\n")

# File processing functions - built-in GenBank parser
def gb_to_files(gbk_path, fna_out, faa_out, ptt_out):
    """Parse GenBank file and extract sequences and annotations."""
    contig_records = []
    protein_records = []
    rows = []  # Location, Strand, Length, PID, Gene, Synonym, Code, COG, Product, Contig
    pid_counter = 0
    rna_counter = 0

    for record in SeqIO.parse(gbk_path, 'genbank'):
        contig_id = record.id
        contig_records.append(SeqRecord(record.seq, id=contig_id, description=""))

        for feat in record.features:
            if feat.type == "CDS":
                start = int(feat.location.start) + 1  # Convert to 1-based
                end = int(feat.location.end)
                strand = "+" if int(feat.location.strand or 1) == 1 else "-"

                pid = (
                    feat.qualifiers.get("protein_id", [None])[0]
                    or feat.qualifiers.get("locus_tag", [None])[0]
                    or f"prot_{pid_counter}"
                )
                pid_counter += 1

                locus_tag = feat.qualifiers.get("locus_tag", ["-"])[0]
                gene = feat.qualifiers.get("gene", ["-"])[0]
                product = feat.qualifiers.get("product", ["Predicted CDS"])[0]

                # Protein sequence
                prot_seq = feat.qualifiers.get("translation", [None])[0]
                if not prot_seq:
                    cds_seq = feat.extract(record.seq)
                    prot_seq = str(cds_seq.translate(to_stop=True))
                prot_seq = prot_seq.replace("*", "")

                protein_records.append(SeqRecord(Seq(prot_seq), id=pid, description=""))

                length = len(prot_seq)
                location = f"{start}..{end}"
                rows.append([location, strand, length, pid, gene, locus_tag, "-", "-", product, contig_id])

            else:
                # Identify ncRNAs
                is_rna_type = feat.type in {"tRNA", "rRNA", "ncRNA", "tmRNA", "misc_RNA", "sRNA", "antisense_RNA"}
                quals = " ".join(
                    sum([[str(v) for v in vs] for vs in feat.qualifiers.values()], [])
                ).lower()
                marked_nucl = "nucl" in quals

                if is_rna_type or marked_nucl:
                    start = int(feat.location.start) + 1  # Convert to 1-based
                    end = int(feat.location.end)
                    strand = "+" if int(feat.location.strand or 1) == 1 else "-"

                    base = (
                        feat.qualifiers.get("locus_tag", [None])[0]
                        or feat.qualifiers.get("gene", [None])[0]
                        or feat.type
                    )
                    rna_id = f"Nucl_{base if base else rna_counter}"
                    rna_counter += 1

                    product = feat.qualifiers.get("product", [feat.type])[0]
                    location = f"{start}..{end}"
                    length = end - start + 1
                    rows.append([location, strand, length, rna_id, "-", "-", "-", "-", product, contig_id])

    # Write outputs
    SeqIO.write(contig_records, fna_out, "fasta")
    SeqIO.write(protein_records, faa_out, "fasta")

    # Write PTT file with header
    total_length = sum(len(r.seq) for r in contig_records)
    first_contig = contig_records[0].id if contig_records else "genome"

    with open(ptt_out, 'w') as f:
        f.write(f"{first_contig} - 1..{total_length}\n")
        f.write(f"{len(rows)} genes\n")
        f.write("Location\tStrand\tLength\tPID\tGene\tSynonym\tCode\tCOG\tProduct\tContig\n")
        for row in rows:
            f.write("\t".join(str(x) for x in row) + "\n")

def _maybe_decompress(filepath, tmp_dir):
    """If filepath is gzipped, decompress to tmp_dir and return new path."""
    import gzip
    if filepath.endswith('.gz'):
        decompressed = os.path.join(tmp_dir, os.path.basename(filepath[:-3]))
        with gzip.open(filepath, 'rb') as f_in, open(decompressed, 'wb') as f_out:
            f_out.write(f_in.read())
        print(f"  Decompressed {os.path.basename(filepath)} -> {os.path.basename(decompressed)}")
        return decompressed
    return filepath


def process_input_files(args, config):
    args.infile = _maybe_decompress(args.infile, args.tmp_dir)

    if args.format == 'gbk':
        gb_to_files(args.infile, config['fna_file'], config['faa_file'], config['ptt_file'])
    else:
        subprocess.run(["cp", args.infile, config['fna_file']], check=True)

        # Prodigal processing
        prodigal_cmd = ["prodigal-gv"] if args.prodigal_mode == "gv" else ["prodigal"]
        if args.prodigal_mode == "meta":
            prodigal_cmd.extend(["-p", "meta"])

        prodigal_out = f"{args.tmp_dir}/prodigal.out"
        with open(args.infile) as infile, open(prodigal_out, 'w') as outfile:
            subprocess.run(prodigal_cmd + ["-f", "gff"], stdin=infile, stdout=outfile, check=True)

        parse_prodigal_into_ptt(prodigal_out, config['ptt_file'], config['fna_file'])
        transform_cds_file(config['ptt_file'])

        tsv_file = f"{config['ptt_file']}.coords"
        ffa_file = f"{args.tmp_dir}/sequence.ffa"

        extract_cds_sequences(config['fna_file'], tsv_file, ffa_file)

        # Translate nucleotide sequences to protein using Biopython
        seq_count = 0
        empty_count = 0
        with open(config['faa_file'], 'w') as f_out:
            for record in SeqIO.parse(ffa_file, "fasta"):
                try:
                    # Use cds=True for complete ORFs from Prodigal
                    protein_seq = str(record.seq.translate(cds=True)).rstrip('*')
                except Exception:
                    # Fallback for non-standard CDS (e.g., partial genes)
                    protein_seq = str(record.seq.translate(to_stop=True))
                if len(protein_seq) == 0:
                    empty_count += 1
                f_out.write(f">{record.id}\n{protein_seq}\n")
                seq_count += 1

        print(f"  Translated {seq_count} CDS to proteins ({empty_count} empty sequences)")

        # Debug: Check first few protein sequences
        if seq_count > 0:
            sample_seqs = []
            for i, record in enumerate(SeqIO.parse(config['faa_file'], "fasta")):
                if i < 3:
                    sample_seqs.append(f"{record.id}: {len(record.seq)} aa")
                else:
                    break
            print(f"  Sample proteins: {', '.join(sample_seqs)}")

def parse_prodigal_into_ptt(infile, outfile, fasta_file):
    total_length = sum(len(seq.seq) for seq in SeqIO.parse(fasta_file, 'fasta'))
    first_seq = next(SeqIO.parse(fasta_file, 'fasta'))
    title = first_seq.description if first_seq.description else first_seq.name

    with open(infile) as infile, open(outfile, "w") as outfile:
        i, content = 0, ""
        for line in infile:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) >= 7 and fields[2] == 'CDS':
                contig, start, end, strand = fields[0], int(fields[3]), int(fields[4]), fields[6]
                i += 1
                length = (end - start + 1) // 3 - 1
                content += f"{start}..{end}\t{strand}\t{length}\tORF{i}\t-\torf{i}\t-\t-\tPredicted ORF by PRODIGAL\t{contig}\n"

        outfile.write(f"{title} - 1..{total_length}\n{i} proteins\nLocation\tStrand\tLength\tPID\tGene\tSynonym\tCode\tCOG\tProduct\tContig\n{content}")

def transform_cds_file(infile):
    tempfile = infile + ".coords"
    with open(infile, 'r') as input_file, open(tempfile, 'w') as outfile:
        for line_no, line in enumerate(input_file, start=1):
            if line_no > 2 and re.match(r'[0-9]+\.\.[0-9]+', line):
                fields = line.split("\t")
                m = re.search(r'([0-9]+)..([0-9]+)', fields[0])
                if m:
                    # Include contig (fields[9]) for multi-contig genomes
                    contig = fields[9].strip() if len(fields) > 9 else ""
                    outfile.write(f"{contig}\t{m.group(1)}\t{m.group(2)}\t{fields[1]}\t{fields[3]}\n")

# Neighborhood building - moving window algorithm
def build_neighborhoods_with_rna_support(ptt_file, args, rna_results=None):
    ptt = pd.read_csv(ptt_file, sep="\t", skiprows=2)
    has_contig = "Contig" in ptt.columns

    # Add new PTT rows for RNA genes detected by MMseqs2
    if rna_results and rna_results.get('new_ptt_rows'):
        print(f"Adding {len(rna_results['new_ptt_rows'])} new RNA gene entries to PTT")
        new_rows_df = pd.DataFrame(rna_results['new_ptt_rows'])
        ptt = pd.concat([ptt, new_rows_df], ignore_index=True)

        # Sort by contig and location to maintain genomic order
        if has_contig:
            ptt = ptt.sort_values(['Contig', 'Location'], key=lambda x: x.map(lambda loc: int(loc.split('..')[0]) if '..' in str(loc) else 0))
        else:
            ptt = ptt.sort_values('Location', key=lambda x: x.map(lambda loc: int(loc.split('..')[0]) if '..' in str(loc) else 0))
        ptt = ptt.reset_index(drop=True)

        # Deduplicate overlapping TA RNA genes (same strand, same contig, ≥50% overlap)
        # Keep the one with highest bitscore
        # Only deduplicate same-strand overlaps (opposite strand = different antisense genes)
        print(f"Deduplicating overlapping TA RNA genes (same strand only)...")
        rna_genes_to_remove = set()
        # Include all TA RNA genes (Type I, III, and VIII)
        rna_genes = [pid for pid in ptt['PID'] if
                     pid.startswith('Nucl_type1_at_') or
                     pid.startswith('Nucl_type3_at_') or
                     pid.startswith('Nucl_type8_')]

        for i, rna1 in enumerate(rna_genes):
            if rna1 in rna_genes_to_remove:
                continue

            row1 = ptt[ptt['PID'] == rna1].iloc[0]
            loc1 = row1['Location']
            strand1 = row1['Strand']
            contig1 = row1['Contig'] if has_contig else None
            start1, end1 = map(int, loc1.split('..'))
            len1 = end1 - start1 + 1
            bitscore1 = rna_results['rna_scores'].get(rna1, (0, 0))[1]

            for rna2 in rna_genes[i+1:]:
                if rna2 in rna_genes_to_remove:
                    continue

                row2 = ptt[ptt['PID'] == rna2].iloc[0]
                loc2 = row2['Location']
                strand2 = row2['Strand']
                contig2 = row2['Contig'] if has_contig else None
                start2, end2 = map(int, loc2.split('..'))
                len2 = end2 - start2 + 1
                bitscore2 = rna_results['rna_scores'].get(rna2, (0, 0))[1]

                # Check if same strand and contig
                if strand1 != strand2:
                    continue
                if has_contig and contig1 != contig2:
                    continue

                # Calculate overlap
                overlap_start = max(start1, start2)
                overlap_end = min(end1, end2)
                overlap_len = max(0, overlap_end - overlap_start + 1)

                # Check if overlap ≥50% of either gene
                overlap_pct1 = overlap_len / len1 if len1 > 0 else 0
                overlap_pct2 = overlap_len / len2 if len2 > 0 else 0

                if overlap_pct1 >= 0.5 or overlap_pct2 >= 0.5:
                    # Mark the one with lower bitscore for removal
                    if bitscore1 >= bitscore2:
                        rna_genes_to_remove.add(rna2)
                    else:
                        rna_genes_to_remove.add(rna1)
                        break  # rna1 is marked, no need to check further

        if rna_genes_to_remove:
            print(f"  Removing {len(rna_genes_to_remove)} duplicate/overlapping RNA genes")
            ptt = ptt[~ptt['PID'].isin(rna_genes_to_remove)].reset_index(drop=True)

            # Also remove from rna_results tracking dicts
            for rna_gene in rna_genes_to_remove:
                rna_results['all_rna_genes'].discard(rna_gene)
                rna_results['rna_toxin_locations'].pop(rna_gene, None)
                rna_results['rna_antitoxin_locations'].pop(rna_gene, None)
                rna_results['type8_genes'].discard(rna_gene)
                rna_results['type3_genes'].discard(rna_gene)
                rna_results['rna_scores'].pop(rna_gene, None)

    rna_gene_mapping = {}  # Maps original gene_id -> gene_id (for compatibility)

    gene_to_neighborhood = {}
    neighborhood_to_genes = {}
    proteins_in_neighborhoods = set()
    neighborhood_counter = 0

    # Load genome sequence for RNA checks (kept for compatibility)
    genome_seq = ""
    try:
        fna_file = f"{args.tmp_dir}/sequence.fna"
        if os.path.exists(fna_file):
            genome_seq = str(next(SeqIO.parse(fna_file, 'fasta')).seq)
    except StopIteration:
        pass

    contigs = ptt['Contig'].unique() if has_contig else [None]

    # Pre-parse all locations once for speed (avoid repeated string parsing in loop)
    locations = ptt['Location'].values
    starts = np.zeros(len(ptt), dtype=np.int64)
    ends = np.zeros(len(ptt), dtype=np.int64)
    for i, loc in enumerate(locations):
        parts = str(loc).split('..')
        if len(parts) == 2:
            starts[i] = int(parts[0])
            ends[i] = int(parts[1])

    pids = ptt['PID'].values
    neighborhood_dist = args.neighborhood_distance
    min_size = 1 if args.search_all else args.min_neighborhood_size
    max_size = len(ptt) + 1 if args.search_all else args.max_neighborhood_size

    for contig in contigs:
        if has_contig:
            contig_mask = (ptt['Contig'] == contig).values
            contig_indices = np.where(contig_mask)[0]
        else:
            contig_indices = np.arange(len(ptt))

        if len(contig_indices) == 0:
            continue

        # Vectorized gap calculation for this contig
        contig_starts = starts[contig_indices]
        contig_ends = ends[contig_indices]
        contig_pids = pids[contig_indices]

        # Gaps between adjacent genes: start[i+1] - end[i]
        if len(contig_indices) > 1:
            gaps = contig_starts[1:] - contig_ends[:-1]
            # Find where gaps exceed threshold (these are segment boundaries)
            break_points = np.where(gaps > neighborhood_dist)[0] + 1
            # Build list of segment boundaries
            all_breaks = np.concatenate(([0], break_points, [len(contig_indices)]))
        else:
            all_breaks = np.array([0, 1])

        # Process each contiguous segment
        for seg_idx in range(len(all_breaks) - 1):
            seg_start = all_breaks[seg_idx]
            seg_end = all_breaks[seg_idx + 1]
            seg_len = seg_end - seg_start

            # Check size criteria
            if seg_len < min_size or seg_len > max_size:
                continue

            # Create neighborhood
            neighborhood_id = f"neighborhood_{neighborhood_counter}"
            neighborhood_counter += 1

            gene_ids = contig_pids[seg_start:seg_end].tolist()
            neighborhood_to_genes[neighborhood_id] = gene_ids

            for gene_id in gene_ids:
                gene_to_neighborhood[gene_id] = neighborhood_id
                if not is_ncrna(gene_id):
                    proteins_in_neighborhoods.add(gene_id)

    print(f"Built {len(neighborhood_to_genes)} neighborhoods with {len(proteins_in_neighborhoods)} proteins")
    return {
        'gene_to_neighborhood': gene_to_neighborhood,
        'neighborhood_to_genes': neighborhood_to_genes,
        'proteins_in_neighborhoods': proteins_in_neighborhoods,
        'ptt': ptt,
        'genome_seq': genome_seq,
        'rna_gene_mapping': rna_gene_mapping  # Track RNA genes that were marked
    }
