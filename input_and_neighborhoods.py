from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import os
import re
import subprocess
from config import MAX_NEIGHBORHOOD_GENES
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

def process_input_files(args, config):
    if args.format == 'gbk':
        gb_to_files(args.input, config['fna_file'], config['faa_file'], config['ptt_file'])
    else:
        subprocess.run(["cp", args.input, config['fna_file']], check=True)

        # Prodigal processing
        prodigal_cmd = ["prodigal-gv"] if args.prodigal_mode == "gv" else ["prodigal"]
        if args.prodigal_mode == "meta":
            prodigal_cmd.extend(["-p", "meta"])

        prodigal_out = f"{args.tmp_dir}/prodigal.out"
        with open(args.input) as infile, open(prodigal_out, 'w') as outfile:
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

    # Mark genes as ncRNA if they're in the validated RNA set
    rna_gene_mapping = {}  # Maps original gene_id -> marked gene_id
    if rna_results and rna_results['validated_rna_genes']:
        print(f"Marking {len(rna_results['validated_rna_genes'])} MMseqs2-detected RNA genes as ncRNA")
        for idx in range(len(ptt)):
            gene_id = ptt.loc[idx, "PID"]
            if gene_id in rna_results['validated_rna_genes']:
                # Mark as ncRNA by prepending 'Nucl' to PID if not already there
                if not gene_id.startswith('Nucl'):
                    marked_id = f"Nucl_{gene_id}"
                    ptt.loc[idx, "PID"] = marked_id
                    rna_gene_mapping[gene_id] = marked_id
                else:
                    rna_gene_mapping[gene_id] = gene_id

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

    for contig in contigs:
        contig_indices = ptt[ptt['Contig'] == contig].index.tolist() if has_contig else list(range(len(ptt)))

        # Moving window algorithm
        i = 0
        while i < len(contig_indices):
            window = []
            window_start_idx = i

            # Add first gene to window
            idx = contig_indices[i]
            window.append(idx)
            prev_loc = ptt.loc[idx, "Location"]
            i += 1

            # Keep adding genes until we hit a gap >200bp or >3 genes (unless search_all)
            max_window = float('inf') if args.search_all else MAX_NEIGHBORHOOD_GENES
            while i < len(contig_indices) and len(window) < max_window:
                idx = contig_indices[i]
                curr_loc = ptt.loc[idx, "Location"]

                # Calculate gap from previous gene
                distance = calculate_distance(prev_loc, curr_loc)

                if distance > args.neighborhood_distance:
                    # Gap too large, stop window
                    break

                # Add to window
                window.append(idx)
                prev_loc = curr_loc
                i += 1

            # If we have 2-3 genes in window, it's a neighborhood
            # (or 1 gene if search_all is enabled)
            if args.search_all or (2 <= len(window) <= MAX_NEIGHBORHOOD_GENES):
                neighborhood_id = f"neighborhood_{neighborhood_counter}"
                neighborhood_counter += 1

                gene_ids = [get_gene_id(ptt, idx, False) for idx in window]
                neighborhood_to_genes[neighborhood_id] = gene_ids.copy()

                for gene_id in gene_ids:
                    gene_to_neighborhood[gene_id] = neighborhood_id
                    if not is_ncrna(gene_id):
                        proteins_in_neighborhoods.add(gene_id)

            # If window only had 1 gene, move to next gene
            # Otherwise we've already moved past the window
            if len(window) == 1:
                # Already incremented i, so continue
                pass

    print(f"Built {len(neighborhood_to_genes)} neighborhoods with {len(proteins_in_neighborhoods)} proteins")
    return {
        'gene_to_neighborhood': gene_to_neighborhood,
        'neighborhood_to_genes': neighborhood_to_genes,
        'proteins_in_neighborhoods': proteins_in_neighborhoods,
        'ptt': ptt,
        'genome_seq': genome_seq,
        'rna_gene_mapping': rna_gene_mapping  # Track RNA genes that were marked
    }
