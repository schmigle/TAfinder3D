from argparse import ArgumentParser, RawDescriptionHelpFormatter
import os

# Configuration constants
MAX_NEIGHBORHOOD_GENES = 5
MAX_RNA_DISTANCE = 150
MAX_NCRNA_DISTANCE = 200

EPILOG = """
Examples:
  # Basic usage with GenBank file
  tafinder -i genome.gbk -f gbk -o results.csv

  # FASTA input with custom database directory
  tafinder -i genome.fna -f fna -o results.csv --db_dir /path/to/db

  # Predict orphans and save neighborhoods
  tafinder -i genome.gbk -o results.csv --predict_orphan --save_unidentified

  # Custom scoring thresholds
  tafinder -i genome.gbk -o results.csv --mmseqs_bitscore 60 --foldseek_bitscore 300

  # Use e-value cutoffs instead of bitscores
  tafinder -i genome.gbk -o results.csv --use_evalue -e 0.001

For more information, visit: https://github.com/schmigle/TAfinder3D
"""

def parse_args():
    parser = ArgumentParser(
        description="TAfinder3D: Toxin-Antitoxin system identification with 3D structure-based searching",
        epilog=EPILOG,
        formatter_class=RawDescriptionHelpFormatter
    )
    parser.add_argument("-i", "--input", type=str, required=True, help="Input genome file")
    parser.add_argument("-f", "--format", type=str, default='gbk', choices=['gbk', 'fna'], help="Input format")
    parser.add_argument("-o", "--output", type=str, default=os.path.join(os.getcwd(), 'result.csv'), help='Output result file')
    parser.add_argument("--use_evalue", action="store_true", help='Use e-value cutoffs instead of bitscore cutoffs')
    parser.add_argument("-e", "--evalue", type=float, default=0.01, help='E-value threshold for all searches (when --use_evalue is set)')
    parser.add_argument("--hmmevalue", type=float, default=0.001, help='HMM E-value threshold (when --use_evalue is set)')
    parser.add_argument("--foldseekevalue", type=float, default=0.001, help='Foldseek E-value threshold (when --use_evalue is set)')
    parser.add_argument("--hmm_bitscore", type=float, default=25, help='HMM bitscore threshold (default: 25)')
    parser.add_argument("--mmseqs_bitscore", type=float, default=50, help='MMseqs2 bitscore threshold (default: 50)')
    parser.add_argument("--foldseek_bitscore", type=float, default=250, help='Foldseek bitscore threshold (default: 250)')
    parser.add_argument("--prodigal_mode", type=str, default="normal", choices=["normal", "meta", "gv"], help='Prodigal mode')
    parser.add_argument("--foldseeksensitivity", type=float, default=9.5, help='Foldseek sensitivity')
    parser.add_argument("--toxin_db", type=str, default='/db/foldseek_tox_db_padded', help='Path to custom toxin FASTA (used for MMseqs2, HMM, and Foldseek databases)')
    parser.add_argument("--antitoxin_db", type=str, default='/db/foldseek_antitox_db_padded', help='Path to custom antitoxin FASTA (used for MMseqs2, HMM, and Foldseek databases)')
    parser.add_argument("--rna_toxin_db", type=str, help='Path to RNA toxin nucleotide database for MMseqs2')
    parser.add_argument("--rna_antitoxin_db", type=str, help='Path to RNA antitoxin nucleotide database for MMseqs2')
    parser.add_argument("-a", "--havalue", type=float, default=0.36, help='Ha-value threshold')
    parser.add_argument("-n", "--minimum_distance", type=int, default=-30, help='Minimum intergenic distances (bp)')
    parser.add_argument("--neighborhood_distance", type=int, default=200, help='Maximum distance between genes in neighborhood (bp)')
    parser.add_argument("-l", "--maximum_length", type=int, default=800, help='Maximum sequence length (a.a.)')
    parser.add_argument("-s", "--extract_seq", action="store_true", help='Extract predicted sequences')
    parser.add_argument("-p", "--predict_orphan", action="store_true", help='Predict orphan toxins and antitoxins')
    parser.add_argument("--save_unidentified", type=str, nargs='?', const='ta_neighborhoods.fna', help='Save unidentified neighborhoods')
    parser.add_argument("-t", "--tmp_dir", type=str, default=os.path.join(os.getcwd(), "tafinder_tmp"), help='Temp directory')
    parser.add_argument("--threads", type=int, default=os.cpu_count(), help='Number of threads')
    parser.add_argument("--no_foldseek", action="store_true", help='Disable Foldseek')
    parser.add_argument("--no_hmm", action="store_true", help='Disable HMM')
    parser.add_argument("--no_mmseqs", action="store_true", help='Disable MMseqs2')
    parser.add_argument("--search_all", action="store_true", help='Search all genes, not just neighborhoods')
    parser.add_argument("--db_dir", type=str, default=os.path.join(os.getcwd(), "tafinder_db"), help='Path to tafinder database directory')
    return parser.parse_args()

def initialize_environment(args):
    script_dir = os.path.dirname(os.path.realpath(__file__))
    os.makedirs(args.tmp_dir, exist_ok=True)

    # Use db_dir from args (defaults to tafinder_db in current directory)
    db_dir = args.db_dir

    # Setup database paths
    return {
        'script_dir': script_dir,
        'fna_file': f"{args.tmp_dir}/sequence.fna",
        'ptt_file': f"{args.tmp_dir}/sequence.ptt",
        'faa_file': f"{args.tmp_dir}/sequence.faa",
        'prostt5_weights': os.path.join(db_dir, "prostt5_weights"),

        # Protein databases (padded for GPU acceleration)
        'mmseqs_tox_db': os.path.join(db_dir, "mmseqs_tox_db_padded"),
        'mmseqs_antitox_db': os.path.join(db_dir, "mmseqs_antitox_db_padded"),
        'hmm_tox_file': os.path.join(db_dir, "tox_domains.hmm"),
        'hmm_antitox_file': os.path.join(db_dir, "antitox_domains.hmm"),
        'foldseek_tox_db': os.path.join(db_dir, "foldseek_tox_db_padded"),
        'foldseek_antitox_db': os.path.join(db_dir, "foldseek_antitox_db_padded"),

        # Type VIII RNA databases (padded for GPU acceleration)
        'mmseqs_type8_tox_db': os.path.join(db_dir, "mmseqs_8T_padded"),
        'mmseqs_type8_antitox_db': os.path.join(db_dir, "mmseqs_8AT_padded"),

        # Type I/III RNA antitoxin database (padded for GPU acceleration)
        'mmseqs_rna_antitox_db': os.path.join(db_dir, "mmseqs_rna_antitox_db_padded"),

        # DNA databases for MMseqs2 (padded for GPU acceleration)
        'mmseqs_dna_tox_db': os.path.join(db_dir, "mmseqs_dna_tox_db_padded"),
        'mmseqs_dna_antitox_db': os.path.join(db_dir, "mmseqs_dna_antitox_db_padded"),

        'use_custom_db': False
    }
