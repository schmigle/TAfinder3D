from argparse import ArgumentParser, RawDescriptionHelpFormatter
import os

# Configuration constants
MAX_RNA_DISTANCE = 150  # For Type I RNA-mRNA pairing
MAX_NCRNA_DISTANCE = 200  # For Type VIII RNA-RNA pairing

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
    parser.add_argument("-i", "--infile", type=str, required=True, help="Input genome file")
    parser.add_argument("-f", "--format", type=str, default='gbk', choices=['gbk', 'fna'], help="Input format")
    parser.add_argument("-o", "--outfile", type=str, default=os.path.join(os.getcwd(), 'result.csv'), help='Output result file')
    parser.add_argument("--mmseqs_evalue", type=float, default=1e-5, help='E-value threshold for MMseqs2 searches (default: 1e-5)')
    parser.add_argument("--hmm_evalue", type=float, default=1e-5, help='HMM E-value threshold (default: 1e-5)')
    parser.add_argument("--foldseek_evalue", type=float, default=1e-5, help='Foldseek E-value threshold (default: 1e-5)')
    parser.add_argument("--hmm_bitscore", type=float, default=25, help='HMM bitscore threshold (default: 25)')
    parser.add_argument("--mmseqs_bitscore", type=float, default=40, help='MMseqs2 bitscore threshold (default: 40)')
    parser.add_argument("--foldseek_bitscore", type=float, default=250, help='Foldseek bitscore threshold (default: 250)')
    parser.add_argument("--prodigal_mode", type=str, default="normal", choices=["normal", "meta", "gv"], help='Prodigal mode')
    parser.add_argument("--mmseqs_sensitivity", type=float, help='MMseqs2 sensitivity (optional, default: auto)')
    parser.add_argument("--foldseek_sensitivity", type=float, help='Foldseek sensitivity (optional, default: auto)')
    parser.add_argument("--toxin_db", type=str, default='/db/foldseek_tox_db_padded', help='Path to custom toxin FASTA (used for MMseqs2, HMM, and Foldseek databases)')
    parser.add_argument("--antitoxin_db", type=str, default='/db/foldseek_antitox_db_padded', help='Path to custom antitoxin FASTA (used for MMseqs2, HMM, and Foldseek databases)')
    parser.add_argument("--rna_toxin_db", type=str, help='Path to RNA toxin nucleotide database for MMseqs2')
    parser.add_argument("--rna_antitoxin_db", type=str, help='Path to RNA antitoxin nucleotide database for MMseqs2')
    parser.add_argument("-a", "--havalue", type=float, default=0.36, help='Ha-value threshold')
    parser.add_argument("-n", "--minimum_distance", type=int, default=-30, help='Minimum intergenic distances (bp)')
    parser.add_argument("--neighborhood_distance", type=int, default=200, help='Maximum distance between genes in neighborhood (bp)')
    parser.add_argument("--max_pairing_distance", type=int, default=150, help='Maximum distance between toxin and antitoxin for pairing (bp, default: 150)')
    parser.add_argument("--min_neighborhood_size", type=int, default=2, help='Minimum number of genes in a neighborhood (default: 2)')
    parser.add_argument("--max_neighborhood_size", type=int, default=5, help='Maximum number of genes in a neighborhood (default: 5)')
    parser.add_argument("-s", "--extract_seq", action="store_true", help='Extract predicted sequences')
    parser.add_argument("-p", "--predict_orphan", action="store_true", help='Predict orphan toxins and antitoxins')
    parser.add_argument("--save_unidentified", type=str, nargs='?', const='ta_neighborhoods.fna', help='Save unidentified neighborhoods')
    parser.add_argument("-t", "--tmp_dir", type=str, default=os.path.join(os.getcwd(), "tafinder_tmp"), help='Temp directory')
    parser.add_argument("--threads", type=int, default=os.cpu_count(), help='Number of threads')
    parser.add_argument("--no_bitscore", action="store_true", help='Disable bitscore filtering (use only e-value cutoffs)')
    parser.add_argument("--no_foldseek", action="store_true", help='Disable Foldseek')
    parser.add_argument("--no_hmm", action="store_true", help='Disable HMM')
    parser.add_argument("--no_mmseqs", action="store_true", help='Disable MMseqs2')
    parser.add_argument("--search_all", action="store_true", help='Search all genes, not just neighborhoods')
    parser.add_argument("--db_dir", type=str, default=os.path.join(os.getcwd(), "tafinder_db"), help='Path to tafinder database directory')
    parser.add_argument("--foldseek_query_db", type=str, help='Path to pre-built Foldseek query database (speeds up analysis by skipping database creation)')
    parser.add_argument("--prostt5_weights", type=str, help='Path to ProstT5 weights directory (default: db_dir/prostt5_weights)')
    return parser.parse_args()

def initialize_environment(args):
    import sys
    import shutil
    script_dir = os.path.dirname(os.path.realpath(__file__))

    # Remove existing tmp directory if it exists
    if os.path.exists(args.tmp_dir):
        shutil.rmtree(args.tmp_dir)

    os.makedirs(args.tmp_dir, exist_ok=True)

    # Detect GPU availability
    gpu_available = False
    try:
        import torch
        gpu_available = torch.cuda.is_available()
        if gpu_available:
            print(f"GPU detected: {torch.cuda.get_device_name(0)}")
        else:
            print("No GPU detected, using CPU")
    except ImportError:
        print("PyTorch not available, using CPU")

    # Determine database directory
    # Priority: 1. User-specified --db_dir, 2. Conda package location, 3. Current directory
    if args.db_dir != os.path.join(os.getcwd(), "tafinder_db"):
        # User explicitly specified a db_dir
        db_dir = args.db_dir
    else:
        # Check for conda-installed database first
        conda_db_dir = os.path.join(sys.prefix, "share", "tafinder3d", "tafinder_db")
        if os.path.exists(conda_db_dir):
            db_dir = conda_db_dir
        else:
            # Fall back to default (current directory)
            db_dir = args.db_dir

    # Determine database suffix based on GPU availability
    db_suffix = "_padded" if gpu_available else ""

    # Setup database paths
    return {
        'script_dir': script_dir,
        'fna_file': f"{args.tmp_dir}/sequence.fna",
        'ptt_file': f"{args.tmp_dir}/sequence.ptt",
        'faa_file': f"{args.tmp_dir}/sequence.faa",
        'prostt5_weights': args.prostt5_weights if args.prostt5_weights else os.path.join(db_dir, "prostt5_weights"),
        'gpu_available': gpu_available,

        # Unified protein databases (automatically use padded versions if GPU available)
        'mmseqs_protein_db': os.path.join(db_dir, f"mmseqs_combined_db{db_suffix}"),
        'foldseek_db': os.path.join(db_dir, f"foldseek_combined_db{db_suffix}"),

        # Legacy separate databases (for backward compatibility)
        'mmseqs_tox_db': os.path.join(db_dir, f"mmseqs_tox_db{db_suffix}"),
        'mmseqs_antitox_db': os.path.join(db_dir, f"mmseqs_antitox_db{db_suffix}"),
        'foldseek_tox_db': os.path.join(db_dir, f"foldseek_tox_db{db_suffix}"),
        'foldseek_antitox_db': os.path.join(db_dir, f"foldseek_antitox_db{db_suffix}"),

        # Unified HMM database (with toxin_/antitoxin_ prefixes)
        'hmm_combined_file': os.path.join(db_dir, "hmm_combined.hmm"),

        # Legacy separate HMM databases (for backward compatibility)
        'hmm_tox_file': os.path.join(db_dir, "tox_domains.hmm"),
        'hmm_antitox_file': os.path.join(db_dir, "antitox_domains.hmm"),

        # Unified nucleotide database (all RNA/DNA with type-specific prefixes)
        # Note: Nucleotide databases cannot be padded, so always use non-padded version
        'mmseqs_nucleotide_db': os.path.join(db_dir, "mmseqs_nucleotide_db"),

        'use_custom_db': False,

        # Family name lookup table
        'family_map': _load_family_map(db_dir),
    }


def _load_family_map(db_dir):
    """Load ta_family_map.tsv and ta_nucleotide_family_map.tsv if they exist.

    Returns a single merged dict of db_id -> protein_name.  The protein and
    nucleotide maps use distinct key prefixes so there are no collisions.
    """
    family_map = {}
    for tsv_name in ("ta_family_map.tsv", "ta_nucleotide_family_map.tsv"):
        map_file = os.path.join(db_dir, tsv_name)
        if not os.path.exists(map_file):
            continue
        count = 0
        with open(map_file) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("db_id"):
                    continue
                parts = line.split("\t", 1)
                if len(parts) == 2:
                    family_map[parts[0]] = parts[1]
                    count += 1
        print(f"Loaded {count} family name mappings from {tsv_name}")
    return family_map
