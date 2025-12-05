# TAfinder3D

TAfinder3D automatically predicts addiction casettes in genomic or metagenomic data. It is an expansion of TAfinder2, which it improves on in the following ways:
* ***Structural search:*** TAfinder3D uses Prostt5 and Foldseek to map input contigs to precomputed structural TA databases. Because TA families often have a small number of members, this can dramatically improve detection rates. Structural search can be disabled with the `--no_foldseek` flag.
* ***Speed-ups:*** TAfinder3D uses MMseqs2 in place of BLAST, supports GPU-accelerated searches, and is parallelizable at each stage.
* ***Metagenomic handling:*** TAfinder3D cleanly parses multi-contig inputs. It provides `prodigal-meta` and `prodigal-gv` for ORF finding, improving reliability in metagenomic and phage queries.
* ***Updated database:*** TAfinder3D includes a few new addiction casettes published between 2023 and 2025; I intend to keep updating it about once a year.

## Getting started

TAfinder3D is provided through my conda channel:

`mamba install -c conda-forge -c bioconda -c schmigle tafinder3d`

It can be run with the following command:

`tafinder -i [INPUT FILE] -f [FILETYPE] -o [OUTPUT FILE] [OPTIONS]`

`-i/--infile`--the genome file you want to annotate. Can be a genbank or fasta file. 

`-f/--format`--format of the genome file. Specify `genbank` or `fna`.

`-o/--outfile`--name of the output file. This will be a `csv` file, though you can name it any string you like.

A complete list of options is provided in the help text, accessible with `tafinder -h`.

## Citations

TAfinder3D is provided under an MIT license. Please additionally cite the following tools when using it:

[MMseqs2](https://www.nature.com/articles/s41592-025-02819-8)

[Foldseek](https://www.nature.com/articles/s41587-023-01773-0)

[Hmmer3](http://eddylab.org/software/hmmer/hmmer.org)

[Prodigal](https://github.com/hyattpd/Prodigal)

[Prodigal-gv](https://github.com/apcamargo/prodigal-gv)

[TAfinder2](https://academic.oup.com/nar/article/52/D1/D784/7332075?login=false)
