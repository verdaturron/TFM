# TFM

SupraWaves RNA-seq pipeline

This project holds a Snakemake workflow that turns raw FASTQ files from dog (MDCK) cells into:

    quality reports (FastQC, MultiQC)

    trimmed reads (fastp)

    transcript counts (Salmon) and gene counts (tximport)

    differential-expression tables (DESeq2)

    plots (volcano, heat-map, PCA)

    optional GO / KEGG enrichment

The code was written for my master thesis “Self-sustained velocity waves and pattern
emergence in tissues: mechanotranscriptomics
in confined epithelia” (Corentin Lacaze, 2025).

1. Quick start

bash

# clone the repo
git clone https://github.com/your-user/TFM.git
cd TFM

# install Snakemake (needs conda or mamba)
mamba create -n snakemake -c conda-forge -c bioconda snakemake=7.32
conda activate snakemake

# run the small test (needs ~5 min, 2 GB RAM)
snakemake -n                    # dry run – shows what will happen
snakemake --use-conda --cores 6 # real run

2. Input files
File / folder	What it is	Notes
fastq_dir/	Pair-end FASTQ files named <sample>_R1.fastq.gz, <sample>_R2.fastq.gz	gzip-compressed
config/samples.tsv	Two columns: sample and condition	one line per biological replicate
config/config.yaml	Settings for the run	edit paths, threads, cut-offs
ref/	Canis lupus familiaris Ensembl r112 cDNA FASTA + GTF	downloaded by the workflow if missing
3. How to run on your data

    Put your FASTQ files in fastq_dir/.

    Fill config/samples.tsv with your sample names and the group they belong to.

    Edit config/config.yaml so paths and cut-offs fit your project.

    Launch:

bash

snakemake --use-conda --cores 8 --rerun-incomplete --printshellcmds

The main results land in results/. A full HTML report is under results/multiqc/.
4. Output overview

    results/fastqc/ – raw read QC

    results/trimmed/ – trimmed FASTQ files

    results/salmon/ – quant.sf files per sample

    results/expression/ – gene_counts.tsv, gene_tpm.tsv

    results/deseq2/ – one TSV per contrast, plus a merged file

    results/plots/ – volcano, heat-map, PCA (PDF + PNG)

    results/enrichment/ – GO / KEGG tables and images (optional)

    results/report/workflow_dag.html – graph of the whole pipeline

5. Software

All tools run inside small conda environments listed in envs/.
Snakemake builds them automatically when you use --use-conda.

Main versions used in the thesis:
Tool	Version
Snakemake	7.32
fastp	0.24
Salmon	1.10.3
DESeq2	1.40
clusterProfiler	4.10
6. Troubleshooting

    Conda solver is slow – try mamba instead of conda.

    File not found – check paths in config.yaml match your folders.

    Few reads mapped – make sure you use the right species reference.

    If the run stops midway, fix the problem and re-launch; Snakemake will resume.

7. Citing

If this pipeline helped your work, please cite:

rust

Lacaze C. (2025)
Self-sustained velocity waves and pattern emergence in tissues: mechano-transcriptomics in confined epithelia.
Master thesis, Universitat oberta de catalunya and Universitat de Barcelona.

8. License

Code: MIT
Documentation and thesis text: CC BY-NC 3.0

Feel free to open an issue or pull request for questions or improvements.
