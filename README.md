# SupraWaves RNA-seq pipeline
![version](https://img.shields.io/badge/version-0.5-blue) ![documentation](https://img.shields.io/badge/documentation-under--construction-red)

## Overview
This project contains a **Snakemake** workflow that processes raw FASTQ files from dog (*MDCK*) cells into:

- ✅ Quality reports (**FastQC**, **MultiQC**)  
- ✂️ Trimmed reads (**fastp**)  
- 🔢 Transcript counts (**Salmon**) and gene counts (**tximport**)  
- 📊 Differential expression tables (**DESeq2**)  
- 📈 Visualizations (volcano plots, heatmaps, PCA)  
- 🧬 *Optional*: GO / KEGG enrichment analysis

---

This workflow was developed as part of my master's thesis:  
**“Self-sustained velocity waves and pattern emergence in tissues: mechanotranscriptomics in confined epithelia”**  
**Corentin Lacaze, 2025**


----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

## 1. Quick start

```bash
# Clone the repo
git clone https://github.com/your-user/TFM.git
cd TFM

# Install Snakemake (requires conda or mamba)
mamba create -n snakemake -c conda-forge -c bioconda snakemake=7.32
conda activate snakemake

# Run the small test (~5 min, 2 GB RAM)
snakemake -n                      # dry run – shows what will happen
snakemake --use-conda -j6        # real run using 6 cores


```

# 2. Input files
| File / Folder        | What it is                                                            | Notes                                                                 |
|----------------------|------------------------------------------------------------------------|-----------------------------------------------------------------------|
| `fastq_dir/`         | Paired-end FASTQ files                                                 | Gzip-compressed, named with `_R1.fastq.gz` and `_R2.fastq.gz` suffixes |
| `config/samples.tsv` | Sample metadata                                                       | Two columns: `sample` and `condition`, one line per biological replicate |
| `config/config.yaml` | Configuration settings for the run                                    | Edit paths, number of threads, cut-off values, etc.                   |
| `ref/`               | Reference files for *Canis lupus familiaris* (Ensembl release 112)    | cDNA FASTA and GTF; auto-downloaded by the workflow if missing        |


# 3. How to run on your data

Put your FASTQ files in fastq_dir/.
Fill config/samples.tsv with your sample names and the group they belong to.
Edit config/config.yaml so paths and cut-offs fit your project.

```bash        
        snakemake --use-conda --cores 8 --rerun-incomplete --printshellcmds
```

The main results land in results/. A full HTML report is under results/multiqc/.



# 4. Output overview

- results/fastqc/ - raw read QC
    
- *results/trimmed/ - trimmed FASTQ files

- results/salmon/ - quant.sf files per sample

- results/expression/ - gene_counts.tsv, gene_tpm.tsv

- results/deseq2/ - one TSV per contrast, plus a merged file

- results/plots/ - volcano, heat-map, PCA (PDF + PNG)

- results/enrichment/ – GO / KEGG tables and images (optional)

- results/report/workflow_dag.html – graph of the whole pipeline


...

    results/report/workflow_dag.svg – graph of the whole pipeline

---

### 🔄 Pipeline DAG

The full structure of the pipeline is shown below (generated automatically by Snakemake):

![Workflow DAG](/docs/workflow_dag.svg)


# 5. Software

All tools run inside small conda environments listed in envs/.
Snakemake builds them automatically when you use --use-conda.

### Main versions used in the thesis

| Tool             | Version |
|------------------|---------|
| Snakemake        | 7.32    |
| fastp            | 0.24    |
| Salmon           | 1.10.3  |
| DESeq2           | 1.40    |
| clusterProfiler  | 4.10    |

# 6. Troubleshooting

- 🐢 Conda solver is slow: try mamba instead of conda.

- ❌ File not found: check paths in config.yaml match your folders.

- 📉 Few reads mapped: make sure you use the right species reference.

- 🔁 If the run stops midway, fix the problem and re-launch; Snakemake will resume.

# 7. Citing

If this pipeline helped your work, please cite:

Lacaze, C. (2025). *Self-sustained velocity waves and pattern emergence in tissues: mechanotranscriptomics in confined epithelia*.  
Master’s thesis, Universitat Oberta de Catalunya and Universitat de Barcelona.

# 8. License

Code: MIT
Documentation and thesis text: CC BY-NC 3.0

Feel free to open an issue or pull request for questions or improvements.
