configfile: "config/config.yaml"

SAMPLES     = config["samples"]
FASTQ_DIR   = config["fastq_dir"]
RESULTS     = config["results_dir"]
THREADS     = config["threads"]

FASTQC  = "fastqc"
FASTP   = "fastp"
MULTIQC = "multiqc"

REF          = config["ref"]
TRANSCRIPTS  = REF["transcript_fasta"]
GTF          = REF["gtf"]
SALMON_INDEX = REF["salmon_index"]

## RÈGLE SOMMAIRE (tout ce qu’on veut en sortie)
import itertools

# dé-duplique et ordonne les conditions déclarées dans le YAML
CONDITIONS = sorted(set(config["conditions"]))

# génère uniquement les combinaisons où cond1 != cond2
PAIRS = [(a, b) for a, b in itertools.combinations(CONDITIONS, 2)]


## Lire la liste des conditions depuis le YAML & générer les paires

CONDITIONS = config["conditions"]
PAIRS = list(itertools.combinations(CONDITIONS, 2))

PAIR_TSVS = [
    f"results/deseq2/deseq2_{c1}_vs_{c2}.tsv"
    for c1, c2 in PAIRS
]

rule all:
    input:
        # 1. FastQC bruts
        expand(f"{RESULTS}/fastqc/{{sample}}_R1_fastqc.html", sample=SAMPLES),
        expand(f"{RESULTS}/fastqc/{{sample}}_R2_fastqc.html", sample=SAMPLES),

        # 2. Lectures trimmed
        expand(f"{RESULTS}/trimmed/{{sample}}_R1_trimmed.fastq.gz", sample=SAMPLES),
        expand(f"{RESULTS}/trimmed/{{sample}}_R2_trimmed.fastq.gz", sample=SAMPLES),

        # 3. FastQC post-trim
        expand(f"{RESULTS}/fastqc_posttrim/{{sample}}_R1_trimmed_fastqc.html", sample=SAMPLES),
        expand(f"{RESULTS}/fastqc_posttrim/{{sample}}_R2_trimmed_fastqc.html", sample=SAMPLES),

        # 4. Quantification Salmon
        expand(f"{RESULTS}/salmon/{{sample}}/quant.sf", sample=SAMPLES),

        # 5. Rapport global
        f"{RESULTS}/multiqc/multiqc_report.html",

        # 6. Expression génique
        f"{RESULTS}/expression/gene_counts.tsv",

        # Les sorties déjà existantes (FastQC, trim, Salmon, tximport, etc.)
        expand(f"{RESULTS}/fastqc/{{sample}}_R1_fastqc.html",sample=SAMPLES), \
        expand(f"{RESULTS}/fastqc/{{sample}}_R2_fastqc.html",sample=SAMPLES), \
        expand(f"{RESULTS}/trimmed/{{sample}}_R1_trimmed.fastq.gz",sample=SAMPLES), \
        expand(f"{RESULTS}/trimmed/{{sample}}_R2_trimmed.fastq.gz",sample=SAMPLES), \
        expand(f"{RESULTS}/fastqc_posttrim/{{sample}}_R1_trimmed_fastqc.html",sample=SAMPLES), \
        expand(f"{RESULTS}/fastqc_posttrim/{{sample}}_R2_trimmed_fastqc.html",sample=SAMPLES), \
        expand(f"{RESULTS}/salmon/{{sample}}/quant.sf",sample=SAMPLES), \
        f"{RESULTS}/multiqc/multiqc_report.html", \
        f"{RESULTS}/expression/gene_counts.tsv",

        # Tous les TSV DE deux-à-deux
        *PAIR_TSVS,

        # Le TSV fusionné
        "results/deseq2/deseq2_all_pairwise.tsv",

        # les volcano-plots PDF
        *[
        f"results/plots/volcano_{c1}_vs_{c2}.pdf"
        for c1, c2 in PAIRS
        ],

        # 7. Volcano plot
        "results/plots/volcano_border_vs_central.pdf",

        # 8. test qualité trim
        expand(f"{RESULTS}/qc_hist/{{sample}}_lenqual.pdf", sample=SAMPLES),

        # 9. Barplot taux d’alignement Salmon

        f"{RESULTS}/qc_hist/salmon_alignment_rates.pdf",

        # 10. Heatmap des gènes DE

        expand("results/plots/heatmap_{c1}_vs_{c2}.pdf", c1=[p[0] for p in PAIRS], c2=[p[1] for p in PAIRS]),
        expand("results/plots/heatmap_{c1}_vs_{c2}.png", c1=[p[0] for p in PAIRS], c2=[p[1] for p in PAIRS]),

        # rapport workflow DAG

        f"{RESULTS}/report/workflow_dag.html",

# ---------- 1. FASTQC sur données brutes ----------

rule fastqc_raw:
    input:
        r1 = lambda w: f"{FASTQ_DIR}/{w.sample}_R1.fastq.gz",
        r2 = lambda w: f"{FASTQ_DIR}/{w.sample}_R2.fastq.gz"
    output:
        r1_html = f"{RESULTS}/fastqc/{{sample}}_R1_fastqc.html",
        r1_zip  = f"{RESULTS}/fastqc/{{sample}}_R1_fastqc.zip",
        r2_html = f"{RESULTS}/fastqc/{{sample}}_R2_fastqc.html",
        r2_zip  = f"{RESULTS}/fastqc/{{sample}}_R2_fastqc.zip"
    threads: THREADS
    shell:
        """
        mkdir -p {RESULTS}/fastqc
        {FASTQC} "{input.r1}" "{input.r2}" --outdir {RESULTS}/fastqc -t {threads}
        """

# ---------- 2. Trim avec fastp ----------

rule fastp:
    input:
        r1 = lambda w: f"{FASTQ_DIR}/{w.sample}_R1.fastq.gz",
        r2 = lambda w: f"{FASTQ_DIR}/{w.sample}_R2.fastq.gz"
    output:
        r1 = f"{RESULTS}/trimmed/{{sample}}_R1_trimmed.fastq.gz",
        r2 = f"{RESULTS}/trimmed/{{sample}}_R2_trimmed.fastq.gz"
    threads: THREADS
    shell:
        """
        mkdir -p {RESULTS}/trimmed
        {FASTP} -i "{input.r1}" -I "{input.r2}" \
                -o "{output.r1}" -O "{output.r2}" \
                -w {threads} --detect_adapter_for_pe
        """

# ---------- 3. FASTQC après trim ----------

rule fastqc_trimmed:
    input:
        r1 = f"{RESULTS}/trimmed/{{sample}}_R1_trimmed.fastq.gz",
        r2 = f"{RESULTS}/trimmed/{{sample}}_R2_trimmed.fastq.gz"
    output:
        r1_html = f"{RESULTS}/fastqc_posttrim/{{sample}}_R1_trimmed_fastqc.html",
        r1_zip  = f"{RESULTS}/fastqc_posttrim/{{sample}}_R1_trimmed_fastqc.zip",
        r2_html = f"{RESULTS}/fastqc_posttrim/{{sample}}_R2_trimmed_fastqc.html",
        r2_zip  = f"{RESULTS}/fastqc_posttrim/{{sample}}_R2_trimmed_fastqc.zip"
    threads: THREADS
    shell:
        """
        mkdir -p {RESULTS}/fastqc_posttrim
        {FASTQC} "{input.r1}" "{input.r2}" --outdir {RESULTS}/fastqc_posttrim -t {threads}
        """

# ---------- Histogrammes qualité / longueur (Before vs After) ----------

rule fastqc_histograms:
    input:
        raw_zip  = f"{RESULTS}/fastqc/{{sample}}_R1_fastqc.zip",
        trim_zip = f"{RESULTS}/fastqc_posttrim/{{sample}}_R1_trimmed_fastqc.zip"
    output:
        pdf = f"{RESULTS}/qc_hist/{{sample}}_lenqual.pdf"
    conda: "envs/qc_plots.yaml"
    script: "scripts/plot_fastqc_hist.R"

# ---------- 4. Rapport global MultiQC ----------

rule multiqc:
    input:
        # zips bruts
        expand(f"{RESULTS}/fastqc/{{sample}}_R1_fastqc.zip", sample=SAMPLES) +
        expand(f"{RESULTS}/fastqc/{{sample}}_R2_fastqc.zip", sample=SAMPLES) +
        # zips trim
        expand(f"{RESULTS}/fastqc_posttrim/{{sample}}_R1_trimmed_fastqc.zip", sample=SAMPLES) +
        expand(f"{RESULTS}/fastqc_posttrim/{{sample}}_R2_trimmed_fastqc.zip", sample=SAMPLES)
    output:
        html = f"{RESULTS}/multiqc/multiqc_report.html"
    threads: 1
    shell:
        """
        mkdir -p {RESULTS}/multiqc
        # On scanne les deux dossiers FastQC
        {MULTIQC} {RESULTS}/fastqc {RESULTS}/fastqc_posttrim -o {RESULTS}/multiqc -f
        """

# ---------- Barplot taux d’alignement Salmon ----------

rule salmon_alignment_barplot:
    input:
        expand(f"{RESULTS}/salmon/{{sample}}/aux_info/meta_info.json",
               sample=SAMPLES)
    output:
        pdf = f"{RESULTS}/qc_hist/salmon_alignment_rates.pdf"
    conda: "envs/salmon_qc.yaml"
    script: "scripts/plot_salmon_alignment.R"

# -------------5 salmon quantification rule -------------

rule download_transcripts:
    output:
        fasta = config["ref"]["transcript_fasta"],
        gtf   = config["ref"]["gtf"]
    params:
        base = "https://ftp.ensembl.org/pub/release-112"
    shell:
        r"""
        mkdir -p ref
        wget -c -O {output.fasta} {params.base}/fasta/canis_lupus_familiaris/cdna/Canis_lupus_familiaris.ROS_Cfam_1.0.cdna.all.fa.gz
        wget -c -O {output.gtf}   {params.base}/gtf/canis_lupus_familiaris/Canis_lupus_familiaris.ROS_Cfam_1.0.112.gtf.gz
        """

rule salmon_index:
    input:
        fasta = TRANSCRIPTS
    output:
        directory(SALMON_INDEX)
    threads: THREADS
    conda: "envs/rnaseq.yaml"
    shell:
        """
        mkdir -p {output}
        salmon index -t {input.fasta} -i {output} -p {threads}
        """

# -------------- 6. Salmon quantification --------------

rule quant:
    input:
        index = SALMON_INDEX,
        r1    = f"{RESULTS}/trimmed/{{sample}}_R1_trimmed.fastq.gz",
        r2    = f"{RESULTS}/trimmed/{{sample}}_R2_trimmed.fastq.gz"
    output:
        quant = f"{RESULTS}/salmon/{{sample}}/quant.sf"
    threads: THREADS
    conda: "envs/salmon.yaml"
    shell:
        """
        mkdir -p $(dirname {output.quant})
        salmon quant \
            -i {input.index} \
            -l A \
            -1 {input.r1} -2 {input.r2} \
            -p {threads} \
            -o $(dirname {output.quant})
        """

# -------------- 7. Agrégation transcript -> gène -----------------

rule tximport:
    input:
        quant = expand(f"{RESULTS}/salmon/{{sample}}/quant.sf", sample=SAMPLES),
        gtf   = GTF
    output:
        counts = f"{RESULTS}/expression/gene_counts.tsv",
        tpm    = f"{RESULTS}/expression/gene_tpm.tsv"
    threads: 1
    conda: "envs/tximport.yaml"
    script: "scripts/tximport.R"

# -------------- 8. deseq2 --------------

rule deseq2_pairwise:
    input:
        counts = "results/expression/gene_counts.tsv",
        meta   = "config/samples.tsv"
    output:
        table = "results/deseq2/deseq2_{cond1}_vs_{cond2}.tsv"
    conda: "envs/deseq2.yaml"
    script: "scripts/deseq2_pairwise.R"

rule merge_pairwise_DE:
    input:
        *PAIR_TSVS
    output:
        merged = "results/deseq2/deseq2_all_pairwise.tsv"
    conda: "envs/deseq2.yaml"
    script: "scripts/merge_pairwise_DE.R"

# -------------- 9. Volcano plot P--------------

rule volcano_pairwise:
    input:
        "results/deseq2/deseq2_{cond1}_vs_{cond2}.tsv"
    output:
        pdf = "results/plots/volcano_{cond1}_vs_{cond2}.pdf",
        png = "results/plots/volcano_{cond1}_vs_{cond2}.png"
    params:
        padj = config["volcano"]["padj_cutoff"],
        lfc  = config["volcano"]["lfc_cutoff"]
    conda: "envs/volcano.yaml"
    script: "scripts/volcano_pairwise.R"


# Heatmap : 30 gènes DE

rule heatmap_pairwise:
    input:
        de      = "results/deseq2/deseq2_{cond1}_vs_{cond2}.tsv",
        counts  = "results/expression/gene_counts.tsv",
        meta    = "config/samples.tsv"
    output:
        pdf = "results/plots/heatmap_{cond1}_vs_{cond2}.pdf",
        png = "results/plots/heatmap_{cond1}_vs_{cond2}.png"
    params:
        top = 30
    conda:
        "envs/heatmap.yaml"
    script:
        "scripts/heatmap_DE.R"

# PCA VST

rule pca_vst:
    input:
        counts = "results/expression/gene_counts.tsv",
        meta   = "config/samples.tsv"
    output:
        pdf = "results/plots/pca_vst.pdf",
        png = "results/plots/pca_vst.png"
    conda: "envs/pca.yaml"
    script: "scripts/pca_vst.R"

# Règle : instantané du DAG Snakemake

rule dag_snapshot:
    output:
        svg  = f"{RESULTS}/report/workflow_dag.svg",
        html = f"{RESULTS}/report/workflow_dag.html"
    conda: "envs/graphviz.yaml"
    shell:
        r"""
        mkdir -p {RESULTS}/report
        snakemake --quiet --dag | dot -Tsvg -o {output.svg}

        # wrapper HTML plein-écran (heredoc non indenté)
cat > {output.html} <<'HTML'
<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>Snakemake DAG</title>
  <style>
    body{{margin:0}}
    svg{{width:100vw;height:100vh}}
  </style>
</head>
<body>
  <object data="workflow_dag.svg" type="image/svg+xml"></object>
</body>
</html>
HTML
        """


# -------------- 10. enrichissemement --------------

rule functional_enrichment:
    input:
        deseq = "results/deseq2/deseq2_results.tsv"
    output:
        dir = directory("results/enrichment")
    conda:
        "envs/clusterprofiler.yaml"
    script:
        "scripts/functional_enrichment.R"

