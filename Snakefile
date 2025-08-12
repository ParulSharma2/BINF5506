# --- VARIABLES (same as your bash script) ---
SRA = "SRR1972739"
REF_ID = "AF086833.2"
RESULTS_FOLDER = "results"
RAW_DIR = f"{RESULTS_FOLDER}/raw"
ALIGNED_DIR = f"{RESULTS_FOLDER}/aligned"
VARIANT_DIR = f"{RESULTS_FOLDER}/variants"
ANNOTATED_DIR = f"{RESULTS_FOLDER}/annotated"
QC_DIR = f"{RESULTS_FOLDER}/qc"
SNPEFF_DIR = f"{RESULTS_FOLDER}/snpEff"
SNPEFF_DATA_DIR = f"{SNPEFF_DIR}/data/reference_db"
SNAKEMAKE_DIR = f"{RESULTS_FOLDER}/snakemake"

# --- S3 SETTINGS ---
BUCKET = "sohail-binf55062"
S3_PREFIX = "ebola"

# --- FINAL TARGETS ---
rule all:
    input:
        f"{SNAKEMAKE_DIR}/.dirs_created",
        f"{RAW_DIR}/reference.fasta",
        f"{RAW_DIR}/{SRA}/{SRA}.sra",
        f"{RAW_DIR}/{SRA}.fastq"

rule create_dirs:
    output:
        marker = f"{SNAKEMAKE_DIR}/.dirs_created"
    shell:
        """
        mkdir -p {RESULTS_FOLDER} {RAW_DIR} {ALIGNED_DIR} {VARIANT_DIR} {ANNOTATED_DIR} {QC_DIR} {SNPEFF_DATA_DIR} {SNAKEMAKE_DIR}
        touch {SNAKEMAKE_DIR}/.dirs_created
       """
# Download reference as FASTA with NCBI Entrez Direct
rule download_reference:
    input:
        f"{SNAKEMAKE_DIR}/.dirs_created"
    output:
        fasta = f"{RAW_DIR}/reference.fasta"
    shell:
        """
        efetch -db nucleotide -id {REF_ID} -format fasta > {output.fasta}
        """

rule download_sra:
    input:
        marker = rules.create_dirs.output.marker
    output:
        sra = f"{RAW_DIR}/{SRA}/{SRA}.sra"
    shell:
        """
        prefetch {SRA} -O {RAW_DIR}
        test -s {output.sra}
        """

rule extract_fastq:
    input:
        sra = rules.download_sra.output.sra
    output:
        fastq = f"{RAW_DIR}/{SRA}.fastq"
    shell:
        """
        fastq-dump -X 10000 {input.sra} -O {RAW_DIR}
        test -s {output.fastq}
        """