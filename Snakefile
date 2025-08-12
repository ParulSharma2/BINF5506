# Parameters and folder shortcuts
SRA = "SRR1972739"            # SRA run to download
REF_ID = "AF086833.2"         # Ebola reference accession
RESULTS_FOLDER = "results"    # base output folder

# Directory shortcuts used in rules (avoids repetition & typos)
RAW_DIR       = f"{RESULTS_FOLDER}/raw"
ALIGNED_DIR   = f"{RESULTS_FOLDER}/aligned"
VARIANT_DIR   = f"{RESULTS_FOLDER}/variants"
ANNOTATED_DIR = f"{RESULTS_FOLDER}/annotated"
QC_DIR        = f"{RESULTS_FOLDER}/qc"
SNPEFF_DIR    = f"{RESULTS_FOLDER}/snpEff"
SNPEFF_DATA_DIR = f"{SNPEFF_DIR}/data/reference_db"
SNAKEMAKE_DIR   = f"{RESULTS_FOLDER}/snakemake"

# Rule all
# The 'all' rule says: “When I run snakemake, these are the final files I want.”
# Snakemake will figure out the dependency chain to produce them.
rule all:
    input:
        f"{RAW_DIR}/reference.fasta",            # downloaded reference
        f"{RAW_DIR}/{SRA}/{SRA}.sra",            # raw SRA container
        f"{RAW_DIR}/{SRA}.fastq",                # extracted FASTQ (subset)
        f"{QC_DIR}/{SRA}_fastqc.html",           # FastQC report
        f"{ALIGNED_DIR}/dedup.bam.bai",          # BAM index (ensures BAM exists + indexed)
        f"{VARIANT_DIR}/filtered_variants.vcf",  # filtered variants
        f"{ANNOTATED_DIR}/annotated_variants.vcf" # snpEff annotated VCF

# Create output folders once and a marker file
rule create_dirs:
    output:
        marker = f"{SNAKEMAKE_DIR}/.dirs_created"
    shell:
        "mkdir -p {RESULTS_FOLDER} {RAW_DIR} {ALIGNED_DIR} {VARIANT_DIR} {ANNOTATED_DIR} {QC_DIR} {SNPEFF_DATA_DIR} {SNAKEMAKE_DIR} && touch {output.marker}"

# Download reference as FASTA with NCBI Entrez Direct
rule download_reference:
    input:
        rules.create_dirs.output.marker   # ensure directories exist
    output:
        f"{RAW_DIR}/reference.fasta"
    shell:
        "efetch -db nucleotide -id {REF_ID} -format fasta > {output}"

# Fetch SRA run (container file) to results/raw/SRA/SRA.sra
rule download_sra:
    input:
        rules.create_dirs.output.marker
    output:
        f"{RAW_DIR}/{SRA}/{SRA}.sra"
    shell:
        "prefetch {SRA} -O {RAW_DIR}"

# Extract a small FASTQ subset (10k spots) for speed
rule extract_sequence:
    input:
        rules.download_sra.output  # the .sra file
    output:
        f"{RAW_DIR}/{SRA}.fastq"
    shell:
        "fastq-dump -X 10000 {input} -O {RAW_DIR}"
