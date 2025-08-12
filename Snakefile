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
