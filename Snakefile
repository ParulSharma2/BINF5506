# Variables (edit if needed)

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


# Final target(s)

rule all:
    input:
        f"{SNAKEMAKE_DIR}/.dirs_created",
        f"{RAW_DIR}/reference.fasta",
        f"{RAW_DIR}/{SRA}/{SRA}.sra",
        f"{RAW_DIR}/{SRA}.fastq",
        f"{QC_DIR}/{SRA}_fastqc.html",
        f"{RAW_DIR}/reference.fasta.fai",
        f"{RAW_DIR}/reference.dict",
        f"{ALIGNED_DIR}/dedup.bam.bai",
        f"{VARIANT_DIR}/raw_variants.vcf",
        f"{VARIANT_DIR}/filtered_variants.vcf",
        f"{SNPEFF_DATA_DIR}/genes.gbk",
        f"{SNPEFF_DIR}/snpEff.config",
        f"{SNPEFF_DIR}/snpEff_reference_db.txt",
        f"{ANNOTATED_DIR}/annotated_variants.vcf",
        f"{SNPEFF_DIR}/snpEff.html"


# Utility: make dirs

rule create_dirs:
    output:
        marker = f"{SNAKEMAKE_DIR}/.dirs_created"
    shell:
        r"""
        mkdir -p {RESULTS_FOLDER} {RAW_DIR} {ALIGNED_DIR} {VARIANT_DIR} {ANNOTATED_DIR} {QC_DIR} {SNPEFF_DATA_DIR} {SNAKEMAKE_DIR}
        touch {output.marker}
        """