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

# Quality check report of the raw reads
rule fastqc_raw:
    input:
        f"{RAW_DIR}/{SRA}.fastq"
    output:
        html=f"{QC_DIR}/{SRA}_fastqc.html"
    shell:
        "fastqc -o {QC_DIR} {input}"

# Create FASTA index (.fai) used by samtools/gatk
rule samtools_faidx:
    input:
        f"{RAW_DIR}/reference.fasta"
    output:
        f"{RAW_DIR}/reference.fasta.fai"
    shell:
        "samtools faidx {input}"

# Build BWA index files (.amb/.ann/.bwt/.pac/.sa)
rule bwa_index:
    input:
        f"{RAW_DIR}/reference.fasta"
    output:
        expand("{ref}.{ext}", ref=f"{RAW_DIR}/reference.fasta", ext=["amb","ann","bwt","pac","sa"])
    shell:
        "bwa index {input}"

# Create a sequence dictionary (.dict) required by GATK
rule gatk_dict:
    input:
        f"{RAW_DIR}/reference.fasta"
    output:
        f"{RAW_DIR}/reference.dict"
    shell:
        "gatk CreateSequenceDictionary -R {input} -O {output}"

# Align reads and read-group info (required by many GATK tools)
rule bwa_mem:
    input:
        ref=f"{RAW_DIR}/reference.fasta",
        fq =f"{RAW_DIR}/{SRA}.fastq"
    output:
        f"{ALIGNED_DIR}/aligned.sam"
    shell:
        r"""bwa mem -R '@RG\tID:1\tLB:lib1\tPL:illumina\tPU:unit1\tSM:sample1' {input.ref} {input.fq} > {output}"""

# Convert SAM to BAM and coordinate-sort
rule sort_bam:
    input:
        f"{ALIGNED_DIR}/aligned.sam"
    output:
        f"{ALIGNED_DIR}/aligned.sorted.bam"
    shell:
        "samtools view -b {input} | samtools sort -o {output}"

# Sanity check (summary) to catch malformed BAMs early
rule validate_bam:
    input:
        f"{ALIGNED_DIR}/aligned.sorted.bam"
    output:
        f"{ALIGNED_DIR}/validated.txt"
    shell:
        "gatk ValidateSamFile -I {input} -MODE SUMMARY > {output}"

# Mark duplicates (PCR/optical dups), produce metrics
rule mark_duplicates:
    input:
        f"{ALIGNED_DIR}/aligned.sorted.bam"
    output:
        bam    = f"{ALIGNED_DIR}/dedup.bam",
        metrics= f"{ALIGNED_DIR}/dup_metrics.txt"
    shell:
        "gatk MarkDuplicates -I {input} -O {output.bam} -M {output.metrics}"

# Make BAM index (.bai) to enables random access for tools like IGV/GATK
rule index_dedup:
    input:
        f"{ALIGNED_DIR}/dedup.bam"
    output:
        f"{ALIGNED_DIR}/dedup.bam.bai"
    shell:
        "samtools index {input}"
