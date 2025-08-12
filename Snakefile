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
       # f"{VARIANT_DIR}/raw_variants.vcf",
       # f"{VARIANT_DIR}/filtered_variants.vcf",
       # f"{SNPEFF_DATA_DIR}/genes.gbk",
       # f"{SNPEFF_DIR}/snpEff.config",
       # f"{SNPEFF_DIR}/snpEff_reference_db.txt",
       # f"{ANNOTATED_DIR}/annotated_variants.vcf",
       # f"{SNPEFF_DIR}/snpEff.html"


# Utility: make dirs

rule create_dirs:
    output:
        marker = f"{SNAKEMAKE_DIR}/.dirs_created"
    shell:
        r"""
        mkdir -p {RESULTS_FOLDER} {RAW_DIR} {ALIGNED_DIR} {VARIANT_DIR} {ANNOTATED_DIR} {QC_DIR} {SNPEFF_DATA_DIR} {SNAKEMAKE_DIR}
        touch {output.marker}
        """
# Data Downloading
# Reference: FASTA file

rule download_reference:
    input:
        rules.create_dirs.output.marker
    output:
        reference_fasta = f"{RAW_DIR}/reference.fasta"
    shell:
        r"""
        echo "Downloading reference genome..."
        efetch -db nucleotide -id {REF_ID} -format fasta > {output.reference_fasta}
        test -s {output.reference_fasta}
        echo "Downloaded reference genome!"
        """


# SRA: download SRA

rule download_sra:
    input:
        rules.create_dirs.output.marker
    output:
        sra = f"{RAW_DIR}/{SRA}/{SRA}.sra"
    shell:
        r"""
        echo "Downloading sequencing data..."
        prefetch {SRA} -O {RAW_DIR}
        test -s {output.sra}
        echo "Downloaded sequencing data!"
        """


# SRA: extract to FASTQ

rule extract_fastq:
    input:
        rules.download_sra.output.sra
    output:
        fastq = f"{RAW_DIR}/{SRA}.fastq"
    shell:
        r"""
        echo "Extracting sequencing data..."
        fastq-dump -X 10000 {input} -O {RAW_DIR}
        test -s {output.fastq}
        echo "Extracted sequencing data!"
        """
# QC: FastQC on FASTQ

rule fastqc_raw:
    input:
        f"{RAW_DIR}/{SRA}.fastq"
    output:
        html = f"{QC_DIR}/{SRA}_fastqc.html",
        zip = f"{QC_DIR}/{SRA}_fastqc.zip"
    shell:
        r"""
        echo "Running FastQC..."
        fastqc -o {QC_DIR} {input}
        test -s {output.html}
        test -s {output.zip}
        """


# Reference indexing (faidx)

rule samtools_faidx:
    input:
        f"{RAW_DIR}/reference.fasta"
    output:
        f"{RAW_DIR}/reference.fasta.fai"
    shell:
        r"""
        echo "Indexing reference with samtools faidx..."
        samtools faidx {input}
        """


# BWA index reference

rule bwa_index:
    input:
        f"{RAW_DIR}/reference.fasta"
    output:
        # BWA creates multiple files; use a sentinel file to avoid listing them all
        touch(f"{RAW_DIR}/.bwa_index_done")
    shell:
        r"""
        echo "Building BWA index..."
        bwa index {input}
        touch {output}
        """


# GATK CreateSequenceDictionary

rule gatk_dict:
    input:
        f"{RAW_DIR}/reference.fasta"
    output:
        f"{RAW_DIR}/reference.dict"
    shell:
        r"""
        echo "Creating sequence dictionary..."
        gatk CreateSequenceDictionary -R {input} -O {output}
        """
# Align with BWA-MEM

rule align_bwa_mem:
    input:
        ref = f"{RAW_DIR}/reference.fasta",
        fastq = f"{RAW_DIR}/{SRA}.fastq",
        bwa_idx_done = f"{RAW_DIR}/.bwa_index_done"
    output:
        sam = f"{ALIGNED_DIR}/aligned.sam"
    threads: 4
    shell:
        r"""
        echo "Aligning reads with BWA-MEM..."
        bwa mem -t {threads} -R '@RG\tID:1\tLB:lib1\tPL:illumina\tPU:unit1\tSM:sample1' {input.ref} {input.fastq} > {output.sam}
        test -s {output.sam}
        """


# Convert to sorted BAM

rule sam_to_sorted_bam:
    input:
        f"{ALIGNED_DIR}/aligned.sam"
    output:
        bam = f"{ALIGNED_DIR}/aligned.sorted.bam"
    threads: 4
    shell:
        r"""
        echo "Converting SAM to sorted BAM..."
        samtools view -@ {threads} -b {input} | samtools sort -@ {threads} -o {output.bam}
        test -s {output.bam}
        """


# Validate BAM (GATK)

rule validate_bam:
    input:
        bam = f"{ALIGNED_DIR}/aligned.sorted.bam"
    output:
        touch(f"{ALIGNED_DIR}/aligned.sorted.bam.validated")
    shell:
        r"""
        echo "Validating BAM..."
        gatk ValidateSamFile -I {input.bam} -MODE SUMMARY
        touch {output}
        """


# Mark duplicates (GATK)

rule mark_duplicates:
    input:
        bam = f"{ALIGNED_DIR}/aligned.sorted.bam",
        ok = f"{ALIGNED_DIR}/aligned.sorted.bam.validated"
    output:
        dedup = f"{ALIGNED_DIR}/dedup.bam",
        metrics = f"{ALIGNED_DIR}/dup_metrics.txt"
    shell:
        r"""
        echo "Marking duplicates..."
        gatk MarkDuplicates -I {input.bam} -O {output.dedup} -M {output.metrics}
        test -s {output.dedup}
        test -s {output.metrics}
        """


# Index deduplicated BAM (bai)

rule index_dedup_bam:
    input:
        f"{ALIGNED_DIR}/dedup.bam"
    output:
        f"{ALIGNED_DIR}/dedup.bam.bai"
    shell:
        r"""
        echo "Indexing deduplicated BAM..."
        samtools index {input}
        """