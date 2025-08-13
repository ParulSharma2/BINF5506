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


# S3 Bucket for data storage
BUCKET = "parul-binf5506"
S3_PREFIX = "ebola"

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
        f"{SNPEFF_DIR}/snpEff.html",
        f"{SNAKEMAKE_DIR}/.s3_upload_done" 

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
# Call variants (HCaller)

rule haplotype_caller:
    input:
        ref = f"{RAW_DIR}/reference.fasta",
        bai = f"{ALIGNED_DIR}/dedup.bam.bai"
    output:
        vcf = f"{VARIANT_DIR}/raw_variants.vcf"
    shell:
        r"""
        echo "Calling variants with GATK HaplotypeCaller..."
        gatk HaplotypeCaller -R {input.ref} -I {ALIGNED_DIR}/dedup.bam -O {output.vcf}
        test -s {output.vcf}
        """


# Filter variants (GATK VFilt)

rule filter_variants:
    input:
        ref = f"{RAW_DIR}/reference.fasta",
        vcf = f"{VARIANT_DIR}/raw_variants.vcf"
    output:
        f"{VARIANT_DIR}/filtered_variants.vcf"
    shell:
        r"""
        echo "Filtering variants..."
        gatk VariantFiltration -R {input.ref} -V {input.vcf} -O {output} \
            --filter-expression "QD < 2.0 || FS > 60.0" --filter-name FILTER
        test -s {output}
        """

# snpEff: download GenBank for custom DB

rule snpeff_fetch_gbk:
    input:
        rules.create_dirs.output.marker
    output:
        gbk = f"{SNPEFF_DATA_DIR}/genes.gbk"
    shell:
        r"""
        echo "Downloading reference GenBank file for snpEff..."
        efetch -db nucleotide -id {REF_ID} -format genbank > {output.gbk}
        test -s {output.gbk}
        echo "Downloaded GenBank file for snpEff!"
        """

# snpEff: write minimal config

rule snpeff_write_config:
    input:
        ref = f"{RAW_DIR}/reference.fasta",
        gbk = f"{SNPEFF_DATA_DIR}/genes.gbk"
    output:
        cfg = f"{SNPEFF_DIR}/snpEff.config"
    shell:
        r"""
        echo "Creating custom snpEff configuration file..."
        mkdir -p {SNPEFF_DIR}
        cat > {output.cfg} <<EOF
# Custom snpEff config for reference_db
reference_db.genome : reference_db
reference_db.fa : {input.ref}
reference_db.genbank : {input.gbk}
EOF
        test -s {output.cfg}
        """

# snpEff: build custom database
rule snpeff_build_db:
    input:
        cfg = f"{SNPEFF_DIR}/snpEff.config"
    output:
        touch(f"{SNPEFF_DIR}/.snpeff_build_done")
    shell:
        r"""
        echo "Building snpEff database..."
        snpEff build -c {input.cfg} -genbank -v -noCheckProtein reference_db
        touch {output}
        """

# snpEff: dump (export)
rule snpeff_dump_db:
    input:
        cfg = f"{SNPEFF_DIR}/snpEff.config",
        built = f"{SNPEFF_DIR}/.snpeff_build_done"
    output:
        dump = f"{SNPEFF_DIR}/snpEff_reference_db.txt"
    shell:
        r"""
        echo "Exporting snpEff database..."
        snpEff dump -c {input.cfg} reference_db > {output.dump}
        test -s {output.dump}
        """

# snpEff: annotate VCF
rule snpeff_annotate:
    input:
        cfg = f"{SNPEFF_DIR}/snpEff.config",
        vcf = f"{VARIANT_DIR}/filtered_variants.vcf",
        built = f"{SNPEFF_DIR}/.snpeff_build_done"
    output:
        anno_vcf = f"{ANNOTATED_DIR}/annotated_variants.vcf",
        html = f"{SNPEFF_DIR}/snpEff.html"
    shell:
        r"""
        echo "Annotating variants with snpEff..."
        snpEff -c {input.cfg} -stats {output.html} reference_db {input.vcf} > {output.anno_vcf}
        test -s {output.anno_vcf}
        test -s {output.html}
        echo "Annotation complete!"
        """

# FINAL STEP: Upload to S3

rule s3_upload:
    input:
        # Depend on end-of-pipeline outputs so this runs last
        rules.fastqc_raw.output.html,
        rules.samtools_faidx.output,
        rules.gatk_dict.output,
        rules.index_dedup_bam.output,
        rules.filter_variants.output,
        rules.snpeff_annotate.output.anno_vcf,
        rules.snpeff_annotate.output.html
    output:
        marker = f"{SNAKEMAKE_DIR}/.s3_upload_done"
    run:
        import os
        import boto3

        bucket = BUCKET
        prefix = S3_PREFIX.strip("/")

        s3 = boto3.client("s3")
        for root, dirs, files in os.walk(RESULTS_FOLDER):
            for file in files:
                local_file = os.path.join(root, file)
                # avoid uploading the marker before it's written
                if os.path.abspath(local_file) == os.path.abspath(output.marker):
                    continue
                rel = os.path.relpath(local_file, RESULTS_FOLDER)
                s3_key = os.path.join(prefix, rel).replace("\\", "/")
                print(f"Uploading {local_file} -> s3://{bucket}/{s3_key}")
                s3.upload_file(local_file, bucket, s3_key)

        with open(output.marker, "w") as f:
            f.write("Upload Complete!\n")
        print("All files uploaded to S3.")