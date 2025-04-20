#!/bin/bash

# Germline Variant Calling and Annotation Pipeline (Simplified)
# Based on GATK4 Best Practices

# Define directories
REF="~/demo/hg38/hg38.fa"
KNOWN_SITES="~/demo/hg38/Homo_sapiens_assembly38.dbsnp138.vcf"
READ1="~/demo/reads/SRR062634_1.filt.fastq.gz"
READ2="~/demo/reads/SRR062634_2.filt.fastq.gz"
OUTDIR="~/demo/results"
FUNCDIR="~/demo/funcotator_dataSources.v1.7.20200521g"

mkdir -p $OUTDIR

# STEP 1: Align with BWA
bwa index $REF
bwa mem -t 4 -R "@RG\tID:sample\tSM:sample\tPL:ILLUMINA" $REF $READ1 $READ2 > $OUTDIR/sample.sam

# STEP 2: Convert, Sort, Mark Duplicates
gatk MarkDuplicatesSpark -I $OUTDIR/sample.sam -O $OUTDIR/sample.sorted.dedup.bam

# STEP 3: BQSR
gatk BaseRecalibrator -I $OUTDIR/sample.sorted.dedup.bam -R $REF --known-sites $KNOWN_SITES -O $OUTDIR/recal_data.table
gatk ApplyBQSR -R $REF -I $OUTDIR/sample.sorted.dedup.bam --bqsr-recal-file $OUTDIR/recal_data.table -O $OUTDIR/sample.bqsr.bam

# STEP 4: Variant Calling
gatk HaplotypeCaller -R $REF -I $OUTDIR/sample.bqsr.bam -O $OUTDIR/raw_variants.vcf

# STEP 5: Select SNPs and INDELs
gatk SelectVariants -R $REF -V $OUTDIR/raw_variants.vcf --select-type SNP -O $OUTDIR/raw_snps.vcf
gatk SelectVariants -R $REF -V $OUTDIR/raw_variants.vcf --select-type INDEL -O $OUTDIR/raw_indels.vcf

# STEP 6: Filter SNPs and INDELs
gatk VariantFiltration -R $REF -V $OUTDIR/raw_snps.vcf -O $OUTDIR/filtered_snps.vcf \
  -filter-name "QD_filter" -filter "QD < 2.0" \
  -filter-name "FS_filter" -filter "FS > 60.0" \
  -filter-name "MQ_filter" -filter "MQ < 40.0"

gatk VariantFiltration -R $REF -V $OUTDIR/raw_indels.vcf -O $OUTDIR/filtered_indels.vcf \
  -filter-name "QD_filter" -filter "QD < 2.0" \
  -filter-name "FS_filter" -filter "FS > 200.0"

# STEP 7: Select passing variants
gatk SelectVariants --exclude-filtered -V $OUTDIR/filtered_snps.vcf -O $OUTDIR/final_snps.vcf
gatk SelectVariants --exclude-filtered -V $OUTDIR/filtered_indels.vcf -O $OUTDIR/final_indels.vcf

# STEP 8: Annotate Variants
gatk Funcotator --variant $OUTDIR/final_snps.vcf --reference $REF --ref-version hg38 \
  --data-sources-path $FUNCDIR --output $OUTDIR/annotated_snps.vcf --output-file-format VCF

gatk Funcotator --variant $OUTDIR/final_indels.vcf --reference $REF --ref-version hg38 \
  --data-sources-path $FUNCDIR --output $OUTDIR/annotated_indels.vcf --output-file-format VCF

# STEP 9: Extract Table
gatk VariantsToTable -V $OUTDIR/annotated_snps.vcf -F AC -F AN -F DP -F AF -F FUNCOTATION -O $OUTDIR/final_output_snps.table

echo "✅ Pipeline Complete: Output -> $OUTDIR/final_output_snps.table"
