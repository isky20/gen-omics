#!/bin/bash

# ============================
# T1D Genomics Analysis Pipeline (FASTQ to GWAS + PRS + HLA + Fine-mapping + Phasing)
# Author: [Your Name]
# ============================

# Set paths
REF=human_g1k_v37.fasta
KNOWN_SITES=known_sites.vcf
THREADS=8
SAMPLES=(sample1 sample2) # list of sample IDs

# Step 1: Alignment and BAM prep
for SAMPLE in "${SAMPLES[@]}"; do
  echo "Processing $SAMPLE..."

  # Align with BWA-MEM
  bwa mem -t $THREADS -R "@RG\tID:$SAMPLE\tSM:$SAMPLE\tPL:ILLUMINA" $REF \
    ${SAMPLE}_R1.fastq.gz ${SAMPLE}_R2.fastq.gz > ${SAMPLE}.sam

  # Convert, sort, and index
  samtools view -bS ${SAMPLE}.sam | samtools sort -o ${SAMPLE}.sorted.bam
  samtools index ${SAMPLE}.sorted.bam

  # Mark duplicates
  gatk MarkDuplicatesSpark -I ${SAMPLE}.sorted.bam -O ${SAMPLE}.dedup.bam \
    -M ${SAMPLE}.metrics.txt
  samtools index ${SAMPLE}.dedup.bam

  # GVCF calling
  gatk HaplotypeCaller -R $REF -I ${SAMPLE}.dedup.bam -O ${SAMPLE}.g.vcf.gz -ERC GVCF

  # HLA Typing with xHLA (optional)
  echo "Running HLA typing for $SAMPLE..."
  xhla --bwa /path/to/bwa \
       --ref /path/to/hla_reference.fa \
       --bam ${SAMPLE}.dedup.bam \
       --sampleID ${SAMPLE} \
       --out hla_typing/${SAMPLE}_HLA

done

# Step 2: Joint Genotyping
GVCF_LIST=$(ls *.g.vcf.gz | sed 's/^/-V /' | tr '\n' ' ')
gatk CombineGVCFs -R $REF $GVCF_LIST -O cohort.g.vcf.gz
gatk GenotypeGVCFs -R $REF -V cohort.g.vcf.gz -O t1d_raw.vcf.gz

# Step 3: Filtering
bcftools filter -i 'QUAL>30 && DP>10' -Oz -o t1d_filtered.vcf.gz t1d_raw.vcf.gz

# Step 4: Annotation (after filtering, with advanced annotations)
vep -i t1d_filtered.vcf.gz -o t1d_annotated.vcf.gz \
  --vcf --cache --offline --assembly GRCh37 \
  --plugin GnomAD,/path/to/gnomad.genomes.r2.1.1.sites.vcf.bgz \
  --plugin CADD,/path/to/whole_genome_SNVs.tsv.gz,/path/to/InDels.tsv.gz \
  --plugin dbNSFP,/path/to/dbNSFP.gz,Polyphen2_HDIV_pred,SIFT_pred,MutationTaster_pred \
  --custom /path/to/clinvar.vcf.gz,ClinVar,vcf,exact,0,CLNSIG \
  --symbol --fork 4 --stats_text

# Step 5: Annotation-based Filtering
bcftools view -i 'INFO/gnomAD_AF<0.01 && INFO/CADD_PHRED>20 && (INFO/Consequence~"missense_variant" || INFO/Consequence~"stop_gained") && INFO/CLNSIG~"Pathogenic"' \
  t1d_annotated.vcf.gz -Oz -o t1d_filtered_final.vcf.gz

# Step 6: Phasing HLA region (chr6:29M-34M) using SHAPEIT5
shapeit5.dup --input t1d_filtered_final.vcf.gz \
  --map /path/to/genetic_map_chr6.txt \
  --region 6:29000000-34000000 \
  --output t1d_chr6_HLA_phased.bcf \
  --thread $THREADS
bcftools index t1d_chr6_HLA_phased.bcf

# Step 7: Convert final filtered VCF to PLINK format
plink --vcf t1d_filtered_final.vcf.gz --make-bed --out t1d_study

# Step 8: Genotype QC
plink --bfile t1d_study --maf 0.01 --geno 0.02 --mind 0.02 --hwe 1e-6 \
  --make-bed --out t1d_study.QC

# Step 9: PCA for population structure
plink --bfile t1d_study.QC --indep-pairwise 50 5 0.2 --out prune_indices
plink --bfile t1d_study.QC --extract prune_indices.prune.in --pca 10 --out t1d_study.PCA

# Step 10: GWAS Association Test
plink --bfile t1d_study.QC --pheno pheno.txt --pheno-name T1D \
  --covar covar.txt --covar-name Sex,Age,PC1,PC2,PC3,PC4,PC5 \
  --logistic hide-covar --ci 0.95 --out T1D_GWAS_results

# Step 11: Fine-mapping HLA region with FINEMAP
plink --bfile t1d_study.QC --chr 6 --from-bp 29000000 --to-bp 34000000 --make-bed --out region_chr6_hla
plink2 --bfile region_chr6_hla --r square --out region_chr6_ld
finemap --sss --in-files region_chr6_input --log --n-causal-snps 5 --out region_chr6_finemap

# Step 12: PRS calculation
plink --bfile t1d_study.QC --score T1D_base_gwas.txt 1 2 3 header --out PRS_scores

echo "Pipeline complete!"
