# 🧬 Variant Discovery Pipeline

- **FastQs**:  
  Raw sequencing reads obtained from the sequencer.  
  **Aim**: Prepare data for quality control and alignment.  
  **Tools**: None (this is your starting point — FASTQ files).

---

- **QC (FastQs)**:  
  Assess quality of sequencing reads (e.g., low-quality bases, adapter contamination).  
  **Aim**: Ensure only good-quality reads are used for downstream analysis.  
  **Tools**: FastQC, MultiQC, Trimmomatic (for trimming), fastp.

---

- **Filtering (Reads)**:  
  Remove low-quality reads, short reads, or reads with too many unknown bases (Ns).  
  **Aim**: Clean reads to improve mapping accuracy.  
  **Tools**: Trimmomatic, fastp.

---

- **Alignment**:  
  Map filtered reads to a reference genome.  
  **Aim**: Generate BAM files showing where each read aligns on the genome.  
  **Tools**: BWA, Bowtie2.

---

- **Post-processing**:  
  Improve the quality of BAM files (e.g., sort, mark duplicates, recalibrate quality scores).  
  **Aim**: Prepare BAM files for accurate variant calling.  
  **Tools**: SAMtools, Picard, GATK BaseRecalibrator.

---

- **Variant Calling (gVCF, VCF)**:  
  Identify SNPs and INDELs from aligned reads.  
  **Aim**: Generate a list of variants for each sample.  
  **Tools**: GATK HaplotypeCaller, FreeBayes, DeepVariant, bcftools call.

---

- **QC (VCF)**:  
  Check quality of called variants (e.g., Ti/Tv ratio, depth, genotype quality).  
  **Aim**: Ensure variants are reliable and biologically meaningful.  
  **Tools**: bcftools stats, vcftools, plot-vcfstats.

---

- **Filtering (Variants)**:  
  Remove low-confidence variants based on quality metrics (e.g., depth, GQ, filter status).  
  **Aim**: Retain only high-quality variant calls.  
  **Tools**: bcftools filter, GATK VariantFiltration.

---

- **Annotation**:  
  Add biological meaning (e.g., what genes are affected, predicted variant effect, disease association).  
  **Aim**: Interpret the variants biologically and clinically.  
  **Tools**: VEP (Variant Effect Predictor), ANNOVAR, SnpEff.

---
