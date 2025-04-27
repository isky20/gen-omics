# 🧬 Germline Variant Calling, QC, Annotation, and Disease-Gene Analysis Pipeline

This pipeline follows the GATK Best Practices with added QC, annotation, and preparation for disease-gene analysis.

---

## 📂 Input Requirements
- FASTQ files (paired-end)
- Reference genome (e.g., hg38)
- Known sites VCF (e.g., dbSNP)
- Funcotator data sources for annotation

---

## 🔹 Pipeline Steps

### 1. Align Reads
- **Tool**: BWA
- **Action**: Index reference and align reads to the genome.
  
### 2. Sort and Mark Duplicates
- **Tool**: GATK MarkDuplicatesSpark
- **Action**: Sort BAM files and mark PCR duplicates.

### 3. Base Quality Score Recalibration (BQSR)
- **Tool**: GATK BaseRecalibrator, ApplyBQSR
- **Action**: Recalibrate base quality scores using known variant sites.

### 4. Variant Calling
- **Tool**: GATK HaplotypeCaller
- **Action**: Call variants and generate a raw VCF.

### 5. VCF Quality Control (QC)
- **Tool**: bcftools stats, plot-vcfstats
- **Action**: Generate variant statistics and visualization plots.

### 6. Variant Separation
- **Tool**: GATK SelectVariants
- **Action**: Separate SNPs and INDELs into different VCF files.

### 7. Variant Filtering
- **Tool**: GATK VariantFiltration
- **Action**: Apply hard filters to SNPs and INDELs based on quality metrics.

### 8. Select High-Confidence Variants
- **Tool**: GATK SelectVariants
- **Action**: Keep only variants that pass filters.

### 9. Variant Annotation
- **Tool**: GATK Funcotator
- **Action**: Annotate variants with functional information.

### 10. Extract Key Variant Information
- **Tool**: GATK VariantsToTable
- **Action**: Extract useful fields like gene names, allele frequencies, and functional annotation.

### 11. Disease-Gene Analysis Preparation
- **Tool**: Custom script / awk
- **Action**: Extract candidate gene list for disease association matching (OMIM, ClinVar, DisGeNET).

---

## 📈 Output
- BAM files after BQSR
- Raw and filtered VCF files
- QC statistics and plots
- Annotated VCF files
- Final table with important variant annotations
- Candidate gene list for downstream disease analysis

---

## ⚡ Example Commands
See the full script inside [`pipeline.sh`](./pipeline.sh) for detailed command usage.

---

## 📚 References
- [GATK Best Practices Workflows](https://gatk.broadinstitute.org/)
- [BWA Manual](http://bio-bwa.sourceforge.net/)
- [bcftools Documentation](http://samtools.github.io/bcftools/)
- [Funcotator User Guide](https://gatk.broadinstitute.org/hc/en-us/articles/360035531132-Funcotator)

---
