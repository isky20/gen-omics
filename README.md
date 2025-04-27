# 🧬 Germline Variant Calling and Annotation Pipeline

This pipeline is a simplified germline variant discovery and annotation workflow based on GATK4 Best Practices.

---

## 📂 Input Requirements
- FASTQ files (paired-end)
- Reference genome (e.g., hg38.fa)
- Known sites VCF for BQSR (e.g., dbSNP v138)
- Funcotator data sources directory for annotation

---

## 🔹 Pipeline Steps

### 1. Align Reads
- **Tool**: BWA
- **Command**:  
  - Index the reference genome.
  - Align paired-end reads to the reference genome using `bwa mem`.

### 2. Sort and Mark Duplicates
- **Tool**: GATK MarkDuplicatesSpark
- **Command**:  
  - Convert SAM to BAM, sort reads, and mark duplicates.

### 3. Base Quality Score Recalibration (BQSR)
- **Tool**: GATK BaseRecalibrator and ApplyBQSR
- **Command**:  
  - Create a recalibration table using known sites.
  - Apply BQSR to the BAM file.

### 4. Variant Calling
- **Tool**: GATK HaplotypeCaller
- **Command**:  
  - Call germline variants (VCF output).

### 5. Select SNPs and INDELs
- **Tool**: GATK SelectVariants
- **Command**:  
  - Separate SNPs and INDELs from the raw variant calls into two VCF files.

### 6. Filter SNPs and INDELs
- **Tool**: GATK VariantFiltration
- **Command**:  
  - Apply hard filters to SNPs and INDELs based on QD, FS, and MQ metrics.

### 7. Select Passing Variants
- **Tool**: GATK SelectVariants
- **Command**:  
  - Select only variants that pass the filters (`PASS` variants).

### 8. Annotate Variants
- **Tool**: GATK Funcotator
- **Command**:  
  - Annotate both SNPs and INDELs using Funcotator with a specified data source.

### 9. Extract Table
- **Tool**: GATK VariantsToTable
- **Command**:  
  - Extract important fields (e.g., AC, AN, DP, AF, FUNCOTATION) from annotated VCFs into a final table.

---

## 📈 Output
- Sorted and duplicate-marked BAM files
- Raw VCF files (variants called)
- Filtered SNP and INDEL VCF files
- Annotated SNP and INDEL VCF files
- Final variant annotation table (`final_output_snps.table`)

---

## 📚 References
- [GATK Best Practices](https://gatk.broadinstitute.org/)
- [BWA Manual](http://bio-bwa.sourceforge.net/)
- [Funcotator Documentation](https://gatk.broadinstitute.org/hc/en-us/articles/360035531132-Funcotator)

---
