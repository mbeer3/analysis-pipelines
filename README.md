# Analysis Pipelines

This repository contains analysis pipelines for multiple genomic assays.

## ATAC-seq Pipeline

### Alignment

- **Input:** Reads
- **Output:** Alignments
- **Software:** `bowtie2-2.2.5`, `samtools-1.9`

### Peak Quantification

- **Input:** Alignments
- **Output:** Peak quantifications
- **Software:** `macs2-2.2.7.1`

### Signal Tracking

- Read-depth signal generated with `bedGraphToBigWig-2.8`

---

## RNA-seq Pipeline

- **Framework:** ENCODE DCC `rna-seq-pipeline`

### Alignment

- **Input:** Reads
- **Output:** Alignments
- **Software:** `STAR_2.5.1b`

### Gene Quantification

- **Input:** Alignments
- **Output:** Gene quantifications
- **Software:** `RSEM v1.2.23`

---

## ChIP-seq Pipeline (Histone)

### Alignment

- **Input:** Reads
- **Output:** Alignments
- **Software:** `bowtie2-2.2.5`, `samtools-1.9`

### Peak Quantification

- **Input:** Alignments
- **Output:** Peak quantifications
- **Software:** `macs2-2.2.7.1`

### Signal Tracking

- Read-depth signal generated with `bedGraphToBigWig-2.8`

---

## TF ChIP-seq Pipeline

### Alignment

- **Input:** Reads
- **Output:** Alignments
- **Software:** `bowtie2-2.2.5`, `samtools-1.9`

### Peak Quantification

- **Input:** Alignments
- **Output:** Peak quantifications
- **Software:** `macs2-2.2.7.1`

### Signal Tracking

- Read-depth signal generated with `bedGraphToBigWig-2.8`

---

## Hi-C Pipeline

### Alignment

- **Input:** Reads
- **Output:** Alignments
- **Software:** `Hi-C-Pro-2-11.4`, `bowtie2-2.4.1`

### Contact Matrix Generation

- **Input:** Alignments
- **Output:** Contact matrix
- **Software:** `juicer-tools-1.22.01`, `samtools-1.9`

---

## scRNA-seq Pipeline

### Count Matrix Generation

- **Input:** Reads
- **Output:** Count matrices
- **Software:** `cellranger-3.1.0`

---

## FACS CRISPR Screen Pipeline

### Count Matrix Generation

- **Input:** Reads
- **Output:** Count matrices
- **Software:** MSKCC custom script; see the FACS CRISPR folder