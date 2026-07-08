## ATAC-seq Pipeline

* **Alignment Stage:**
* **Input:** reads
* **Output:** alignments
* **Software:** bowtie2-2.2.5, samtools-1.9


* **Peak Quantification Stage:**
* **Input:** alignments
* **Output:** peak quantifications
* **Software:** macs2-2.2.7.1


* **Signal Tracking:**
* Read depth signal generated with bedGraphToBigWig-2.8



---

## RNA-seq Pipeline

* **Framework:** ENCODE DCC rna-seq-pipeline
* **Alignment Stage:**
* **Input:** reads
* **Output:** alignments
* **Software:** STAR_2.5.1b


* **Gene Quantification Stage:**
* **Input:** alignments
* **Output:** gene quantifications
* **Software:** RSEM v1.2.23



---

## ChIP-seq Pipeline (Histone)

* **Alignment Stage:**
* **Input:** reads
* **Output:** alignments
* **Software:** bowtie2-2.2.5, samtools-1.9


* **Peak Quantification Stage:**
* **Input:** alignments
* **Output:** peak quantifications
* **Software:** macs2-2.2.7.1


* **Signal Tracking:**
* Read depth signal generated with bedGraphToBigWig-2.8



---

## TF ChIP-seq Pipeline

* **Alignment Stage:**
* **Input:** reads
* **Output:** alignments
* **Software:** bowtie2-2.2.5, samtools-1.9


* **Peak Quantification Stage:**
* **Input:** alignments
* **Output:** peak quantifications
* **Software:** macs2-2.2.7.1


* **Signal Tracking:**
* Read depth signal generated with bedGraphToBigWig-2.8



---

## Hi-C Pipeline

* **Alignment Stage:**
* **Input:** reads
* **Output:** alignments
* **Software:** Hi-C-Pro-2-11.4, bowtie2-2.4.1


* **Contact Matrix Stage:**
* **Input:** alignments
* **Output:** contact matrix
* **Software:** juicer-tools-1.22.01, samtools-1.9



---

## scRNA-seq Pipeline

* **Count Matrices Stage:**
* **Input:** reads
* **Output:** count matrices
* **Software:** cellranger-3.1.0



---

## FACS CRISPR Screen Pipeline

* **Count Matrices Stage:**
* **Input:** reads
* **Output:** count matrices
* **Software:** MSKCC Custom Script (see FACS CRISPR folder)
