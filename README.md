ATAC-seq pipeline:
	 input: reads	output: alignments	software:	bowtie2-2.2.5, samtools-1.9
	 input: alignments	output: peak quantifications	software:	macs2-2.2.7.1

RNA-seq pipeline:
ENCODE DCC rna-seq-pipeline
       input: reads	output: alignments	software:	STAR_2.5.1b
       input: alignments	output: gene quantifications	software:	RSEM v1.2.23

ChIP-seq pipeline (histone):
	 input: reads	output: alignments	software:	bowtie2-2.2.5, samtools-1.9
	 input: alignments	output: peak quantifications	software:	macs2-2.2.7.1

TF ChIP-seq pipeline:
	 input: reads	output: alignments	software:	bowtie2-2.2.5, samtools-1.9
	 input: alignments	output: peak quantifications	software:	macs2-2.2.7.1

Hi-C pipeline:

scRNA-seq pipeline:
	  input: reads	output: count matrices	software:	cellranger-3.1.0

FACS CRISPR screen pipeline:
	  input: reads	output: count matrices	software:	cellranger-3.1.0
	