# Codon base editing validation

Scripts to process the validation NGS data for the codon base editing target sites.

Main scripts including two steps:
1. `step1_samstack.sh`: generated stacked sequences files from the mapped SAM file.
2. `step2_visualization.sh`: generate Seqlogo graphs for the validation sites.

Running the scripts requires:
* [BWA](http://bio-bwa.sourceforge.net/): align FASTQ files.
* [samtools](https://www.htslib.org/): process SAM/BAM files.
* [cutadapt](https://cutadapt.readthedocs.io/en/stable/): remove adapters from FASTQ files.
* [weblogo](http://weblogo.threeplusone.com/): generate sequence logo graphs for the sites.
* python with pandas and biopython packages.
* R with ggplot2, tidyr, and grid packages.

One info.txt files should be prepared for the analysis with following columns:
1. number: id of data
2. type: experiment or control
3. gene: gene
4. lysine: number of lysine in the protein to be edited
5. sgRNA: sgRNA sequence
6. a.first.ref: first A position in the reference file
7. R1: file path of FASTQ R1 file
8. R2: file path of FASTQ R2 file

Reference FASTA files together with the BWA indexes should also be provided.
