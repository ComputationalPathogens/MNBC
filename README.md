# MNBC

The MNBC (minimizer-based Naive Bayes classifier) read binner

*********************************************************************************************************  
<b>Prerequisites:</b>  
1. Please download and install JDK (Java Development Kit) from https://www.oracle.com/ca-en/java/technologies/downloads  
2. Please download the MNBC.jar, eclipse-collections-11.1.0.jar, eclipse-collections-api-11.1.0.jar files from this repository  
*********************************************************************************************************  

The 'example' folder includes a demo, which is described below to demonstrate the usage of the binner:

<b>Problem description:</b>  
The 'reads.fasta' file contains ten short read sequences to be classified. Five reads, whose headers start with SRR227300, were sequenced from the E. coli O104:H4 strain. The other five reads, whose headers start with SRR032501, from the Yersinia rohdei ATCC_43380 strain. The reference database contains two complete genomes obtained from NCBI RefSeq; GCF_022869985.1 belongs to the E. coli O104:H4 strain, and GCF_000834455.1 belongs to the Yersinia rohdei YRA strain. From the result file 'classify_result.txt', it can be seen that all the ten reads were correctly classified.

<b>Tool usage:</b>  
1. Run GenomeTaxidFinder with the following command to find the taxonomic IDs of the reference genomes:  
<b>java -cp MNBC.jar -Xmx1G GenomeTaxidFinder RefSeq_genomes/ assembly_summary_refseq.txt taxdmp/nodes.dmp metainfo.txt</b>  
RefSeq_genomes: Input directory containing reference sequence files downloaded from RefSeq (*.fna.gz)  
assembly_summary_refseq.txt: Downloaded from  RefSeq https://ftp.ncbi.nlm.nih.gov/genomes/.vol2/refseq/  
nodes.dmp: Taxonomy file in the taxdmp folder (taxdmp.zip downloaded from https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/)  
metainfo.txt: Output file

2. Run NBCBuild with the following command to build the reference database:  
<b>java -cp MNBC.jar:eclipse-collections-api-11.1.0.jar:eclipse-collections-11.1.0.jar -Xmx1G NBCBuild -k 15 -c 2 -f 300000 -i RefSeq_genomes/ -o db/</b>  
(The following help menu will be displayed by using '-h')  
-k:	K-mer length  
-c:	Number of threads  
-i:	Input directory containing the reference sequence files downloaded from  RefSeq (*.fna.gz or *.fna)
-o: Existing output directory which is the reference database
-f (optional): The minimum threshold on the sequence length (an integer not smaller than 0, the default value is 0). Sequences with lengths below this threshold will be removed, and all plasmids will be removed.
-b (optional): Log file of the previous abnormally killed run (.out file in Slurm)

3. Run NBCClassify with the following command to classify the reads against the database:  
<b>java -cp MNBC.jar:eclipse-collections-api-11.1.0.jar:eclipse-collections-11.1.0.jar -Xmx2G NBCClassify -k 15 -c 2 -d db/ -m metainfo.txt -o result.txt -t 1 reads.fasta</b>  
(The following help menu will be displayed by using '-h')  
-k: K-mer length  
-c: Number of threads  
-d: Input database directory  
-m:	Input metainfo file  
-o:	Final classification file  
-t:	Type of reads (Paired-end: 2, Single-end: 1). Paired-end reads have two .fasta/.fastq files (can be gzipped) following; single-end reads have one (can be gzipped).  
-p (optional): Penalty for non-existent k-mers (default -2000)
