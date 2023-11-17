# MNBC

The MNBC (minimizer-based Naive Bayes classifier) read binner

*********************************************************************************************************  
<b>Prerequisites:</b>  
1. Please download and install JDK (Java Development Kit) from https://www.oracle.com/ca-en/java/technologies/downloads  
2. Please download the MNBC.jar, eclipse-collections-11.1.0.jar, eclipse-collections-api-11.1.0.jar files from this repository  
*********************************************************************************************************  

The 'example' folder includes a demo, which is described below to demonstrate how to use the binner:

<b>Problem description:</b>  
The 'reads.fasta' file contains ten short read sequences to be classified. Five reads, whose headers start with SRR227300, were sequenced from the E. coli O104:H4 strain. The other five reads, whose headers start with SRR032501, from the Yersinia rohdei ATCC_43380 strain. The reference database contains two complete genomes obtained from NCBI RefSeq; GCF_022869985.1 belongs to the E. coli O104:H4 strain, and GCF_000834455.1 belongs to the Yersinia rohdei YRA strain. From the result file 'classify_result.txt', it can be seen that all the ten reads were correctly classified.

<b>Tool usage:</b>  
1. Run MNBC_taxonomy with the following command to find the taxonomic IDs of the reference genomes:  
<b>java -cp MNBC.jar -Xmx1G MNBC_taxonomy -i RefSeq_genomes/ -a assembly_summary_refseq.txt -n taxdmp/nodes.dmp -o taxonomy.txt</b>  
(The following help menu will be displayed by using '-h')  
-a:	Assembly summary file downloaded from NCBI (e.g. assembly_summary_refseq.txt from https://ftp.ncbi.nlm.nih.gov/genomes/refseq/))  
-n:	Taxonomy nodes.dmp file downoaded from NCBI (e.g. taxdmp.zip from https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/)  
-i:	Input directory containing the (gzipped) files of reference sequences in the database (e.g. *.fna.gz downloaded from RefSeq)  
-o:	Output taxonomy file for the database

2. Run MNBC_build with the following command to build the reference database:  
<b>java -cp MNBC.jar:eclipse-collections-api-11.1.0.jar:eclipse-collections-11.1.0.jar -Xmx1G MNBC_build -k 15 -c 2 -f 300000 -i RefSeq_genomes/ -o db/</b>  
(The following help menu will be displayed by using '-h')  
-k:	K-mer length  
-c:	Number of threads  
-i:	Input directory containing the (gzipped) files of reference sequences (e.g. *.fna.gz downloaded from RefSeq)  
-o: Existing Output database directory  
-f (optional): The minimum filtering threshold on the sequence length (an integer not smaller than 0, the default value is 0). Sequences with lengths below this threshold will be ignored, and all plasmids will be ignored.  
-b (optional): Log file of the previous prematurely killed run (i.e. Slurm .out file)

3. Run MNBC_classify with the following command to classify the reads against the database:  
<b>java -cp MNBC.jar:eclipse-collections-api-11.1.0.jar:eclipse-collections-11.1.0.jar -Xmx1G MNBC_classify -k 15 -c 2 -d db/ -m metainfo.txt -o result.txt -t 1 reads.fasta</b>  
(The following help menu will be displayed by using '-h')  
-k: K-mer length  
-c: Number of threads  
-d: Input database directory  
-m:	Input taxonomy file  
-o:	Output classification file  
-t:	Type of reads (paired-end: 2, single-end: 1). Paired-end reads have two following (gzipped) .fasta/.fastq files. Single-end reads have one following (gzipped) .fasta/.fastq file.  
-p (optional): Penalty for non-existent k-mers (default -2000)
