# MNBC

The MNBC (multithreaded Minimizer-based Naive Bayes Classifier) read classifier

*********************************************************************************************************  
<b>Prerequisites:</b>  
1. Please download and install Java JDK (version >= 17.0.4) from https://www.oracle.com/ca-en/java/technologies/downloads.  
For example, for the Linux operating system, download the 'jdk-17_linux-x64_bin.tar.gz' file (https://download.oracle.com/java/17/latest/jdk-17_linux-x64_bin.tar.gz) to the '/home' directory, and decompress it with the command <b>tar -xzvf jdk-17_linux-x64_bin.tar.gz</b>. Then a new folder 'jdk-17.0.10' appears in the '/home' directory.  
2. Please download MNBC.jar and the example folder from this repository  
*********************************************************************************************************  

The 'example' folder includes a small demo, which is described below to demonstrate how to use the tool. All input files used in this demo and output files produced by the following commands are included in the folder.  

<b>Problem description:</b>  
The 'reads.fasta' file contains ten short-read sequences to be classified. Five reads, whose headers start with SRR227300, were sequenced from the E. coli O104:H4 strain. The other five reads, whose headers start with SRR032501, from the Yersinia rohdei ATCC_43380 strain. The reference database contains two complete genomes obtained from RefSeq: GCF_022869985.1 belongs to the E. coli O104:H4 strain, and GCF_000834455.1 belongs to the Yersinia rohdei YRA strain. From the result file 'result.txt', it can be seen that all ten reads were correctly classified.

<b>Tool usage:</b>  
(Please change the following path '/home/jdk-17.0.10/bin/java' accordingly if the folder 'jdk-17.0.10' is in another directory other than '/home')  
1. Run the following command in a terminal window to generate the taxonomy file of the reference database:  
<b>/home/jdk-17.0.10/bin/java -cp MNBC.jar -Xmx1G MNBC taxonomy -i RefSeq_genomes/ -a assembly_summary_refseq.txt -n taxdmp/nodes.dmp -o taxonomy.txt</b>  
(The following help menu displays by using '-h')  
-a:	Assembly summary file downloaded from NCBI (e.g. assembly_summary_refseq.txt from https://ftp.ncbi.nlm.nih.gov/genomes/refseq/))  
-n:	Taxonomy nodes.dmp file downoaded from NCBI (Please decompress the file 'taxdmp.zip', then the folder 'taxdmp' appears. The file 'taxdmp.zip' is downloaded from https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/)  
-i:	Input directory containing the (gzipped) files of reference sequences in the database (e.g. GCF_000009045.1_ASM904v1_genomic.fna.gz is a reference genome sequence file downloaded from RefSeq)  
-o:	Output taxonomy file for the database

2. Run the following command in a terminal window to build the database:  
<b>/home/jdk-17.0.10/bin/java -cp MNBC.jar -Xmx1G MNBC build -k 15 -c 2 -f 300000 -i RefSeq_genomes/ -o db/</b>  
(The following help menu displays by using '-h')  
-k:	K-mer length  
-c:	Number of threads  
-i:	Input directory containing the (gzipped) files of reference sequences (e.g. GCF_000009045.1_ASM904v1_genomic.fna.gz is a reference genome sequence file downloaded from RefSeq)  
-o: Existing output database directory (please first make this directory if it doesn't already exist)  
-f (optional): Filtering threshold on the sequence length (an integer >= 0). Chromosomes with lengths below this threshold are ignored as well as all plasmids. The default value is 0 (i.e. all chromosomes are retained).  
-b (optional): Log file of the previous prematurely killed run (i.e. .out file in Slurm). This allows breakpoint resumption after the previous run exits abnormally.

3. Run the following command in a terminal window to classify the reads against the database:  
<b>/home/jdk-17.0.10/bin/java -cp MNBC.jar -Xmx1G MNBC classify -k 15 -c 2 -d db/ -m taxonomy.txt -o result.txt -t 1 reads.fasta</b>  
(The following help menu displays by using '-h')  
-k: K-mer length  
-c: Number of threads  
-d: Input database directory  
-m:	Input taxonomy file  
-o:	Output classification file  
-t:	Type of reads (paired-end: 2, single-end: 1). Paired-end reads have two following (gzipped) .fasta/.fastq files. Single-end reads have one following (gzipped) .fasta/.fastq file.  
-p (optional): Penalty for absent minimizers (default -2000)  
-e (optional): Threshold on the difference between adjacent scores (default 1500)
