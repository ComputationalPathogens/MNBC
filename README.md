# MNBC

The MNBC (multithreaded Minimizer-based Naive Bayes Classifier) read classifier

*********************************************************************************************************  
<b>Prerequisites:</b>  
1. Please download and install Java JDK (version >= 17.0.4) from https://www.oracle.com/ca-en/java/technologies/downloads.  
For the Linux operating system, using the command "<b>wget https://download.oracle.com/java/17/latest/jdk-17_linux-x64_bin.tar.gz</b>" to download the 'jdk-17_linux-x64_bin.tar.gz' file to your home directory, and decompress it with the command "<b>tar -xzvf jdk-17_linux-x64_bin.tar.gz</b>". Then a new folder 'jdk-17.0.10' appears.  
2. Please download this repository using the command "<b>git clone https://github.com/ComputationalPathogens/MNBC</b>", then a new folder 'MNBC' appears.  
In the 'example' subfolder, decompress the file 'taxdmp.zip' using the command <b>"unzip taxdmp.zip"</b>.   
*********************************************************************************************************  

The 'example' folder includes a small demo, which is described below to demonstrate how to use the tool. All input files used in this demo and output files produced in the following 3 steps are included in this folder, so the 3 commands can be directly run in a terminal window.  

<b>Problem description:</b>  
The 'reads.fasta' file contains ten short-read sequences to be classified. Five reads, whose headers start with SRR227300, were sequenced from the E. coli O104:H4 strain. The other five reads, whose headers start with SRR032501, from the Yersinia rohdei ATCC_43380 strain. The reference database (i.e. the 'RefSeq_genomes' folder) contains two complete genomes obtained from RefSeq: GCF_022869985.1 belongs to the E. coli O104:H4 strain, and GCF_000834455.1 belongs to the Yersinia rohdei YRA strain.  

<b>Tool usage (3 steps):</b>  
(Please first open a terminal window, and change to the folder 'MNBC' by using the command "</b>cd MNBC</b>")  
(Please change the path '\~/jdk-17.0.10/bin/java' in the following 3 commands accordingly, if the folder 'jdk-17.0.10' is in another directory other than your home directory)  
(Alternatively you can first add Java to Linux PATH variable as follows: Add this line "export PATH=~/jdk-17.0.10/bin/java:${PATH}" to the end of the hidden '.bashrc' file in your home directory, run the command "<b>source ~/.bashrc</b>", then the following 3 commands can be simplified to "<b>java -cp MNBC.jar -Xmx1G MNBC ...</b>")  
1. Copy the following command to a terminal window to generate the taxonomy file of the reference database:  
<b>~/jdk-17.0.10/bin/java -cp MNBC.jar -Xmx1G MNBC taxonomy -i example/RefSeq_genomes/ -a example/assembly_summary_refseq.txt -n example/nodes.dmp -o example/taxonomy.txt</b>  
(The following help menu displays by using '-h')  
-a:	Assembly summary file downloaded from NCBI (e.g. assembly_summary_refseq.txt downloaded from https://ftp.ncbi.nlm.nih.gov/genomes/refseq/))  
-n:	Taxonomy nodes.dmp file downoaded from NCBI (The file 'taxdmp.zip' is downloaded from https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/)  
-i:	Input directory containing the (gzipped) files of reference sequences in the database (e.g. GCF_000834455.1_ASM83445v1_genomic.fna.gz is a reference genome sequence file downloaded from RefSeq)  
-o:	Output taxonomy file for the database

2. Copy the following command to a terminal window to build the database:  
<b>~/jdk-17.0.10/bin/java -cp MNBC.jar -Xmx1G MNBC build -k 15 -c 2 -f 300000 -i example/RefSeq_genomes/ -o example/db/</b>  
(The following help menu displays by using '-h')  
-k:	K-mer length  
-c:	Number of threads  
-i:	Input directory containing the (gzipped) files of reference sequences (e.g. GCF_000834455.1_ASM83445v1_genomic.fna.gz is a reference genome sequence file downloaded from RefSeq)  
-o: Existing output database directory (please first make this directory if it doesn't already exist)  
-f (optional): Filtering threshold on the sequence length (an integer >= 0). Chromosomes with lengths below this threshold are ignored as well as all plasmids. The default value is 0 (i.e. all chromosomes are retained).  
-b (optional): Log file of the previous prematurely killed run (i.e. .out file in Slurm). This allows breakpoint resumption after the previous run exits abnormally.

3. Copy the following command to a terminal window to classify the reads against the database:  
<b>~/jdk-17.0.10/bin/java -cp MNBC.jar -Xmx1G MNBC classify -k 15 -c 2 -d example/db/ -m example/taxonomy.txt -o example/result.txt -t 1 example/reads.fasta</b>  
(The following help menu displays by using '-h')  
-k: K-mer length  
-c: Number of threads  
-d: Input database directory  
-m:	Input taxonomy file  
-o:	Output classification file  
-t:	Type of reads (paired-end: 2, single-end: 1). Paired-end reads have two following (gzipped) .fasta/.fastq files. Single-end reads have one following (gzipped) .fasta/.fastq file.  
-p (optional): Penalty for absent minimizers (default -2000)  
-e (optional): Threshold on the difference between adjacent scores (default 1500)

<b>Format of final classification file:</b>  
In the final tab-delimited classification file 'result.txt' produced by the last command, the 1st row contains column headers, and each subsequent row gives the classification for a read.  

Read	Genome	Species	Genus	Family	Order	Class	Phylum	Superkingdom  
SRR227300.1.1	GCF_022869985.1	562	561	543	91347	1236	1224	2  
SRR227300.2.1	GCF_022869985.1	562	561	543	91347	1236	1224	2  
SRR227300.4.1	GCF_022869985.1	562	561	543	91347	1236	1224	2  
SRR227300.3.1	GCF_022869985.1	562	561	543	91347	1236	1224	2  
SRR227300.5.1	GCF_022869985.1	562	561	543	91347	1236	1224	2  
SRR032501.2.2	GCF_000834455.1	29485	629	1903411	91347	1236	1224	2  
SRR032501.1.2	GCF_000834455.1	29485	629	1903411	91347	1236	1224	2  
SRR032501.3.2	GCF_000834455.1	29485	629	1903411	91347	1236	1224	2  
SRR032501.4.2	GCF_000834455.1	29485	629	1903411	91347	1236	1224	2  
SRR032501.5.2	GCF_000834455.1	29485	629	1903411	91347	1236	1224	2  

The 1st column is the read ID, the 2nd column is the genome ID assigned to the read, and the next 7 columns are the assigned taxon ID numbers from the species level to the domain level. It can be seen that all 10 reads were assigned the correct species-level taxon IDs (i.e. E. coli -- 562, Yersinia rohdei -- 29485).