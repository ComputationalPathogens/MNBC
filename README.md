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

<b>Binner usage:</b>  
1. Run GenomeTaxidFinder with the following command to find the taxonomic IDs of the reference genomes:  
<b>java -cp MA-NBC.jar -Xmx2G GenomeTaxidFinder RefSeq_genomes/ assembly_summary_refseq.txt taxdmp/nodes.dmp metainfo.txt<b>  
RefSeq_genomes: Input directory corresponding to the reference database, containing sequence files downloaded from RefSeq (*.fna.gz)  
assembly_summary_refseq.txt: Downloaded from  RefSeq https://ftp.ncbi.nlm.nih.gov/genomes/.vol2/refseq/  
nodes.dmp: Taxonomy file in the taxdmp folder (taxdmp.zip downloaded from https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/)  
metainfo.txt: Output file
