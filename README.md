# MNBC

The MNBC (multithreaded Minimizer-based Naive Bayes Classifier) read classifier

*********************************************************************************************************  
## Install Java
Please download and install Java JDK (version >= 17.0.4) from https://www.oracle.com/ca-en/java/technologies/downloads. For example, use the following Linux command to download the 'jdk-17_linux-x64_bin.tar.gz' file.  
````
wget https://download.oracle.com/java/17/archive/jdk-17.0.12_linux-x64_bin.tar.gz
````
Decompress it with the following command.  
````
tar -xzvf jdk-17_linux-x64_bin.tar.gz
````
Then a new folder 'jdk-17.0.12' appears.<br/>

(Alternatively, you can use the command "<b>mamba create -n java -c conda-forge openjdk</b>" to install Java JDK in mamba/conda, then run the command "<b>mamba activate java</b>" to activate the Java environment)  
## Install MNBC
Please download this repository using the following command, then a new folder 'MNBC' appears.  
````
git clone https://github.com/ComputationalPathogens/MNBC
````
Change to the 'example' subfolder using the following command.  
````
cd MNBC/example
````
Decompress the file 'taxdmp.zip' using the following command, then the file 'nodes.dmp' appears.  
````
unzip taxdmp.zip
````
*********************************************************************************************************  

## Run MNBC on the example
The 'example' folder includes a small demo, which is described below to demonstrate how to use the tool. All input files used and output files produced in the following 3 steps are included in this folder, so the commands can be directly run in a terminal window.  

<b>Problem description:</b>  
The 'reads.fasta' file contains 15 short-read sequences to be classified. The 5 reads whose headers start with SRR227300 were sequenced from the E. coli O104:H4 strain. The 5 reads whose headers start with SRR032501 were sequenced from the Yersinia rohdei ATCC_43380 strain. The 5 reads whose headers start with SRR095845 were sequenced from the Arabidopsis thaliana WalhaesB4 strain. The reference database (i.e. the 'RefSeq_genomes' folder) contains two complete genomes obtained from RefSeq: GCF_022869985.1 belongs to the E. coli O104:H4 strain, and GCF_000834455.1 belongs to the Yersinia rohdei YRA strain.  

<b>Tool usage (3 steps):</b>  
Change to the 'MNBC' folder using the following command.  
````
cd ..
````

<b>(Optional) add Java to Linux PATH variable:</b>  
Replace <parent_directory> in the following line with the actual path, then append this line to the end of the '.bashrc' file in your home directory.  
````
export PATH=<parent_directory>/jdk-17.0.12/bin:${PATH}
````
Either open a new terminal window, or run the following command in the current window.
````
source ~/.bashrc
````
Then the commands in the following 3 steps can be simplified to "<b>java -cp MNBC.jar -Xmx1G MNBC ...</b>".

(If you installed Java JDK in mamba/conda, then the 3 commands can be also simplified to "<b>java -cp MNBC.jar -Xmx1G MNBC ...</b>")  

<b>Step 1</b>:  
Run the following command to generate the taxonomy file of the reference database:  
````
../jdk-17.0.12/bin/java -cp MNBC.jar -Xmx1G MNBC taxonomy -i example/RefSeq_genomes/ -a example/assembly_summary_refseq.txt -n example/nodes.dmp -o example/taxonomy.txt
````
(The following help menu displays by using ```-h```)  
```-a```:	Assembly summary file downloaded from NCBI (e.g. assembly_summary_refseq.txt downloaded from https://ftp.ncbi.nlm.nih.gov/genomes/refseq/))  
```-n```:	Taxonomy nodes.dmp file downoaded from NCBI (the file 'taxdmp.zip' is downloaded from https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/)  
```-i```:	Input directory containing the (gzipped) files of reference sequences in the database (e.g. GCF_000834455.1_ASM83445v1_genomic.fna.gz is a reference genome sequence file downloaded from RefSeq)  
```-o```:	Output taxonomy file for the database

<b>Step 2</b>:  
Run the following command to build the database:  
````
../jdk-17.0.12/bin/java -cp MNBC.jar -Xmx1G MNBC build -k 15 -c 2 -f 300000 -i example/RefSeq_genomes/ -o example/db/
````
(The following help menu displays by using ```-h```)  
```-k```:	K-mer length  
```-c```:	Number of threads  
```-i```:	Input directory containing the (gzipped) files of reference sequences (e.g. GCF_000834455.1_ASM83445v1_genomic.fna.gz is a reference genome sequence file downloaded from RefSeq)  
```-o```: Existing output database directory (please first make this directory if it doesn't already exist)  
```-f (optional)```: Filtering threshold on the sequence length (an integer >= 0). Chromosomes with lengths below this threshold are ignored as well as all plasmids. The default value is 0 (i.e. all chromosomes are retained).  
```-b (optional)```: Log file of the previous prematurely killed run (i.e. .out file in Slurm). This allows breakpoint resumption after the previous run exits abnormally.

<b>Step 3</b>:  
Run the following command to classify the reads against the database:  
````
../jdk-17.0.12/bin/java -cp MNBC.jar -Xmx1G MNBC classify -k 15 -c 2 -d example/db/ -m example/taxonomy.txt -o example/result.txt -t 1 example/reads.fasta
````
(The following help menu displays by using ```-h```)  
```-k```: K-mer length  
```-c```: Number of threads  
```-d```: Input database directory  
```-m```:	Input taxonomy file  
```-o```:	Output classification file  
```-t```:	Type of reads (paired-end: 2, single-end: 1). Paired-end reads have two following (gzipped) .fasta/.fastq files. Single-end reads have one following (gzipped) .fasta/.fastq file.  
```-u (optional)```: Filtering threshold on the ratio of common read-genome minimizers over all read minimizers (default 0.35). When set to 0, only reference genomes sharing minimizers with the read are considered, so that the read will be left unclassified if no genome has common minimizers with it.  
```-p (optional)```: Penalty for absent minimizers (default -2000)  
```-e (optional)```: Threshold on the difference between adjacent scores (default 1500)

<b>Tips</b>:
When using a large reference database (e.g. one containing all prokaryotic and viral complete genomes with OK taxonomy check status from RefSeq), increase the memory amount that MNBC can use in Steps 2 and 3 by adjusting the '-Xmx' parameter (e.g. -Xmx200G), and also increase the number of CPU cores by adjusting the '-c' parameter (e.g. -c 100) to accelerate.

## Format of the taxonomy file
In the tab-delimited taxonomy file 'taxonomy.txt' produced in Step 1, the 1st row contains column headers, and each subsequent row gives the taxonomy information for a reference genome in the database:
````
Assembly ID	taxid.species	taxid.genus	taxid.family	taxid.order	taxid.class	taxid.phylum	taxid.superkingdom	Organism name
GCF_000834455.1	29485	629	1903411	91347	1236	1224	2	Yersinia rohdei strain=YRA
GCF_022869985.1	562	561	543	91347	1236	1224	2	Escherichia coli O104:H4 strain=PartB-EcoliO104H4-RM8376
````
The 1st column is the genome assembly ID, the 2nd to 8th columns are the taxon ID numbers of the genome from the species level to the domain level, and the 9th column is the string name (species and strain) of the genome.

## Format of the classification file
In the tab-delimited classification file 'result.txt' produced in Step 3, the 1st row contains column headers, and each subsequent row gives the classification for a read:  
````
Read	Species	Genus	Family	Order	Class	Phylum	Domain	CandidateGenomes  
SRR227300.1.1	562	561	543	91347	1236	1224	2	GCF_022869985.1  
SRR227300.2.1	562	561	543	91347	1236	1224	2	GCF_022869985.1  
SRR227300.26.1	562	561	543	91347	1236	1224	2	GCF_022869985.1  
SRR227300.27.1	562	561	543	91347	1236	1224	2	GCF_022869985.1  
SRR227300.28.1	562	561	543	91347	1236	1224	2	GCF_022869985.1  
SRR032501.1.2	29485	629	1903411	91347	1236	1224	2	GCF_000834455.1  
SRR032501.2.2	29485	629	1903411	91347	1236	1224	2	GCF_000834455.1  
SRR032501.4.2	29485	629	1903411	91347	1236	1224	2	GCF_000834455.1  
SRR032501.3.2	29485	629	1903411	91347	1236	1224	2	GCF_000834455.1  
SRR095845.1744	unclassified  
SRR095845.1746	unclassified  
SRR032501.5.2	29485	629	1903411	91347	1236	1224	2	GCF_000834455.1  
SRR095845.1747	unclassified  
SRR095845.1745	unclassified  
SRR095845.1748	unclassified  
````
The 1st column is the read ID, the 2nd to 8th columns are assigned taxon ID numbers from the species level to the domain level, and the last column is the IDs of all candidate reference genomes determining the predicted species. It can be seen that all 15 reads were correctly classified (i.e. assigned the correct species-level taxon IDs, or labelled as unclassified).