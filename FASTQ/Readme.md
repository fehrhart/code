# FASTQ data from mRNA sequencing - from FASTQ to raw reads
## Setup computing environment
STAR runs in Linux environment. For every set of paired-end sequencing data, it's best to have about 4 GB of storage for the raw data. You'll also need around 30 GB for the indexed reference genome. The space for storing STAR alignment results depends on how you set up the alignment process. For instance, if you save the BAM files, you'll need additional space about equal to the raw data size.
For running STAR effectively, you should have at least 32 GB of RAM, but I recommend at least 64 GB for smoother operation. More CPU cores mean faster processing. 
In summary, as an example, if you're handling 100 samples, you need about 1 TB of hard drive space, 64 GB of RAM, and as many CPU cores as you can get for the best performance.
Small samples you can run on your laptop, but for larger experiments a good computation server is needed. 

1. Setup in Maastricht DSRI environment: Register account at DSRI server: https://dsri.maastrichtuniversity.nl/register/ and use "ubuntu with webinterface" as tool. 
2. Setup using mobaXterm: mobaXterm is a software that connects PC/laptop with a server. It has a free version that is sufficient for analysis. You need to create an SSH session, have an IP address and username/password for access to a Linux server

## Analysis steps
### Upload data
Create a folder on the server for your data and upload your data there. From SURF e.g. with a link to download: wget https://surfdrive.surf.nl/files/index.php/s/DATALINK/download and unzip. 
Also, create an output folder for the results.

### Trimming (not provided)
Trimming is removing adapter sequences. This step is usually already done by the sequencing centre. If needed, the software for this is "trimmomatic".

### Install STAR
STAR software is for for RNA-seq aligning. Get latest STAR source from releases, a manual with actual installation instruction can be found there: https://github.com/alexdobin/STAR, Reference: STAR: ultrafast universal RNA-seq aligner - [PMID:23104886](https://pubmed.ncbi.nlm.nih.gov/23104886/)

### Download Genecode
GENCODE - Human Release 45 (gencodegenes.org) https://www.gencodegenes.org/human/

### Genome Indexing
Genome indexing is needed to get the count matrix. The software prepares genome references from fasta and annotation files. https://github.com/fehrhart/code/blob/master/FASTQ/STAR_GenomeIndexing.sh
Variables:
* SBATCH = these are server specifications, check with local server admins.
* runThreadN 16 = number of CPUs on server. 16 was used for a standalone linux server, on Maastricht DSRI up to 128 can be used.
* runMode genomeGenerate = runs the genome indexing
* genomeDir = path to output folder on the server (make sure to create an output folder before!)
* genomeFastaFiler = path to fasta files folder
* sjdbGTFfile = path to annotation folder
* sjdbOverhang 150 = depends on sequencing technique: mRNA generally average size of transcripts 150 bp
The result is an index file in the output folder.

### Counting
https://github.com/fehrhart/code/blob/master/FASTQ/STAR_count.sh 
This will read the files, first $R1 then $R2 and then combine the results - check the names of the files accordingly.
Variables:
* quantMode GeneCounts = count gene mode
* genomeDir = output of the previous indexing step
* outSAMtype None = unless you do want the SAM file (remember to have extra space)
* outFileNamePrefix = takes the file name and adds out and puts it in output directory
* readFilesCommand zcat = tells STAR the fastq input files are in gz format
* other variables - for paired end mRNA they are ok
* outTmDir = needs a temporary output file, not needed, STAR can generate one itself but better to define one if there are big files, temp will by default safe on home directory, which may be too small
When done, the output is SampleNameReadsPerGene.out.tab. The first 4 lines about all genes in sample. The first result column is unstranded Rnaseq, the others are R1 and R2. Column 2 is result.

### Log file 
Gives information on the quality of sequencing. Uniquely mapped reads over 60% are acceptable, the higher the better. % of reads unmapped should not be higher than 40/50%.
