# FASTQ data from mRNA sequencing - from FASTQ to raw reads
## Resources and information
* GENCODE - Human Release 45 (gencodegenes.org) https://www.gencodegenes.org/human/
* STAR software for RNA-seq aligning https://github.com/alexdobin/STAR, Reference: STAR: ultrafast universal RNA-seq aligner - [PMID:23104886](https://pubmed.ncbi.nlm.nih.gov/23104886/)

## Setup computing environment
STAR runs in Linux environment. For every set of paired-end sequencing data, it's best to have about 4 GB of storage for the raw data. You'll also need around 30 GB for the indexed reference genome. The space for storing STAR alignment results depends on how you set up the alignment process. For instance, if you save the BAM files, you'll need additional space about equal to the raw data size.
For running STAR effectively, you should have at least 32 GB of RAM, but I recommend at least 64 GB for smoother operation. More CPU cores mean faster processing. 
In summary, as an example, if you're handling 100 samples, you need about 1 TB of hard drive space, 64 GB of RAM, and as many CPU cores as you can get for the best performance.

1. Setup in Maastricht DSRI environment: Register account at DSRI server: https://dsri.maastrichtuniversity.nl/register/ and use "ubuntu with webinterface" as tool. 
2. Setup using mobaXterm: mobaXterm is a software that connects PC/laptop with a server. It has a free version that is sufficient for analysis. You need to create an SSH session, have an IP address and username/password for access to a Linux server

## Install STAR
Get latest STAR source from releases, a manual with actual installation instruction can be found there: https://github.com/alexdobin/STAR
