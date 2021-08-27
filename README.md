# Analysis of rRNA data from BSF Metatrascriptomic Study
This readme illustrates the workflow and pipeline for the analysis of rRNA data from BSF metatrascriptomic.
## Introduction
The Black Soldier Fly (Hermetia illucens; BSF) is a useful tool in valorising organic biomass and other biodegradable wastes. In this study, the BSF larvae were bred under different diets selected based on increasing lignocellulose content. These diets were: processed chicken feed (CF), chicken manure (CM), Brewer’s spent grain (BSG), and Water Hyacinth (WH). An additional diet Feed Mix (FM), consisting of the four diets in equal proportions was also incorporated. The different metatranscriptomes were sequenced using the PCR-cDNA approach on the ONT MinION platform. While the work, using ONT, _aimed to identify and functionally characterise lignocellulosic biomass-degrading microbes_, the mRNA enrichment protocol still retained some rRNAs, which were filtered out using SortMeRNA (_Kopylova et al., 2012_).This is an explorative study. We aim to explore what we can glean from the data.
# Getting started
## Creating and activating conda environment
First a working conda environment was created.
````
conda create transcriptomics
conda activate transcriptomics
````
## The pipeline.
![image](https://i.imgur.com/Syr0V2T.png)
## Workflow
### Step 1. Analysing Sequence Quality with FastQC
"FastQC aims to provide a simple way to do some quality control checks on raw sequence data coming from high throughput sequencing pipelines. It provides a modular set of analyses which you can use to give a quick impression of whether your data has any problems of which you should be aware before doing any further analysis."

The first step before processing any samples is to analyze the quality of the data. Within the `fastq` file is quality information that refers to the accuracy (% confidence) of each base call. **FastQC** looks at different aspects of the sample sequences to determine any irregularies or features that make affect your results (adapter contamination, sequence duplication levels, etc.)
#### command
````
#This script loads fastqc from hpc.
#Through a loop it runs the fastqc trough the sample and redirects the output to a directory named as fastqc.
module load fastqc/0.11.9

for file in ./data/*.fq;
do
        fastqc $file -o fastqc1
done
#Loads multiqc from hpc and opens fastqc directory
module load multiqc/1.4
cd ./fastqc1/
#Runs multiqc
multiqc ./
````
#### Output in fastqc directory
````
./fastqc
-merg_bc07_rna_fastqc.html
-merg_bc08_rna_fastqc.html
-merg_bc09_rna_fastqc.html
-merg_bc10_rna_fastqc.html
-merg_bc11_rna_fastqc.html
-multiqc_report.html
````
The quality of the reads were good.
#### Adapter content
![image](https://i.imgur.com/H9pxbJu.png)
#### Sequence counts
![image](https://i.imgur.com/eQtSG3U.png)
#### Sequence duplication levels
![image](https://i.imgur.com/58tJC5t.png)
#### Overrepresented sequences
![image](https://i.imgur.com/CtEkMsk.png)

### Step 2. Filtering and sorting rRNA Sequences with SortMeRNA
#### Description

"SortMeRNA is a program tool for filtering, mapping and OTU-picking NGS reads in metatranscriptomic and metagenomic data. The core algorithm is based on approximate seeds and allows for fast and sensitive analyses of nucleotide sequences. The main application of SortMeRNA is filtering ribosomal RNA from metatranscriptomic data."

Once we have removed low quality sequences and remove any adapter contamination, we can then proceed to an additional (and optional) step to remove rRNA sequences from the samples. If your samples were not prepared with an rRNA depletion protocol before library preparation, it is reccomended to run this step to computational remove any rRNA sequence contiamation that may otheriwse take up a majority of the aligned sequences.
#### Installation
````
conda install sortmerna
````
## Preparation
Before we can run the `sortmerna` command, we must first download  the eukaryotic, archeal and bacterial rRNA databases. The `sortmerna_db/` folder will be the location that we will keep the files necessary to run **SortMeRNA**. These databases only need to be created once, so any future RNAseq experiements can use these files.
````
 # Download the sortmerna package (~2min) into sortmerna_db folder
    wget -P sortmerna_db https://github.com/biocore/sortmerna/archive/2.1b.zip

    # Decompress folder 
    unzip sortmerna_db/2.1b.zip -d sortmerna_db

    # Move the database into the correct folder
    mv sortmerna_db/sortmerna-2.1b/rRNA_databases/ sortmerna_db/

    # Remove unnecessary folders
    rm sortmerna_db/2.1b.zip
    rm -r sortmerna_db/sortmerna-2.1b
 
````
#### sorting
````
#This script takes in our fastq files and sorts the rRNA based on ribosomal subunits.

mkdir bac-16S
mkdir bac-23S
mkdir arc-16S
mkdir arc-23S
mkdir euk-18S
mkdir euk-28S
mkdir rfam-5S
mkdir rfam-5.8S
#Each loop generates mapped rRNA reads ,unmapped rRNA reads with respect to each sample and a log file.
#Removes the kvdb file to allows next run.(sortmerna guidlines)
for file in ./data/*.fq;
do
        name=$(basename ${file} .fq)
        echo $name
        sortmerna --ref ./databases/rRNA_databases/silva-bac-16s-id90.fasta  --reads $file  --aligned ./aligned/bac-16S-mapped-$name\
         --other ./aligned/unmapped-bac-16S-$name \
         --fastx ./aligned/unmapped-bac-16S-$name --threads 7
                rm -rfv /home/oscarmwaura/sortmerna/run/kvdb
                mv ./aligned/bac-16S* bac-16S
done

# Takes in the unmapped reads from bac-16srRna sortmerna run out put and runs to filter from next database.
for file in ./aligned/unmapped-bac-16S*;
do
        name=$(basename ${file} .fq)
        prefix="unmapped-bac-16S-"
        newname=${name#"$prefix"}
        echo "${newname}"
        sortmerna --ref ./databases/rRNA_databases/silva-bac-23s-id98.fasta --reads $file --aligned ./aligned/bac-23S-mapped-$newname\
         --other ./aligned/unmapped-bac-23S-$newname --fastx ./aligned/unmapped-bac-23S-$newname.fq --threads 7
               rm -rfv /home/oscarmwaura/sortmerna/run/kvdb
                mv ./aligned/bac-23* bac-23S
done

# Takes in the unmapped reads from bac-23srRna sortmerna run out put and runs to filter from next database.
for file in ./aligned/unmapped-bac-23*;
do
        name=$(basename ${file} .fq)
        prefix="unmapped-bac-23S-"
        newname=${name#"$prefix"}
        echo "${newname}"
        sortmerna --ref ./databases/rRNA_databases/silva-arc-16s-id95.fasta --reads $file --aligned ./aligned/arc-16S-mapped-$newname\
         --other ./aligned/unmapped-arc-16S-$newname --fastx ./aligned/unmapped-arc-16S-$newname.fq --threads 7
               rm -rfv /home/oscarmwaura/sortmerna/run/kvdb
                mv ./aligned/arc-16S* arc-16S
done

# Takes in the unmapped reads from arc-16SrRna sortmerna run out put and runs to filter from next database.
for file in ./aligned/unmapped-arc-16S*;
do
        name=$(basename ${file} .fq)
        prefix="unmapped-arc-16S-"
        newname=${name#"$prefix"}
        echo "${newname}"
        sortmerna --ref ./databases/rRNA_databases/silva-arc-23s-id98.fasta --reads $file --aligned ./aligned/arc-23S-mapped-$newname\
         --other ./aligned/unmapped-arc-23S-$newname --fastx ./aligned/unmapped-arc-23S-$newname.fq --threads 7
               rm -rfv /home/oscarmwaura/sortmerna/run/kvdb
                mv ./aligned/arc-23S* arc-23S
done

# Takes in the unmapped reads from arc-23SrRna sortmerna run out put and runs to filter from next database.
for file in ./aligned/unmapped-arc-23S*;
do
        name=$(basename ${file} .fq)
        prefix="unmapped-arc-23S-"
        newname=${name#"$prefix"}
        echo "${newname}"
        sortmerna --ref ./databases/rRNA_databases/silva-euk-18s-id95.fasta --reads $file --aligned ./aligned/euk-18S-mapped-$newname\
         --other ./aligned/unmapped-euk-18S-$newname --fastx ./aligned/unmapped-euk-18S-$newname.fq --threads 7
               rm -rfv /home/oscarmwaura/sortmerna/run/kvdb
                mv ./aligned/euk-18S* euk-18S
done

# Takes in the unmapped reads from euk-18SrRna sortmerna run out put and runs to filter from next database.
for file in ./aligned/unmapped-euk-18S*;
do
        name=$(basename ${file} .fq)
        prefix="unmapped-euk-18S-"
        newname=${name#"$prefix"}
        echo "${newname}"
        sortmerna --ref ./databases/rRNA_databases/silva-euk-28s-id98.fasta --reads $file --aligned ./aligned/euk-28S-mapped-$newname\
         --other ./aligned/unmapped-euk-28S-$newname --fastx ./aligned/unmapped-euk-28S-$newname.fq --threads 7
               rm -rfv /home/oscarmwaura/sortmerna/run/kvdb
                mv ./aligned/euk-28S* euk-28S
done

# Takes in the unmapped reads from euk-28srRna sortmerna run out put and runs to filter from next database.
for file in ./aligned/unmapped-euk-28S*;
do
        name=$(basename ${file} .fq)
        prefix="unmapped-euk-28S-"
        newname=${name#"$prefix"}
        echo "${newname}"
        sortmerna --ref ./databases/rRNA_databases/rfam-5s-database-id98.fasta --reads $file --aligned ./aligned/rfam-5S-mapped-$newname\
         --other ./aligned/unmapped-rfam-5S-$newname --fastx ./aligned/unmapped-rfam-5S-$newname.fq --threads 7
               rm -rfv /home/oscarmwaura/sortmerna/run/kvdb
                mv ./aligned/rfam-5S* rfam-5S
 done

# Takes in the unmapped reads from rfam-5SrRna sortmerna run out put and runs to filter from next database.
for file in ./aligned/unmapped-rfam-5S*;
do
        name=$(basename ${file} .fq)
        prefix="unmapped-rfam-5S-"
        newname=${name#"$prefix"}
        echo "${newname}"
        sortmerna --ref ./databases/rRNA_databases/rfam-5.8s-database-id98.fasta --reads $file --aligned ./aligned/rfam-5.8S-mapped-$newname\
         --other ./aligned/unmapped-rfam-5.8S-$newname --fastx ./aligned/unmapped-rfam-5.8S-$newname.fq --threads 7
               rm -rfv /home/oscarmwaura/sortmerna/run/kvdb
                mv ./aligned/rfam-5.8S* rfam-5.8S
done
````
#### Output
````
./Transcriptomics/
-mapped reads ./bac-16s/  bac-16S-mapped-merg_bc08_rna.fq <- sequences with rRNA
-Unmapped reads ./untarget_reads/ bac-16S-unmapped-merg_bc08_rna.fq  <- sequences without any rRNA
-logs .aligned/ bac-16S-mapped-merg_bc08_rna.log  <- log from SortMeRNA analysis
````
#### Example output of one of the mapped files(bac-16S)
![image](https://i.imgur.com/tVpV8WI.png)

## Taxonomic classification
### qiime2
#### Importing data
