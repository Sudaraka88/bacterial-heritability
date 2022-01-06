# Code, data and figures for bacterial-heritability (https://www.biorxiv.org/content/10.1101/2021.10.04.462983v1)


![Circos](/figures/circ.png)

## Guide for running the code
This code is prepared to run on Linux based systems, but should be portable to any OS with minimum effort...

### Preparing the fasta multiple sequence alignment (MSA)
- Download the read sequence data in **Accession_lane.csv** (use [sratoolkit](https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit))
- Download and setup [snippy](https://github.com/tseemann/snippy)
- Prepare the MSA with the [ATCC700669](https://www.ncbi.nlm.nih.gov/assembly/GCA_000026665.1) reference genome

### Preparing phenotype data
- Maela swab data from [Chewapreecha et al. 2014](https://www.nature.com/articles/ng.2895) can be obtained from the authours
    - Use **pheno_AMR.R** to format and save antimicrobial resistance data from raw files
    - Use **pheno_mic.R** to format and save minimum inhibitory concentration (MIC) data from raw files
    - First run **pheno_pathGen.R** to generate carriage paths from the raw files (helper functions are available in **fn**)
    - Then run **pheno_carDur.R** to extract carriage episodes with durations
    - Then run **pheno_carCount.R** to extract further details about carriage episodes/counts (optional)

Final outputs will be written to a folder named **Out**, temporary data will be written to **TempData**. 
- Prepare two text files **cd_isolates** and **mic_isolates**, each containing the subset of isolates linked with each phenotype

### Preparing data for pangenome analysis (optional)
- Download and setup [panaroo](https://github.com/gtonkinhill/panaroo)
    - Use the downloaded read sequence data with and the annotated reference genome to prepare a gene presence/absence CSV file in roary format

### Partioning and convering fasta genotype data into rds format
- Extract the sequences in **cd_isolates** and **mic_isolates** into separate fasta files (use [seqtk subseq](https://github.com/lh3/seqtk))

For analysis using **R**, it is convenient to convert these MSA fasta files into native **rds** format

    - Run **geno_extractSNP.R** to filter, extract SNPs and save the resulting matrix as an RDS file
        - Alternative filtering ethods are available by changing the output folder name
    - Then run **geno_numCodeSNPs.R** to perform allele frequency coding (AFC) on the previously generated R matrices

## File description
### Accession Lanes
- Accession_lanes.csv   - Genotype data

### Helper Files
- init.R			- Information on installing the required packages
- *reference_genome.rds*	- SPN23F reference genome

### Genotype Preparation
- extractSNPs.sh		- Use with geno_extractSNPs.R
- extractSEQstats.sh	- Extract statistics from fasta files
- createFreqFile.sh		- Extract sequences from fasta files
- geno_extractSNPs.R	- Extract SNPs from core genome alignments
- geno_getAccGenDat.R	- Extract accessory genome data from panaroo output
- geno_numCodeSNPs.R	- Numerically code {A,C,G,T,N} SNP and accessory genes
- tableC.cpp			- Quick ‘table()’ implementation using Rcpp
- Rcppf.cpp			- Helper function for geno_extractSNPs.R

### Phenotype Preparation
- pheno_AMR.R			- Extract anti microbial resistance data
- pheno_carCount.R		- Extract carriage count data
- pheno_carDur.R		- Extract carriage duration data
- pheno_mic.R			- Extract minimum inhibitory concentration data
- pheno_pathGen.R		- Generate paths for carriage duration estimation
- **fn** 				- Contains helper functions for phenotype preparation

### Analysis
- gls_GWAS.R			- Perform GWAS on all phenotypes (num-coded)
- ldsc_h2.R 			- Estimate heritability using LDSC model (num-coded)
- PRED_enet.R			- Predict phenotype with elastic net model (num-coded)
- PRED_MLR.R			- Predict phenotype with linear regression (num-coded) - less accurate
- fec_SNPs.R			- Isolate large effect sized SNPs using clumping
- summarise_preds.R		- Summarise prediction outcomes
- map2gene.R			- Search for identified gene in the ref. Genome
- GWAS_gaps             - Run the GAP test using lme4qtl (less accurate)
- GWAS_snps.R           - Run the SNP test using lme4qtl
- GWAS_gapsnp_combo.R   - Combine the outputs from gap/snp testing pipeline

### Phylogeny and Clustering
- build_trees.sh		- Script for iqtree2 and ClonalFrameML
- **Trees**				- Contains output generated from built_trees.sh
- fastbaps_clusts.R		- Generate fastbaps clusters from core genome alignments

### pyseer
- pyseer_tsttrnPheno.R	- Break phenotype into test/train portions
- pyseer_mkPheno.R		- Prepare phenotype for pyseer analysis
- pyseer_mkCov.R		- Prepare covariates for pyseer analysis
- pyseer_GWAS.sh		- Script to perform GWAS using pyseer
- pyseer_loso.sh		- Script to perform loso analysis using pyseer
- pyseer_xfcv.sh		- Script to perform xfcv analysis using pyseer
- pyseer_predict.sh		- Script to prepare prediction models for pyseer
- pyseer_VCF2bialle.R	- Code the VCF file to bi-allelic for pyseer analyses
- pyseer_com_mae_r2.R	- Summarise prediction outcomes from pyseer
- pyseer_GAP.sh         - Run the GAP test using pyseer (preferred option)

### treeWAS
- treeWAS_cd.R		- perform treeWAS analysis for carriage duration
- treeWAS_mic.R		- perform treeWAS analysis for MIC phenotypes

### simulations
- **AssociationTesting**    - Codes used for the association testing simulation (detailed versions available in **Analysis**)
- **Heritability**          - Codes used for the heritability simulation (detailed versions available in **Analysis**)

### figures
All paper figures in high quality format

