# Code, data and figures for bacterial-heritability (https://www.biorxiv.org/content/10.1101/2021.10.04.462983v1)


![Circos](/figures/circ.png)

## Setting up required packages and files
This code is prepared to run on Linux based systems, but should be portable to any OS with minimum effort. Please refer to **init.R** for details on installing the required packages.

### Preparing the fasta multiple sequence alignment (MSA)
- Download the read sequence data in **Accession_lane.csv** (use [sratoolkit](https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit))
- Download and setup [snippy](https://github.com/tseemann/snippy)
- Prepare the MSA with the [ATCC700669](https://www.ncbi.nlm.nih.gov/assembly/GCA_000026665.1) reference genome
- Save a VCF copy of fasta MSA for pyseer analyses (can use [snp-sites](http://sanger-pathogens.github.io/snp-sites/))

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
- Extract the sequences in **cd_isolates** and **mic_isolates** into separate fasta files (can use [seqtk subseq](https://github.com/lh3/seqtk))
- Save a VCF copy of fasta MSAs for pyseer analyses (can use [snp-sites](http://sanger-pathogens.github.io/snp-sites/))

For analysis using **R**, it is convenient to convert these MSA fasta files into native **rds** format
    - Run **geno_extractSNPs.R** to filter, extract SNPs and save the resulting matrix
    - Run **geno_extractGAPs.R** to filter, extract gaps and save the resulting matrix
        - Alternative filtering methods are available by changing the output folder name
    - Then run **geno_numCodeSNPs.R** to perform allele frequency coding (AFC) on the previously generated R matrices
    - For pangenome analysis (optional), use **geno_getAccGenDat.R** to prepare gene count matrix

### Phylogeny and Clustering
- Install and setup [iqtree](http://www.iqtree.org/doc/Quickstart)
- Install and setup [fastbaps](https://github.com/gtonkinhill/fastbaps)
- Install and setup [ClonalFrameML] (https://github.com/xavierdidelot/ClonalFrameML)

To generate the phylogeny from the fasta MSA, use **build_trees.sh**. It will also perform the ClonalFrameML inference of bacterial microevolution recommended for [treeWAS](https://github.com/caitiecollins/treeWAS)

- Run **fastbaps_clusts.R** to perform Bayesian population clustering of isolates in the paritioned MSA fasta files using the baps algorithm


## Genome-wide association analysis
### Perform gap/snp testing
- For population structure correction, generate phylogeny based distance and similarity matrices using pyseer (see [here](https://pyseer.readthedocs.io/en/master/usage.html#population-structure) for details)
    - The gap/snp testing pipeline only requires the phylogeny similarity matrix in default (tsv) format
- Gap testing can be performed using two approaches
    - **GWAS_gaps.R** uses [lme4qtl](https://github.com/variani/lme4qtl) to fit an approximate linear model on the extracted gap matrix - directly compatible with **geno_extractGAPs.R** output
    - **pyseer_GAP.sh** uses [pyseer](https://pyseer.readthedocs.io/en/master/) to fit the exact [FaST-LMM model](https://github.com/fastlmm/FaST-LMM) model to gap data (recommended)
        - Use **pyseer_GAP2VCF.R** to convert the output to VCF format as required by pyseer (more information below)
- SNP testing fits the **lme4qtl** approximate linear model on allele frequency coded data after omitting gaps
    - Run **GWAS_snps.R** for snp testing (directly compatible with the **geno_numCodeSNPs.R** output)
- Finally, run **GWAS_gapsnp_combo.R** to combine the gap/snp statistics (options to use Stouffer's method or max)
    - This code currently accepts the gap testing output from **pyseer_GAP.sh**, to use the output from **GWAS_gaps.R**, simple read the file in using **readRDS()**
- For permutation testing, run **runPermi.R**, this file is a combination of **pyseer_GAP.sh**, **GWAS_snps.R** and **GWAS_gapsnp_combo.R**
    

### Perform major allele testing
- Major allele testing can be performed using multiple methods
    - Run **gls_GWAS.R** to fit the approximate linear model in lme4qtl to either AFC or major allele coded (MAC) SNPs - requires the output from **geno_extractSNPs.R** and/or **geno_numCodeSNPs.R**
    - Run **pyseer_GWAS.sh** to fit the FaST-LMM model
        - Use **pyseer_mkPheno.R, pyseer_mkCov.R** to prepare the phenotype and covariates
        - Use **pyseer_VCF2bialle.R** to force gaps to be considered as major alleles (Warning! Might lead to spurious associations)
    - Run **treeWAS_cd/mic.R** to perform treeWAS analysis
        - Requires treeWAS to be setup and installed, refer to **init.R**
        - The use of the ClonalFrameML phylogeny (see above) is recommended for recombinant species such as *S. pneumoniae*

## Heritability analysis
- Run **ldsc_h2.R** to use the LDSC model for heritability estimation and partioning
    - Requires LD scores to be computed using **ldsc_compute_r2.R**
    - These two scripts are compatible with AFC output generated from **geno_numCodeSNPs.R** and/or **geno_getAccGenDat.R**
- Run **pyseer_GWAS.sh** to print an LMM based heritability estimate (see [here](https://pyseer.readthedocs.io/en/master/usage.html#mixed-model-fast-lmm) for details)
- See below for elastic net based heritability estimation


## Phenotype prediction
### Allele frequency coding
- Run **PRED_enet.R** with **geno_numCodeSNPs.R** and/or **geno_getAccGenDat.R** outputs
    - Uses GLMNET to fit a whole-genome elastic net model similar to pyseer
    - Code supports partioning for 10-fold (xfcv) and leave-one-strain-out (LOSO) cross validation for benchmarking
- Run **PRED_MLR.R** with **geno_numCodeSNPs.R** and/or **geno_getAccGenDat.R** outputs (less accurate and not recommended)
    - Fits a multiple linear regression model to benchmark prediction accuracy
- Finally, use **summarise_preds.R** to summarise prediction outcomes, compute summary statistics and generate plots

### Major allele coding
- Pyseer uses a whole-genome elastic net model based on glmnet for phenotype prediction, use **pyseer_predict.sh**
    - Compatible with generated VCF files and outputs from **pyseer_mkPheno.R, pyseer_mkCov.R**
- For controlled benchmarking, first use **pyseer_txttrnPheno.R** to partition the phenotype data set as desired
    - Code supports partioning for 10-fold (xfcv) and leave-one-strain-out (LOSO) cross validation
    - LOSO partioning requires clusters, compatible with baps clusters from **fastbaps_clusts.R**
- Then **pyseer_loso/xfcv.sh** can be used to perform 10-fold or LOSO benchmarking
- Finally, use **pyseer_comp_mae_r2.R** to summarise prediction outcomes, compute summary statistics and generate plots

Both whole genome elastic net models (major allele and allele frequency coded) output a heritability estimate as described in the main text.

## Quick File description
### Accession Lanes
- Accession_lanes.csv   - Genotype data

### Helper Files
- init.R			- Information on installing the required packages
- *reference_genome.rds*	- SPN23F reference genome

### Genotype Preparation
- extractGAPs.sh		- Use with geno_extractGAPs.R
- extractSNPs.sh		- Use with geno_extractSNPs.R
- extractSEQstats.sh	- Extract statistics from fasta files
- createFreqFile.sh		- Extract sequences from fasta files
- geno_extractSNPs.R	- Extract SNPs from core genome alignments - save both SNP and major allele coded matrix
- geno_getAccGenDat.R	- Extract accessory genome data from panaroo output
- geno_numCodeSNPs.R	- AFC coding of SNP matrix
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
- gls_GWAS.R			- Perform GWAS on all phenotypes (AFC)
- ldsc_compute_r2.R     - Calculate LD scores and weights necessary to use **ldsc_h2.R**
- ldsc_h2.R 			- Estimate heritability using LDSC model (AFC)
- PRED_enet.R			- Predict phenotype with elastic net model (AFC)
- PRED_MLR.R			- Predict phenotype with linear regression (AFC) - less accurate
- fec_SNPs.R			- Isolate large effect sized SNPs using clumping
- summarise_preds.R		- Summarise prediction outcomes
- map2gene.R			- Search for identified gene in the ref. Genome
- GWAS_gaps.R           - Run the GAP test using lme4qtl (less accurate)
- GWAS_snps.R           - Run the SNP test using lme4qtl
- GWAS_gapsnp_combo.R   - Combine the outputs from gap/snp testing pipeline
- runPermi.R            - Code for permutation testing

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
- pyseer_GAP2VCF.R
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

