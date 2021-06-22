init.R			- Information on installing the required packages

extractSNPs.sh		- Use with geno_extractSNPs.R
extractSEQstats.sh	- Extract statistics from fasta files
createFreqFile.sh		- Extract sequences from fasta files
geno_extractSNPs.R	- Extract SNPs from core genome alignments
geno_getAccGenDat.R	- Extract accessory genome data from panaroo output
geno_numCodeSNPs.R	- Numerically code {A,C,G,T,N} SNP and accessory genes
tableC.cpp			- Quick ‘table()’ implementation using Rcpp
Rcppf.cpp			- Helper function for geno_extractSNPs.R

fastbaps_clusts.R		- Generate fastbaps clusters from core genome alignments

pheno_AMR.R			- Extract anti microbial resistance data
pheno_carCount.R		- Extract carriage count data
pheno_carDur.R		- Extract carriage duration data
pheno_mic.R			- Extract minimum inhibitory concentration data
pheno_pathGen.R		- Generate paths for carriage duration estimation
fn 				- Contains helper functions for phenotype preparation

gls_GWAS.R			- Perform GWAS on all phenotypes (num-coded)
ldsc_h2.R			- Estimate heritability using LDSC model (num-coded)
PRED_enet.R			- Predict phenotype with elastic net model (num-coded)
PRED_MLR.R			- Predict phenotype with linear regression (num-coded)
fec_SNPs.R			- Isolate large effect sized SNPs using clumping
summarise_preds.R		- Summarise prediction outcomes

map2gene.R			- Search for identified gene in the ref. Genome
build_trees.sh		- Script for iqtree2 and ClonalFrameML
Trees				- Contains output generated from built_trees.sh
reference_genome.rds	- SPN23F reference genome 

pyseer_tsttrnPheno.R	- Break phenotype into test/train portions
pyseer_mkPheno.R		- Prepare phenotype for pyseer analysis
pyseer_mkCov.R		- Prepare covariates for pyseer analysis
pyseer_GWAS.sh		- Script to perform GWAS using pyseer
pyseer_loso.sh		- Script to perform loso analysis using pyseer
pyseer_xfcv.sh		- Script to perform xfcv analysis using pyseer
pyseer_predict.sh		- Script to prepare prediction models for pyseer
pyseer_VCF2bialle.R	- Code the VCF file to bi-allelic for pyseer analyses
pyseer_com_mae_r2.R	- Summarise prediction outcomes from pyseer

treeWAS_cd.R		- perform treeWAS analysis for carriage duration
treeWAS_mic.R		- perform treeWAS analysis for MIC phenotypes




