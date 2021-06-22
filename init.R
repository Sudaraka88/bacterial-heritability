############# These R packages are required, run to install if missing ###########
list.of.packages <- c("ggplot2", "Rcpp", "gridExtra", "data.table", "foreach", 
                      "doParallel", "MASS", "glmnet", "devtools", "ggtree", 
                      "phytools", "DECIPHER", "dplyr", "colorspace", "ggpubr", 
                      "ape", "vcfR", "stringr", "msm", "minqa", "zoo")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

################## These packages are required from github #####################
library(devtools)
# Install fastbaps
install_github("gtonkinhill/fastbaps")

# Install lme4qtl: More information https://github.com/variani/lme4qtl
install_github("variani/lme4qtl")
install_github("variani/matlm")
install_github("variani/wlm")
install_github("variani/qq")

# Install treeWAS: More information https://github.com/caitiecollins/treeWAS
install_github("caitiecollins/treeWAS")

################# These packages should be installed manually ##################
# pyseer: More information https://pyseer.readthedocs.io/en/master/installation.html
# iqtree2: More information http://www.iqtree.org/
# ClonalFrameML: https://github.com/xavierdidelot/ClonalFrameML

########## These packages were used for SNP extraction and filtering ###########
# samtools: http://www.htslib.org/
# snp-sites: https://github.com/sanger-pathogens/snp-sites
# vcftools: https://vcftools.github.io/index.html
