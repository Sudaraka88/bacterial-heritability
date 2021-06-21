# These R packages are needed, run to install if missing...
list.of.packages <- c("ggplot2", "Rcpp", "gridExtra", "data.table", "foreach", "doParallel", "MASS", "glmnet", "devtools",
                      "ggtree", "phytools", "DECIPHER", "dplyr", "colorspace", "ggpubr", "ape")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# Install lme4qtl: More information https://github.com/variani/lme4qtl
library(devtools)
install_github("variani/lme4qtl")
install_github("variani/matlm")
install_github("variani/wlm")
install_github("variani/qq")

# Install pyseer: More information https://pyseer.readthedocs.io/en/master/installation.html

# Install iqtree2: More information http://www.iqtree.org/

# Install ClonalFrameML: https://github.com/xavierdidelot/ClonalFrameML

# Install treeWAS: More information https://github.com/caitiecollins/treeWAS
install_github("caitiecollins/treeWAS")

