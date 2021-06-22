# This code onverts a categorically coded {A,C,G,T,N} SNP to the numerically coded version
if(rstudioapi::isAvailable()) setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # WORKING DIRECTORY
# Checked 20210621
library(Rcpp)

# Switch datasets here
fasta_ = "cd_isolates"
# fasta_ = "mic_isolates"

folder = "FASTA"
# This code can numerically code a categorical coded SNP based on allele frequencies
sourceCpp("tableC.cpp")
num_code_SNP = function(x,n){
  num_SNP = rep(0,n)
  tbl = sort(tableC(x)/n, decreasing = T)
  Nloc = which(names(tbl)%in%"N")
  if(length(Nloc)>0) { # there are gaps
    if(length(tbl) == 2){ # This is the case with exactly 1 allele and gaps
      F_MA = 1 - (tbl[1]) # This is the frequency set to the major allele
    } else {
      F_MA = 1 - (tbl[1] + tbl[Nloc]) # This is the frequency set to the major allele
    }
    
    tbl = tbl[-c(1,which(names(tbl)%in%"N"))] # Exclude major allele and gap from table
    if(length(tbl) > 0) {# Other minor alleles exist
      FA = c(F_MA, tbl/sum(tbl)*(1-F_MA))
    } else {
      FA = F_MA
    }
  } else { # No gaps present
    F_MA = 1 - tbl[1] # Set major allele frequency
    tbl = tbl[-1] # Exclude major allele, there must be minor alleles
    FA = c(F_MA, tbl/sum(tbl)*(1-F_MA))
  }
  for(i in 1:length(FA)){
    num_SNP[which(x%in%names(FA[i]))] = unname(FA[i]) 
  }
  return(as.numeric(scale(num_SNP)))
}

SNP = readRDS(file.path(folder, paste(fasta_, "_mx.rds", sep = "")))
n = length(SNP[,2])
num_SNP = matrix(nrow = nrow(SNP), ncol = ncol(SNP))
for(i in 1:ncol(SNP)){
  num_SNP[,i] = num_code_SNP(SNP[,i], n)
}
# Check and remove NAs
ind = which(is.na(num_SNP))
c = floor((ind-1) / nrow(num_SNP)) + 1
rm_idx = unique(c) # These SNPs contain no variation or have poor quality
# Edit pos
pos = readRDS(file.path(folder, paste(fasta_, "_pos.rds", sep = "")))
pos = pos[-rm_idx]
num_SNP = num_SNP[,-rm_idx]
rownames(num_SNP) = rownames(SNP)
colnames(num_SNP) = as.character(pos)

paste("Any more NAs?", any(is.na(num_SNP)))
saveRDS(num_SNP, file.path(folder, paste(fasta_, "_mx_NC.rds", sep = "")))
saveRDS(pos, file.path(folder, paste(fasta_, "_pos_NC.rds", sep = "")))
