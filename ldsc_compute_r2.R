if(rstudioapi::isAvailable()) setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # WORKING DIRECTORY
# Checked 20220105
print("Start compute_r2.R for MIC")
library(foreach)
#library(ccaPP)
library(doParallel)
# library(treeWAS)
library(Rcpp)
registerDoParallel(cores = 12L)
sourceCpp("tableC.cpp")
getSNP = function(X){
  u = unique(X)
  if(length(u) == 1) return(0) else return(length(u))
  # 0 if no SNP, else # of SNPs returned
}

filter_alt = function(mxCC){ # This function returns the mx directly
  print("Perform filtering...")
  t_fs = Sys.time()
  uqmx = unique(mxCC, MARGIN = 2) # This only works with NCmx, much faster! No correlation filter applied here
  # Drop non-snps
  tbls = apply(unname(uqmx), 2, function(x) length(tableC(x))) 
  rmidx = which(tbls == 1)
  if(length(rmidx) > 0) uqmx = uqmx[,-rmidx] 
  print(paste("Filtering time:" , Sys.time() - t_fs))
  return(uqmx)
}

compute_r2 = function(mxCC, retained_snps, nSEQs){
  mxCC = unname(mxCC)
  print("Commence r2 computation")
  t = Sys.time()
  # r2scores = foreach(i = retained_snps, .combine = "c") %dopar% { # Pearson's
  #   r2 = as.numeric(cor(mxCC[,i],mxCC)^2)
  #   sum(r2 - (1-r2)/(nSEQs-2)) # LDSC bias adjustment 
  # }
  
  # r2scores = foreach(snp = retained_snps, .combine = "c") %dopar% { # Spearman's
  #   rsp = unname(apply(mxCC, 2, function(x) corSpearman(mxCC[,snp], x)^2)) # spearman's correlation 
  #   sum(rsp*(1 + (1-rsp^2)/(2*(nSEQs-3)))) # Olkin Pratt bias adjustment
  # }
  
  r2scores = foreach(i = retained_snps, .combine = "c") %dopar% { # RANKIT 
    r2 = as.numeric(cor(mxCC[,i],mxCC)^2) # Pearson's correlation 
    # sum(r2*(1 + (1-r2^2)/(2*(nSEQs-3)))) # Olkin Pratt bias adjustment
    sum(r2 - (1-r2)/(nSEQs-2)) # LDSC bias adjustment 
  }
  
  print("r2 computation complete")
  Sys.time() - t
  
  gc()
  # we can use 1/r2s as weights - omit below and return r2s
  
  # print("Commence weight computation")
  # mxCC = mxCC[, retained_snps]
  # print(str(mxCC))
  # t = Sys.time()
  # r2weights = foreach(i = 1:ncol(mxCC), .combine = "c") %dopar% {
  #   r2 = as.numeric(cor(mxCC[,i],mxCC)^2)
  #   sum(r2 - (1-r2)/(nSEQs-2))
  # }
  # Sys.time() - t
  # print("Weight computation complete")
  # gc()
  
  # return(list(r2scores = r2scores, r2weights = r2weights))
  return(list(r2scores = r2scores, retained_snps = retained_snps))
}

# Code to compute sum(r^2) for the SNP mx
# For MIC
# We need an R^2 + W computation without accessory genes
# We need an R^2 + W computation with accessory genes

print("Reading mxCC [mic]...")
# mxCC = readRDS("../../Accessory_Genome/DATASET2/Pensar2019/acc_mic_core_mx_NC.rds")
# mxCC = readRDS("acc_mic_core_mx_NC.rds") # place file here for spartan
# mxCC = readRDS("acc_mic_core_mx_NC3.rds") # place file here for spartan
# mxCC = readRDS("combi_mic_core_mx_NC.rds") # place file here for spartan
mxCC = readRDS("acc_mic_core_pys_mx_NC.rds") # For RANKIT, we should transform

# Let's quickly do the filtering
nSEQs = dim(mxCC)[1]
mxCC_filt = filter_alt(mxCC)
nSNPs = ncol(mxCC_filt)
print(paste("Successful! Found", nSNPs, "SNPs and", nSEQs, "SEQs"))
retained_snps = sort(which(colnames(mxCC) %in% colnames(mxCC_filt)))
print(str(retained_snps))

print("Performing rankit transform")
mxCC_r = apply(mxCC, 2, function(x) rank(x - 0.5)/nrow(mxCC))

rm(mxCC_filt, mxCC)

gc()

# op = compute_r2(mxCC, retained_snps, nSEQs)
op = compute_r2(mxCC_r, retained_snps, nSEQs)
saveRDS(op, "r2s_mic_pys_rkit.rds")
print("op saved")
gc()


################################
rm(op, mxCC_r, nSEQs, nSNPs, retained_snps) # REMOVE PREV VARS
################################
print("Reading mxCC [cdacute]...")
# For cd/acute
# mxCC = readRDS("../../Accessory_Genome/DATASET2/Pensar2019/acc_cd.acute_core_mx_NC.rds")
# mxCC = readRDS("acc_cd.acute_core_mx_NC.rds") # place file here for spartan
# mxCC = readRDS("acc_cd.acute_core_mx_NC3.rds") # place file here for spartan
# mxCC = readRDS("combi_cd.acute_core_mx_NC.rds") # place file here for spartan
mxCC = readRDS("cd_core_uqepi_pys_mx_NC.rds") # place file here for spartan
# Let's quickly do the filtering
nSEQs = dim(mxCC)[1]
mxCC_filt = filter_alt(mxCC)
nSNPs = ncol(mxCC_filt)
print(paste("Successful! Found", nSNPs, "SNPs and", nSEQs, "SEQs"))
retained_snps = sort(which(colnames(mxCC) %in% colnames(mxCC_filt)))
print(str(retained_snps))

print("Performing rankit transform")
mxCC_r = apply(mxCC, 2, function(x) rank(x - 0.5)/nrow(mxCC))

rm(mxCC_filt, mxCC)

gc()

# op = compute_r2(mxCC, retained_snps, nSEQs)
op = compute_r2(mxCC_r, retained_snps, nSEQs)
saveRDS(op, "r2s_cdacute_pys_rkit.rds")
print("op saved")
