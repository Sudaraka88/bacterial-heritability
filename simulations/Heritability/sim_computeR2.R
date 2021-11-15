if(rstudioapi::isAvailable()) setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # WORKING DIRECTORY

library(foreach)
library(doParallel)
library(Rcpp)
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
  
  r2scores = foreach(i = retained_snps, .combine = "c") %dopar% { # RANKIT 
    r2 = as.numeric(cor(mxCC[,i],mxCC)^2) # Pearson's correlation 
    # sum(r2*(1 + (1-r2^2)/(2*(nSEQs-3)))) # Olkin Pratt bias adjustment
    sum(r2 - (1-r2)/(nSEQs-2)) # LDSC bias adjustment 
  }
  
  print("r2 computation complete")
  Sys.time() - t
  
  gc()
  
  return(list(r2scores = r2scores, retained_snps = retained_snps))
}

registerDoParallel(cores = 4L)
sourceCpp("../Genotype2/tableC.cpp")

folder_ = dir(pattern = "^[:h:][:0:]")

# These are identical, no point in repeating!
# scanning step
for(folder in folder_){
  fasta_ = paste(folder, "/simulations/genSim/genSim", sep = "")
  if(file.exists(paste(fasta_, "_r2s.rds", sep = ""))){
    op = readRDS(paste(fasta_, "_r2s.rds", sep = ""))
  } else {
    print("No pre-computed instance!")
  }
}


for(folder in folder_){
  print(paste("Now processing", folder))
  fasta_ = paste(folder, "/simulations/genSim/genSim", sep = "")
  if(!file.exists(paste(fasta_, "_r2s.rds", sep = ""))){
    if(exists("op")){
      saveRDS(op, paste(fasta_, "_r2s.rds", sep = ""))
      print("op copied")
    } else {
      getSNP = function(X){
        u = unique(X)
        if(length(u) == 1) return(0) else return(length(u))
        # 0 if no SNP, else # of SNPs returned
      }
      
      print("Reading mxAFC...")
      mxCC = readRDS(paste(fasta_, "_AFC.rds", sep = "")) # For RANKIT, we should transform
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
      op = compute_r2(mxCC_r, retained_snps, nSEQs)
      saveRDS(op, paste(fasta_, "_r2s.rds", sep = ""))
      print("op saved")
    }
  } else {
    print("op exists")
  }
}
