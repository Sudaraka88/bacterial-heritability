if(rstudioapi::isAvailable()) setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # WORKING DIRECTORY
library(foreach)
library(doParallel)
# library(ggplot2)
# library(gridExtra)
library(Rcpp)
library(RcppArmadillo)
# library(fastDummies)

sourceCpp("tableC.cpp")
registerDoParallel(cores = 8L)
# library(numbers)
rename_simmx = function(nm){
  x = unlist(strsplit(nm,"_"))
  return(paste(x[1],"_",x[2],"#",x[3],sep = ""))
}
checkAlignment = function(nm1,nm2){
  print(paste("Pass alignment test?", all(nm1 == nm2)))
}
dist_mx_cleaner = function(dmx, pheno_order){
  print("Cleaning dmx")
  rnms = rownames(dmx)
  rnms = unname(sapply(rnms, function(x) rename_simmx(x)))
  rownames(dmx) = rnms
  colnames(dmx) = rnms
  dmx = as.matrix(dmx)
  # Now we need to arrange the sim_mx in order
  matches = c()
  for(i in pheno_order) matches[i] = which(rnms %in% i)
  checkAlignment(rownames(dmx)[matches], pheno_order)
  dmx = dmx[matches, matches]
  return(dmx)
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
acc_gene_cleaner = function(acc_genes){
  # acc_genes = acc_genes[,grep("group", colnames(acc_genes))] # Dropping non-group genes - not a good approach 20210514
  # Let's try the 95% rule instead
  gv = apply(acc_genes, 2, getVariation) # These show no variation
  acc_genes = acc_genes[,-which(gv == 0)]
  pop_test = apply(acc_genes, 2, function(x) sort(table(x)/nrow(acc_genes), decreasing = T)[1] < 0.95)
  acc_genes = acc_genes[,pop_test]
  acc_genes = apply(acc_genes, 2, function(x) x - mean(x)) # normalise
  return(unique.matrix(acc_genes))
}
getVariation = function(X){
  u = unique(X)
  if(length(u) == 1) return(0) else return(length(u))
  # 0 if no SNP, else # of SNPs returned
}

# LDSC model
perform_ldsc = function(y_, r2scores, w, N, M, mxCC, pcax){ #, fec = FALSE){
  
  t = Sys.time()
  
  out_df = data.frame()
  for(j in 1:ncol(pcax)){
    
    X_ = as.matrix(cbind(1, pcax[,1:j]))
   
    # PCA SS Explained
    av = anova(lm(y~., data.frame(y = y_, pcax[,1:j])))
    rss = av[nrow(av), 2]
    tss = sum(av[,2])
    ssexp = (tss - rss)/tss
  
    chi2 = foreach(i = 1:ncol(mxCC), .combine = "c") %dopar% {
      mdlx_ = fastLm(X = cbind(mxCC[,i] , X_), y = y_)
     
      (mdlx_$coefficients[1]/mdlx_$stderr[1])^2
    }
    
    ldsc = lm.wfit(x = as.matrix(cbind(1, r2scores)), y = chi2, w = w)
    h2 = unname(coef(ldsc)[2])*M/(N*(1-1/N))
    print(paste("ssexp = ", ssexp, ", h2 = ", h2, sep = ""))
    
    out_df = rbind(out_df, c(j, rss, tss, ssexp, h2)) 
  }
  colnames(out_df) = c("PCAs", "RSS", "TSS", "SSexp", "h2")
  print(Sys.time() - t)
  return(out_df)
  
}

out_full_full = data.frame()
folder_ = dir(pattern = "^[:h:][:0:]")

TT = Sys.time()
for(folder in folder_){
  if(file.exists(paste(folder, "/ldsc.rds", sep = ""))){
    print(paste("File", paste(folder, "/ldsc.rds", sep = ""), "exists"))
    out_full = readRDS(paste(folder, "/ldsc.rds", sep = ""))
  } else {
    t0 = Sys.time()
    out_full = data.frame()
    r2 = readRDS(paste(folder, "/simulations/genSim/genSim_r2s.rds", sep = "")) # These are the r2 outputs
    mxCC = readRDS(paste(folder, "/simulations/genSim/genSim_AFC.rds", sep = ""))
    dmx = readRDS(paste(folder, "/simulations/genSim/genSim_phylo_dist.rds", sep = ""))
    pca = prcomp(dmx, scale. = T, rank. = 10) # use first 20
    pcax = pca$x
    if(length(r2$retained_snps) > 0){
      mxCC = mxCC[,r2$retained_snps]
    } else {
      mxCC = filter_alt(mxCC)
    }
    r2scores = r2$r2scores
    names(r2scores) = colnames(mxCC)
    
    
    w = 1/sapply(r2$r2scores, function(x) max(c(1, x)))
  
    N = nrow(mxCC)
    M = ncol(mxCC)
    
    i_ = length(dir(paste(folder, "/simulations/phenSim/", sep = ""))) - 1
    print(paste("Found", i_+1, "phenos"))
    for(i in 0:i_){
      print(paste("Pheno", i))
      if(file.exists(paste(folder, "/simulations/phenSim/",i,"/ldsc_out.rds", sep = ""))){
        print(paste("File", paste(folder, "/simulations/phenSim/",i,"/ldsc_out.rds", sep = ""), "exists"))
        out_df = readRDS(paste(folder, "/simulations/phenSim/",i,"/ldsc_out.rds", sep = ""))
      } else{
        y = read.table(paste(folder, "/simulations/phenSim/",i,"/phenSim.phen", sep = ""))$V3
      
        
        out_df = perform_ldsc(y, r2scores, w, N, M, mxCC, pcax)
        out_df = cbind(out_df, phen = i)
        
        saveRDS(out_df, paste(folder, "/simulations/phenSim/",i,"/ldsc_out.rds", sep = ""))
      }
      
      out_full = rbind(out_full, out_df)
     
    }
    out_full = cbind(out_full, folder)
    saveRDS(out_full, paste(folder, "/ldsc.rds", sep = ""))
    print(Sys.time() - t0)
  }
  
  out_full_full = rbind(out_full_full, out_full)

}
saveRDS(out_full_full, "ldsc_out.rds")
print(Sys.time() - TT)
