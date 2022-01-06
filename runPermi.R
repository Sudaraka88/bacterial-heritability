if(rstudioapi::isAvailable()) setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # WORKING DIRECTORY
# Code for permutation testing and analysis 
# This code is largely a combination of GWAS_snps.R, pyseer_GAP.sh and GWAS_gapsnp_combo.R - Please refer to original codes for detailed comments
# Checked 20220105

library(lme4qtl)
library(Rcpp)
library(doParallel)
library(foreach)
library(ggplot2)
library(RcppArmadillo)
registerDoParallel(cores = 8L)

sourceCpp("tableC.cpp")
checkAlignment = function(nm1,nm2){
  print(paste("Pass alignment test?", all(nm1 == nm2)))
}
rename_simmx = function(nm){
  x = unlist(strsplit(nm,"_"))
  return(paste(x[1],"_",x[2],"_",x[3],sep = ""))
}
sim_mx_cleaner = function(sim_mx, pheno_order){ # fn to clean phylogenetic_similarity mx and re-arrange
  print("Cleaning sim_mx")
  rnms = rownames(sim_mx)
  if(all(pheno_order %in% rnms)) rnms = unname(sapply(rnms, function(x) rename_simmx(x)))
  rownames(sim_mx) = rnms
  colnames(sim_mx) = rnms
  sim_mx = as.matrix(sim_mx)
  # Now we need to arrange the sim_mx in order
  matches = c()
  for(i in pheno_order) matches[i] = which(rnms %in% i)
  checkAlignment(rownames(sim_mx)[matches], pheno_order)
  sim_mx = sim_mx[matches, matches]
  return(sim_mx)
}
getPhenoOrder = function(premutation, N){
  set.seed(premutation, kind = "Mersenne-Twister") # Fix permutation order
  return(sample(N,N, replace = F))
}
refmt = function(nm){
  if(length(grep("#", nm))){
    temp = unlist(strsplit(nm, split = "#"))
    return(paste(temp[1], "_", temp[2], sep = ""))
  } else {
    return(nm)
  }
}

phen_name = "cef.mic" # Give cef.mic/pen.mic - not run for CD because of no hits
if(phen_name == "cd"){
  gap_freq = readRDS("cd_core_uqepi_gap_freq.rds")
} else {
  gap_freq = readRDS("acc_mic_core_gap_freq.rds") # Distribution of gaps - not provided, can be generated with provided code
}

pheno = readRDS(paste("../Accessory_Genome/DATASET2/pen.mic_cef.mic_pheno_reduced.rds"))

pheno$sampleID = sapply(pheno$sampleID, refmt)
n = nrow(pheno)

ph = data.frame(samples = pheno$sampleID, y = log10(pheno$Ceftriaxone.MIC))
st = 1
nperm = 500

t0 = Sys.time()
foreach(perm = st:nperm) %dopar% { # Perform permutation testing
  fldr = paste("p", perm, sep = "")
  if(!file.exists(fldr)){
    dir.create(fldr) # save output to new folder
  }
  idxord = getPhenoOrder(perm, n)
  if(!file.exists(file.path(fldr, "perm.order"))){
    write.table(idxord, file.path(fldr, "perm.order"), quote = F, row.names = F, col.names = F) # write permutation order
  }
  ph_ = data.frame(samples = ph$samples, y = ph$y[idxord])
  # mxCC_ = mxCC[idxord,]
  
  # gap test
  if(!file.exists(file.path(fldr, paste("pheno_", phen_name, sep = "")))){
    write.table(ph_, file = file.path(fldr, paste("pheno_", phen_name, sep = "")), quote = F, row.names = F, sep = "\t")
  } 
  if(!file.exists(file.path(fldr, "pyseer_GAP.sh"))) file.copy("pyseer_GAP.sh", file.path(fldr, "pyseer_GAP.sh"))
  
  if(!file.exists(file.path(fldr, paste("res_", phen_name, "_lmm_ps.txt", sep = "")))){
    setwd(fldr)
    system("./pyseer_GAP.sh", intern = TRUE)
    setwd("../../PermutationTesting/")
  } # else { # won't work with parallel implementation
  #   print("Gap test data available from previous run")
  # }
  ###############################################################################
  # Reading in gap data
  # print("Reading gap data")
  gap_out =  read.csv(file.path(fldr, paste("res_", phen_name, "_lmm_ps.txt", sep = "")), sep = "\t")
  gap_out$variant = unname(sapply(gap_out$variant, function(x) unlist(strsplit(x, "_"))[2] ))
  gf = gap_freq[which(as.numeric(names(gap_freq)) %in% gap_out$variant)]
  if(length(which(!gap_out$variant %in% names(gf))) > 0) {
    gap_out = gap_out[-which(!gap_out$variant %in% names(gf)), ]
  }
  gap_out = cbind(gap_out, gf = gf)
  gap_out = gap_out[,-c(2,3,8)]
  colnames(gap_out) = c("pos","p", "beta","se", "h2", "gf")
  gap_out = apply(gap_out, 2, as.numeric)
  rm_idx = which(is.na(gap_out[,1]) | is.na(gap_out[,2]) | is.na(gap_out[,3]) | is.na(gap_out[,4]) | is.na(gap_out[,5]))
  if(length(rm_idx) > 0) gap_out = gap_out[-rm_idx, ]
  gap_out = as.data.frame(gap_out)
  
  ###############################################################################
  # snp test
  if(!(file.exists(file.path(fldr, paste("snpgwas_", phen_name, sep = ""))) | file.exists(file.path("spartan/",fldr, paste("snpgwas_", phen_name, sep = ""))))){
    if(!file.exists(file.path(fldr, paste("mod_", phen_name, ".rds", sep = "")))){
      sim_mx_ = sim_mx_cleaner(sim_mx = sim_mx, pheno_order = ph_$samples)
      mod = lme4qtl::relmatLmer(y ~ (1|samples), ph_, relmat = list(samples = sim_mx_))
      V <- lme4qtl::varcov(mod, idvar = "samples")
      decomp <- wlm::decompose_varcov(V, method = "evd", output = "all")
      W <- decomp$transform
      saveRDS(object = list(mod = mod, V = V, W = W), file = file.path(fldr, paste("mod_", phen_name, ".rds", sep = "")))
    } else {
      W = readRDS(file.path(fldr, paste("mod_", phen_name, ".rds", sep = "")))$W
    }
    
    y_t <- crossprod(W, ph_$y)
    # t = Sys.time()
    snp_out = data.frame()
    # snp_out = foreach(i = 1:ncol(mxCC), .combine='rbind') %dopar% {
    # i = 1
    for(i in 1:ncol(mxCC)){
      tbl = sort(tableC(mxCC[,i])/n, decreasing = T)
      gap_locs = which(names(tbl) == "N" |names(tbl) == "-")
      gf = 0
      if(length(gap_locs > 0)) {
        gf = sum(tbl[gap_locs])
        tbl = tbl[-gap_locs]
      }
      rm_idx = which(mxCC[,i] == "N" | mxCC[,i] == "-")
      # perform transformation
      if(length(rm_idx) > 0){
        gap_count = length(rm_idx)
        # AFC
        num_SNP = rep(0,n)
        for(j in 1:length(tbl)){
          num_SNP[which(mxCC[,i]%in%names(tbl[j]))] = unname(tbl[j])
        }
        num_SNP = num_SNP[-rm_idx]
        num_SNP_s = as.numeric(scale(num_SNP))
        # In rare cases, allele frequencies can be even across populations, resulting in 1/af distributions to all alleles
        if(any(is.nan(num_SNP_s))) {
          num_SNP = rep(0,n)
          k = 0
          for(j in 1:length(tbl)){
            num_SNP[which(mxCC[,i]%in%names(tbl[j]))] = k
            k = k+1
          }
          num_SNP = num_SNP[-rm_idx]
          num_SNP_s = scale(num_SNP)
        }
        
        X = unname(crossprod(W[-rm_idx, -rm_idx], as.numeric(scale(num_SNP))))
        
        mdlx = fastLmPure(X = cbind(X, 1), y_t[-rm_idx])
        
        
        # mdlx2 = lm.fit(x = cbind(X,1), y = y_t[-rm_idx])
        
        # mdlx_ = lm('y~.', data.frame(y = y_t[-rm_idx], x = X))
        
      } else {
        gap_count = 0
        num_SNP = rep(0,n)
        for(j in 1:length(tbl)){
          num_SNP[which(mxCC[,i]%in%names(tbl[j]))] = unname(tbl[j])
        }
        num_SNP_s = as.numeric(scale(num_SNP))
        # In rare cases, allele frequencies can be even across populations, resulting in 1/af distributions to all alleles
        if(any(is.nan(num_SNP_s))) {
          k = 0
          for(j in 1:length(tbl)){
            num_SNP[which(mxCC[,i]%in%names(tbl[j]))] = k
            k = k+1
          }
          num_SNP_s = scale(num_SNP)
        }
        
        X = unname(crossprod(W, as.numeric(num_SNP_s)))
        mdlx = fastLmPure(X = cbind(X, 1), y_t)
        # mdlx = lm('y~.', data.frame(y = y_t, x = X))
      }
      # unname(c(summary(mdlx_)$coefficients[scol,c(1,2,4)], tbl[1], gf , (n-gap_count)))
      snp_out = rbind(snp_out, unname(c(mdlx$coefficients[1], mdlx$stderr[1], 2*pt(abs(mdlx$coefficients[1]/mdlx$stderr[1]), mdlx$df.residual, lower.tail = FALSE), # p-val
                                        tbl[1], gf , (n-gap_count))))
    }
    
    snp_out = data.frame(snp_out)
    snp_out = cbind(as.numeric(colnames(mxCC)), snp_out)
    
    colnames(snp_out) = c("pos", "beta", "se", "p","MAF", "gf", "samples")
    # print(paste("SNP test complete. Elapsed time:", Sys.time() - t))
    saveRDS(snp_out, file.path(fldr, paste("snpgwas_", phen_name, sep = "")))
    
    # } else {
    #   # print("SNP test data available from previous run")
    docombi = T
    if(file.exists(file.path(fldr, paste("snpgwas_", phen_name, sep = "")))) {
      snp_out = readRDS(file.path(fldr, paste("snpgwas_", phen_name, sep = "")))
      snp_out$pos = as.numeric(snp_out$pos)
    } else if(file.exists(file.path("spartan/", fldr, paste("snpgwas_", phen_name, sep = "")))){
      snp_out = readRDS(file.path("spartan/", fldr, paste("snpgwas_", phen_name, sep = "")))
      snp_out$pos = as.numeric(snp_out$pos)
    } else {
      docombi = F
    }
    #   
  }
  ###############################################################################
  # combi & max
  if(docombi & !file.exists(file.path(fldr, paste("gwas_", phen_name, sep = "")))){
    
    gap_out = cbind(gap_out, z = gap_out$beta/gap_out$se)
    snp_out = cbind(snp_out, z = snp_out$beta/snp_out$se)
    
    combi = match(gap_out$pos, snp_out$pos)
    combi_idx = gap_out$pos %in% snp_out$pos
    
    paste("combi check passed?", all(gap_out$pos[combi_idx] == snp_out$pos[combi[combi_idx]]))
    
    gwas_full = data.frame(pos = gap_out$pos[combi_idx],
                           z.combi = (gap_out$z[combi_idx] + snp_out$z[combi[combi_idx]])/sqrt(2),
                           p.combi = apply(cbind(gap_out$p[combi_idx], snp_out$p[combi[combi_idx]]), MARGIN = 1, function(x) survcomp::combine.test(p = x,method = "z.transform")),
                           Test.combi = rep("combi", length(which(combi_idx))),
                           z.max = apply(cbind(gap_out$z[combi_idx], snp_out$z[combi[combi_idx]]), 1, max),
                           p.max = apply(cbind(gap_out$p[combi_idx], snp_out$p[combi[combi_idx]]), 1, min),
                           Test.max = apply(cbind(gap_out$z[combi_idx], snp_out$z[combi[combi_idx]]), 1, function(x) ifelse(x[1] >= x[2], 'gap', 'snp')),
                           gf = gap_out$gf[combi_idx])
    
    gap_only_idx = seq(1,nrow(gap_out))[-which(gap_out$pos %in% gwas_full$pos)]
    snp_only_idx = seq(1,nrow(snp_out))[-which(snp_out$pos %in% gwas_full$pos)]
    
    gwas_full = rbind(gwas_full,
                      data.frame(pos = gap_out$pos[gap_only_idx],
                                 z.combi = gap_out$z[gap_only_idx],
                                 p.combi = gap_out$p[gap_only_idx],
                                 Test.combi = rep("gap", length(gap_only_idx)),
                                 z.max = gap_out$z[gap_only_idx],
                                 p.max = gap_out$p[gap_only_idx],
                                 Test.max = rep("gap", length(gap_only_idx)),
                                 gf = gap_out$gf[gap_only_idx]),
                      data.frame(pos = snp_out$pos[snp_only_idx],
                                 z.combi = snp_out$z[snp_only_idx],
                                 p.combi = snp_out$p[snp_only_idx],
                                 Test.combi = rep("snp", length(snp_only_idx)),
                                 z.max = snp_out$z[snp_only_idx],
                                 p.max = snp_out$p[snp_only_idx],
                                 Test.max = rep("snp", length(snp_only_idx)),
                                 gf = snp_out$gf[snp_only_idx]))
    
    
    gwas_full = gwas_full[order(gwas_full$pos),]
    gwas_full$pos = as.numeric(gwas_full$pos)
    gwas_full = gwas_full[order(gwas_full$pos), ]
    saveRDS(gwas_full, file.path(fldr, paste("gwas_", phen_name, sep = "")))
  }
  # } else {
  #   # gwas_full = readRDS(file.path(fldr, paste("gwas_", phen_name, sep = "")))
  # }
  # 
  # # Make plots
  # # print("Preparing outputs, existing will be overwritten...")
  # # plt1 = makeplots(phen_name, as.numeric(gwas_full$pos), gwas_full$p.combi, gwas_full$gf, gwas_full$Test.combi, "combi")
  # # ggsave(plot = plt1, device = "png",  filename = file.path(fldr, paste(phen_name, "_combi.png", sep = "")), width = 15, height = 12)
  # # 
  # # plt2 = makeplots(phen_name, as.numeric(gwas_full$pos), gwas_full$p.max, gwas_full$gf, gwas_full$Test.max, "max")
  # # ggsave(plot = plt1, device = "png",  filename = file.path(fldr, paste(phen_name, "_max.png", sep = "")), width = 15, height = 12)
  ###############################################################################
} 

print(paste("Pipeline complete. Elapsed time:", Sys.time() - t0))