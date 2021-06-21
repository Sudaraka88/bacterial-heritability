# see: https://www.mv.helsinki.fi/home/mjxpirin/GWAS_course/material/GWAS2.html for a description of GWAS statistics
if(rstudioapi::isAvailable()) setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # WORKING DIRECTORY
library(foreach)
library(doParallel)
library(ggplot2)
library(gridExtra)
# library(numbers)
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
# LDSC model output viewer
view_op = function(ph, ldsc, df, M, N, nbin = 0, alpha_ = 0.75){
  # Prepare function for pub 20210527
  df = df[order(df$LDscrs),] 
  if(nbin > 0){
    # Let's try to bin this
    df_plt = c()
    for(i in 1:nbin){
      idx_l = (i-1)*ceiling(M/nbin)+1;
      idx_h = min(i*ceiling(M/nbin), M)
      df_plt = rbind(df_plt, c(mean(df$chi2[idx_l:idx_h]), 
                               mean(df$LDscrs[idx_l:idx_h]),
                               max(df$weights[idx_l:idx_h])))
      if(idx_h == M) break
    }
    
    df_plt = as.data.frame(df_plt)
    colnames(df_plt) = names(df)
  } else {
    df_plt = df
  }
  dfln = data.frame(intercept = unname(coef(ldsc)[1]) ,
                    slope = unname(coef(ldsc)[2]))
  
  gp = ggplot() + geom_point(data = df_plt, aes(x=LDscrs, y = chi2, color = weights), size = 2, shape = 19, alpha = alpha_) + 
    scale_color_gradientn(colours = rainbow(7)) +
    geom_abline(data = dfln, aes(intercept = intercept, slope = slope))
  gp = gp + ggtitle(paste("Heritability:", round(unname(coef(ldsc)[2])*M/(N*(1-1/N)), 2))) + #,
    # ", nSNP = ", M, ", nSEQ =", N, sep = "")) +
    xlab("LD Score Bin") + ylab(expression("Mean"~""*chi^2*"")) + theme(text = element_text(size=15)) 
  return(gp)
}
# LDSC model
perform_ldsc = function(ph, acc = FALSE, fec = FALSE){
  t = Sys.time()
  if(ph=="cd" | ph == "acute"){
    if(acc){
      r2 = readRDS("computeR2/r2s_cdacute_acc.rds") # These are the r2 outputs
      core_snps = readRDS("computeR2/acc_cd.acute_core_mx_NC.rds") # core_genome
      acc_genes = readRDS("computeR2/cd_acute_accessory_genome.rds")
      acc_genes = acc_gene_cleaner(acc_genes)
      mxCC = cbind(core_snps, acc_genes)
    } else{
      r2 = readRDS("computeR2/r2s_cdacute.rds") # These are the r2 outputs
      mxCC = readRDS("computeR2/acc_cd.acute_core_mx_NC.rds")
    }
    r2scores = r2$r2scores
    r2w = r2$r2weights
    pheno = readRDS("../Clustering/Results/FBCLUSTER_GLS/mapped_pheno_Wclust.rds") # load pheno
    # pos_original = readRDS("../Accessory_Genome/DATASET2/Pensar2019/acc_cd.acute_core_pos_NC.rds") #  load all pos in mxCC
    passoc = readRDS("RUN03/OUT/cd_passoc_gls_nocov.rds") # This is the GLS output, only to find retained pos
    if(!acc){
      # acc_genes = grep("group", passoc$tab$predictor) # idxes of acc_genes
      acc_genes = which(is.na(as.numeric(passoc$tab$predictor))) # Accessory gene names cannot be converted to numeric, catch them using NAs during coercion
      pos_retained = as.numeric(passoc$tab$predictor[-acc_genes]) # Retained SNPs
    } else {
      pos_retained = passoc$tab$predictor
    }
    names(r2scores) = pos_retained
    matches = which(colnames(mxCC) %in% pos_retained)
    mxCC = mxCC[,matches] # This is the reduced mxCC
  } else if(ph=="cef.mic" | ph == "pen.mic"){
    if(fec){
      if(ph == "cef.mic"){
        fecs = readRDS("FECs_cef.mic.rds")
      } else {
        fecs = readRDS("FECs_cef.mic.rds")
      }
    }
    if(acc){
      r2 = readRDS("computeR2/r2s_mic_acc.rds") # These are the r2 scores
      core_snps = readRDS("computeR2/acc_mic_core_mx_NC.rds") # core_genome
      acc_genes = readRDS("computeR2/pen.mic_cef.mic_accessory_genome.rds")
      acc_genes = acc_gene_cleaner(acc_genes)
      mxCC = cbind(core_snps, acc_genes)
    } else {
      r2 = readRDS("computeR2/r2s_mic.rds") # These are the r2 scores
      mxCC = readRDS("computeR2/acc_mic_core_mx_NC.rds")
    }
    r2scores = r2$r2scores
    r2w = r2$r2weights
    pheno = readRDS("../Clustering/Results/FBCLUSTER_GLS/mapped_pheno_Wclust_MIC.rds") # load pheno
    # pos_original = readRDS("../Accessory_Genome/DATASET2/Pensar2019/acc_mic_core_pos_NC.rds")
    passoc = readRDS("RUN03/OUT/pen.mic_passoc_gls_nocov.rds")
    if(!acc){
      # acc_genes = grep("group", passoc$tab$predictor) # idxes of acc_genes
      acc_genes = which(is.na(as.numeric(passoc$tab$predictor))) # Accessory gene names cannot be converted to numeric, catch them using NAs during coercion
      pos_retained = as.numeric(passoc$tab$predictor[-acc_genes]) # Retained SNPs
    } else {
      pos_retained = passoc$tab$predictor
    }
    names(r2scores) = pos_retained
    matches = which(colnames(mxCC) %in% pos_retained)
    mxCC = mxCC[,matches] # This is the reduced mxCC
  }
  w = 1/(r2w)
  w[which(w <= 0)] = 1
  N = nrow(pheno)
  M = ncol(mxCC)
  if(ph == "cd") {
    pheno_data = data.frame(ph = log10(pheno$carriage_duration))
  } else {
    if(fec){
     if(ph == "cef.mic") {
        pheno_data = data.frame(ph = log10(pheno$Ceftriaxone.MIC), fec = fecs) # Put the fixed effect covariates to reduce chi^2 values
      } else if(ph == "pen.mic") {
        pheno_data = data.frame(ph = log10(pheno$Penicillin.MIC), fec = fecs)
      }
    } else{
     if(ph == "cef.mic") {
        pheno_data = data.frame(ph = log10(pheno$Ceftriaxone.MIC))
      } else if(ph == "pen.mic") {
        pheno_data = data.frame(ph = log10(pheno$Penicillin.MIC))
      } 
    }
  }
  
  
  registerDoParallel(cores = 8L)
  chi2 = foreach(i = 1:M, .combine = "c") %dopar% {
    df_temp = data.frame(mxCC[,i], pheno_data)
    mdlx = lm(formula = 'ph ~ .', data = df_temp)
    summary(mdlx)$coeff[2,3]^2
  }
  df = data.frame(chi2 = chi2, LDscrs = r2scores, weights = w)
  ldsc = lm('chi2~LDscrs', df, weights = weights)
  print(Sys.time() - t)
  return(list(df=df, ldsc=ldsc, M=M, N=N))
}
############# cd ####################
ldsc_cd_acc = perform_ldsc("cd", acc = T)
plt_cd_acc = view_op("Log Carriage Duration (days)", ldsc = ldsc_cd_acc$ldsc, df = ldsc_cd_acc$df, M = ldsc_cd_acc$M, N = ldsc_cd_acc$N, nbin = 250)

ldsc_cd = perform_ldsc("cd", acc = F)
plt_cd = view_op("Log Carriage Duration (days)", ldsc = ldsc_cd$ldsc, df = ldsc_cd$df, M = ldsc_cd$M, N = ldsc_cd$N, nbin = 250)
ggsave("PLOTS/ldsc_cd.pdf", arrangeGrob(plt_cd, plt_cd_acc, nrow = 1), width = 12)

############## cef.mic #############
ldsc_cef.mic_acc = perform_ldsc("cef.mic", acc = T, fec = T)
plt_cef.mic_acc = view_op("Ceftriaxone MIC", ldsc = ldsc_cef.mic_acc$ldsc, df = ldsc_cef.mic_acc$df, M = ldsc_cef.mic_acc$M, N = ldsc_cef.mic_acc$N, nbin = 100)

ldsc_cef.mic = perform_ldsc("cef.mic", acc = F, fec = T)
plt_cef.mic = view_op("Ceftriaxone MIC", ldsc = ldsc_cef.mic$ldsc, df = ldsc_cef.mic$df, M = ldsc_cef.mic$M, N = ldsc_cef.mic$N, nbin = 100)
ggsave("PLOTS/ldsc_cefmic_wFec.pdf", arrangeGrob(plt_cef.mic, plt_cef.mic_acc, nrow = 1), width = 12)

ldsc_cef.mic_acc = perform_ldsc("cef.mic", acc = T, fec = F)
plt_cef.mic_acc = view_op("Ceftriaxone MIC", ldsc = ldsc_cef.mic_acc$ldsc, df = ldsc_cef.mic_acc$df, M = ldsc_cef.mic_acc$M, N = ldsc_cef.mic_acc$N, nbin = 100)

ldsc_cef.mic = perform_ldsc("cef.mic", acc = F, fec = F)
plt_cef.mic = view_op("Ceftriaxone MIC", ldsc = ldsc_cef.mic$ldsc, df = ldsc_cef.mic$df, M = ldsc_cef.mic$M, N = ldsc_cef.mic$N, nbin = 100)
ggsave("PLOTS/ldsc_cefmic.pdf", arrangeGrob(plt_cef.mic, plt_cef.mic_acc, nrow = 1), width = 12)


############## pen.mic #############
ldsc_pen.mic_acc = perform_ldsc("pen.mic", acc = T, fec = T)
plt_pen.mic_acc = view_op("Penicillin MIC", ldsc = ldsc_pen.mic_acc$ldsc, df = ldsc_pen.mic_acc$df, M = ldsc_pen.mic_acc$M, N = ldsc_pen.mic_acc$N, nbin = 100)

ldsc_pen.mic = perform_ldsc("pen.mic", acc = F, fec = T)
plt_pen.mic = view_op("Penicillin MIC", ldsc = ldsc_pen.mic$ldsc, df = ldsc_pen.mic$df, M = ldsc_pen.mic$M, N = ldsc_pen.mic$N, nbin = 100)
ggsave("PLOTS/ldsc_penmic_wFec.pdf", arrangeGrob(plt_pen.mic, plt_pen.mic_acc, nrow = 1), width = 12)

ldsc_pen.mic_acc = perform_ldsc("pen.mic", acc = T, fec = F)
plt_pen.mic_acc = view_op("Penicillin MIC", ldsc = ldsc_pen.mic_acc$ldsc, df = ldsc_pen.mic_acc$df, M = ldsc_pen.mic_acc$M, N = ldsc_pen.mic_acc$N, nbin = 100)

ldsc_pen.mic = perform_ldsc("pen.mic", acc = F, fec = F)
plt_pen.mic = view_op("Penicillin MIC", ldsc = ldsc_pen.mic$ldsc, df = ldsc_pen.mic$df, M = ldsc_pen.mic$M, N = ldsc_pen.mic$N, nbin = 100)
ggsave("PLOTS/ldsc_penmic.pdf", arrangeGrob(plt_pen.mic, plt_pen.mic_acc, nrow = 1), width = 12)
