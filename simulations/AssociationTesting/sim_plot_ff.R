# this is the analysis of false values
if(rstudioapi::isAvailable()) setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # WORKING DIRECTORY
library(ggplot2)
gap_freq = readRDS("cd_core_uqepi_gap_freq.rds")
tt = data.frame()

fp_gapsnp = c()
fp_ma = c()
n_f_gapsnp = c()
n_f_ma = c()

dirs_ = c("RUN03", "RUN04", "RUN05")
for(dirs in dirs_){
  ph_ = dir(dirs)[grep(pattern = "^phe", dir(dirs))]
  if(length(grep(ph_, pattern = "pkl")) > 0) ph_ = ph_[!(ph_ %in%  ph_[grep(ph_, pattern = "pkl")])] # remove pkl
  
  for (ph in ph_){
    tt = rbind(tt, readRDS(paste(dirs, "/OUT/max/", ph, "_tt.rds", sep = "")))
    gwas = readRDS(paste(dirs, "/OUT/max/", ph, "_combi_gwas.rds", sep = ""))
    gwas = gwas[-which(gwas$pos %in% tt$variant), ]
    n_f_gapsnp = c(n_f_gapsnp, nrow(gwas))
    fp_gapsnp = c(fp_gapsnp, length(which(-log10(gwas$p) > -log10(0.05/100000))))
    
    resfiles = dir(dirs, pattern = "res_")
    resfiles = resfiles[grep(ph, resfiles)]
    lmm_res = read.csv(file.path(dirs, resfiles[grep(paste(ph,"lmm",sep = "_"), resfiles)]), sep = "\t")
    lmm_res$variant = unname(sapply(lmm_res$variant, function(x) unlist(strsplit(x, "_"))[2] ))
    lmm_res = lmm_res[-which(lmm_res$variant %in% tt$variant), ]
    n_f_ma = c(n_f_ma, nrow(lmm_res))
    
    fp_ma = c(fp_ma, length(which(-log10(lmm_res$lrt.pvalue) > -log10(0.05/100000)))) 
  }
}
plot(fp_ma, fp_gapsnp) + abline(0, 1)
