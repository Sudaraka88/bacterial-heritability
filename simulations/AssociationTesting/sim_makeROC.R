# this is the analysis of true values
if(rstudioapi::isAvailable()) setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # WORKING DIRECTORY
library(ggplot2)
library(pROC)

makeTWROC = function(tw_res){
  t_null = sort(tw_res$terminal$corr.sim, decreasing = T)
  t_sim = tw_res$terminal$corr.dat
  roc_t = c()
  for(i in t_null) roc_t = c(roc_t, length(which(t_sim > i)))
}

gap_freq = readRDS("../Accessory_Genome/DATASET2/gapDat/cd_core_uqepi_gap_freq.rds")

dirs_ = c("RUN03", "RUN04","RUN05")
gwas = data.frame()

# remove treeWAS from this
if(file.exists("ROC_gwas.rds")){
  gwas = readRDS("ROC_gwas.rds")
} else{
  for(dirs in dirs_){
    ph_ = dir(dirs)[grep(pattern = "^phe", dir(dirs))]
    if(length(grep(ph_, pattern = "pkl")) > 0) ph_ = ph_[!(ph_ %in%  ph_[grep(ph_, pattern = "pkl")])] # remove pkl
    
    resfiles_ = dir(dirs, pattern = "res_")
    for (ph in ph_){
      gwas_full_combi = readRDS(paste(dirs, "/OUT/combi/", ph, "_combi_gwas.rds", sep = ""))[,c(1,3)]
      gwas_full_max = readRDS(paste(dirs, "/OUT/max/", ph, "_combi_gwas.rds", sep = ""))[,c(1,3)]
      
      resfiles = resfiles_[grep(ph, resfiles_)]
      lmm_res = read.csv(file.path(dirs, resfiles[grep(paste(ph,"lmm",sep = "_"), resfiles)]), sep = "\t")
      lmm_res$variant = unname(sapply(lmm_res$variant, function(x) unlist(strsplit(x, "_"))[2] ))
      
    
      tt = readRDS(paste(dirs, "/OUT/combi/", ph, "_tt.rds", sep = "")) # max run
      
      lmm_res_ = lmm_res[which(lmm_res$variant %in% gwas_full_combi$pos),]
      gwas_full_combi_ = gwas_full_combi[which(gwas_full_combi$pos %in% lmm_res_$variant),]
      gwas_full_max_ = gwas_full_max[which(gwas_full_max$pos %in% lmm_res_$variant),]
      tt_ = tt[which(tt$pos %in% lmm_res_$variant),]
  
      tf = rep(FALSE, nrow(gwas_full_combi_))
      tf[which(gwas_full_combi_$pos %in% tt_$variant)] = TRUE
      
    
      gwas = rbind(gwas, cbind(p_max = -log10(gwas_full_max_$p),  p_combi = -log10(gwas_full_combi_$p), p_ma = -log10(lmm_res_$lrt.pvalue), tf))
    }
  }
  
  saveRDS(gwas, "ROC_gwas.rds")
}


# Let's make the treeWAS datasets now
dirs_ = c("RUN03", "RUN04","RUN05")
tw_op = data.frame()
for(dirs in dirs_){
  ph_ = dir(dirs)[grep(pattern = "^phe", dir(dirs))]
  if(length(grep(ph_, pattern = "pkl")) > 0) ph_ = ph_[!(ph_ %in%  ph_[grep(ph_, pattern = "pkl")])] # remove pkl
  for (ph in ph_){
    tw_res = readRDS(file.path("treeWAS_simulation", dirs, paste("tw_sig_", ph, ".rds", sep = "")))
     tt = readRDS(paste(dirs, "/OUT/combi/", ph, "_tt.rds", sep = ""))
    tf = rep(FALSE, length(tw_res$terminal$p.vals))
    tf[which(names(tw_res$terminal$corr.dat) %in% tt$variant)] = TRUE
     tw_op = rbind(tw_op,
                  cbind(ter = abs(tw_res$terminal$corr.dat), sim = abs(tw_res$simultaneous$corr.dat), sub = abs(tw_res$subsequent$corr.dat), tf))
    
  }
}


r_ter = roc(response = tw_op$tf, predictor = tw_op$ter, ci = T)
r_sim = roc(response = tw_op$tf, predictor = tw_op$sim, ci = T)
r_sub = roc(response = tw_op$tf, predictor = tw_op$sub, ci = T)


r1 = roc(response = gwas$tf, predictor = gwas$p_max, ci = T)
r2 = roc(response = gwas$tf, predictor = gwas$p_combi, ci = T)
r3 = roc(response = gwas$tf, predictor = gwas$p_ma, ci = T)

r_tw = c()

lenlen = min(length(r_ter$sensitivities), length(r_sim$sensitivities), length(r_sub$sensitivities))
r_tw$sensitivities = apply(cbind(r_ter$sensitivities[1:lenlen], r_sim$sensitivities[1:lenlen], r_sub$sensitivities[1:lenlen]), 1, max)
r_tw$specificities = r_sim$specificities[1:lenlen]

par(mar=c(1,1,1,1))
png("roc7.png", res = 200, height = 1500, width = 1500)

plot(x = 0:1, y = 0:1, xlim = c(1,0), ylim = c(0,1), type = "n", xlab = "Specificity", ylab = "Sensitivity",
     axes = T, cex.lab = 1.5, cex.axis = 1.2) +
  lines(x = r1$specificities, y = r1$sensitivities, col = "#F60000", lwd = 4) +
  lines(x = r2$specificities, y = r2$sensitivities, col = "#FF8C00", lwd = 4) +
  lines(x = r3$specificities, y = r3$sensitivities, col = "#FFEE00", lwd = 4) +
  lines(x = r_tw$specificities, y = r_tw$sensitivities, col = "#4815AA", lwd = 4) +
  abline(1,-1, col = "#6A6C6D")
legend("bottomright", legend=c("max: 0.91 ±0.01",
                               "combi: 0.89 ±0.01",
                               "MA: 0.85 ±0.01",
                               "treeWAS: 0.83 ±0.01"),
       col=c("#F60000","#FF8C00", "#FFEE00", "#4815AA"), lwd=4, title = "AUC", cex = 1.2) # "#4DE94C", "#3783FF", 
dev.off()


AUC_tw =round(sum(r_tw$specificities [1:lenlen]*diff(c(0, 1 - r_tw$sensitivities [1:lenlen]))),2)


round(c(r1$ci), 3)
round(c(r2$ci), 3)
round(c(r3$ci), 3)
round(c(r_ter$ci), 3)
round(c(r_sim$ci), 3)
round(c(r_sub$ci), 3)


library(caret)

confusionMatrix(data = as.factor(gwas$enet), reference = as.factor(gwas$tf), positive = "1")
confusionMatrix(data = as.factor(gwas$treeWAS), reference = as.factor(gwas$tf), positive = "1")
confusionMatrix(data = factor(gwas$p_max > -log10(0.05/1e5), labels = c("0", "1")), reference = as.factor(gwas$tf), positive = "1")
confusionMatrix(data = factor(gwas$p_combi > -log10(0.05/1e5), labels = c("0", "1")), reference = as.factor(gwas$tf), positive = "1")
confusionMatrix(data = factor(gwas$p_ma > -log10(0.05/1e5), labels = c("0", "1")), reference = as.factor(gwas$tf), positive = "1")
