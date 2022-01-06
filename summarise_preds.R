# This code is to summarise prediction results for each phenotype
if(rstudioapi::isAvailable()) setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # WORKING DIRECTORY
# Checked 20220105

library(ggplot2)
rm(list=ls())
viewop = function(bla, nm = "out") {
  print(paste(nm, ": ", round(mean(bla),3), " (", round(sd(bla),3), "), max:", max(bla), ", min:", min(bla), sep = ""))
}

ph_lab = c("Log carriage duration", "Log ceftriaxone MIC", "Log penicillin MIC")
ph_ = c("cd", "cef.mic", "pen.mic") # "cd", "pen.mic"
# fec = T
method_ = c("loso", "xfcv") # "loso""xfcv"

## For enet
fldr = "OUT"
a = 0

cnt = 0
for(ph in ph_){
  cnt = cnt + 1
  for(method in method_){
    mae = c()
    mse = c()
    cor = c()
    keep = c()
    tru = c()
    pred = c()
    clust = c()
    files = dir(paste(file.path(fldr,method)), pattern = paste(ph, "_", method, "_a_", a, "_", sep = "" ))
    
    for(i in files){
      op = readRDS(file.path(fldr, method, i))
      keep = c(keep, op$N)
      mae = c(mae, mean(abs(op$truth-op$preds)))
      cor = c(cor, op$cor)
      mse = c(mse, op$mse)
      tru = c(tru, op$truth)
      pred = c(pred, op$preds)
      clust = c(clust, rep(as.numeric(unlist(strsplit(i, "enet_"))[2]), length(op$truth)))
    }
    viewop(keep, "keep")
    viewop(mse, "mse")
    viewop(cor, "cor")
    
    
    
    if(method == "xfcv"){
      df_plt = data.frame(tru = tru, pred = pred, Step = as.factor(clust))
    } else {
      df_plt = data.frame(tru = tru, pred = pred, Cluster = as.factor(clust))
    }
    
    if(length(which(abs(df_plt$pred) > 100)) > 0) df_plt = df_plt[-which(abs(df_plt$pred) > 100), ]
    coef = lm(df_plt$pred ~ df_plt$tru)
    if(method == "xfcv") {
      plt = ggplot(df_plt) + geom_point(aes(x = tru, y = pred, color=Step), size = 0.7) + 
        geom_abline(slope = coef$coefficients[2], intercept = coef$coefficients[1]) + xlab(ph_lab[cnt]) + ylab("Prediction") + ylim(range(df_plt$pred))
    } else {
      plt = ggplot(df_plt) + geom_point(aes(x = tru, y = pred, color=Cluster), size = 0.7) + 
        geom_abline(slope = coef$coefficients[2], intercept = coef$coefficients[1]) + xlab(ph_lab[cnt]) + ylab("Prediction") + ylim(range(df_plt$pred))
    }
    ggsave(paste("PLOTS/", ph, "_", method,".png", sep = ""), plt, width = 7, height = 4)
  }
}
#############################################################################################################################################################
# # For MLR method (not recommended)
############################################################################################################################################################# 
ph = "pen.mic"
method = "loso"
files = dir(paste(file.path("RUN02/OUT",method)), pattern = ph )
fec = T
# Read in the pheno data as well

if(ph == "cd"){
  pheno = readRDS("cd_pheno.rds")
  pheno_data = data.frame(ids = pheno$sampleID, y = log10(pheno$carriage_duration), carried = pheno$carried)
} else if(ph == "cef.mic"){
  pheno = readRDS("mic_pheno.rds")
  pheno_data = data.frame(ids = pheno$sampleID, y = log10(pheno$Ceftriaxone.MIC), acute = as.numeric(pheno$acute=="No"),
                          cat = as.numeric(pheno$category=="Infant"))
} else if(ph == "cef.mic"){
  pheno = readRDS("mic_pheno.rds")
  pheno_data = data.frame(ids = pheno$sampleID, y = log10(pheno$Penicillin.MIC), acute = as.numeric(pheno$acute=="No"),
                          cat = as.numeric(pheno$category=="Infant"))
}

mae_test = c()
mse_test = c()
mae_train = c()
mse_train = c()
cor_test = c()
cor_train = c()
keep = c()
n_train = c()
n_test = c()

for(i in files){
  op = readRDS(file.path("OUT", method, i))
  if(fec){
    keep = c(keep, length(op$full_cov$train$keep))
    mae_test = c(mae_test, op$full_cov$predict$mae_test)
    cor_test = c(cor_test, op$full_cov$predict$cor_test)
    if(ph == "cd") mse_test = c(mse_test, mean((pheno_data$y[as.numeric(names(op$full_cov$predict$predict_test))] - op$full_cov$predict$predict_test)^2))
    else mse_test = c(mse_test, mean((pheno_data$y[match(names(op$full_cov$predict$predict_test), pheno_data$ids)] - op$full_cov$predict$predict_test)^2))

    mae_train = c(mae_train, op$full_cov$predict$mae_train)
    cor_train = c(cor_train, op$full_cov$predict$cor_train)
    if(ph == "cd") mse_train = c(mse_train, mean((pheno_data$y[as.numeric(names(op$full_cov$predict$predict_train))] - op$full_cov$predict$predict_train)^2))
    else mse_train = c(mse_train, mean((pheno_data$y[match(names(op$full_cov$predict$predict_train), pheno_data$ids)] - op$full_cov$predict$predict_train)^2))

  } else{
    keep = c(keep, length(op$nocov$train$keep))
    mae_test = c(mae_test, op$nocov$predict$mae_test)
    cor_test = c(cor_test, op$nocov$predict$cor_test)
    n_test = c(n_test, length(op$nocov$predict$predict_test))
    if(ph == "cd") mse_test = c(mse_test, mean((pheno_data$y[as.numeric(names(op$nocov$predict$predict_test))] - op$nocov$predict$predict_test)^2))
    else mse_test = c(mse_test, mean((pheno_data$y[match(names(op$nocov$predict$predict_test), pheno_data$ids)] - op$nocov$predict$predict_test)^2))
    # mse_test = c(mse_test, mean((pheno_data$y[as.numeric(names(op$nocov$predict$predict_test))] - op$nocov$predict$predict_test)^2))

    mae_train = c(mae_train, op$nocov$predict$mae_train)
    cor_train = c(cor_train, op$nocov$predict$cor_train)
    n_train = c(n_train, length(op$nocov$predict$predict_train))
    # mse_train = c(mse_train, mean((pheno_data$y[as.numeric(names(op$nocov$predict$predict_train))] - op$nocov$predict$predict_train)^2))
    if(ph == "cd") mse_train = c(mse_train, mean((pheno_data$y[as.numeric(names(op$nocov$predict$predict_train))] - op$nocov$predict$predict_train)^2))
    else mse_train = c(mse_train, mean((pheno_data$y[match(names(op$nocov$predict$predict_train), pheno_data$ids)] - op$nocov$predict$predict_train)^2))
  }
}

viewop(mse_test, "MSE (test)")
viewop(cor_test, "cor (test)")

viewop(mae_test, "MAE (test)")
viewop(cor_test, "cor (test)")
viewop(keep, 'keep')
viewop(mae_train, "MAE (train)")
viewop(cor_train, "cor (train)")





