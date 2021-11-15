if(rstudioapi::isAvailable()) setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # WORKING DIRECTORY
# library(ggplot2)
# library(lme4qtl)
# library(wlm)
# library(gridExtra)
# library(Rcpp)
# library(MASS)
library(glmnet)
library(doParallel)
registerDoParallel(cores = 6L)
# Rcpp::sourceCpp("tableC.cpp")

geth2est = function(snps, y, alpha, w, nfolds = 10){
  print("Estimating h^2 with glmnet")
  glmfit_cv = cv.glmnet(x = snps, y = y, alpha = alpha, weights = w, family = "gaussian", parallel = TRUE, nfolds = nfolds)
  h2 = glmfit_cv$glmnet.fit$dev.ratio[which(glmfit_cv$glmnet.fit$lambda == glmfit_cv$lambda.min)]
  print(paste("h2 =", h2))
  return(h2)
}

alpha = 0

TT = Sys.time()
out_full_full = data.frame()
folder_ = dir(pattern = "^[:h:][:0:]")

for(folder in folder_){
  if(file.exists(paste(folder, "/enet.rds", sep = ""))){
    print(paste("File", paste(folder, "/enet.rds", sep = ""), "exists"))
    out_full = readRDS(paste(folder, "/enet.rds", sep = ""))
  } else {
    t0 = Sys.time()
    out_full = data.frame()
    mxCC = readRDS(paste(folder, "/simulations/genSim/genSim_AFC.rds", sep = ""))

    N = nrow(mxCC)
    M = ncol(mxCC)
    i_ = length(dir(paste(folder, "/simulations/phenSim/", sep = ""))) - 1
    print(paste("Found", i_+1, "phenos"))
    
    for(i in 0:i_){
      print(paste("Pheno", i))
      if(file.exists(paste(folder, "/simulations/phenSim/",i,"/enet_out.rds", sep = ""))){
        print(paste("File", paste(folder, "/simulations/phenSim/",i,"/enet_out.rds", sep = ""), "exists"))
        out_df = readRDS(paste(folder, "/simulations/phenSim/",i,"/enet_out.rds", sep = ""))
      } else{
        y = read.table(paste(folder, "/simulations/phenSim/",i,"/phenSim.phen", sep = ""))$V3
        out_df = geth2est(snps = mxCC, y = y, alpha = alpha, w = rep(1,nrow(mxCC)))
        out_df = cbind(h2 = out_df, phen = i)
        
        saveRDS(out_df, paste(folder, "/simulations/phenSim/",i,"/enet_out.rds", sep = ""))
      }
      
      out_full = rbind(out_full, out_df)
    }
    out_full = cbind(out_full, folder)
    saveRDS(out_full, paste(folder, "/enet.rds", sep = ""))
    print(Sys.time() - t0)
  }
  
  out_full_full = rbind(out_full_full, out_full)
  
}
saveRDS(out_full_full, "enet_out.rds")
print(Sys.time() - TT)

