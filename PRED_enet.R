if(rstudioapi::isAvailable()) setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # WORKING DIRECTORY
# This is a combined run with core and accessory genome data 20210430
library(ggplot2)
library(lme4qtl)
library(wlm)
library(gridExtra)
library(Rcpp)
library(MASS)
library(glmnet)
library(doParallel)
sourceCpp("../Genotype2/tableC.cpp")
validation_method_ = c("loso", "xfcv")
ph_ = c("cd","pen.mic", "cef.mic")#,"cef.mic", "cd")
getVariation = function(X){
  u = unique(X)
  if(length(u) == 1) return(0) else return(length(u))
  # 0 if no gene variations, else # of variations returned
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
gwas_pipeline = function(path, formi = NULL, pheno_data = NULL, sim_mx = NULL, mxCC = NULL, overwrite = FALSE, W = NULL){
  if(file.exists(path) && !overwrite){
    print("Loading saved model")
    passoc_gls = readRDS(path)
  } else {
    if(is.null(formi) | is.null(pheno_data) | is.null(sim_mx) | is.null(mxCC)) stop("No GLS at path or overwrite requested, data must be provided for GWAS.")
    print("Initiating GWAS pipeline")
    if(is.null(W)){
      print("Performing population structure correction")
      t = Sys.time()
      mod = lme4qtl::relmatLmer(y ~ (1|ids), pheno_data, relmat = list(ids = sim_mx))
      print(lme4qtl::VarProp(mod))
      V <- lme4qtl::varcov(mod, idvar = "ids")
      decomp <- wlm::decompose_varcov(V, method = "evd", output = "all")
      W <- decomp$transform
      print(paste("Completed: elapsed", Sys.time() - t))
    } 
    print("Population structure matrix provided")
    passoc_gls <- matlm::matlm(formi, pheno_data, pred = mxCC, ids = ids, transform = W, batch_size = 1000, verbose = 2)
    saveRDS(passoc_gls, path)
  }
  return(list(W = W, passoc_gls = passoc_gls))
}
create_formula = function(pheno_data, type = "nocov"){
  # Return the formula and the data frame for training
  cnms = colnames(pheno_data)
  # 1st column must always be y, if not delete it and warn
  if(cnms[1] != "y") { stop("First column of pheno_data must be y")}
  acc_vars = grep("a_", cnms)
  fec_vars = grep("f_", cnms)
  others = seq(1:ncol(pheno_data))
  if(length(acc_vars)>0 & length(fec_vars)>0) {
    others = others[-c(acc_vars, fec_vars)]
  } else if(length(acc_vars) > 0 & length(fec_vars)==0) {
    others = others[-acc_vars]
  } else if(length(acc_vars) == 0 & length(fec_vars) > 0) {
    others = others[-fec_vars]
  } else {
    warning("Pheno data contains no accessory or fixed effect covariates")
  }
  if(type=="nocov"){
    return(list(formi = DF2formula(pheno_data[others]), pheno_data = pheno_data[others]))
  } else if(type == "acc_cov"){
    return(list(formi = DF2formula(pheno_data[c(others, acc_vars)]), pheno_data = pheno_data[c(others, acc_vars)]))
  } else{
    return(list(formi = DF2formula(pheno_data), pheno_data = pheno_data))
  }
}
sim_mx_cleaner = function(sim_mx, pheno_order){
  print("Cleaning sim_mx")
  rnms = rownames(sim_mx)
  rnms = unname(sapply(rnms, function(x) rename_simmx(x)))
  rownames(sim_mx) = rnms
  colnames(sim_mx) = rnms
  sim_mx = as.matrix(sim_mx)
  # Now we need to arrange the sim_mx in order
  matches = c()
  for(i in pheno_order) matches[i] = which(rnms %in% i)
  checkAlignment(rownames(sim_mx)[matches], pheno_order)
  return(sim_mx[matches, matches])
}
rename_simmx = function(nm){
  x = unlist(strsplit(nm,"_"))
  return(paste(x[1],"_",x[2],"#",x[3],sep = ""))
}
checkAlignment = function(nm1,nm2){
  print(paste("Pass alignment test?", all(nm1 == nm2)))
}
subset_snps = function(passoc_gls, mxCC, cor_thresh = 0.5, retain_thresh = 0.75){
  t = Sys.time()
  print(paste("Initiating clumping pipeline at", t))
  # Let's clump and select a set of SNPs based on p-value order
  snporder = order(-log10(passoc_gls$tab$pval), decreasing = T) # snporder[1] is the most important SNP
  keep = snporder[1]
  NN = round(retain_thresh*length(snporder))
  for(i in 2:NN){
    print(i)
    idx = snporder[i]
    if( all(abs(stats::cor(mxCC[,idx], mxCC[,keep])) < cor_thresh) ) {
      keep = c(keep, idx)
    }
  }
  saveRDS(keep, paste("OUT/keep_", ph, ".rds", sep = ""))
  print(paste("Retained:", length(keep), "variants" ))
  return(keep)
}
alpha_ = c(0, 0.01)

registerDoParallel(2)
# ph = ph_[1]
# validation_method = validation_method_[1]
for(ph in ph_){
  t_ph = Sys.time() # Timer for phenotype
  for(validation_method in validation_method_){
    t_val = Sys.time() # Timer for validation
    # ph = "pen.mic" # "cd", "cef.mic"
    # validation_method = "xfcv" # "loso"
    
    print(paste("Now Processing", ph, "with", validation_method))
    if(ph == "cd"){
      print("Reading phenotype: ../Accessory_Genome/DATASET2/cd_acute_pheno_reduced.rds")
      pheno = readRDS("../Accessory_Genome/DATASET2/cd_acute_pheno_reduced.rds")
      w_vec = 1/table(pheno$cluster)
      w = rep(0, nrow(pheno))
      for(p in 1:length(w_vec)) w[which(pheno$cluster == names(w_vec)[p])] = unname(w_vec[p])
      
      print("Reading acc_genes: ../Accessory_Genome/DATASET2/cd_acute_accessory_genome.rds")
      acc_genes = readRDS("../Accessory_Genome/DATASET2/cd_acute_accessory_genome.rds")
      # checkAlignment(pheno$sampleID, rownames(acc_genes))
      acc_genes = acc_gene_cleaner(acc_genes)
      
      print("Reading core_snps: ../Accessory_Genome/DATASET2/Pensar2019/acc_cd.acute_core_mx_NC.rds")
      core_snps = readRDS("../Accessory_Genome/DATASET2/Pensar2019/acc_cd.acute_core_mx_NC.rds")
      # checkAlignment(pheno$sampleID, rownames(core_snps))
      mxCC = cbind(core_snps, acc_genes)
      rm(core_snps, acc_genes); gc()
      
      mxCC = filter_alt(mxCC)
      # Let's make the phenotype table
      pheno_data = data.frame(ids = pheno$sampleID, y = log10(pheno$carriage_duration), carried = pheno$carried, weights = w)
      # formi = create_formula(pheno_data[,-1])$formi
      # keep_idx = readRDS("OUT/keep_cd.rds")
      print("Reading similarity mx: ../pySEER/DATASET2/phylosim_cdacute.tsv")
      sim_mx = read.table("../pySEER/DATASET2/phylosim_cdacute.tsv") # These row names/col names need to be modified
    } else if(ph == "cef.mic" | ph == "pen.mic"){
      pheno = readRDS("../Accessory_Genome/DATASET2/pen.mic_cef.mic_pheno_reduced.rds")
      
      acc_genes = readRDS("../Accessory_Genome/DATASET2/pen.mic_cef.mic_accessory_genome.rds")
      # checkAlignment(pheno$sampleID, rownames(acc_genes))
      acc_genes = acc_gene_cleaner(acc_genes)
      
      core_snps = readRDS("../Accessory_Genome/DATASET2/Pensar2019/acc_mic_core_mx_NC.rds")
      # checkAlignment(pheno$sampleID, rownames(core_snps))
      mxCC = cbind(core_snps, acc_genes)
      rm(core_snps, acc_genes); gc()
      
      mxCC = filter_alt(mxCC)
      # Let's make the phenotype table
      if(ph == "cef.mic") {
        # fec_cv = readRDS("FECs_cef.mic.rds")
        # shape_snps = as.numeric(colnames(fec_cv))
        # colnames(fec_cv) = paste("f_", colnames(fec_cv), sep = "")
        pheno_data = data.frame(ids = pheno$sampleID, y = log10(pheno$Ceftriaxone.MIC), acute = as.numeric(pheno$acute=="No"), 
                                cat = as.numeric(pheno$category=="Infant")) 
        # keep_idx = readRDS("OUT/keep_cef.mic.rds")
        # lab_pos = list(nocov = c(-2,-6,-10,-14,-18,-2,-6,-10,-14,-2,-2,-2),
        #                fecov = c(-2,-2, -2))
      } else if(ph == "pen.mic") {
        # fec_cv = readRDS("FECs_pen.mic.rds")
        # shape_snps = as.numeric(colnames(fec_cv))
        # colnames(fec_cv) = paste("f_", colnames(fec_cv), sep = "")
        pheno_data = data.frame(ids = pheno$sampleID, y = log10(pheno$Penicillin.MIC), acute = as.numeric(pheno$acute=="No"), 
                                cat = as.numeric(pheno$category=="Infant"))
        # keep_idx = readRDS("OUT/keep_pen.mic.rds")
        # lab_pos = list(nocov = c(-2,-6,-10, -14, -2,-6,-10, -14),
        #                fecov = c(-2, -2))
      }
      # formi = list(nocov = create_formula(pheno_data[,-1], type = "nocov")$formi,
      #              allcov = create_formula(pheno_data[,-1], type = "all_cov")$formi)
      
      print("Reading similarity mx: ../pySEER/DATASET2/phylosim_mic.tsv")
      sim_mx = read.table("../pySEER/DATASET2/phylosim_mic.tsv") # These row names/col names need to be modified
    }
    # # Let's introduce a clumping pipeline here
    outpath = "OUT"; if(!file.exists(outpath)) dir.create(outpath) # Outputs
    
    if(file.exists(file.path("OUT", paste(ph, "_passoc_gls.rds", sep = "")))){
      passoc_gls = readRDS(file.path("OUT", paste(ph, "_passoc_gls.rds", sep = "")))
    } else {
      formi = create_formula(pheno_data[,-1])$formi
      sim_mx = sim_mx_cleaner(sim_mx, pheno_data$ids)
      gwas_out = gwas_pipeline(path =  file.path("OUT", paste(ph, "_passoc_gls.rds", sep = "")),
                               formi = formi, pheno_data = pheno_data, sim_mx = sim_mx, mxCC = mxCC)
      passoc_gls = gwas_out$passoc_gls
    }
    # 
    # snporder = order(-log10(passoc_gls$tab$pval), decreasing = T)
    # keep_idx = snporder[1:round(0.75*length(snporder))]
    # 
    if(file.exists(paste("OUT/keep_", ph, ".rds", sep = ""))){
      keep_idx = readRDS(paste("OUT/keep_", ph, ".rds", sep = ""))
    } else {
      keep_idx = subset_snps(passoc_gls, mxCC, cor_thresh = 0.5, retain_thresh = 0.75) # This is very slow!
    }
    
    
    # t_temp = Sys.time()
    # for (i in 0:20) {
    #   assign(paste("fit", i, sep=""), cv.glmnet(mxCC, pheno_data$y, type.measure="mse", 
    #                                             alpha=i/100,family="gaussian", parallel = T))
    #   print(paste(i, "done, T = ", Sys.time() - t_temp))
    # }
    
    #   plotpath = "PLOTS"; if(!file.exists(plotpath)) dir.create(plotpath) # Plots
    if(validation_method == "xfcv"){
      print("Commencing xfcv")
      # Code for 10-fold cross validation
      outpath = file.path(outpath,"xfcv"); if(!file.exists(outpath)) dir.create(file.path(outpath))
      # plotpath = file.path(plotpath, "xfcv"); if(!file.exists(plotpath)) dir.create(file.path(plotpath))
      
      set.seed(23)
      tp = 0.1
      x = 1
      y = round(nrow(pheno)*tp)
      testing_ord = sample(nrow(pheno), nrow(pheno))
      testing_idx_ = list()
      for(i in 1:round(1/tp)){
        testing_idx_[[i]] = testing_ord[x:y]
        x = y + 1
        y = round(y + (nrow(pheno)*tp))
        if(x > length(testing_ord)) break
        if(y > length(testing_ord) | i == ((round(1/tp))-1)) y = length(testing_ord)
        
      }
    } else if(validation_method == "loso"){
      print("Commencing loso")
      # Code for loso
      outpath = file.path(outpath,"loso"); if(!file.exists(outpath)) dir.create(file.path(outpath))
      # plotpath = file.path(plotpath,"loso"); if(!file.exists(plotpath)) dir.create(file.path(plotpath))
      
      clusts = sort(unique(pheno$cluster))
      print(paste("Identified", length(clusts), "Clusters"))
      testing_idx_ = list()
      idx = 1
      for(i in clusts){
        testing_idx_[[idx]] = which(pheno$cluster == i)
        idx = idx+1
      }
    }
    print(paste("Check for correct training/testing sets passed?", all(sort(unlist(testing_idx_)) == seq(1:nrow(pheno)))))
    #   
    #   layout_mx = rbind(c(3,3),c(1,2)) # This is for arrangegrob
    #   
    xxx = 1 # counter
    op = list()
    for(testing_idx in testing_idx_){
      t_val_step = Sys.time()
      # if(!file.exists(file.path(outpath, paste(ph, "_", validation_method, "_s_", xxx, sep = "")))){
      print(paste("Working on step", xxx))
      # tp = 0.5 # training portion
      tp = length(testing_idx)/nrow(pheno)
      training_idx = seq(1, nrow(pheno))
      # training_idx = sort(sample(nrow(pheno), round(nrow(pheno)*tp)))
      training_idx = training_idx[-testing_idx]
      # pheno_train = pheno[training_idx,]
      # We should do a re-filtering of the training SNPs, we can't predict using unseen data!
      mxCC_train = mxCC[training_idx,keep_idx]
      
      # mxCC_train = filter_alt(mxCC_train)   # Filtering Pipeline
      # n_acc_ = length(grep("group", colnames(mxCC_train)))
      # n_core_ = ncol(mxCC_train) - n_acc_
      # print(paste("After filtering, nSEQ = ", nrow(mxCC_train), ", nSNP_core = ", n_core_, ", nGenes_acc = ", n_acc_, " in training dataset", sep = ""))
      
      # t = Sys.time()
      # print(paste("Initiating GLS pipelines at", t))
      
      pheno_data_train = pheno_data[training_idx, ]
      pheno_data_test = pheno_data[testing_idx, ]
      # 10 fold CV using glmnet
      for (alpha in alpha_){
        gc()
        # alpha = 1
        # glmnet_cv = cv.glmnet(x = mxCC_train, y = pheno_data_train$y, weights = pheno_data_train$weights, alpha = alpha, family = "gaussian", parallel = TRUE)
        glmnet_cv = cv.glmnet(x = mxCC_train, y = pheno_data_train$y, alpha = alpha, family = "gaussian", parallel = TRUE)
        
        
        preds = predict(object = glmnet_cv, newx = mxCC[testing_idx, keep_idx], s = glmnet_cv$lambda.min)
        best_idx = which(glmnet_cv$lambda == glmnet_cv$lambda.min)
        save_op = list(glmnet = glmnet_cv, N = glmnet_cv$nzero[best_idx], preds = preds, truth = pheno_data_test$y, cor = cor(preds, pheno_data_test$y), mse = mean((preds - pheno_data_test$y)^2))
        saveRDS(save_op, file.path(outpath, paste(ph, "_", validation_method, "_a_",alpha, "_enet_", xxx, sep = "")))
        print(paste(validation_method,"step=", xxx,"for ph",ph,", alpha=",alpha,"mse=",save_op$mse,"cor=",save_op$cor,"N=",save_op$N))
        rm(glmnet_cv)
      }
      
      # plot(preds~pheno_data_test$y)
      
      # plot(fit.cv.glmnet)
      
      # best_idx = which(fit.cv.glmnet$lambda == fit.cv.glmnet$lambda.min)
      # print(paste(ph, ": MSE: ", fit.cv.glmnet$cvm[best_idx], 
      #             "SD: ", fit.cv.glmnet$cvsd[best_idx],
      #             "n_preds:", fit.cv.glmnet$nzero[best_idx]))
      
      # sim_mx_train = sim_mx_cleaner(sim_mx, pheno_data_train$ids)
      
      # mod <- relmatLmer(y ~ (1|ids), pheno_data_train, relmat = list(ids = sim_mx_train))
      # print(VarProp(mod))
      # V <- varcov(mod, idvar = "ids")
      # # Matrix::image(V[1:20, 1:20], main = "Estimated V (with artifacts)") # some artifacts close to zero due to limited numeric precision
      # # Matrix::image(V_thr[1:20, 1:20], main = "Estimated V (with artifacts removed)")
      # decomp <- decompose_varcov(V, method = "evd", output = "all")
      # W <- decomp$transform
      
      
      # analysis = "nocov"
      # formi_nocov = create_formula(pheno_data_train[,-1], analysis) # first column is ids
      # 
      # train_nocov = trainer(formi_nocov, mxCC_train, W)
      # predict_nocov = predictor(train_nocov$keep, formi_nocov$pheno_data, create_formula(pheno_data_test[,-1])$pheno_data, training_idx, testing_idx, mxCC, W)
      # 
      # predplt_nocov = mkplt_pred(formi_nocov$pheno_data, create_formula(pheno_data_test[,-1], analysis)$pheno_data, predict_nocov, tp)
      # gwas_nocov = mkplt_gwas(passoc_gls =train_nocov$passoc_gls, pheno_type = ph, sim_mx_type = paste("pred:", length(train_nocov$keep), "of", ncol(mxCC_train), "SNPs used"),
      #                         lab_pos = NULL,  shape_snps =   colnames(mxCC_train)[train_nocov$keep])
      # gwas_nocov
      # nocov_out = list(formi = formi_nocov, train = train_nocov, predict = predict_nocov, pred_plt = predplt_nocov, gwas_plot = gwas_nocov)
      # ggsave(file.path(plotpath, paste(ph, "_nocov_", validation_method, "_s_", xxx, ".png", sep = "")),
      #        arrangeGrob(predplt_nocov$plt_train, predplt_nocov$plt_test, gwas_nocov, layout_matrix = layout_mx), width = 12)
      
      # analysis = "acc_cov" # This is no longer a separate part!
      # formi_acc_cov = create_formula(pheno_data_train[,-1], analysis) # first column is ids
      # train_acc_cov = trainer(formi_acc_cov, mxCC_train, W)
      # predict_acc_cov = predictor(train_acc_cov$keep, formi_acc_cov$pheno_data, create_formula(pheno_data_test[,-1], analysis)$pheno_data, training_idx, testing_idx, mxCC)
      # predplt_acc_cov = mkplt_pred(formi_acc_cov$pheno_data, create_formula(pheno_data_test[,-1], analysis)$pheno_data, predict_acc_cov, tp)
      # # predplt_acc_cov$plt_test
      # gwas_acc_cov = mkplt_gwas(passoc_gls = train_acc_cov$passoc_gls, pheno_type = ph,
      #                           sim_mx_type =  paste("pred:", length(train_acc_cov$keep), "of", ncol(mxCC_train), "SNPs used"), color_snps = colnames(mxCC_train)[train_acc_cov$keep])
      # # gwas_acc_cov
      # acc_cov_out = list(formi = formi_acc_cov, train = train_acc_cov, predict = predict_acc_cov, pred_plt = predplt_acc_cov, gwas_plot = gwas_acc_cov)
      # ggsave(file.path(plotpath, paste(ph, "_acc_cov_", validation_method, "_s_", xxx, ".png", sep = "")),
      #        arrangeGrob(predplt_acc_cov$plt_train, predplt_acc_cov$plt_test, gwas_acc_cov, layout_matrix = layout_mx), width = 12)
      
      
      # if(ph == "cef.mic" | ph == "pen.mic") {
      #   analysis = "all_cov"
      #   formi_all_cov = create_formula(pheno_data_train[,-1], analysis) # first column is ids
      #   train_all_cov = trainer(formi_all_cov, mxCC_train, W)
      #   predict_all_cov = predictor(train_all_cov$keep, formi_all_cov$pheno_data, create_formula(pheno_data_test[,-1], analysis)$pheno_data, training_idx, testing_idx, mxCC, W)
      #   predplt_all_cov = mkplt_pred(formi_all_cov$pheno_data, create_formula(pheno_data_test[,-1], analysis)$pheno_data, predict_all_cov, tp)
      #   predplt_all_cov$plt_test
      #   gwas_all_cov = mkplt_gwas(passoc_gls = train_all_cov$passoc_gls, pheno_type = ph, sim_mx_type =  paste("pred:", length(train_all_cov$keep), "of", ncol(mxCC_train), "SNPs used"),
      #                             lab_pos = NULL, shape_snps = colnames(mxCC_train)[train_all_cov$keep])
      #   gwas_all_cov
      #   full_cov_out = list(formi = formi_all_cov, train = train_all_cov, predict = predict_all_cov, pred_plt = predplt_all_cov, gwas_plot = gwas_all_cov)
      #   save_op = list(nocov = nocov_out, full_cov = full_cov_out)
      #   ggsave(file.path(plotpath, paste(ph, "_all_cov_", validation_method, "_s_", xxx, ".png", sep = "")),
      #          arrangeGrob(predplt_all_cov$plt_train, predplt_all_cov$plt_test, gwas_all_cov, layout_matrix = layout_mx), width = 12)
      # } else {
      #   save_op = list(nocov = nocov_out)
      # }
      
      # } else {
      #   print(paste("File:", file.path(outpath, paste(ph, "_", validation_method, "_s_", xxx, sep = "")),"exists, skipping..."))
      # }
      xxx = xxx + 1
      print(paste("Elapsed time for step:", Sys.time() - t_val_step))
    }
    print(paste("Elapsed time for validation:", Sys.time() - t_val))
  }
  print(paste("Elapsed time for phenotype:", Sys.time() - t_ph))
}
