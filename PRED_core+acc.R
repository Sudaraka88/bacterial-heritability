if(rstudioapi::isAvailable()) setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # WORKING DIRECTORY
# This is a combined run with core and accessory genome data 20210430
library(ggplot2)
library(lme4qtl)
library(wlm)
library(gridExtra)
library(Rcpp)
library(MASS)
sourceCpp("../Genotype2/tableC.cpp")
# Gene mapping
source("map2gene.R")
gene_mapper = function(x){
  gene_list = c()
  snp_list = c()
  starts = c()
  ends = c()
  for(i in x){
    tmp = map2gene(i, gene_data)
    if(tmp$dist_from_gene < 1000) {
      gene_list = c(gene_list, tmp$gd$gene[1])
      starts = c(starts, tmp$gd$start[1])
      ends = c(ends, tmp$gd$end[1])
      snp_list = c(snp_list, i)
    }
  }
  gene_summary = data.frame(snp = snp_list,
                            gene = gene_list,
                            starts = starts,
                            ends = ends)
  
  # remove duplicated entries
  gene_summary = gene_summary[!duplicated(gene_summary$gene),]
  print(noquote(paste0(gene_summary$gene, ",", sep = "", collapse = NULL)))
  return(gene_summary)
}
mkplt_pred = function(train_df, test_df, pred_out, tp){
  df_train = data.frame(Observation = train_df$y, Prediction = pred_out$predict_train)
  line = lm('Prediction ~ Observation', df_train)$coefficients
  plt_train = ggplot(df_train) + geom_point(aes(x = Observation, y = Prediction)) + geom_abline(slope = line[2], intercept = line[1]) + 
    ggtitle(paste("Training: ", round((1-tp)*100,1),  "%, MAE: ", round(pred_out$mae_train, 2), ", cor: ", round(pred_out$cor_train,2), sep = "" )) + 
    theme(text = element_text(size=15)) 
  
  df_test = data.frame(Observation = test_df$y, Prediction = pred_out$predict_test)
  line = lm('Prediction ~ Observation', df_test)$coefficients
  plt_test = ggplot(df_test) + geom_point(aes(x = Observation, y = Prediction)) + geom_abline(slope = line[2], intercept = line[1]) + 
    ggtitle(paste("Testing: ", round(tp*100,1), "%, MAE: ", round(pred_out$mae_test, 2), ", cor: ", round(pred_out$cor_test,2), sep = "")) + 
    theme(text = element_text(size=15)) 
  
  return(list(plt_train=plt_train, plt_test=plt_test))
}
mkplt_gwas = function(passoc_gls, lab_pos = 0, pheno_type, sim_mx_type, shape_snps = NULL, offset = 2350000){
  # Increase font size and remove titles for pub 20210527
  # Plotting code
  lgpvs = -log10(passoc_gls$tab$pval)
  bonf_lines = data.frame(intercepts = c(-log10(0.05/length(passoc_gls$tab$predictor)),
                                         -log10(0.01/length(passoc_gls$tab$predictor))),
                          slopes = c(0,0),
                          Bonf = c("0.05","0.01"))
  
  # Core/acc
  coreacc = rep(TRUE, length(passoc_gls$tab$predictor))
  coreacc[which(is.na(as.numeric(passoc_gls$tab$predictor)))] = FALSE # 20210521 - isolate acc_genes using NA technique
  df.plot_lmm = data.frame(pos = as.numeric(passoc_gls$tab$predictor[coreacc]), pval = lgpvs[coreacc])
  df.plot_lmm = rbind(df.plot_lmm, data.frame(pos = seq(offset, by = 100, length = length(which(!coreacc))), pval = lgpvs[!coreacc]))
  df.plot_lmm = cbind(df.plot_lmm, Type = as.factor(c(rep("core", length(which(coreacc))),rep("acc", length(which(!coreacc))))))
  
  if(!is.null(shape_snps)){
    covs = rep(0, length(lgpvs))
    covs[which(as.numeric(passoc_gls$tab$predictor[coreacc]) %in% shape_snps)] = 1
    df.plot_lmm = cbind(df.plot_lmm, fec = as.factor(covs))
  }
  
  
  gw_sig_idx = which(df.plot_lmm$pval[coreacc] > bonf_lines$intercepts[2])
  gene_summary = gene_mapper(as.numeric(df.plot_lmm$pos[gw_sig_idx]))
  # The plotting needs to be done now
  if(is.null(lab_pos)){
    ymin = 0
  } else{
    if(length(lab_pos) == 1) {
      if(lab_pos == 0) {
        if(nrow(gene_summary) > 0){
          lab_pos = rep(-2, nrow(gene_summary))
          print("Warning! lab_pos must be set") 
        } else {
          lab_pos = c()
          print("No associated genes were identified!")
        }
      }
    }
    if(length(lab_pos) > 0) {
      ymin = min(lab_pos)-2
    } else {
      ymin = 0
    }
  }
  ymax = max(c(lgpvs, bonf_lines$intercepts))
  if(is.null(shape_snps)){
    plt_mic = ggplot(df.plot_lmm) + geom_point(aes(x = pos, y = pval, col = Type))
  } else {
    plt_mic = ggplot() + geom_point(data = df.plot_lmm[which(df.plot_lmm$fec==0),], mapping = aes(x = pos, y = pval, col = Type)) +
      geom_point(df.plot_lmm[which(df.plot_lmm$fec==1),], mapping = aes(x = pos, y = pval),col = "red", shape = 24, fill = "red", size = 2)
  }
  plt_mic = plt_mic + geom_abline(data = bonf_lines, aes(intercept = intercepts, slope = slopes, linetype = Bonf))
  if(!is.null(lab_pos)){
    if(length(nrow(gene_summary) > 0)){
      plt_mic = plt_mic + annotate("text", x =  gene_summary$starts, y = lab_pos, label = gene_summary$gene, angle = 90) 
    } 
  }
  plt_mic = plt_mic + xlab("Basepair position") + ylab("-10log(p-value)") + ylim(ymin, ymax) + theme(text = element_text(size=15)) 
  # ggtitle(paste(pheno_type, ", lme4qtl, ", sim_mx_type, sep = ""))
  plt_mic
  return(plt_mic)
}
checkAlignment = function(nm1,nm2){
  print(paste("Pass alignment test?", all(nm1 == nm2)))
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
rename_simmx = function(nm){
  x = unlist(strsplit(nm,"_"))
  return(paste(x[1],"_",x[2],"#",x[3],sep = ""))
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
getVariation = function(X){
  u = unique(X)
  if(length(u) == 1) return(0) else return(length(u))
  # 0 if no gene variations, else # of variations returned
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
subset_snps = function(passoc_gls, pheno_data_train, mxCC_train, cor_thresh = 0.2, retain_thresh = 0.3, AIC_thresh = -1e5){
  # This is a temporary hack to overwrite cor_thresh and retain_thresh 20210524
  # cor_thresh = 0.4
  # retain_thresh = 0.5
  ############################################################################
  t = Sys.time()
  print(paste("Initiating LS pipelines at", t))
  # Let's clump and select a set of SNPs based on p-value order to fit a simple Linear Model
  snporder = order(-log10(passoc_gls$tab$pval), decreasing = T) # snporder[1] is the most important SNP
  keep = keep_temp = snporder[1]
  # cor_thresh = 0.2 # correlation
  # NN = round(0.3*length(snporder)) # Roughly the top 30% of SNPs only
  NN = round(retain_thresh*length(snporder)) # Roughly the top 30% of SNPs only
  # lm_temp = lm('y~.', data.frame(pheno_data_train[,-1], X = mxCC_train[,keep]))
  lm_temp = lm(y~., data.frame(pheno_data_train, X = mxCC_train[,keep]))
  AIC0 = AIC(lm_temp)
  for(i in 2:NN){
    # i = i + 1
    idx = snporder[i]
    # cor(mxCC_clump[,idx], mxCC_clump[,rem])^2
    # any(cor(mxCC_clump[,idx], mxCC_clump[,rem])^2>thresh)
    if( all(abs(cor(mxCC_train[,idx], mxCC_train[,keep])) < cor_thresh) ) {
      keep_temp = c(keep, idx)
      # lm_temp = lm('y~.', data.frame(pheno_data_train[,-1], X = mxCC_train[,keep_temp]))
      lm_temp = lm(y~., data.frame(pheno_data_train, X = mxCC_train[,keep_temp]))
      AIC_t = AIC(lm_temp)
      if(AIC_t < AIC0){ # model improved
        print(paste("Retained", length(keep_temp), "SNPs with AIC:", AIC_t))
        keep = keep_temp
        AIC0 = AIC_t
        if(AIC0 < AIC_thresh) break
      }
    }
    # if(i >= length(snporder)) {
    #   print(paste("Retained", length(which(markers==TRUE)), "SNPs"))
    #   break
    # }
  }
  return(list(keep = keep, lm = lm_temp))
}
trainer = function(formi, mxCC_train, W){
  passoc_gls <- matlm::matlm(formi$formi, cbind(ids = rownames(formi$pheno_data), formi$pheno_data), pred = mxCC_train, ids = ids, 
                             transform = W, batch_size = 1000, verbose = 2, stats_full = TRUE)
  print("done!")
  subset = subset_snps(passoc_gls = passoc_gls, pheno_data_train = formi$pheno_data, mxCC_train = mxCC_train)
  return(list(keep = subset$keep, passoc_gls = passoc_gls))
}
predictor = function(keep, pheno_data_train, pheno_data_test, training_idx, testing_idx, mxCC, W){
  # predict_nocov = predictor(train_nocov$keep, formi_nocov$pheno_data, create_formula(pheno_data_test[,-1])$pheno_data, training_idx, testing_idx, mxCC)
  train_df = data.frame(pheno_data_train, mxCC[training_idx,keep])
  # gls_reduced = lm.gls(y~., train_df, W)
  # predict_train = predict.lm(gls_reduced, train_df)
  
  lm_reduced = lm(y~., train_df)
  predict_train = predict(lm_reduced, train_df)
  mae_train = mean(abs(predict_train - train_df$y))
  cor_train = cor(predict_train, train_df$y)
  
  test_df = data.frame(pheno_data_test, mxCC[testing_idx,keep])
  # predict_test = predict.lm(gls_reduced, test_df)
  
  predict_test =predict(lm_reduced, test_df)
  mae_test = mean(abs(predict_test - test_df$y))
  cor_test = cor(predict_test, test_df$y)
  return(list(lm_reduced = lm_reduced, predict_train = predict_train, mae_train = mae_train, cor_train = cor_train,
              predict_test = predict_test, mae_test = mae_test, cor_test = cor_test))
  # return(list(lm_reduced = gls_reduced, predict_train = predict_train, mae_train = mae_train, cor_train = cor_train,
  #             predict_test = predict_test, mae_test = mae_test, cor_test = cor_test))
}

validation_method_ = c("xfcv") #, "loso")
ph_ = c("cef.mic") # "pen.mic", "cd"

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
      
      print("Reading acc_genes: ../Accessory_Genome/DATASET2/cd_acute_accessory_genome.rds")
      acc_genes = readRDS("../Accessory_Genome/DATASET2/cd_acute_accessory_genome.rds")
      checkAlignment(pheno$sampleID, rownames(acc_genes))
      acc_genes = acc_gene_cleaner(acc_genes)
      
      print("Reading core_snps: ../Accessory_Genome/DATASET2/Pensar2019/acc_cd.acute_core_mx_NC.rds")
      core_snps = readRDS("../Accessory_Genome/DATASET2/Pensar2019/acc_cd.acute_core_mx_NC.rds")
      checkAlignment(pheno$sampleID, rownames(core_snps))
      mxCC = cbind(core_snps, acc_genes)
      rm(core_snps, acc_genes); gc()
      
      mxCC = filter_alt(mxCC)
      # Let's make the phenotype table
      pheno_data = data.frame(ids = pheno$sampleID, y = log10(pheno$carriage_duration), carried = pheno$carried)
      formi = create_formula(pheno_data[,-1])$formi
      
      print("Reading similarity mx: ../pySEER/DATASET2/phylosim_cdacute.tsv")
      sim_mx = read.table("../pySEER/DATASET2/phylosim_cdacute.tsv") # These row names/col names need to be modified
    } else if(ph == "cef.mic" | ph == "pen.mic"){
      pheno = readRDS("../Accessory_Genome/DATASET2/pen.mic_cef.mic_pheno_reduced.rds")
      
      acc_genes = readRDS("../Accessory_Genome/DATASET2/pen.mic_cef.mic_accessory_genome.rds")
      checkAlignment(pheno$sampleID, rownames(acc_genes))
      acc_genes = acc_gene_cleaner(acc_genes)
      
      core_snps = readRDS("../Accessory_Genome/DATASET2/Pensar2019/acc_mic_core_mx_NC.rds")
      checkAlignment(pheno$sampleID, rownames(core_snps))
      mxCC = cbind(core_snps, acc_genes)
      rm(core_snps, acc_genes); gc()
      
      mxCC = filter_alt(mxCC)
      # Let's make the phenotype table
      if(ph == "cef.mic") {
        fec_cv = readRDS("FECs_cef.mic.rds")
        shape_snps = as.numeric(colnames(fec_cv))
        colnames(fec_cv) = paste("f_", colnames(fec_cv), sep = "")
        pheno_data = data.frame(ids = pheno$sampleID, y = log10(pheno$Ceftriaxone.MIC), acute = as.numeric(pheno$acute=="No"), 
                                cat = as.numeric(pheno$category=="Infant"), fec_cv) 
        lab_pos = list(nocov = c(-2,-6,-10,-14,-18,-2,-6,-10,-14,-2,-2,-2),
                       fecov = c(-2,-2, -2))
      } else if(ph == "pen.mic") {
        fec_cv = readRDS("FECs_pen.mic.rds")
        shape_snps = as.numeric(colnames(fec_cv))
        colnames(fec_cv) = paste("f_", colnames(fec_cv), sep = "")
        pheno_data = data.frame(ids = pheno$sampleID, y = log10(pheno$Penicillin.MIC), acute = as.numeric(pheno$acute=="No"), 
                                cat = as.numeric(pheno$category=="Infant"), fec_cv)
        lab_pos = list(nocov = c(-2,-6,-10, -14, -2,-6,-10, -14),
                       fecov = c(-2, -2))
      }
      formi = list(nocov = create_formula(pheno_data[,-1], type = "nocov")$formi,
                   allcov = create_formula(pheno_data[,-1], type = "all_cov")$formi)
      
      print("Reading similarity mx: ../pySEER/DATASET2/phylosim_mic.tsv")
      sim_mx = read.table("../pySEER/DATASET2/phylosim_mic.tsv") # These row names/col names need to be modified
    }
    
    n_acc = length(grep("group", colnames(mxCC)))
    n_core = ncol(mxCC) - n_acc
    print(paste("After filtering, nSEQ = ", nrow(mxCC), ", nSNP_core = ", n_core, ", nGenes_acc = ", n_acc, " in full dataset", sep = ""))
    outpath = "OUT"; if(!file.exists(outpath)) dir.create(outpath) # Outputs
    plotpath = "PLOTS"; if(!file.exists(plotpath)) dir.create(plotpath) # Plots
    if(validation_method == "xfcv"){
      print("Commencing xfcv")
      # Code for 10-fold cross validation
      outpath = file.path(outpath,"xfcv"); if(!file.exists(outpath)) dir.create(file.path(outpath))
      plotpath = file.path(plotpath, "xfcv"); if(!file.exists(plotpath)) dir.create(file.path(plotpath))
    
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
      plotpath = file.path(plotpath,"loso"); if(!file.exists(plotpath)) dir.create(file.path(plotpath))
      
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
    
    layout_mx = rbind(c(3,3),c(1,2)) # This is for arrangegrob
    
    xxx = 1 # counter
    for(testing_idx in testing_idx_){
      t_val_step = Sys.time()
      if(!file.exists(file.path(outpath, paste(ph, "_", validation_method, "_s_", xxx, sep = "")))){
        print(paste("Working on step", xxx))
        # tp = 0.5 # training portion
        tp = length(testing_idx)/nrow(pheno)
        training_idx = seq(1, nrow(pheno))
        # training_idx = sort(sample(nrow(pheno), round(nrow(pheno)*tp)))
        training_idx = training_idx[-testing_idx]
        # pheno_train = pheno[training_idx,]
        # We should do a re-filtering of the training SNPs, we can't predict using unseen data!
        mxCC_train = mxCC[training_idx,]
        mxCC_train = filter_alt(mxCC_train)   # Filtering Pipeline 
        n_acc_ = length(grep("group", colnames(mxCC_train)))
        n_core_ = ncol(mxCC_train) - n_acc
        print(paste("After filtering, nSEQ = ", nrow(mxCC_train), ", nSNP_core = ", n_core_, ", nGenes_acc = ", n_acc_, " in training dataset", sep = ""))
        
        # t = Sys.time()
        # print(paste("Initiating GLS pipelines at", t))
        
        pheno_data_train = pheno_data[training_idx, ]
        pheno_data_test = pheno_data[testing_idx, ]
        sim_mx_train = sim_mx_cleaner(sim_mx, pheno_data_train$ids)
        
        mod <- relmatLmer(y ~ (1|ids), pheno_data_train, relmat = list(ids = sim_mx_train))
        print(VarProp(mod))
        V <- varcov(mod, idvar = "ids")
        # Matrix::image(V[1:20, 1:20], main = "Estimated V (with artifacts)") # some artifacts close to zero due to limited numeric precision
        # Matrix::image(V_thr[1:20, 1:20], main = "Estimated V (with artifacts removed)")
        decomp <- decompose_varcov(V, method = "evd", output = "all")
        W <- decomp$transform
        
        
        analysis = "nocov"
        formi_nocov = create_formula(pheno_data_train[,-1], analysis) # first column is ids
        
        train_nocov = trainer(formi_nocov, mxCC_train, W)
        predict_nocov = predictor(train_nocov$keep, formi_nocov$pheno_data, create_formula(pheno_data_test[,-1])$pheno_data, training_idx, testing_idx, mxCC, W)
        
        predplt_nocov = mkplt_pred(formi_nocov$pheno_data, create_formula(pheno_data_test[,-1], analysis)$pheno_data, predict_nocov, tp)
        gwas_nocov = mkplt_gwas(passoc_gls =train_nocov$passoc_gls, pheno_type = ph, sim_mx_type = paste("pred:", length(train_nocov$keep), "of", ncol(mxCC_train), "SNPs used"),
                                lab_pos = NULL,  shape_snps =   colnames(mxCC_train)[train_nocov$keep])
        gwas_nocov
        nocov_out = list(formi = formi_nocov, train = train_nocov, predict = predict_nocov, pred_plt = predplt_nocov, gwas_plot = gwas_nocov)
        ggsave(file.path(plotpath, paste(ph, "_nocov_", validation_method, "_s_", xxx, ".png", sep = "")), 
               arrangeGrob(predplt_nocov$plt_train, predplt_nocov$plt_test, gwas_nocov, layout_matrix = layout_mx), width = 12)
        
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
        
        
        if(ph == "cef.mic" | ph == "pen.mic") {
          analysis = "all_cov"
          formi_all_cov = create_formula(pheno_data_train[,-1], analysis) # first column is ids
          train_all_cov = trainer(formi_all_cov, mxCC_train, W)
          predict_all_cov = predictor(train_all_cov$keep, formi_all_cov$pheno_data, create_formula(pheno_data_test[,-1], analysis)$pheno_data, training_idx, testing_idx, mxCC, W)
          predplt_all_cov = mkplt_pred(formi_all_cov$pheno_data, create_formula(pheno_data_test[,-1], analysis)$pheno_data, predict_all_cov, tp)
          predplt_all_cov$plt_test
          gwas_all_cov = mkplt_gwas(passoc_gls = train_all_cov$passoc_gls, pheno_type = ph, sim_mx_type =  paste("pred:", length(train_all_cov$keep), "of", ncol(mxCC_train), "SNPs used"), 
                                    lab_pos = NULL, shape_snps = colnames(mxCC_train)[train_all_cov$keep])
          gwas_all_cov
          full_cov_out = list(formi = formi_all_cov, train = train_all_cov, predict = predict_all_cov, pred_plt = predplt_all_cov, gwas_plot = gwas_all_cov)
          save_op = list(nocov = nocov_out, full_cov = full_cov_out)
          ggsave(file.path(plotpath, paste(ph, "_all_cov_", validation_method, "_s_", xxx, ".png", sep = "")), 
                 arrangeGrob(predplt_all_cov$plt_train, predplt_all_cov$plt_test, gwas_all_cov, layout_matrix = layout_mx), width = 12)
        } else {
          save_op = list(nocov = nocov_out)
        }
        saveRDS(save_op, file.path(outpath, paste(ph, "_", validation_method, "_s_", xxx, sep = "")))
      } else {
        print(paste("File:", file.path(outpath, paste(ph, "_", validation_method, "_s_", xxx, sep = "")),"exists, skipping..."))
      }
      xxx = xxx + 1
      print(paste("Elapsed time for step:", Sys.time() - t_val_step))
    }
    print(paste("Elapsed time for validation:", Sys.time() - t_val))
  }
  print(paste("Elapsed time for phenotype:", Sys.time() - t_ph))
}
