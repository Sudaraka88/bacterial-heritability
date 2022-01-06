if(rstudioapi::isAvailable()) setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # WORKING DIRECTORY
# Checked 20220105

library(lme4qtl)
library(wlm)
library(Rcpp)
library(ggplot2)
library(gridExtra)
sourceCpp("tableC.cpp")

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
# Generate a GWAS Manhattan plot 
mkplt_gwas = function(passoc_gls, lab_pos = 0, pheno_type, sim_mx_type, shape_snps = NULL, offset = 2350000){
  # Plotting code
  lgpvs = -log10(passoc_gls$tab$pval)
  bonf_lines = data.frame(intercepts = c(-log10(0.05/length(passoc_gls$tab$predictor)),
                                         -log10(0.01/length(passoc_gls$tab$predictor))),
                          slopes = c(0,0),
                          Bonf = c("0.05","0.01"))
  
  # Core/acc
  coreacc = rep(TRUE, length(passoc_gls$tab$predictor))
  coreacc[which(is.na(as.numeric(passoc_gls$tab$predictor)))] = FALSE # isolate acc_genes using NA technique
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
  plt_mic = plt_mic + xlab("Basepair position") + ylab(expression("-log"[10]*"(p-value)")) + ylim(ymin, ymax) + theme(text = element_text(size=23)) 
    # ggtitle(paste(pheno_type, ", lme4qtl, ", sim_mx_type, sep = ""))
  plt_mic
  return(plt_mic)
}
checkAlignment = function(nm1,nm2){
  print(paste("Pass alignment test?", all(nm1 == nm2)))
}
# Drop non-snps and perform basic filtering
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
getVariation = function(X){
  u = unique(X)
  if(length(u) == 1) return(0) else return(length(u))
  # 0 if no gene variations, else # of variations returned
}
# format formula for lme4qtl
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
acc_gene_cleaner = function(acc_genes){
  # 95% rule 
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
# rename and arrange the phylo_sim matrix
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
  sim_mx = sim_mx[matches, matches]
  return(sim_mx)
  # Let's do a quick normalisation - this fails due to non-positive-definite nature
  # temp = matrix(nrow = nrow(sim_mx), ncol = ncol(sim_mx))
  # for(i in 1:nrow(sim_mx)){
  #   temp[i,] = sim_mx[i,]/sim_mx[i,i]
  #   temp[,i] = temp[i,]
  # }
  # rownames(temp) = rownames(sim_mx)
  # colnames(temp) = colnames(sim_mx)
  # return(temp)
}
# perform GWAS
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
      mod = lme4qtl::relmatLmer(y ~ (1|ids), pheno_data, relmat = list(ids = sim_mx)) # Full model
      print(lme4qtl::VarProp(mod))
      V <- lme4qtl::varcov(mod, idvar = "ids")
      decomp <- wlm::decompose_varcov(V, method = "evd", output = "all")
      W <- decomp$transform
      print(paste("Completed: elapsed", Sys.time() - t))
    } 
    print("Population structure matrix provided")
    passoc_gls <- matlm::matlm(formi, pheno_data, pred = mxCC, ids = ids, transform = W, batch_size = 1000, verbose = 2) # Approximate linear model
    saveRDS(passoc_gls, path) # save output 
  }
  return(list(W = W, passoc_gls = passoc_gls))
}

OUTPATH = "OUT"; dir.create(OUTPATH) # outputs will be saved to OUT
PLOTPATH = "PLOTS"; dir.create(PLOTPATH) # plots will be saved to PLOTS

ph_ = c("cd","cef.mic","pen.mic") # subset of phenotypes to analyse
# ph = "pen.mic"
for(ph in ph_){
  T0 = Sys.time()
  print(paste("Processing", ph))
  if(ph == "cd"){
    print("Reading phenotype: cd_pheno.rds")
    pheno = readRDS("cd_pheno.rds")
    
    print("Reading acc_genes: cd_acc.rds")
    acc_genes = readRDS("cd_acc.rds")
    checkAlignment(pheno$sampleID, rownames(acc_genes))
    acc_genes = acc_gene_cleaner(acc_genes)
    
    print("Reading core_snps: cd_core_mx_NC.rds")
    core_snps = readRDS("cd_core_mx_NC.rds")
    checkAlignment(pheno$sampleID, rownames(core_snps))
    mxCC = cbind(core_snps, acc_genes)
    rm(core_snps, acc_genes); gc()
    
    mxCC = filter_alt(mxCC)
    # Let's make the phenotype table
    pheno_data = data.frame(ids = pheno$sampleID, y = log10(pheno$carriage_duration), carried = pheno$carried)
    formi = create_formula(pheno_data[,-1])$formi # create formula compatible with lme4qtl
    
    print("Reading similarity mx: phylosim_cdacute.tsv")
    sim_mx = read.table("phylosim_cdacute.tsv") # These row names/col names need to be modified
    sim_mx = sim_mx_cleaner(sim_mx, pheno$sampleID)
  } else if(ph == "cef.mic" | ph == "pen.mic"){
    pheno = readRDS("mic_pheno.rds")
    
    acc_genes = readRDS("mic_acc.rds")
    checkAlignment(pheno$sampleID, rownames(acc_genes))
    acc_genes = acc_gene_cleaner(acc_genes)
    
    core_snps = readRDS("mic_core_mx_NC.rds")
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
      lab_pos = list(nocov = c(-2,-6,-10,-14,-18,-2,-6,-10,-14,-2,-2,-2), # arbitrary positions for identified Gene labels
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
    
    print("Reading similarity mx: phylosim_mic.tsv")
    sim_mx = read.table("phylosim_mic.tsv") # These row names/col names need to be modified
    sim_mx = sim_mx_cleaner(sim_mx, pheno$sampleID)
  }
  
  n_acc = length(grep("group", colnames(mxCC)))
  n_core = ncol(mxCC) - n_acc
  print(paste("After filtering, nSEQ = ", nrow(mxCC), ", nSNP_core = ", n_core, ", nGenes_acc = ", n_acc, sep = ""))
  
  # GWAS pipeline
  if(ph == "cd"){
    gwas_out_nocov = gwas_pipeline(path = file.path(OUTPATH, paste(ph, "_passoc_gls_nocov.rds", sep = "")), 
                                   formi = formi, pheno_data = pheno_data, sim_mx = sim_mx, mxCC = mxCC)
    # plt = mkplt_gwas(gwas_out_nocov$passoc_gls, lab_pos = c(0) , ph, "phylosim")
    plt = mkplt_gwas(gwas_out_nocov$passoc_gls, lab_pos = NULL , ph, "phylosim")
    plt
    ggsave(file.path(PLOTPATH, paste(ph, "_mhp_plot_nocov.png", sep = "")), width = 14, height = 9)
  } else if(ph == "cef.mic" | ph == "pen.mic"){
    # 2 runs neede for mic phenotypes
    gwas_out_nocov = gwas_pipeline(path = file.path(OUTPATH, paste(ph, "_passoc_gls_nocov.rds", sep = "")),
                                   formi = formi[[1]], pheno_data = pheno_data, sim_mx = sim_mx, mxCC = mxCC)
    # plt_nocov = mkplt_gwas(gwas_out_nocov$passoc_gls, lab_pos = lab_pos[[1]] , ph,  paste("phylosim,", length(shape_snps) ,"FEC cov"), shape_snps = shape_snps)
    plt_nocov = mkplt_gwas(gwas_out_nocov$passoc_gls, lab_pos = NULL , ph,  paste("phylosim,", length(shape_snps) ,"FEC cov"), shape_snps = shape_snps)
    plt_nocov
    ggsave(file.path(PLOTPATH, paste(ph, "_mhp_plot_nocov_.png", sep = "")), width = 14, height = 9)
    
    gwas_out_fecov = gwas_pipeline(path =  file.path(OUTPATH, paste(ph, "_passoc_gls_fecov.rds", sep = "")), 
                                   formi = formi[[2]], pheno_data = pheno_data, sim_mx = sim_mx, mxCC = mxCC, W = gwas_out_nocov$W)
    
    # plt_fecov = mkplt_gwas(gwas_out_fecov$passoc_gls, lab_pos = lab_pos[[2]] , ph, "phylosim")
    plt_fecov = mkplt_gwas(gwas_out_fecov$passoc_gls, lab_pos = NULL , ph, "phylosim")
    plt_fecov
    ggsave(file.path(PLOTPATH, paste(ph, "_mhp_plot_fecov.png", sep = "")), width = 14, height = 9)
    
    ggsave(file.path(PLOTPATH, paste(ph, "_mhp.png", sep = "")), arrangeGrob(plt_nocov, plt_fecov), height = 15, width = 12)
  }
  print(paste("Pipeline Time Elapsed:", Sys.time() - T0 ))
}
