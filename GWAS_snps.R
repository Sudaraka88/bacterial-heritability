if(rstudioapi::isAvailable()) setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # WORKING DIRECTORY
# Try to perform GWAS using 2df analysis
library(lme4qtl)
library(wlm)
library(Rcpp)
library(ggplot2)
library(gridExtra)
library(lmtest)
library(fastDummies)
library(foreach)
library(doParallel)
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
    if(tmp$dist_from_gene < 2) {
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
  
  
  gw_sig_idx = which(df.plot_lmm$pval[coreacc] > bonf_lines$intercepts[1])
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
    plt_mic = ggplot(df.plot_lmm) + geom_point(aes(x = pos, y = pval, col = Type))  + scale_x_continuous(breaks = c(0, 500000, 1000000, 1500000, 2000000), name = "Basepair position", labels = scales::comma)
  } else {
    plt_mic = ggplot() + geom_point(data = df.plot_lmm[which(df.plot_lmm$fec==0),], mapping = aes(x = pos, y = pval, col = Type)) +
      geom_point(df.plot_lmm[which(df.plot_lmm$fec==1),], mapping = aes(x = pos, y = pval),col = "black", shape = 24, fill = "black", size = 3) +
      scale_x_continuous(breaks = c(0, 500000, 1000000, 1500000, 2000000), name = "Basepair position", labels = scales::comma)
  }
  plt_mic = plt_mic + geom_abline(data = bonf_lines, aes(intercept = intercepts, slope = slopes, linetype = Bonf))
  if(!is.null(lab_pos)){
    if(length(nrow(gene_summary) > 0)){
      plt_mic = plt_mic + annotate("text", x =  gene_summary$starts, y = lab_pos, label = gene_summary$gene, angle = 90) 
    } 
  }
  plt_mic = plt_mic + ylab(expression("-log"[10]*"(p-value)")) + ylim(ymin, ymax) + theme(text = element_text(size=16)) 
  # ggtitle(paste(pheno_type, ", lme4qtl, ", sim_mx_type, sep = ""))
  plt_mic
  return(plt_mic)
}
checkAlignment = function(nm1,nm2){
  print(paste("Pass alignment test?", all(nm1 == nm2)))
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
getVariation = function(X){
  u = unique(X)
  if(length(u) == 1) return(0) else return(length(u))
  # 0 if no gene variations, else # of variations returned
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
acc_gene_cleaner = function(acc_genes){
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
  sim_mx = sim_mx[matches, matches]
  return(sim_mx)
 
}
gwas_pipeline = function(ph, formi, pheno_data, sim_mx, mxCC){
  if(ph == "cd") scol = 3 else if(ph == "cef.mic" | ph == "pen.mic") scol = 3 else print("Unknown phenotype")
  print("Performing population structure correction")
  t = Sys.time()
  if(file.exists(paste("OUT/mod_", ph, ".rds", sep = ""))){
    mod_f = readRDS(paste("OUT/mod_", ph, ".rds", sep = ""))
    W = mod_f$W
    rm(mod_f)
  } else {
    if(ph == "cd"){
      mod = lme4qtl::relmatLmer(y ~ (1|ids) + (1|hostID), pheno_data, relmat = list(ids = sim_mx))
    } else {
      mod = lme4qtl::relmatLmer(y ~ (1|ids), pheno_data, relmat = list(ids = sim_mx))
    }
    
    print(lme4qtl::VarProp(mod_f$mod))
    V <- lme4qtl::varcov(mod, idvar = "ids")
    decomp <- wlm::decompose_varcov(V, method = "evd", output = "all")
    W <- decomp$transform
    saveRDS(list(mod = mod, V = V, W = W), paste("OUT/mod_", ph, ".rds", sep = ""))
  }
  print(paste("Completed: elapsed", Sys.time() - t))
  
  n = nrow(mxCC)
  y <- model.extract(model.frame(formi, pheno_data), "response")
  C <- model.matrix(formi, pheno_data)
  
  y_t <- crossprod(W, y)
  C_t <- crossprod(W, C)
  
  
  ######################### BELOW IS THE GWAS PIPELINE #########################
  print("Initiating GWAS pipeline")
  pos = as.numeric(colnames(mxCC))
  t = Sys.time()
  # mdl0 = lm('y~.', data.frame(y = y_t, c = C_t)) # only needed if df ncol (X) > 1
  doParallel::registerDoParallel(cores = 4L)
  if(!file.exists(paste("OUT/snpgwas_", ph, "_afc.rds", sep = ""))){
    # t = Sys.time()
    out = foreach(i = 1:ncol(mxCC), .combine='rbind') %dopar% {
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
        mdlx = lm('y~.', data.frame(y = y_t[-rm_idx], c = C_t[-rm_idx,], x = X))
      } else {
        gap_count = 0
        # y_t = y_t_d
        # C_t = C_t_d
        # binary code major allele
        # num_SNP = as.numeric(mxCC[,i] == names(tbl)[1])
        
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
        mdlx = lm('y~.', data.frame(y = y_t, c = C_t, x = X))
      }
      
      
      unname(c(summary(mdlx)$coefficients[scol,c(1,2,4)], tbl[1], gf , (n-gap_count)))
    }
    # Sys.time() - t
    out = data.frame(out)
    out = cbind(pos, out)
    
    colnames(out) = c("pos", "beta", "se", "p","MAF", "gf", "samples")
    saveRDS(out, paste("OUT/snpgwas_", ph, "_afc.rds", sep = ""))
  } else {
    out = readRDS(paste("OUT/snpgwas_", ph, "_afc.rds", sep = ""))
  }
  
  print(paste("Completed: elapsed", Sys.time() - t))
  
  bonf_lines = data.frame(intercepts = c(-log10(0.05/ncol(mxCC)),
                                         -log10(0.01/ncol(mxCC))),
                          slopes = c(0,0),
                          Bonf = c("0.05","0.01"))
  
  ymin = 0
  ymax = max(c(-log10(out$p), bonf_lines$intercepts))
  
  plt_dat = data.frame(pos = out$pos, p = -log10(out$p), gap_freq = out$gf)
  
  gw_sig_idx = which(plt_dat[,2] > bonf_lines$intercepts[1])
  gene_summary = gene_mapper(as.numeric(plt_dat$pos[gw_sig_idx]))
  
  plt = ggplot(plt_dat) + geom_point(aes(x = pos, y = p, col = gap_freq), shape = 20, size = 1) + ylim(ymin, ymax) + geom_abline(data = bonf_lines, aes(intercept = intercepts, slope = slopes, linetype = Bonf)) +
    scale_x_continuous(breaks = c(0, 500000, 1000000, 1500000, 2000000), name = "Basepair position", labels = scales::comma) +
    ylab(expression("-log"[10]*"(p-value)"))
  ggsave(paste("PLOTS/snpgwas_", ph, "_afc.png", sep = ""), plt, width = 12)
  ######################### ABOVE IS THE GWAS PIPELINE #########################
  
}

OUTPATH = "OUT"; dir.create(OUTPATH)
PLOTPATH = "PLOTS"; dir.create(PLOTPATH)

ph_ = c("cd","cef.mic", "pen.mic") # "cd", "cef.mic" 
# ph = "cd"
for(ph in ph_){
  gc()
  T0 = Sys.time()
  print(paste("Processing", ph))
  if(ph == "cd"){
    print("Reading phenotype: pheno_cd_uqepi_last.rds")
    pheno = readRDS("pheno_cd_uqepi_last.rds") # This is the full phenotype
    
    print("Reading snp_dat: cd_core_uqepi_mx.rds")
    
    mxCC = readRDS("cd_core_uqepi_mx.rds")
    mx_idxord = c()
    for(i in pheno$sampleID) mx_idxord = c(mx_idxord, which(rownames(mxCC) %in% i))
    mxCC = mxCC[mx_idxord, ]
    checkAlignment(pheno$sampleID, rownames(mxCC))
    
   
    pheno_data = data.frame(ids = pheno$sampleID, y = log10(pheno$carriage_duration), dummy = rep(1,nrow(pheno)), hostID = pheno$subject)
    formi = create_formula(pheno_data[,-c(1,4)])$formi
    
    print("Reading similarity mx: phylosim_cdacute.tsv")
    sim_mx = read.table("phylosim_cdacute.tsv") # These row names/col names need to be modified
    
    sim_mx = sim_mx_cleaner(sim_mx, pheno$sampleID)
  } else if(ph == "cef.mic" | ph == "pen.mic"){
    pheno = readRDS("pen.mic_cef.mic_pheno_reduced.rds")
    
    # new numerical coding 20210706
    print("Reading snp_dat: acc_mic_core_mx.rds")
    mxCC = readRDS("acc_mic_core_mx.rds")
      checkAlignment(pheno$sampleID, rownames(mxCC))
   
    
    # Let's make the phenotype table
    if(ph == "cef.mic") {
      pheno_data = data.frame(ids = pheno$sampleID, y = log10(pheno$Ceftriaxone.MIC), dummy = rep(1,nrow(pheno)))
      
    } else if(ph == "pen.mic") {
      pheno_data = data.frame(ids = pheno$sampleID, y = log10(pheno$Penicillin.MIC), dummy = rep(1,nrow(pheno)))
      
    }
   
    formi = create_formula(pheno_data[,-1])$formi
    
    print("Reading similarity mx: phylosim_mic.tsv")
    sim_mx = read.table("phylosim_mic.tsv") # These row names/col names need to be modified
   
    sim_mx = sim_mx_cleaner(sim_mx, pheno$sampleID)
  }
  
  n_acc = length(grep("group", colnames(mxCC)))
  n_core = ncol(mxCC) - n_acc
  print(paste("After filtering, nSEQ = ", nrow(mxCC), ", nSNP_core = ", n_core, ", nGenes_acc = ", n_acc, sep = ""))
  
  gwas_pipeline(ph, formi, pheno_data, sim_mx, mxCC)
  print(paste("Pipeline Time Elapsed:", Sys.time() - T0 ))
}
