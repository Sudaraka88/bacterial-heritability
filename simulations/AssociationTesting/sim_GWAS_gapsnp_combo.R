if(rstudioapi::isAvailable()) setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # WORKING DIRECTORY
library(survcomp)
library(ggplot2)
# Attempt to combine z-scores and get a single MHP
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
  hits = nrow(gene_summary)
  print(sort(table(gene_summary$gene), decreasing = T))
  gene_summary = gene_summary[!duplicated(gene_summary$gene),]
  print(paste(hits, "of", length(x),"hits in", nrow(gene_summary), "genes identified"))
  # print(noquote(paste0(gene_summary$gene, ",", sep = "", collapse = NULL)))
  return(gene_summary)
}
rename_simmx = function(nm){
  x = unlist(strsplit(nm,"_"))
  return(paste(x[1],"_",x[2],"#",x[3],sep = ""))
}
usecombi = T
noplots = T

ph_ = dir()[grep(pattern = "^phe", dir())]
gap_freq = readRDS("cd_core_uqepi_gap_freq.rds")
for (ph in ph_){
  pheno = read.table(paste(ph, sep = ""), header = T) # This is the full phenotype
  pheno$samples = sapply(pheno$samples, rename_simmx)
  
  gap_out =  read.csv(paste("res_",ph,"_gap_lmm_ps.txt", sep = ""), sep = "\t")
  gap_out$variant = unname(sapply(gap_out$variant, function(x) unlist(strsplit(x, "_"))[2] ))
  gf = gap_freq[which(as.numeric(names(gap_freq)) %in% gap_out$variant)]
  if(length(which(!gap_out$variant %in% names(gf))) > 0) {
    gap_out = gap_out[-which(!gap_out$variant %in% names(gf)), ]
  }
  gap_out = cbind(gap_out, gf = gf)
  gap_out = gap_out[,-c(2,3,8)]
  colnames(gap_out) = c("pos","p", "beta","se", "h2", "gf")
  gap_out$pos = as.numeric(gap_out$pos)
  
  snp_out = readRDS(file.path("OUT", paste("snpgwas_", ph, "_afc.rds", sep = "")))
  gap_out = cbind(gap_out, z = gap_out$beta/gap_out$se)
  snp_out = cbind(snp_out, z = snp_out$beta/snp_out$se)
  
  combi = match(gap_out$pos, snp_out$pos)
  combi_idx = gap_out$pos %in% snp_out$pos
  
  paste("combi check passed?", all(gap_out$pos[combi_idx] == snp_out$pos[combi[combi_idx]]))
  if(usecombi){
    gwas_full = data.frame(pos = gap_out$pos[combi_idx], z = (gap_out$z[combi_idx] + snp_out$z[combi[combi_idx]])/sqrt(2),
                           p = apply(cbind(gap_out$p[combi_idx], snp_out$p[combi[combi_idx]]), MARGIN = 1, function(x) survcomp::combine.test(p = x,method = "z.transform")),
                           gf = gap_out$gf[combi_idx], Test = rep("combi", length(which(combi_idx))))
  } else {
    gwas_full = data.frame(pos = gap_out$pos[combi_idx], z = (gap_out$z[combi_idx] + snp_out$z[combi[combi_idx]])/sqrt(2),
                           p = apply(cbind(gap_out$p[combi_idx], snp_out$p[combi[combi_idx]]), 1, min),
                           gf = gap_out$gf[combi_idx], 
                           Test = apply(cbind(gap_out$p[combi_idx], snp_out$p[combi[combi_idx]]), 1, function(x) ifelse(x[1] >= x[2], 'gap', 'snp')))
  }
  
  gap_only_idx = seq(1,nrow(gap_out))[-which(gap_out$pos %in% gwas_full$pos)]
  snp_only_idx = seq(1,nrow(snp_out))[-which(snp_out$pos %in% gwas_full$pos)]
  
  gwas_full = rbind(gwas_full,
                    data.frame(pos = gap_out$pos[gap_only_idx], z = gap_out$z[gap_only_idx],
                               p = gap_out$p[gap_only_idx], gf = gap_out$gf[gap_only_idx], Test = rep("gap", length(gap_only_idx))),
                    data.frame(pos = snp_out$pos[snp_only_idx], z = snp_out$z[snp_only_idx],
                               p = snp_out$p[snp_only_idx], gf = snp_out$gf[snp_only_idx], Test = rep("snp", length(snp_only_idx)))
  )
  
  
  gwas_full = gwas_full[order(gwas_full$pos),]
  
  plt_dat = data.frame(pos = gwas_full$pos, p = -log10(gwas_full$p), Gap_frq = gwas_full$gf, Test = as.factor(gwas_full$Test))
  
  # let's segment plt_dat based on gf
  wgap_idx = sort(c(which(plt_dat$Gap_frq > (nrow(pheno)-10)/nrow(pheno) & plt_dat$Test == "gap"),
                    which(plt_dat$Gap_frq < 10/nrow(pheno) & plt_dat$Test == "gap"))) # these are low gaps and weak tests
  strong_idx = which(!(seq(1, nrow(plt_dat)) %in% wgap_idx))
  
  rm_gap_idx =  sort(c(which(plt_dat$Gap_frq < 10/nrow(pheno) & plt_dat$Test != "snp"), 
                       which(plt_dat$Gap_frq > (nrow(pheno)-10)/nrow(pheno) & plt_dat$Test != "snp"))) # these are poor quality data, must be removed for bonf. level + qqplot
  kp_gap_idx = which(!(seq(1, nrow(plt_dat)) %in% rm_gap_idx))
  
  test_lines = data.frame(intercepts = c(-log10(0.05/nrow(gwas_full)),
                                         -log10(0.01/nrow(gwas_full))),
                          slopes = c(0,0),
                          FWER = c("0.05","0.01"))
  
  ymin = 0
  
  ymax = max(c(-log10(gwas_full$p), test_lines$intercepts))
  
  gw_sig_idx = which(plt_dat[,2] > test_lines$intercepts[1])

  print(table(plt_dat$Test[gw_sig_idx]))
 
  plt_dat_weak = plt_dat[rm_gap_idx[which(rm_gap_idx %in% which(plt_dat$p > test_lines$intercepts[2]))], ]
  
  plt_dat_strong = plt_dat[which(!(seq(1,nrow(plt_dat)) %in% rm_gap_idx[which(rm_gap_idx %in% which(plt_dat$p > test_lines$intercepts[2]))])), ]
  
  saveRDS(plt_dat, file = paste("OUT/combi/", ph, "_combi_plt.rds", sep = ""))
  saveRDS(gwas_full, file = paste("OUT/combi/", ph, "_combi_gwas.rds", sep = ""))
  
  if(!noplots){
    plt = ggplot() + geom_point(data = plt_dat_strong, aes(x = pos, y = p, col = Gap_frq, shape = Test), size = 4) +
      geom_point(data = plt_dat_weak, aes(x = pos, y = p, shape = Test), size = 4, col = "black") +
      ylim(ymin, ymax) +
      geom_abline(data = test_lines, aes(intercept = intercepts, slope = slopes, linetype = FWER)) +
      scale_x_continuous(breaks = c(0, 500000, 1000000, 1500000, 2000000), name = "Basepair position", labels = scales::comma) +
      ylab(expression("-log"[10]*"(p-value)")) + theme(text = element_text(size=24)) + guides( color = guide_colorbar(order = 1), fill = guide_legend(order = 1)) 
    # plt
    ggsave(paste("PLOTS/", ph, ".png", sep = ""), plt, width = 15, height = 7)
    # ggsave(paste("PLOTS/combigwas_", ph, ".pdf", sep = ""), plt, width = 15, height = 7)
    
    qqplt = gridExtra::arrangeGrob(qq::qq_plot(gwas_full$p[kp_gap_idx]))
    ggsave(qqplt, device = "png", filename = paste("PLOTS/qq_", ph, ".png", sep = ""), width = 7, height = 7)
    # ggsave(qqplt, device = "pdf", filename = paste("PLOTS/qq_", ph, ".pdf", sep = ""), width = 7, height = 7)
    
    plt_full = gridExtra::arrangeGrob(plt, qqplt, ncol = 1)
    ggsave(plt_full, device = "png", filename = paste("PLOTS/", ph, "_gwas.png", sep = ""), width = 15, height = 10)
  }
  
  pltpath = "RUN03" # full test
  resfiles = dir(pltpath, pattern = "res_")
  resfiles = resfiles[grep(ph, resfiles)]
  lmm_res = read.csv(file.path(pltpath, resfiles[grep(paste(ph,"lmm",sep = "_"), resfiles)]), sep = "\t")
  lmm_res$variant = unname(sapply(lmm_res$variant, function(x) unlist(strsplit(x, "_"))[2] ))
  
  lmm_res_ = lmm_res[which(lmm_res$variant %in% gwas_full$pos),]
  gwas_full_ = gwas_full[which(gwas_full$pos %in% lmm_res_$variant),]
  flmm_combi_df = cbind(lmm_res_, gwas_full_)
  
  if(!noplots){
    plt_cmp = ggplot(data = flmm_combi_df) + geom_point(aes(x = -log10(lrt.pvalue), y = -log10(p), col = gf)) +
      geom_abline(slope = 1, intercept = 0) +
      facet_wrap(~Test, ncol = 2) + xlab("-log10(FaSTLMM p-vals)") + ylab("-log10(Combi p-vals)") +
      theme(text = element_text(size=24)) 
    ggsave(plt_cmp, device = "png", filename = paste("PLOTS/", ph, "_breakdown.png", sep = ""), width = 15, height = 10)
  }
  
  target = read.csv(paste("targets",unlist(strsplit(ph, "y"))[2], sep = ""), header = F)
  # target$V1 = colnames(mxCC)[as.numeric(target$V1)]
  u = c()
  uu = c(5, 10, 15, 20, 25)
  x = 1
  for(ui in 1:nrow(target)) {
    u = c(u, uu[x])
    if(x < 5){
      x = x + 1
    } else {
      x = 1
    }
  }
  
  t1 = cbind(lmm_res_[which(lmm_res_$variant %in% target$V1),])
  t2 = cbind(gwas_full_[which(gwas_full_$pos %in% target$V1),])
  tt = cbind(t1, t2, u)
  if(!noplots){
    plt_cmp2 = ggplot(tt, aes(x = -log10(lrt.pvalue), y = -log10(p))) + geom_point(aes(col = gf, shape = Test), size = 3) +
      geom_abline(slope = 1, intercept = 0) + geom_abline(slope = 0, intercept = 6.431431) + geom_vline(xintercept = 6.431431) + 
      xlab("-log10 p-vals MA") + ylab("-log10 p-vals gap/snp") +
      theme(text = element_text(size=24)) 
    ggsave(plt_cmp2, device = "png", filename = paste("PLOTS/", ph, "_cmp.png", sep = ""), width = 15, height = 10)
  }
  saveRDS(tt, file = paste("OUT/combi/", ph, "_tt.rds", sep = ""))
}
