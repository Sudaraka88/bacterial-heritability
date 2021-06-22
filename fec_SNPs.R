if(rstudioapi::isAvailable()) setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # WORKING DIRECTORY
# Checked 20210621
# Build a least set of SNPs from cef.mic/pen.mic highly associated SNPs to be used as fixed effect covariates in LMM
# library(treeWAS)
library(ggplot2)
library(Rcpp)
getSNP = function(X){
  u = unique(X)
  if(length(u) == 1) return(0) else return(length(u))
  # 0 if no SNP, else # of SNPs returned
}
sourceCpp("tableC.cpp")
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
# Core Only - not provided
# mxCC = readRDS("mapped_mx_NC_MIC.rds") # load SNPs
# pheno = readRDS("pen.mic_cef.mic_pheno_reduced.rds")
# mxCC = readRDS("mapped_mx_NC_MIC.rds") # load SNPs
# pheno = readRDS("mapped_pheno_Wclust_MIC.rds")

# Core + Accessory
mxCC = readRDS("acc_mic_core_mx_NC.rds") # not provided 
pheno = readRDS("mic_pheno.rds") # not provided
mxCC = filter_alt(mxCC)

nSNPs = ncol(mxCC)
nSEQs = nrow(mxCC)
print(paste("Successful! Found", nSNPs, "SNPs and", nSEQs, "SEQs"))

bonf_lines = data.frame(intercepts = c(-log10(0.05/ncol(mxCC)), -log10(0.01/ncol(mxCC))),
                        slopes = c(0,0),
                        gw_bonf = c("0.05","0.01"))

# pen/cef decision
passoc_gls = readRDS("OUT/cef.mic_passoc_gls_nocov.rds") # cef.mic
# passoc_gls = readRDS("OUT/passoc_gls_pen.mic.rds") # pen.mic
idxes = which(passoc_gls$tab$predictor %in% colnames(mxCC))
print(paste("Alignment check?", all(passoc_gls$tab$predictor[idxes] == colnames(mxCC))))
lgpvs = -log10(passoc_gls$tab$pval)[idxes]


# Let's start with the gene summary to find the starts/ends
# gene_summary = data.frame(snp = snp_list,
#                           gene = gene_list,
#                           starts = starts,
#                           ends = ends)
# # Manually, start = 284650, 340532
gw_sig_idx = which(lgpvs > bonf_lines$intercepts[2])
x = as.numeric(colnames(mxCC))[gw_sig_idx]
candidate_x = x[which(x >= 284650 & x <= 340532)]
mxCC_clump = mxCC[,which(as.numeric(colnames(mxCC)) %in% candidate_x)]
candidate_pvs = lgpvs[which(as.numeric(colnames(mxCC)) %in% candidate_x)]
names(candidate_pvs) = colnames(mxCC_clump)

snporder = order(candidate_pvs, decreasing = T) # snporder[1] is the most important SNP
keep = snporder[1]
thresh = 0.25
NN = length(snporder)
for(i in 2:NN){
  # i = i + 1
  idx = snporder[i]
  # cor(mxCC_clump[,idx], mxCC_clump[,rem])^2
  # any(cor(mxCC_clump[,idx], mxCC_clump[,rem])^2>thresh)
  if( all(abs(cor(mxCC_clump[,idx], mxCC_clump[,keep])) < thresh) ) {
    keep = c(keep, idx)
  }
  # if(i >= length(snporder)) {
  #   print(paste("Retained", length(which(markers==TRUE)), "SNPs"))
  #   break
  # }
}
print(length(keep))
mxCC_ret = mxCC_clump[,keep]
cor_ret = cor(mxCC_ret)
diag(cor_ret) = 0
print(max(cor_ret))

# The plotting can to be done now if needed

# covs = rep(0, length(lgpvs))
# covs[which(colnames(mxCC) %in% colnames(mxCC_ret))] = 1
# 
# df.plot_lmm = data.frame(pos = as.numeric(colnames(mxCC)), pval = lgpvs, ret = as.factor(covs))
# gw_sig_idx = which(df.plot_lmm$pval > bonf_lines$intercepts[2])
# x = df.plot_lmm$pos[gw_sig_idx]
# 

# Gene mapping
# source("map2gene.R")
# 
# gene_list = c()
# snp_list = c()
# starts = c()
# ends = c()
# for(i in x){
#   tmp = map2gene(i, gene_data)
#   if(tmp$dist_from_gene < 1000) {
#     gene_list = c(gene_list, tmp$gd$gene[1])
#     starts = c(starts, tmp$gd$start[1])
#     ends = c(ends, tmp$gd$end[1])
#     snp_list = c(snp_list, i)
#   }
# }
# gene_summary = data.frame(snp = snp_list,
#                           gene = gene_list,
#                           starts = starts,
#                           ends = ends)
# 
# # remove duplicated entries
# gene_summary = gene_summary[!duplicated(gene_summary$gene),]
# noquote(paste0(gene_summary$gene, ",", sep = "", collapse = NULL))
# # lab_pos = c(-10,-20,-30,-40,-50,-60,-50,-40,-30,-20,-10,-10, -20,-10) # pen.mic
# 
# # lab_pos = c(-10,-10,-20,-30,-40,-50,-60,-70,-80,-10,-20,-30, -40,-50,-10,-10,-10,-10,-20,-10,-20,-10,-10,-20, -10,-10,-10) # cef.mic
# 
# # lab_pos = c(-10,-10,-20,-30,-40,-50,-60,-70,-80,-10,-20,-30, -40,-50,-10,-10,-10,-10,-20,-10,-20,-10,-10,-10, -20,-10,-10)
# 
# lab_pos = c(-10,-20,-30,-40,-50,-60,-10,-20,-30,-10,-10,-10)/2
# # lab_pos = c(-2,-6,-10,-14,-2,-6, -10, -14)
# ymin = min(lab_pos) - 2
# 
# ymax = max(c(lgpvs, bonf_lines$intercepts))
# plt.mic = ggplot(df.plot_lmm) + geom_point(aes(x = pos, y = pval, col = ret), alpha = 0.5) + 
#   geom_abline(data = bonf_lines, aes(intercept = intercepts, slope = slopes, linetype = gw_bonf))  + 
#   annotate("text", x =  gene_summary$starts, y = lab_pos, label = gene_summary$gene, angle = 90) +
#   xlab("SNP position") + ylab("-log10(pval): LMM") + ylim(ymin, ymax) + 
#   ggtitle(paste("Pheno: cef.mic, lme4qtl, MI", sep = ""))
# 
# plt.mic
saveRDS(mxCC_ret, "FECs_cef.mic.rds")
# saveRDS(mxCC_ret, "FECs_pen.mic.rds")


