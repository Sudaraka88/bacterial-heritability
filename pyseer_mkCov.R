# This file creates the covariate information required for pySEER
if(rstudioapi::isAvailable()) setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # WORKING DIRECTORY

refmt = function(nm){
  temp = unlist(strsplit(unlist(strsplit(nm, split = "_S_"))[1], split = "_"))
  return(paste(temp[1], "_", temp[2], "#", temp[3], sep = ""))
}

# pheno = readRDS("../../Clustering/Results/FBCLUSTER_GLS/mapped_pheno_Wclust_MIC.rds"); fasta_order_seqs = read.table("seqs.mic", as.is = T, comment.char = "") # for MIC
pheno = readRDS("../../Clustering/Results/FBCLUSTER_GLS/mapped_pheno_Wclust.rds"); fasta_order_seqs = read.table("seqs.cdacute", as.is = T, comment.char = "") # for cd/acute

fasta_order_fixed_seqs = apply(fasta_order_seqs, 1, refmt)
idxord = unname(sapply(fasta_order_fixed_seqs, function(x) which(pheno$sampleID %in% x)))

# ## FOR MIC UNCOMMENT BELOW ###
# cluster = paste("clust", pheno$cluster[idxord], sep ="")
# category = pheno$category[idxord]
# acute = pheno$acute[idxord]
# cef.mic_df = data.frame(samples = fasta_order_seqs[,1], cluster = cluster, category = category, acute = acute)
# # # #
# write.table(cef.mic_df, file = "cov.mic", quote = FALSE, row.names = FALSE, sep = "\t")
# # make lineage file for mic
# cluster = pheno$cluster[idxord]
# mic_lineage_df = data.frame(fasta_order_seqs[,1], cluster)
# write.table(mic_lineage_df, file = "fastbaps_mic.txt", quote = FALSE, row.names = FALSE, sep = "\t", col.names = FALSE)

################################################################################################
## FOR CD/ACUTE UNCOMMENT BELOW ###
cluster = paste("clust", pheno$cluster[idxord], sep ="")
acute = ifelse(pheno$acute[idxord], "yes", "no")
cd_df = data.frame(samples = fasta_order_seqs[,1], cluster = cluster, acute = acute)
#
write.table(cd_df, file = "cov.cdacute", quote = FALSE, row.names = FALSE, sep = "\t")
# make lineage file for cd/acute
cluster = pheno$cluster[idxord]
cd_lineage_df = data.frame(fasta_order_seqs[,1], paste("BAPS_",cluster,sep = ""))
write.table(cd_lineage_df, file = "fastbaps_cd.acute.txt", quote = FALSE, row.names = FALSE, sep = "\t", col.names = FALSE)
################################################################################################

