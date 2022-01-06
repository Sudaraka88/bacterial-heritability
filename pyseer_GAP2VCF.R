# This code converts the gap data to VCF format
if(rstudioapi::isAvailable()) setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # WORKING DIRECTORY
# Checked 20220105
library(vcfR)
refmt = function(nm){
  if(length(grep("#", nm))){
    temp = unlist(strsplit(nm, split = "#"))
    return(paste(temp[1], "_", temp[2], sep = ""))
  } else {
    return(nm)
  }
}

######################################################################################################################################
pos = readRDS("../../Accessory_Genome/DATASET2/gapDat/acc_mic_core_pos.rds") # New filtering 20210802
tempVCF = read.vcfR("../DATASET2/snps_mic.vcf") # for MIC dataset
gap_dat = readRDS("../../Accessory_Genome/DATASET2/gapDat/acc_mic_core_gap.rds")

# pos = readRDS("../../Accessory_Genome/DATASET2/gapDat/cd_core_uqepi_pos.rds") # New filtering 20210802
# tempVCF = read.vcfR("../../Accessory_Genome/DATASET2/cd_core_uqepi.fasta.vcf") # for CD/acute dataset
# gap_dat = readRDS("../../Accessory_Genome/DATASET2/gapDat/cd_core_uqepi_gap.rds")

rownames(gap_dat) = unname(sapply(rownames(gap_dat), refmt))
gap_dat = t(gap_dat)
rownames(gap_dat) = NULL
gap_dat = -1*(gap_dat - 1)
# gap_dat = apply(gap_dat, 2, as.character)
gap_dat = cbind(FORMAT = rep("GT", nrow(gap_dat)), gap_dat)

vcf_pos = as.numeric(unname(tempVCF@fix[,2]))
retain_idx = vcf_pos %in% pos # retained core genome SNPs

new_fix = tempVCF@fix[retain_idx,]
tempVCF@meta[2]
new_fix[,5] = rep('*',nrow(new_fix))
tempVCF@fix = new_fix
gc()

tempVCF@gt = gap_dat
gc()

write.vcf(tempVCF, "gap_mic.vcf.gz"); write.table(colnames(tempVCF@gt), "seqs.mic", quote = F, row.names = F, col.names = F)
# write.vcf(tempVCF, "gap_cd.vcf.gz"); write.table(colnames(tempVCF@gt), "seqs.cd", quote = F, row.names = F, col.names = F)
