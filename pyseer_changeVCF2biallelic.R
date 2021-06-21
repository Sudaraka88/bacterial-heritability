# This code adds the accessory genome data in binary coding to the end of the VCF
# Change vcf file to bi-allelic type
if(rstudioapi::isAvailable()) setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # WORKING DIRECTORY

library(vcfR)
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

acc_gene_bincoder = function(acc_genes){
  acc_genes_temp = matrix(0, nrow = nrow(acc_genes), ncol = ncol(acc_genes))
  for(x in 1:ncol(acc_genes)){
    acc_genes_temp[which(acc_genes[,x] == max(acc_genes[,x])),x] = 1
  }
  acc_genes_temp = unique.matrix(acc_genes_temp)
}
######################################################################################################################################
# None of the conventional filtering techniques seem to work, best bet is to manually filter and recode the VCF file
# Let's stick with the original VCF file

# tempVCF = read.vcfR("maela3K_MIC_outlierRmd_snps.vcf") # for MIC dataset
# tempVCF = read.vcfR("Old/maela3K_MIC_outlierRmd_wDashes_snps.vcf") # this is the new MIC snps file with dashes in FASTA

# tempVCF = read.vcfR("Old/maela3K_CD_snps.vcf") # for CD/acute dataset

# Let's hack the vcf file and see what happens...

# Quickly replace * with ., this doesn't make much sense, stop using
# ALT_ = unname(sapply(tempVCF@fix[,5], function(x) gsub('[*]',".",x)))
# tempVCF@fix[,5] = ALT_
# write.vcf(tempVCF, "Old/maela3K_CD_gapsFxd.vcf.gz")

# Replace the alternative alleles with gaps 
# [This loses too much power! Fixed below]
# tempVCF@fix[,5] = rep('*',length(tempVCF@fix[,5]))
# This function recreates the ALT allele 
replace_alt = function(alt){
  # in the matrix, 0 marks major allele, other numbers are in splt order
  # {*,A,T}:= 1 - gap, 2 - A, 3 - T, 0 will be the major allele
  splt = unlist(strsplit(alt,","))
  nAll = length(splt)
  gap_idx = which(splt%in%"*") # if gap is present get idx
  r_list = NA
  # rp: indicate whether the VCF file needs modification
  if(length(gap_idx) == 0) gap_idx = NA # no gap idx, replace with NA
  if(length(splt) == 1){
    r_list = list(r_alt = splt, gap_idx = gap_idx, nAll = nAll, alt = alt, rp = 0)
  } else if(length(splt) == 2){ # There are two alt alleles
    if(!is.na(gap_idx)){ # one of the alternative alleles is gap, this is allowed
      r_list = list(r_alt = paste(splt[c(gap_idx, seq(1,2)[-gap_idx])],sep = "", collapse = ","), gap_idx = gap_idx, nAll = nAll, alt = alt, rp = 1)
    } else { # The alternative allele is also a nuc, can't have that
      r_list = list(r_alt = splt[1], gap_idx = gap_idx, nAll = nAll, alt = alt, rp = 1)
    }
  } else if(length(splt) == 3){ # There are three alt alleles
    if(!is.na(gap_idx)){ # one of the alternative alleles is gap, this is allowed
      r_list = list(r_alt = paste(splt[c(gap_idx, sample(seq(1,3)[-gap_idx], 1))],sep = "", collapse = ","), 
                    gap_idx = gap_idx, nAll = nAll, alt = alt, rp = 1)
    } else { # The alternative allele is also a nuc, can't have that
      r_list = list(r_alt = splt[1], gap_idx = gap_idx, nAll = nAll, alt = alt, rp = 1)
    }
  } else if(length(splt) == 4){ # There are three alt alleles
    if(!is.na(gap_idx)){ # one of the alternative alleles is gap, this is allowed
      r_list = list(r_alt = paste(splt[c(gap_idx, sample(seq(1,4)[-gap_idx], 1))],sep = "", collapse = ","), 
                    gap_idx = gap_idx, nAll = nAll, alt = alt, rp = 1)
    } else { # The alternative allele is also a nuc, can't have that
      r_list = list(r_alt = splt[1], gap_idx = gap_idx, nAll = nAll, alt = alt, rp = 1)
    }
  } else if(length(splt) == 5){ # There are three alt alleles
    if(!is.na(gap_idx)){ # one of the alternative alleles is gap, this is allowed
      r_list = list(r_alt = paste(splt[c(gap_idx, sample(seq(1,5)[-gap_idx], 1), gap_idx)],sep = "", collapse = ","), 
                    gap_idx = gap_idx, nAll = nAll, alt = alt, rp = 1)
    } else { # The alternative allele is also a nuc, can't have that
      r_list = list(r_alt = splt[1], gap_idx = gap_idx, nAll = nAll, alt = alt, rp = 1)
    }
  } else {
    print(paste("Unknown input", alt))
  }
  return(r_list)
  
}

######################################################################################################################################
# pos = readRDS("../../Accessory_Genome/DATASET2/Pensar2019/acc_cd.acute_core_pos_NC.rds") # The two positions after filtering
# tempVCF = read.vcfR("../DATASET2/snps_cdacute.vcf") # for CD/acute dataset
# acc_genes = readRDS("../../Accessory_Genome/DATASET2/cd_acute_accessory_genome.rds")

pos = readRDS("../../Accessory_Genome/DATASET2/Pensar2019/acc_mic_core_pos_NC.rds") # The two positions after filtering
tempVCF = read.vcfR("../DATASET2/snps_mic.vcf") # for MIC dataset
acc_genes = readRDS("../../Accessory_Genome/DATASET2/pen.mic_cef.mic_accessory_genome.rds")

rnms = unname(sapply(rownames(acc_genes), function(x) sub("#", "_", x)))
acc_genes = acc_gene_bincoder(acc_gene_cleaner(acc_genes)) # This can be added to gt now
rownames(acc_genes) = rnms
n_acc = ncol(acc_genes)
# Code to map pos in vcf to picked pos in <pos>
vcf_pos = as.numeric(unname(tempVCF@fix[,2]))
retain_idx = vcf_pos %in% pos # retained core genome SNPs

new_fix = tempVCF@fix[retain_idx,]
# Create a fix_acc matrix to append to new_fix
# "CHROM" "POS"  "ID" "REF"  "ALT" "QUAL"  "FILTER" "INFO"
offset = 2350000 # offset by the same amount

tempVCF@meta[2]
set.seed(100)
acc_fix = data.frame("CHROM" = rep(1, n_acc), 
                     "POS" = seq(offset, by = 100, length = n_acc),
                     "ID" = rep(NA, n_acc), 
                     "REF" = sample(c("A","C","G","T"), n_acc, replace = T),
                     "ALT" = rep('*', n_acc),
                     "QUAL" = rep(NA, n_acc), 
                     "FILTER" = rep(NA, n_acc), 
                     "INFO" = rep(NA, n_acc))

new_fix = as.matrix(rbind(new_fix, acc_fix)) # This seems to be OK
rm(acc_fix)

new_gt = tempVCF@gt[retain_idx,]
# Should also do a sequence mapping with acc_data
cnms = colnames(new_gt) # These are the sequences appearing in VCF, only a subset in acc_gene
match = sapply(rownames(acc_genes), function(x) which(x == cnms))
new_gt = cbind("FORMAT" = new_gt[,1], new_gt[,match])

acc_gt = data.frame("FORMAT" = rep("GT", n_acc),
                    t(acc_genes))
colnames(acc_gt) = c("FORMAT", rnms)

new_gt = as.matrix(rbind(new_gt, acc_gt))
rm(acc_gt)
# Try force back these values to tempVCF
tempVCF@fix = new_fix
gc()
tempVCF@gt = new_gt
gc()
rm(new_gt, new_fix)
gc()

# Now let's see if we can binary code this matrix
tempVCF@fix[,5] = rep('*',length(tempVCF@fix[,5]))
gc()

# nAll = rep(0, nrow(tempVCF@fix))
# gap_idx = nAll
# for(i in 1:length(tempVCF@fix[,5])){
#   temp = replace_alt(unname(tempVCF@fix[i,5]))
#   if(temp$rp) {
#     tempVCF@fix[i,5] = temp$r_alt # replacement needed
#     nAll[i] = temp$nAll # Can use this as a marker, if 0, nothing needs to been done to the specific row in gtmx
#     gap_idx[i] = temp$gap_idx
#   }
# }


ncols = ncol(tempVCF@gt)
nrows = nrow(tempVCF@gt)

head(tempVCF@gt[,2])

# This is way too slow! can we do a matrix wide operation instead?
# recodeSNP = function(snp){
#   snp[which(snp>1)] = 1
#   return(snp)
# }


# for(i in 1:nrows){
#   tempVCF@gt[i,2:ncols] = recodeSNP(tempVCF@gt[i,2:ncols])
# }


cnm = colnames(tempVCF@gt)
rnm = row.names(tempVCF@gt)

gtmx = tempVCF@gt
colnames(gtmx) = NULL
rownames(gtmx) = NULL
col1 = gtmx[,1]
gtmx = gtmx[,-1]
gtmx = as.numeric(gtmx)
gtmx[gtmx > 1] = 1
gtmx = matrix(c(col1, gtmx), ncol = length(cnm))
colnames(gtmx) = cnm
gc()
# saveRDS(gtmx, "maela3K_MIC_outlierRmd_gtmx.rds")
# saveRDS(gtmx, "maela3K_CD_gtmx.rds")
tempVCF@gt = gtmx
gc()
# write.vcf(tempVCF, "biallelicsnps_cdacute.vcf.gz"); write.table(rnms, "seqs.cdacute", quote = F, row.names = F, col.names = F)
write.vcf(tempVCF, "biallelicsnps_mic.vcf.gz"); write.table(rnms, "seqs.mic", quote = F, row.names = F, col.names = F)
