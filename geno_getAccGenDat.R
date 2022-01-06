if(rstudioapi::isAvailable()) setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # WORKING DIRECTORY
# Get Accessory Genome Data, along with aligned phenotype files
# Checked 20210621

library(data.table)
library(stringr)

# gene_presence_absence_roary.csv is not provided, it can be generated using panaroo [https://github.com/gtonkinhill/panaroo]

acc_dat = data.frame(fread("gene_presence_absence_roary.csv", select = c(1)))

# Can we identify accessory genome data from this?
summary(acc_dat$`No. isolates`)

# # 
# # # 2677 from bash: head -1 gene_presence_absence_roary.csv | sed 's/[^,]//g' | wc -c
for(i in 2:2677) {
 temp = fread("gene_presence_absence_roary.csv", select = c(i))
 acc_dat = cbind(acc_dat, temp)
 print(i)
}
#
# saveRDS(acc_dat, "accessory_genome_data.rds")

acc_dat = readRDS("accessory_genome_data.rds")
acc_dat_reduced = data.frame(acc_dat[,15:ncol(acc_dat)])
rownames(acc_dat_reduced) = acc_dat$Gene
colnames(acc_dat_reduced) = colnames(acc_dat[15:ncol(acc_dat_reduced)])
# 

# refmt = function(vals, colnm) return(str_count(vals, colnm))
# for(i in 1:ncol(acc_dat_reduced)) acc_dat_reduced[,i] = refmt(acc_dat_reduced[,i], colnames(acc_dat_reduced)[i])
# saveRDS(acc_dat_reduced, "acc_dat_reduced.rds")
################################################################################

acc_dat_reduced = readRDS("acc_dat_reduced.rds")

############################## cd/acute ########################################
# Let's make accessory genome datasets for each analysis
# pheno = readRDS("raw_pheno_cd_core.rds") # core only
pheno = readRDS("raw_pheno_cd_coreacc.rds") # core + acc

pheno_reduced = c()
accessory_genome = c()
for(i in 1:nrow(pheno)){
  tmp = which(colnames(acc_dat_reduced) %in% pheno$sampleID[i])
  # It appears that some sequences have been dropped? Create reduced phenotype and genotype datasets
  if(length(tmp) != 1) {
    print(paste("Isolate", i, ":", pheno$sampleID[i], "not found, dropping! Length(tmp) = " , length(tmp)))
  } else {
    pheno_reduced = rbind(pheno_reduced, pheno[i,])
    accessory_genome = rbind(accessory_genome, acc_dat_reduced[,tmp])
  }
}

rownames(accessory_genome) = pheno_reduced$sampleID
colnames(accessory_genome) = rownames(acc_dat_reduced)
write.table(pheno_reduced$sampleID, "seqs_cd.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
saveRDS(accessory_genome, "cd_acc.rds")
saveRDS(pheno_reduced, "cd_pheno.rds")
############################## cd/acute ########################################

# Let's make accessory genome datasets for each analysis
########################## cef.mic/pen.mic #####################################
# pheno = readRDS("raw_pheno_mic_core.rds") # core only
pheno = readRDS("raw_pheno_mic_coreacc.rds") # core + acc

pheno_reduced = c()
accessory_genome = c()
for(i in 1:nrow(pheno)){
  tmp = which(colnames(acc_dat_reduced) %in% pheno$sampleID[i])
  # It appears that some sequences have been dropped? Create reduced phenotype and genotype datasets
  if(length(tmp) != 1) {
    print(paste("Isolate", i, ":", pheno$sampleID[i], "not found, dropping! Length(tmp) = " , length(tmp)))
  } else {
    pheno_reduced = rbind(pheno_reduced, pheno[i,])
    accessory_genome = rbind(accessory_genome, acc_dat_reduced[,tmp])
  }
}

rownames(accessory_genome) = pheno_reduced$sampleID
colnames(accessory_genome) = rownames(acc_dat_reduced)
write.table(pheno_reduced$sampleID, "seqs_mic.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
saveRDS(accessory_genome, "mic_acc.rds")
saveRDS(pheno_reduced, "pheno_mic.rds")

