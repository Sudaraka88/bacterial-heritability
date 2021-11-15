if(rstudioapi::isAvailable()) setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # WORKING DIRECTORY

# out = data.frame(X = "chr - pn1 1 0 2221315 chr1")
gene_data_red = gene_data[!is.na(gene_data$protein_id), ]

e1 = gene_data_red$end[1]
for(i in 2:nrow(gene_data_red)) if(gene_data_red$start[i] < e1) {print(i); gene_data_red$start[i] = e1 + 10 ;e1 = gene_data_red$end[i]} else e1 = gene_data_red$end[i]

e1 = gene_data_red$end[1]
for(i in 2:nrow(gene_data_red)) if(gene_data_red$start[i] < e1) {print(i); e1 = gene_data_red$end[i]} else e1 = gene_data_red$end[i]

out = as.data.frame(cbind("pn1", gene_data_red$start, gene_data_red$end))
out2 = data.frame("pn1", "1", out$V2[1], "0")
colnames(out2) = c("V1", "V2", "V3", "V4")
for(i in 2:nrow(out)){
  temp = rbind(c("pn1", out$V3[i-1], out$V2[i], 0), cbind(out[i,], 1))
  colnames(temp) = c("V1", "V2", "V3", "V4")
  out2  = rbind(out2, temp)
}
out2 = as.data.frame(out2)
# colnames(bla) = "X"
# 
# out = rbind(out, bla)

# write.table(out, file="circos-0.69-9/data/karyotype/genes.pneumo.txt", quote = F, row.names = F, col.names = F)
write.table(out2, file="data/karyotype/genes.pneumo.txt", quote = F, row.names = F, col.names = F)

cefmic = readRDS("../GWAS+PRED+LDSC/OUT/cef.mic_combi_gwas.rds")
cefmic$p = -log10(cefmic$p)
cefmic = cefmic[which(cefmic$p > 6.194522), ]

cefmic2 = as.data.frame(cbind("pn1", cefmic$pos-2, cefmic$pos+2, sapply(cefmic$Test, function(x) ifelse(x=="gap",1,2) )))
cefmic3 = data.frame("pn1", "1", cefmic2$V2[1], "0")
colnames(cefmic3) = c("V1", "V2", "V3", "V4")
for(i in 2:nrow(cefmic2)){
  temp = rbind(c("pn1", cefmic2$V3[i-1], cefmic2$V2[i], 0), cbind(cefmic2[i,]))
  colnames(temp) = c("V1", "V2", "V3", "V4")
  cefmic3  = rbind(cefmic3, temp)
}
cefmic3 = rbind(cefmic3, c("pn1", cefmic3$V3[nrow(cefmic3)], 2221315, 0))
# write.table(cbind(rep("pn1", nrow(cefmic)), cefmic$pos, cefmic$pos), file="data/karyotype/snps.pneumo_cefmic.txt", quote = F, row.names = F, col.names = F)
write.table(cefmic3, file="data/karyotype/snps.pneumo_cefmic.txt", quote = F, row.names = F, col.names = F)

penmic = readRDS("../GWAS+PRED+LDSC/OUT/pen.mic_combi_gwas.rds")
penmic$p = -log10(penmic$p)
penmic = penmic[which(penmic$p > 6.194522), ]

penmic2 = as.data.frame(cbind("pn1", penmic$pos-2, penmic$pos+2, sapply(penmic$Test, function(x) ifelse(x=="gap",1,2) )))
penmic3 = data.frame("pn1", "1", penmic2$V2[1], "0")
colnames(penmic3) = c("V1", "V2", "V3", "V4")
for(i in 2:nrow(penmic2)){
  temp = rbind(c("pn1", penmic2$V3[i-1], penmic2$V2[i], 0), cbind(penmic2[i,]))
  colnames(temp) = c("V1", "V2", "V3", "V4")
  penmic3  = rbind(penmic3, temp)
}
penmic3 = rbind(penmic3, c("pn1", penmic3$V3[nrow(penmic3)], 2221315, 0))
# write.table(cbind(rep("pn1", nrow(cefmic)), penmic$pos, penmic$pos), file="data/karyotype/snps.pneumo_penmic.txt", quote = F, row.names = F, col.names = F)
write.table(penmic3, file="data/karyotype/snps.pneumo_penmic.txt", quote = F, row.names = F, col.names = F)

gaps = readRDS("../Accessory_Genome/DATASET2/pyseer/acc_mic_core_gap_freq.rds")
gap_2 = data.frame()
x = 1
for(i in 1:13462){
  gap_2 = rbind(gap_2, as.numeric(c(names(gaps[x]), names(gaps[x+10]), mean(gaps[x:(x+10)]))))
  x = x + 10
}

write.table(cbind(rep("pn1", nrow(gap_2)), gap_2[,1], gap_2[,2], gap_2[,3]), file="data/karyotype/gaps.pneumo.txt", quote = F, row.names = F, col.names = F)

snps = readRDS("../Accessory_Genome/DATASET2/pyseer/acc_mic_core_ma_freq.rds")
maf_2 = data.frame()
x = 1
for(i in 1:13462){
  maf_2 = rbind(maf_2, as.numeric(c(names(snps[x]), names(snps[x+10]), mean(snps[x:(x+10)]))))
  x = x + 10
}

write.table(cbind(rep("pn1", nrow(maf_2)), maf_2[,1], maf_2[,2], maf_2[,3]), file="data/karyotype/maf.pneumo.txt", quote = F, row.names = F, col.names = F)
