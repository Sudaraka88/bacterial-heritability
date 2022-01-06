# This file creates the phenotype information required for pySEER
if(rstudioapi::isAvailable()) setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # WORKING DIRECTORY
# Checked 20220105

refmt = function(nm){
  temp = unlist(strsplit(unlist(strsplit(nm, split = "_S_"))[1], split = "_"))
  return(paste(temp[1], "_", temp[2], "#", temp[3], sep = ""))
}

pheno = readRDS("mic_pheno.rds"); fasta_order_seqs = read.table("seqs.mic", as.is = T, comment.char = "") # for MIC
# pheno = readRDS("cd_pheno.rds"); fasta_order_seqs = read.table("seqs.cdacute", as.is = T, comment.char = "") # for cd/acute

fasta_order_fixed_seqs = apply(fasta_order_seqs, 1, refmt)
idxord = unname(sapply(fasta_order_fixed_seqs, function(x) which(pheno$sampleID %in% x)))


### FOR MIC UNCOMMENT BELOW ###
# # cef.mic
cef.mic = log10(pheno$Ceftriaxone.MIC[idxord])
cef.mic_df = data.frame(samples = fasta_order_seqs[,1], cef.mic = cef.mic)

write.table(cef.mic_df, file = "pheno_cef.mic", quote = FALSE, row.names = FALSE, sep = "\t")

# pen.mic
pen.mic = log10(pheno$Penicillin.MIC[idxord])
pen.mic_df = data.frame(samples = fasta_order_seqs[,1], pen.mic = pen.mic)

write.table(pen.mic_df, file = "pheno_pen.mic", quote = FALSE, row.names = FALSE, sep = "\t")
# #################################################################################################


### FOR CD UNCOMMENT BELOW ###
# # cd
# cd = log10(pheno$carriage_duration[idxord])
# cd_df = data.frame(samples = fasta_order_seqs[,1], cd = cd)
# write.table(cd_df, col.names = T, file = "pheno_cd", quote = FALSE, row.names = FALSE, sep = "\t")
#  
# # # acute
# acute = pheno$acute[idxord]
# acute_df = data.frame(samples = fasta_order_seqs[,1], acute = acute)
# # 
# write.table(acute_df, file = "pheno_acute", quote = FALSE, row.names = FALSE, sep = "\t")
#################################################################################################