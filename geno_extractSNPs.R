if(rstudioapi::isAvailable()) setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # WORKING DIRECTORY
# Checked 20220105
library(Rcpp) # required for C code interfacing
# This code is used to extract SNPs from fasta files
fasta_ = "cd_isolates"; pheno = readRDS("pheno_cd.rds") # fasta alignment and phenotype not provided
# fasta_ = "mic_isolates"; pheno = readRDS("pheno_mic.rds") # fasta alignment and phenotype not provided
fasta = paste(fasta_,".fasta",sep = "") 

folder = "Method2" # path to output folder

view_filt_hist = function(freq_table_filt){ # quickly view the ACGTN distribution by parsing the filtered freq_table
  par(mfrow = c(2,3))
  hist(freq_table_filt$A)
  hist(freq_table_filt$C)
  hist(freq_table_filt$G)
  hist(freq_table_filt$T)
  hist(freq_table_filt$N)
  plot(freq_table_filt$POS, pch = '.')
}

# generate the frq file if not there already
if(!file.exists(paste(fasta,".frq",sep=""))) system(paste("sh createFreqFile.sh", fasta)) # do not create if exists 
freq_dat = read.table(paste(fasta,".frq",sep=""), header = FALSE, fill = TRUE, col.names = c("CHR", "POS", "N_ALL", "N_CHR", "A1", "A2", "A3", "A4", "A5"), stringsAsFactors = FALSE) # header format
# strip unnecessary data and reformat
freq_dat = freq_dat[-1,]
freq_dat = freq_dat[,c(-1,-4)]
freq_dat[,1] = as.numeric(freq_dat[,1])
freq_dat[,2] = as.numeric(freq_dat[,2])
rownames(freq_dat) = NULL 

t_snps = dim(freq_dat)[1] # total number of SNPs

Rcpp::sourceCpp("Rcppf.cpp")
freq_out = matrix(rep(0, 5*t_snps), nrow = t_snps) # update with freq data
makeFreqTable2(freq_dat$N_ALL, unname(as.matrix(freq_dat[3:7])), t_snps, freq_out) # Order A, C, G, T, N 
# Save to nice data.frame
freq_table = data.frame(freq_dat$POS, freq_dat$N_ALL, freq_out)
names(freq_table) = c("POS", "N_ALL", "A", "C", "G", "T", "N")
saveRDS(freq_table, paste(fasta_, "_freq_table.rds", sep = "")) # save freq_table


############################################ SNP FILTERING ################################################################
MAF12out_ng = apply(freq_out[,1:4], 1, function(x) sort(x, decreasing = T)[1:2]) # Get major allele freq (non-gap)
MAF12out_wg = apply(freq_out[,1:5], 1, function(x) sort(x, decreasing = T)[1:2]) # Get major allele freq (with-gap)
if(folder == "Method1"){
  # Filtering -
  rm_list1 = which(freq_table$N > freq_table$A & freq_table$N > freq_table$C &
                     freq_table$N > freq_table$G & freq_table$N > freq_table$T) # Remove all where N is the major allele
  rm_list2 = which(MAF12out_ng[1,] > 0.99) # Remove major allele > 0.99
  rm_list3 = which(MAF12out_ng[2,] < 0.01) # Remove minor allele < 0.01
  rm_list4 = which(freq_out[,5]  > 0.05) # Remove gap_freq > 0.05
  rm_list = union(rm_list1, union(rm_list2, union(rm_list3, rm_list4)))
  freq_table_filt = freq_table[-rm_list,]
  view_filt_hist(freq_table_filt)
} else if(folder == "Method2"){
  # Filtering - Pensar 2019
  dir.create(folder)
  rm_list1 = which(freq_out[,5] > 0.15)  # gap frequency filter
  rm_list2 = which(MAF12out_ng[2,]<0.01)
  rm_list = union(rm_list1, rm_list2)
  freq_table_filt = freq_table[-rm_list,]
  view_filt_hist(freq_table_filt)
} else if(folder == "Method3") {
  # Filtering - Include more SNPs
  rm_list1 = which(freq_table$N > freq_table$A & freq_table$N > freq_table$C &
                     freq_table$N > freq_table$G & freq_table$N > freq_table$T) # Remove all where N is the major allele
  rm_list2 = which(MAF12out_wg[1,] >= 0.99) # Remove major allele freq > 0.99
  # rm_list3 = which(MAF12out_ng[2,] < 0.01) # Remove minor allele freq < 0.01
  rm_list3 = which(freq_out[,5]  > 0.01) # Remove gap_freq > 0.2
  rm_list = union(rm_list1, union(rm_list2, rm_list3))
  freq_table_filt = freq_table[-rm_list,]
  view_filt_hist(freq_table_filt)
}
###########################################################################################################################
saveRDS(freq_table_filt, file.path(folder,paste(fasta_, "_freq_table.rds", sep = "")))
pos = freq_table_filt$POS
saveRDS(pos, file.path(folder, paste(fasta_, "_pos.rds", sep = "")))
names = pheno$sampleID # This is the correct order

library(foreach)
library(doParallel)
cores = parallel::detectCores()
registerDoParallel(cores-2)
# opts = list(progress = function(n) setTxtProgressBar(txtProgressBar(max = length(names), style = 3), n)) 
t = Sys.time()
dna_seq = foreach(i = 1:length(names), .combine = "rbind", .verbose = TRUE) %dopar% { # t < 7m
  dna = unlist(strsplit(system(paste("samtools faidx " ,fasta , names[i], " | perl -pe 'chomp unless /^>/' | tail -1", collapse = ""), intern = TRUE), split = ""))
  stringr::str_to_upper(dna[pos])
}
Sys.time() - t
rownames(dna_seq) = names
colnames(dna_seq) = as.character(pos)


# order needs to be fixed, quick hack

saveRDS(dna_seq,file.path(folder, paste(fasta_, "_mx.rds", sep = "")))


