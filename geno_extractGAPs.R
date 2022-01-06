# This code is for the new AccessoryGenome
library(foreach)
library(doParallel)
# make a binary gap allele at each loci 20210715
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # WORKING DIRECTORY
fasta_ = "cd_isolates"; pheno = readRDS("pheno_cd.rds") # fasta alignment and phenotype not provided
# fasta_ = "mic_isolates"; pheno = readRDS("pheno_mic.rds") # fasta alignment and phenotype not provided
folder = "gapDat"

view_filt_hist = function(freq_table_filt){
  par(mfrow = c(2,3))
  hist(freq_table_filt$A)
  hist(freq_table_filt$C)
  hist(freq_table_filt$G)
  hist(freq_table_filt$T)
  hist(freq_table_filt$N)
  plot(freq_table_filt$POS, pch = '.')
}

# generate the frq file if not there already
if(!file.exists(paste(fasta,".frq",sep=""))) system(paste("sh createFreqFile.sh", fasta))
freq_dat = read.table(paste(fasta,".frq",sep=""), header = FALSE, fill = TRUE, col.names = c("CHR", "POS", "N_ALL", "N_CHR", "A1", "A2", "A3", "A4", "A5"), stringsAsFactors = FALSE)
# strip unnecessary data and reformat
freq_dat = freq_dat[-1,]
freq_dat = freq_dat[,c(-1,-4)]
freq_dat[,1] = as.numeric(freq_dat[,1])
freq_dat[,2] = as.numeric(freq_dat[,2])
rownames(freq_dat) = NULL

t_snps = dim(freq_dat)[1] # total number of SNPs
Rcpp::sourceCpp("../../Genotype2/Rcppf.cpp")
freq_out = matrix(rep(0, 5*t_snps), nrow = t_snps)
makeFreqTable2(freq_dat$N_ALL, unname(as.matrix(freq_dat[3:7])), t_snps, freq_out) # Order A, C, G, T, N
# Save to nice data.frame
freq_table = data.frame(freq_dat$POS, freq_dat$N_ALL, freq_out)
names(freq_table) = c("POS", "N_ALL", "A", "C", "G", "T", "N")
saveRDS(freq_table, paste(fasta_, "_freq_table.rds", sep = ""))


############################################ SNP FILTERING ################################################################
MAF12out_ng = apply(freq_out[,1:4], 1, function(x) sort(x, decreasing = T)[1:2]) # Get major allele freq (non-gap)
# MAF12out_wg = apply(freq_out[,1:5], 1, function(x) sort(x, decreasing = T)[1:2]) # Get major allele freq (with-gap)
# symmetric filtering
rm_list1 = which(freq_table$N > (nrow(pheno)-10)/nrow(pheno)) # Remove all loci less than 10 allele count
rm_list2 = which(freq_table$N < 10/nrow(pheno)) # Remove all with 10 gap count
rm_list = unique(sort(union(rm_list1, rm_list2)))
# all gaps
# rm_list = which(freq_table$N == 0)
freq_table_filt = freq_table[-rm_list,]
view_filt_hist(freq_table_filt)
gap_freq = freq_table_filt$N
names(gap_freq) = as.character(freq_table_filt$POS)
saveRDS(gap_freq, file.path(folder, paste(fasta_, "_gap_freq.rds", sep = "")))
###########################################################################################################################
saveRDS(freq_table_filt, file.path(folder,paste(fasta_, "_freq_table_gapDat.rds", sep = "")))
pos = freq_table_filt$POS
saveRDS(pos, file.path(folder, paste(fasta_, "_pos.rds", sep = "")))
# write.table(pos, file.path(folder,"acc_cd.acute_pos.txt"), row.names = FALSE, col.names = FALSE) # This is a temporary writing place for positions
#################################################################################
# Extract SNPs - This requires a lot of memory, better to load pos and run fresh
#################################################################################
# R Attempt
# pos = readRDS(file.path(folder, "acc_cd.acute_pos.rds")) # Extract positions from the temporary writing place
# if(!file.exists(paste(fasta,".fai",sep=""))) system(paste("samtools faidx", fasta)) # Index the fasta file
# fastainfo = read.delim(paste(fasta,".fai",sep = ""), header = FALSE, fill = TRUE, stringsAsFactors = FALSE)
# names = fastainfo[,1] # These shouldn't be used, order will be different

names = pheno$sampleID # This is the correct order


# cores = parallel::detectCores()
registerDoParallel(cores = 6L)
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

gap_seq = matrix(nrow = nrow(dna_seq), ncol = ncol(dna_seq))

for(i in 1:ncol(dna_seq)) gap_seq[,i] = as.numeric(dna_seq[,i] == "N" | dna_seq[,i] == "-")
rownames(gap_seq) = names
colnames(gap_seq) = as.character(pos)

saveRDS(gap_seq, file.path(folder, paste(fasta_, "_gap.rds", sep = "")))