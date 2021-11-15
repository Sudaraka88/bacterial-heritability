if(rstudioapi::isAvailable()) setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # WORKING DIRECTORY

# we need to code the simulated FASTA files into Numerically/Binary coded SNPs
# Also, are there gaps?

library(foreach)
library(doParallel)
library(Rcpp)
view_filt_hist = function(freq_table_filt){
  par(mfrow = c(2,3))
  hist(freq_table_filt$A)
  hist(freq_table_filt$C)
  hist(freq_table_filt$G)
  hist(freq_table_filt$T)
  hist(freq_table_filt$N)
  plot(freq_table_filt$POS, pch = '.')
}
Rcpp::sourceCpp("../Genotype2/Rcppf.cpp")
sourceCpp("../Genotype2/tableC.cpp")
num_code_SNP = function(x,n){
  num_SNP = rep(0,n)
  tbl = sort(tableC(x)/n, decreasing = T)
  for(i in 1:length(tbl)){
    num_SNP[which(x%in%names(tbl[i]))] = unname(tbl[i])
  }
  num_SNP = as.numeric(scale(num_SNP))
  # In rare cases, allele frequencies can be even across sequences, resulting in 1/af distributions to all alleles
  if(any(is.nan(num_SNP))) {
    k = 0
    for(i in 1:length(tbl)){
      num_SNP[which(x%in%names(tbl[i]))] = k
      k = k+1
    }
    num_SNP = scale(num_SNP)
    if(!any(is.nan(num_SNP))){
      print(paste("Caught equally distributed variant with", k, "levels"))
    } else {
      print(paste("Caught equally distributed variant with", k, "levels, failed to correct"))
    }
  }
  return(num_SNP) 
}

folder_ = dir(pattern = "^[:h:][:0:]")
# These are identical, no point in repeating!
# scanning step
for(folder in folder_){
  fasta_ = paste(folder, "/simulations/genSim/genSim", sep = "")
  if(file.exists(paste(fasta_, "_posAFC.rds", sep = ""))){
    freq_table = readRDS(paste(fasta_, "_freq_table.rds", sep = ""))
    pos = readRDS(paste(fasta_, "_pos.rds", sep = ""))
    freq_table_filt = readRDS(paste(fasta_, "_freq_table_filt.rds", sep = ""))
    gap_freq = readRDS(paste(fasta_, "_gap_freq.rds", sep = ""))
    mxCC = readRDS(paste(fasta_, "_mx.rds", sep = ""))
    mxBC = readRDS(paste(fasta_, "_MAC.rds", sep = ""))
    mxNC = readRDS(paste(fasta_, "_AFC.rds", sep = ""))
    pos_afc = readRDS(paste(fasta_, "_posAFC.rds", sep = ""))
    print("Pre-computed instance found")
    break
  } else {
    print("No pre-computed instance!")
  }
}




for(folder in folder_){
  print(paste("Analysing folder", folder))
  fasta_ = paste(folder, "/simulations/genSim/genSim", sep = "")
  if(!file.exists(paste(fasta_, "_posAFC.rds", sep = ""))){ # Check for file in path
    if(exists("mxNC")){
      print(paste("Copying to", folder))
      saveRDS(freq_table, paste(fasta_, "_freq_table.rds", sep = ""))
      saveRDS(pos, paste(fasta_, "_pos.rds", sep = ""))
      saveRDS(freq_table_filt, paste(fasta_, "_freq_table_filt.rds", sep = ""))
      saveRDS(gap_freq, paste(fasta_, "_gap_freq.rds", sep = ""))
      saveRDS(mxCC, paste(fasta_, "_mx.rds", sep = ""))
      saveRDS(mxBC, paste(fasta_, "_MAC.rds", sep = ""))
      saveRDS(mxNC, paste(fasta_, "_AFC.rds", sep = ""))
      saveRDS(pos_afc, paste(fasta_, "_posAFC.rds", sep = ""))
    } else {
      y = read.table(paste(folder, "/simulations/phenSim/0/phenSim.phen", sep = ""))$V3
      
      fasta = paste(fasta_,".fasta",sep = "")
      
      if(file.exists(paste(fasta_, "_freq_table.rds", sep = ""))){
        freq_table = readRDS(paste(fasta_, "_freq_table.rds", sep = ""))
      } else{
        if(!file.exists(paste(fasta,".frq",sep=""))) system(paste("sh createFreqFile.sh", fasta))
        Sys.sleep(5)
        freq_dat = read.table(paste(fasta,".frq",sep=""), header = FALSE, fill = TRUE, col.names = c("CHR", "POS", "N_ALL", "N_CHR", "A1", "A2", "A3", "A4", "A5"), stringsAsFactors = FALSE)
        # strip unnecessary data and reformat
        freq_dat = freq_dat[-1,]
        freq_dat = freq_dat[,c(-1,-4)]
        freq_dat[,1] = as.numeric(freq_dat[,1])
        freq_dat[,2] = as.numeric(freq_dat[,2])
        rownames(freq_dat) = NULL
        
        t_snps = dim(freq_dat)[1] # total number of SNPs
        
        freq_out = matrix(rep(0, 5*t_snps), nrow = t_snps)
        makeFreqTable2(freq_dat$N_ALL, unname(as.matrix(freq_dat[3:7])), t_snps, freq_out) # Order A, C, G, T, N
        # Save to nice data.frame
        freq_table = data.frame(freq_dat$POS, freq_dat$N_ALL, freq_out)
        names(freq_table) = c("POS", "N_ALL", "A", "C", "G", "T", "N")
        saveRDS(freq_table, paste(fasta_, "_freq_table.rds", sep = ""))
      }
      
      MAF12out_ng = apply(freq_table[,3:6], 1, function(x) sort(x, decreasing = T)[1:2]) # Get major allele freq (non-gap)
      pheno = data.frame(sampleID = 0:(length(y)-1), y = y)
      
      snp_Seq = seq(1,nrow(freq_table))
      # There are no gaps, this filter does not work yet
      # rm_list1 = which(freq_table$N < 10/nrow(pheno))  # gap frequency filter
      # rm_list2 = which(freq_table$N > (nrow(pheno)-10)/nrow(pheno))  # gap frequency filter
      # snpseql1 = snp_Seq[-which(snp_Seq %in% sort(unique(union(rm_list1, rm_list2))))]
      rm_list3 = which(MAF12out_ng[2,] < 10/nrow(pheno))
      snpseq2 = snp_Seq[-which(snp_Seq %in% rm_list3)]
      # snpseq = sort(unique(union(snpseql1, snpseql2)))
      # freq_table_filt = freq_table[snpseq,]
      freq_table_filt = freq_table[snpseq2,]
      view_filt_hist(freq_table_filt)
      
      pos = freq_table_filt$POS
      saveRDS(pos, paste(fasta_, "_pos.rds", sep = ""))
      saveRDS(freq_table_filt, paste(fasta_, "_freq_table_filt.rds", sep = ""))
      
      gap_freq = freq_table_filt$N
      names(gap_freq) = as.character(freq_table_filt$POS)
      saveRDS(gap_freq, paste(fasta_, "_gap_freq.rds", sep = ""))
      
      names = pheno$sampleID # This is the correct order
      registerDoParallel(cores = 6L)
      # opts = list(progress = function(n) setTxtProgressBar(txtProgressBar(max = length(names), style = 3), n)) 
      t = Sys.time()
      mxCC = foreach(i = 1:length(names), .combine = "rbind", .verbose = FALSE) %dopar% { # t < 7m
        dna = unlist(strsplit(system(paste("samtools faidx " ,fasta , names[i], " | perl -pe 'chomp unless /^>/' | tail -1", collapse = ""), intern = TRUE), split = ""))
        stringr::str_to_upper(dna[pos])
      }
      Sys.time() - t
      rownames(mxCC) = names
      colnames(mxCC) = as.character(pos)
      saveRDS(mxCC, paste(fasta_, "_mx.rds", sep = ""))
      
      
      # Major allele coding
      print("Major Allele Coding")
      mxCC = readRDS(paste(fasta_, "_mx.rds", sep = ""))
      mxBC = matrix(nrow = nrow(mxCC), ncol = ncol(mxCC))
      Rcpp::sourceCpp("../Clustering/tableC.cpp")
      for(i in 1:ncol(mxCC)){
        MA = names(sort(tableC(mxCC[,i]), decreasing = T)[1])
        mxBC[,i] = as.numeric(mxCC[,i] == MA)
      }
      colnames(mxBC) = colnames(mxCC)
      rownames(mxBC) = rownames(mxCC)
      saveRDS(mxBC, paste(fasta_, "_MAC.rds", sep = ""))
      
      # Allele frequency coding
      print("Allele Frequency Coding")
      mxCC = readRDS(paste(fasta_, "_mx.rds", sep = ""))
      n = length(mxCC[,2])
      mxNC = matrix(nrow = nrow(mxCC), ncol = ncol(mxCC))
      for(i in 1:ncol(mxCC)){
        mxNC[,i] = num_code_SNP(mxCC[,i], n)
      }
      ind = which(is.na(mxNC))
      if(length(ind) > 0) {
        c = floor((ind-1) / nrow(mxNC)) + 1
        rm_idx = unique(c) # These SNPs contain no variation or have poor quality
        # Edit pos
        pos = readRDS(paste(fasta_, "_pos.rds", sep = ""))
        pos_afc = pos[-rm_idx]
        mxNC = mxNC[,-rm_idx]
        rownames(mxNC) = rownames(mxCC)
        colnames(mxNC) = as.character(pos_afc)
        paste("Any more NAs?", any(is.na(mxNC)))
      } else {
        pos_afc = readRDS(paste(fasta_, "_pos.rds", sep = ""))
      }
      colnames(mxNC) = as.character(pos_afc)
      rownames(mxNC) = rownames(mxCC)
      saveRDS(mxNC, paste(fasta_, "_AFC.rds", sep = ""))
      saveRDS(pos_afc, paste(fasta_, "_posAFC.rds", sep = ""))
    }
  } else {
    print(paste("Previously analysed:", folder))
  }
  
}

