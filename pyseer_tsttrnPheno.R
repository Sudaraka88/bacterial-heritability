# Break pheno in to train/test sets for xfcv and loso
if(rstudioapi::isAvailable()) setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # WORKING DIRECTORY
# Checked 20220105

write_pheno = function(samples, y, testing_idx, ph_nm, fldr){
  n = length(testing_idx_)
  for(i in 1:n){
    train_idx = rep(TRUE, nrow(pheno))
    train_idx[testing_idx[[i]]] = FALSE
    train = data.frame(samples = samples[train_idx], ph = y[train_idx])
    test = data.frame(samples = samples[!train_idx], ph = y[!train_idx])
    
    write.table(x = train, file.path(fldr, paste(ph_nm, i, "trn", sep = "_")), quote = F, row.names = F, sep = "\t")
    write.table(x = test, file.path(fldr, paste(ph_nm, i, "tst", sep = "_")), quote = F, row.names = F, sep = "\t")
  }
}

ph_ = c("cd", "cef.mic", "pen.mic")
for(ph in ph_){
  # get phenotype
  if(ph == "cd"){
    pheno = readRDS("cd_pheno.rds")
    pheno_vals = log10(pheno$carriage_duration)
  } else {
    pheno = readRDS("mic_pheno.rds")
    if(ph == "cef.mic") pheno_vals = log10(pheno$Ceftriaxone.MIC)
    if(ph == "pen.mic") pheno_vals = log10(pheno$Penicillin.MIC)
  }
  pheno$sampleID = gsub("#", "_", pheno$sampleID) # pyseer requirement
  
  print("Commencing xfcv")
  set.seed(23) # hard code seed for reproducibility
  tp = 0.1
  x = 1
  y = round(nrow(pheno)*tp)
  testing_ord = sample(nrow(pheno), nrow(pheno))
  testing_idx_ = list()
  for(i in 1:round(1/tp)){
    testing_idx_[[i]] = testing_ord[x:y]
    x = y + 1
    y = round(y + (nrow(pheno)*tp))
    if(x > length(testing_ord)) break
    if(y > length(testing_ord) | i == ((round(1/tp))-1)) y = length(testing_ord)
  }
  print(paste("Check for correct training/testing sets passed?", all(sort(unlist(testing_idx_)) == seq(1:nrow(pheno)))))
  write_pheno(samples = pheno$sampleID, y = pheno_vals, testing_idx = testing_idx_, ph_nm = ph, fldr = "XFCV")
  
  print("Commencing loso")
  clusts = sort(unique(pheno$cluster)) # Fastbaps clusters should be provided (saved to pheno for convenience)
  print(paste("Identified", length(clusts), "Clusters"))
  testing_idx_ = list()
  idx = 1
  for(i in clusts){
    testing_idx_[[idx]] = which(pheno$cluster == i)
    idx = idx+1
  }
  print(paste("Check for correct training/testing sets passed?", all(sort(unlist(testing_idx_)) == seq(1:nrow(pheno)))))
  write_pheno(samples = pheno$sampleID, y = pheno_vals, testing_idx = testing_idx_, ph_nm = ph, fldr = "LOSO")
}
