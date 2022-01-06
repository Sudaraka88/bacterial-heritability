if(rstudioapi::isAvailable()) setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # WORKING DIRECTORY
# Checked 20220205
library(fastbaps)


makefastbapsclust = function(fasta.file.name){
  t = Sys.time()
  # fasta.file.name = "cd+acute_acc+core_isolates.fasta"
  sparse.data <- import_fasta_sparse_nt(paste(fasta.file.name, ".fasta", sep = ""))
  saveRDS(sparse.data, paste(fasta.file.name, "_sparse_data_unoptimised.rds", sep = ""))
  print(paste("Step 1", Sys.time() - t ))
  
  sparse.data <- optimise_prior(sparse.data, type = "optimise.symmetric")
  saveRDS(sparse.data, paste(fasta.file.name, "_data_optimised_symmetric.rds", sep = ""))
  print(paste("Step 2", Sys.time() - t ))
  
  baps.hc <- fast_baps(sparse.data, n.cores = 8L)
  saveRDS(baps.hc, paste(fasta.file.name, "_baps_hc.rds", sep = ""))
  print(paste("Step 3", Sys.time() - t ))
  
  best.partition <- best_baps_partition(sparse.data, baps.hc)
  saveRDS(best.partition, paste(fasta.file.name, "_best.partition.rds", sep = ""))
  print(paste("Step 4", Sys.time() - t ))
}

makefastbapsclust("cd_isolates") # not provided
makefastbapsclust("mic_isolates") # not provided
