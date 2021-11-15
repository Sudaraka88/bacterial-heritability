if(rstudioapi::isAvailable()) setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # WORKING DIRECTORY
library(ape)
alignMx = function(mx, order){
  matches = c()
  for(i in order) matches[i] = which(rownames(mx) %in% i)
  checkAlignment(rownames(mx)[matches], order)
  return(mx[matches, matches])
}
checkAlignment = function(nm1,nm2){
  print(paste("Pass alignment test?", all(nm1 == nm2)))
}
folder_ = dir(pattern = "^[:h:][:0:]")
for(folder in folder_){
  print(paste("Folder:", folder))
  # folder = "h0.3"
  fasta_ = paste(folder, "/simulations/genSim/genSim", sep = "")
  phylo = read.tree(paste(folder, "/simulations/genSim/phylogeny.nwk", sep = "")) # 0 is labelled as zero, change first
  y = read.table(paste(folder, "/simulations/phenSim/0/phenSim.phen", sep = ""))$V3
  
  phylo$tip.label[which(phylo$tip.label == "zero")] = "0"
  
  D = cophenetic.phylo(phylo)
  
  pheno = data.frame(sampleID = 0:(length(y)-1), y = y)
  
  D = alignMx(D, as.character(pheno$sampleID))
  saveRDS(D, paste(fasta_, "_phylo_dist.rds", sep = ""))
  
  node_depths = node.depth.edgelength(phylo)[(length(phylo$tip.label)+1):(dim(phylo$edge)[1]+1)]
  names(node_depths) = as.character((length(phylo$tip.label)+1):(dim(phylo$edge)[1]+1))
  
  sMx = matrix(nrow = length(phylo$tip.label), ncol = length(phylo$tip.label))
  for(i in 1:length(phylo$tip.label)){
    for(j in 1:length(phylo$tip.label)){
      mrca_node = getMRCA(phylo, c(phylo$tip.label[i], phylo$tip.label[j]))
      sMx[i,j] = sMx[j,i] = node_depths[which(names(node_depths) %in% mrca_node)]
    }
  }
  colnames(sMx) = rownames(sMx) = phylo$tip.label
  sMx = alignMx(sMx, as.character(pheno$sampleID))
  saveRDS(sMx, paste(fasta_, "_phylo_sim.rds", sep = ""))
  
}
