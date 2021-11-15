if(rstudioapi::isAvailable()) setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # WORKING DIRECTORY
library(lme4qtl)

geth2est = function(pheno_data, sim_mx){
  mod = relmatLmer(y ~ (1|ids), pheno_data, relmat = list(ids = sim_mx), calc.derivs = FALSE)
  h2 = VarProp(mod)$prop[1]
  print(paste("h2 = ", h2))
  return(h2)
}

TT = Sys.time()
out_full_full = data.frame()
folder_ = dir(pattern = "^[:h:][:0:]")[8:9]

for(folder in folder_){
  if(file.exists(paste(folder, "/lmm.rds", sep = ""))){
    print(paste("File", paste(folder, "/lmm.rds", sep = ""), "exists"))
    out_full = readRDS(paste(folder, "/lmm.rds", sep = ""))
  } else {
    t0 = Sys.time()
    out_full = data.frame()
     sim_mx = readRDS(paste(folder, "/simulations/genSim/genSim_phylo_sim.rds", sep = ""))
    i_ = length(dir(paste(folder, "/simulations/phenSim/", sep = ""))) - 1
    print(paste("Found", i_+1, "phenos"))
    for(i in 0:i_){
      print(paste("Pheno", i))
      if(file.exists(paste(folder, "/simulations/phenSim/",i,"/lmm_out.rds", sep = ""))){
        print(paste("File", paste(folder, "/simulations/phenSim/",i,"/lmm_out.rds", sep = ""), "exists"))
        out_df = readRDS(paste(folder, "/simulations/phenSim/",i,"/lmm_out.rds", sep = ""))
      } else{
        y = read.table(paste(folder, "/simulations/phenSim/",i,"/phenSim.phen", sep = ""))$V3
        pheno_data = data.frame(ids = as.character(0:999), y = y)
        
        out_df = geth2est(pheno_data, sim_mx)
        out_df = cbind(h2 = out_df, phen = i)
        
        saveRDS(out_df, paste(folder, "/simulations/phenSim/",i,"/lmm_out.rds", sep = ""))
      }
      
      out_full = rbind(out_full, out_df)
    }
    out_full = cbind(out_full, folder)
    saveRDS(out_full, paste(folder, "/lmm.rds", sep = ""))
    print(Sys.time() - t0)
  }
  
  out_full_full = rbind(out_full_full, out_full)
  
}
saveRDS(out_full_full, "lmm_out.rds")
print(Sys.time() - TT)
