if(rstudioapi::isAvailable()) setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # WORKING DIRECTORY
# Checked 20210621
# This code is intended to make a separate phenotype file for mic data

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # WORKING DIRECTORY
seq_data <- dplyr::as_tibble(read.delim("RAW/sequence_metadata.txt",header=T,stringsAsFactors = F))
AMR_data <- dplyr::as_tibble(read.csv("RAW/ARI_AMR_data.csv",header = T, stringsAsFactors = F ))
AMR_data <- AMR_data[-c(10,4986,4401,819,4481,638,5229,5231,5229),] # Duplicate rows
penmicidx = which(!is.na(AMR_data$Penicillin.MIC..mg.L.))
cefmicidx = which(!is.na(AMR_data$Ceftriaxone.MIC..mg.L.)) # These two are the same


carry_data = read.csv("TempData/carryOp_datefix.csv")

AMR_data = AMR_data[penmicidx,]
# Let's also drop the cultured in error samples 
AMR_data = AMR_data[-which(AMR_data$Cultured.in.error == "Yes"),] # 3504 left

AMR_mapper = data.frame(codenum = AMR_data$codenum, 
                        category = toupper(AMR_data$category),
                        specdate = as.Date(AMR_data$specdate, "%d/%m/%Y"),
                        serotype = AMR_data$serotype)

SEQ_mapper = data.frame(codenum = seq_data$codenum, 
                        category = toupper(seq_data$category), 
                        specdate = as.Date(seq_data$specdate, "%d-%b-%y"), 
                        serotype = seq_data$serotype)

mapped_idxes = apply(AMR_mapper, 1, function(x) row.match(x, SEQ_mapper))

AMR_data$lanes = seq_data$lane[mapped_idxes]
AMR_data = AMR_data[which(!is.na(AMR_data$lanes)),]

mapped_pheno_AMR = data.frame(sampleID = AMR_data$lanes, 
                              category = AMR_data$category,
                              serotype = AMR_data$serotype,
                              Penicillin = AMR_data$Penicillin,
                              Ceftriaxone = AMR_data$Ceftriaxone,
                              Chloramphenicol = AMR_data$Chloramphenicol,
                              Clindamycin = AMR_data$Clindamycin,
                              Erythromycin = AMR_data$Erythromycin,
                              Sulpha.trimethoprim = AMR_data$Sulpha.trimethoprim,
                              Tetracycline = AMR_data$Tetracycline,
                              Penicillin.MIC = AMR_data$Penicillin.MIC..mg.L.,
                              Ceftriaxone.MIC = AMR_data$Ceftriaxone.MIC..mg.L.,
                              acute = AMR_data$acute)

# Let's try to map the carried variable
sep = function(ln){
  if(is.na(ln)) return(rep(NA,9))
  tmp = unlist(strsplit(ln, ','))
  return(c(tmp, rep(NA, (9-length(tmp)))))
}

match = function(y){
  tmp = apply(all_lanes, 2, function(x) which(y == x))
  if(length(tmp) == 0) return(NA) else return(unname(tmp[which(sapply(tmp, length) == 1)]))
  
}

all_lanes = data.frame(matrix(rep(NA,4392*9), ncol = 9))
for(i in 1:4392){
  all_lanes[i,] = sep(carry_data$lanes[i])
}
# Let's do a column by colum search now
bla = unname(unlist(sapply(AMR_data$lanes, function(y) match(y))))


write.table(mapped_pheno_AMR, "Out/MIC_phenoData.csv", row.names = FALSE)
