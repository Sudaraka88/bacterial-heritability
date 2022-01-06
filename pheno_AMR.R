if(rstudioapi::isAvailable()) setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # WORKING DIRECTORY
# Checked 20210621

# Run this code to generate the anti-microbial resistance phenotypes
# Load Packages
library(dplyr) 
library(msm)
require(minqa)
require(data.table)
library(foreach)
library(doParallel)
# Try to map lanes to routine and immunology data
# Files in RAW folder can be obtained from the authors of Chewapreecha et al. 2014 (https://www.nature.com/articles/ng.2895)
routine_data <- dplyr::as_tibble(read.csv("RAW/20151026_maela_rout_sweep_data.csv",header = T, stringsAsFactors = F))
routine_data <- routine_data[-c(409,800,1112,1396,1531,1629,1636,1818,2036,2214,2410,3072,3086,3273,4490,4673,4718,4921,5066,6386,6388,7592,7837,7950,8324,8488,8550,8566),] # These sequences are duplicates

immunology_data <- dplyr::as_tibble(read.csv("RAW/20151026_maela_imm_who_data.csv",header = T, stringsAsFactors = F ))
immunology_data <- immunology_data[-c(10,332,911,929,1142,1757,4103,5678,5680,5796,6340,6695),] # These sequences are duplicates

AMR_data <- dplyr::as_tibble(read.csv("RAW/ARI_AMR_data.csv",header = T, stringsAsFactors = F ))
AMR_data <- AMR_data[-c(10,4986,4401,819,4481,638,5229,5231,5229),] # These sequences are duplicates

# Performing basic filtering
seq_data <- dplyr::as_data_frame(read.delim("RAW/sequence_metadata.txt",header=T,stringsAsFactors = F))
seq_data = mutate(seq_data,specdate=as.Date(seq_data$specdate, "%d-%b-%y"))
immunology_data = mutate(immunology_data,specdate=as.Date(immunology_data$specdate,"%d/%m/%Y"))
routine_data = mutate(routine_data,specdate=as.Date(routine_data$specdate,"%d/%m/%Y"))
AMR_data = mutate(AMR_data,specdate=as.Date(AMR_data$specdate,"%d/%m/%Y"))

rout_idx = list()
imun_idx = list()
amr_idx = list()
no_match = c()
amr_no_match = c()
immunology_data$lane = NA
routine_data$lane = NA
AMR_data$lane = NA

for(i in 1:dim(seq_data)[1]){
  rout_idx[[i]] = intersect(intersect(which(seq_data$specdate[i] == routine_data$specdate), which(toupper(seq_data$category[i]) == toupper(routine_data$category))),
                                (which(seq_data$codenum[i]==routine_data$codenum))) # routine
  imun_idx[[i]] = intersect(intersect(which(seq_data$specdate[i] == immunology_data$specdate), which(toupper(seq_data$category[i]) == toupper(immunology_data$category))),
                            (which(seq_data$codenum[i]==immunology_data$codenum))) # immunology
  amr_idx[[i]] = intersect(intersect(intersect(which(seq_data$specdate[i] == AMR_data$specdate), which(toupper(seq_data$category[i]) == toupper(AMR_data$category))),
                           (which(seq_data$codenum[i]==AMR_data$codenum))), which(seq_data$serotype[i] == AMR_data$serotype))
  if(length(amr_idx[[i]]) == 0) amr_no_match = c(amr_no_match, i)
  else{
    if(is.na(AMR_data$lane[amr_idx[[i]]])) AMR_data$lane[amr_idx[[i]]] = seq_data$lane[i]
    else AMR_data$lane[amr_idx[[i]]] = paste(AMR_data$lane[amr_idx[[i]]], seq_data$lane[i], sep = "; ")
  }
  
  if(length(imun_idx[[i]]) == 1){
    if(length(rout_idx[[i]]) == 1) print("Both rout and immuno seem to match")
    else {
      if(is.na(immunology_data$lane[imun_idx[[i]]])) immunology_data$lane[imun_idx[[i]]] = seq_data$lane[i]
      else immunology_data$lane[imun_idx[[i]]] = paste(immunology_data$lane[imun_idx[[i]]], seq_data$lane[i], sep = "; ")
    }
  } 
  else if(length(rout_idx[[i]]) == 1) {
    if(is.na(routine_data$lane[rout_idx[[i]]])) routine_data$lane[rout_idx[[i]]] = seq_data$lane[i]
    else routine_data$lane[rout_idx[[i]]] = paste(routine_data$lane[rout_idx[[i]]], seq_data$lane[i], sep = "; ")
  }
  else {print("No match!"); no_match = c(no_match,i)}
}
    
AMR_phenoData = data.table(sampleID = AMR_data$lane, category = AMR_data$category, serotype = AMR_data$serotype, Penicillin = AMR_data$Penicillin, Ceftriaxone = AMR_data$Ceftriaxone, Chloramphenicol = AMR_data$Chloramphenicol,
                           Clindamycin = AMR_data$Clindamycin, Erythromycin = AMR_data$Erythromycin, Sulpha.trimethoprim = AMR_data$Sulpha.trimethoprim, 
                           Tetracycline = AMR_data$Tetracycline, Penicillin.MIC = AMR_data$Penicillin.MIC..mg.L., Ceftriaxone.MIC = AMR_data$Ceftriaxone.MIC..mg.L.)

AMR_phenoData = AMR_phenoData[!is.na(AMR_phenoData$sampleID),]

ceftriaxoneBlankList = which(AMR_phenoData$Ceftriaxone=="")
for(i in ceftriaxoneBlankList){
  if(AMR_phenoData$Penicillin[i] == "SENSITIVE") AMR_phenoData$Ceftriaxone[i] = "SENSITIVE"
}

no_data_list = c()
for(i in 1:3137){
  if(AMR_phenoData$Penicillin[i] == ""){
    if(AMR_phenoData$Ceftriaxone[i] == ""){
      if(AMR_phenoData$Clindamycin[i] == ""){
        if(AMR_phenoData$Erythromycin[i] == ""){
          if(AMR_phenoData$Sulpha.trimethoprim[i] == ""){
            if(AMR_phenoData$Tetracycline[i] == ""){
              no_data_list = c(no_data_list, i)
            }
          }
        }
      }
    }
  }
}
AMR_phenoData = AMR_phenoData[-no_data_list,]
track = unname(sapply(AMR_phenoData$sampleID, function(x) grepl(";",x)))
idxes = which(track==TRUE)
for(i in idxes){
  temp = strsplit(AMR_phenoData$sampleID[i], '; ')[[1]]
  AMR_phenoData = rbind(AMR_phenoData, AMR_phenoData[i,])
  AMR_phenoData$sampleID[i] = temp[1]
  AMR_phenoData$sampleID[dim(AMR_phenoData)[1]] = temp[2]
}

#ComboSero
L_15BC = union(which(AMR_phenoData$serotype == "15C"), which(AMR_phenoData$serotype == "15B"))
AMR_phenoData$serotype[L_15BC] = "15B/C"
L_6AC = union(which(AMR_phenoData$serotype == "6A"), which(AMR_phenoData$serotype == "6C"))
AMR_phenoData$serotype[L_6AC] = "6A/C"

#Lane Check
MSAlanes = read.delim("RAW/maela3K_lanes.txt", stringsAsFactors = F)
track = sapply(AMR_phenoData$sampleID, function(x) grepl(paste(">",x,"_S_pneumoniae_Spanish_23F.dna",sep=""), MSAlanes))
plot(track)
write.csv(AMR_phenoData, 'Out/AMR_phenoData.csv', row.names = FALSE)
