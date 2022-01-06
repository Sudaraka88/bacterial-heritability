if(rstudioapi::isAvailable()) setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # WORKING DIRECTORY
# Checked 20220105

# Load Packages
library(dplyr) 
library(msm)
require(minqa)
require(data.table)
library(foreach)
library(doParallel)
# Files in RAW folder can be obtained from the authors of Chewapreecha et al. 2014 (https://www.nature.com/articles/ng.2895)


# Set Working Directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # WORKING DIRECTORY
# Read observation data in
routine_data <- dplyr::as_tibble(read.csv("RAW/20151026_maela_rout_sweep_data.csv",header = T, stringsAsFactors = F))
#DupRows: 410 801 1113 1397 1531 1630 1637 1818 2036 2215 2411 3073 3087 3274 4491 4674 4719 4922 5067 6387 6389 7593 7838 7950 8324 8489 8551 8567
#RemRows: 409 800 1112 1396 1531 1629 1636 1818 2036 2214 2410 3072 3086 3273 4490 4673 4718 4921 5066 6386 6388 7592 7837 7950 8324 8488 8550 8566
routine_data <- routine_data[-c(409,800,1112,1396,1531,1629,1636,1818,2036,2214,2410,3072,3086,3273,4490,4673,4718,4921,5066,6386,6388,7592,7837,7950,8324,8488,8550,8566),]

immunology_data <- dplyr::as_tibble(read.csv("RAW/20151026_maela_imm_who_data.csv",header = T, stringsAsFactors = F ))
#DupRows: 11 332 912 930 1143 1757 4103 5679 5681 5797 6341 6696
#RemRows: 10 332 911 929 1142 1757 4103 5678 5680 5796 6340 6695
immunology_data <- immunology_data[-c(10,332,911,929,1142,1757,4103,5678,5680,5796,6340,6695),]

source('fn/fn_viewOutputs.R')
#T = table(immunology_data$whoserotype1); viewDistro(T[-1])
#T = table(routine_data$whoserotype1); viewDistro(T[-1])

# Make Serotype list
obsSeroList = unique(na.omit(unlist(c(immunology_data[,10:12], routine_data[,c(10:13,15:17)]))))
# 6A/6C are combined (as cannot be separated by sweep). NTs are treated separately
# 15B/C spontaneously intraconvert, so removed and added later
rm_var = c("","6A","6C","15B","15C","NT")
for (i in 1:length(rm_var)) obsSeroList = obsSeroList[obsSeroList != rm_var[i]]
obsSeroList <- c(obsSeroList, "15B/C")

# Added the acute column 20201012
# Format data
source('fn/fn_comboSero.R')
# Immunology Data
immunology_data %<>% 
  filter(category == "Infant") %>%
  mutate(collection = "immunology") %>%
  mutate(specdate=as.Date(specdate,"%d/%m/%Y"),sampleday=as.numeric(specdate-min(specdate))) %>% 
  group_by(codenum) %>%
  #distinct(age_d, .keep_all = TRUE) %>% # Remove samples obtained on the same day # THIS WON'T DO!
  rowwise() %>%
  mutate(serotype=paste(comboSero(c(whoserotype1, whoserotype2, whoserotype3)), collapse = ",")) %>%
  ungroup() %>%
  arrange(codenum, sampleday) %>%
  dplyr::select(codenum, collection, age_d, specdate, sampleday, whopnc, serotype, acute) %>%
  distinct()

# Routine Data
routine_data %<>%
  mutate(collection = "routine") %>%
  mutate(specdate=as.Date(specdate,"%d/%m/%Y"),sampleday=as.numeric(specdate-min(specdate))) %>% 
  group_by(codenum) %>%
  #distinct(age_d, .keep_all = TRUE) %>% # Remove samples obtained on the same day
  rowwise() %>%
  mutate(serotype=paste(comboSero(c(whoserotype1,whoserotype2,whoserotype3),c(sweepserotype1,sweepserotype2,sweepserotype3,sweepserotype4)),collapse=",")) %>%
  ungroup() %>%
  group_by(codenum, sampleday) %>%
  filter(n() == 1 | whopnc != "") %>% # Identifies cases with two observations, and takes better of two
  ungroup() %>%
  arrange(codenum, sampleday) %>%
  dplyr::select(codenum, collection, age_d, specdate, sampleday, whopnc, serotype, acute) %>%
  distinct()

all_obs <- bind_rows(immunology_data, routine_data) # Combine cohorts

# Let's convert acute to a 1/0 variable
all_obs$acute[which(all_obs$acute != "Yes")] = 0
all_obs$acute[which(all_obs$acute == "Yes")] = 1

# Normalize Times
time_dev <- sd(all_obs$sampleday)
all_obs <- mutate(all_obs, normday=sampleday/time_dev)

# Previous Carry
prev_carry = rep(0,nrow(all_obs))
prv_id = ""
for (i in 1:nrow(all_obs)){
  id = all_obs$codenum[i]
  if (prv_id!=id) carry=0 # Has ID changed? New ID start with carry 0 
  if(all_obs$serotype[i]!="") carry=1 # Does this entry show carry?
  prev_carry[i]=carry
  prv_id = id
}
all_obs %<>% mutate(carried=prev_carry) %>% group_by(codenum) %>% filter(n()>1) # Remove entries with just one obs.

# Create data for NT serotype
source('fn/fn_createDummy2S.R')
NT_obs_with_state <- filter(all_obs, whopnc != "") %>% # whopnc needs to be done to confirm NT
  group_by(codenum) %>%
  filter(n() > 1) %>% # Need more than one observation
  rowwise() %>% mutate(state=createDummy2S("NT",serotype))

# HMM Models

source('fn/fn_fitMSM.R')
# setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]) #not to overload your computer
registerDoParallel(cl)
 
modOut_noNT <- foreach (i=1:length(obsSeroList), .combine = 'cbind', .packages = c("dplyr","msm","minqa")) %dopar% {
  sero <- obsSeroList[i]
  obs_with_state <- rowwise(all_obs) %>% mutate(state=createDummy2S(sero, serotype))
  temp <- fitMSM(obs_with_state)
}
stopCluster(cl)
name_list = as.character()
for (i in obsSeroList) name_list = c(name_list, paste(i,'M'),paste(i, 'T'),paste(i, 'V'))
names(modOut_noNT) = name_list

if(!file.exists("TempData")) dir.create("TempData") # Output paths will be saved to a folder named TempData
saveRDS(modOut_noNT,file="TempData/modOut_noNT.Rdata")

 
# Run HMM for NT - DO NOT RE-RUN BAD MODEL
# temp <- fitMSM(NT_obs_with_state)
# modOut_wNT <- cbind(modOut_noNT, temp)
# e21 observation is elevated - high FN count

# Re-run HMM for NT with e21 = 0.12 ~~~ Why 0.12? ~~~~
temp = fitMSM(NT_obs_with_state, emI = rbind(c(0,0),c(0.12,0)), fp=c(3))
modOut_wNT = cbind(modOut_noNT, temp)
name_list = c(name_list, 'NT M','NT T','NT V')
names(modOut_wNT) = name_list[1:length(modOut_wNT)]
saveRDS(modOut_wNT,file="TempData/modOut_wNT.Rdata")
# Observing Outputs
goodSeroList = c('19F','23F','6A/C','6B','14','NT','15B/C') # Why not 15B/C?
e21 = e21_Out(modOut_wNT,goodSeroList); viewDistro(e21[goodSeroList],-0.03); abline(h=e21[goodSeroList])
t21 = t21_Out(modOut_wNT,goodSeroList); viewDistro(t21[goodSeroList],-0.03); abline(h=t21[goodSeroList])
c21 = counts21_Out(modOut_wNT,goodSeroList); viewDistro(c21[goodSeroList],-30); abline(h=c21[goodSeroList])
hess = h_Out(modOut_wNT,goodSeroList); paste(goodSeroList, hess)


# Create Viterbi Paths for bad serotypes using 19F
modOut_wNT = readRDS('TempData/modOut_wNT.Rdata') # Can run from here using temp data
modelRef = modOut_wNT$`19F M` # Reference model
badSeroList = obsSeroList[!(obsSeroList%in%goodSeroList)]
foreach(i=1:length(badSeroList)) %do% {
  modelRef$data = modOut_wNT[[paste(badSeroList[i], 'M')]]$data
  modOut_wNT[[paste(badSeroList[i], 'V')]] = viterbi.msm(modelRef)
}

# Save paths, models and tables separately
paths = list()
models = list()
tables = list()
for (i in c(obsSeroList,'NT')){
  paths[[i]] = modOut_wNT[[paste(i,'V')]]
  models[[i]] = modOut_wNT[[paste(i,'M')]]
  tables[[i]] = modOut_wNT[[paste(i,'T')]]
}
if(!file.exists("Out")) dir.create("Out") # Output paths will be saved to a folder named Out
saveRDS(paths,file="Out/paths.Rdata")
saveRDS(models,file="Out/models.Rdata")
saveRDS(tables,file="Out/tables.Rdata")
