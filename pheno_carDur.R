if(rstudioapi::isAvailable()) setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # WORKING DIRECTORY
# Checked 20210621

# requires paths - Run code to generate paths.Rdata first (pheno_pathGen.R)
# Read sequence metadata, format and save as Rdata
library(dplyr)
library(data.table)

seq_data <- dplyr::as_tibble(read.delim("RAW/sequence_metadata.txt",header=T,stringsAsFactors = F))
seq_data%<>%
  mutate(specdate=as.Date(seq_data$specdate,"%d-%b-%y")) %>%
  filter(category=="INFANT") %>%
  dplyr::select(lane, specdate, codenum, serotype) %>%
  arrange(codenum, specdate) %>%
  filter(codenum %in% all_obs$codenum)

paths = readRDS('Out/paths.Rdata')
source('fn/fn_getPheno.R')

pathcopy <- paths
carry_output = data.frame()

for(serotype in names(pathcopy)){
  # If the obsreved value is 2, then the fitted value must be 2 (no false positives). 
  # Why don't we just do this to start with? Also, is this a good idea?
  pathcopy[[serotype]][pathcopy[[serotype]][,3]==2, 4] = 2
  
  # Adding lag/lead fitted states (below and above)
  pathcopy_new <- as.data.table(pathcopy[[serotype]])
  pathcopy_new[, ('lag') := shift(fitted, 1), by = subject]
  pathcopy_new[, ('lead') := shift(fitted, 1, type = 'lead'), by = subject]
  
  # Find start and end points in one go! WOW! I didn't think of this?
  start_points <- which(pathcopy_new$fitted==2 & (pathcopy_new$lag != 2 | is.na(pathcopy_new$lag)))
  end_points <- which(pathcopy_new$fitted==2 & (pathcopy_new$lead != 2 | is.na(pathcopy_new$lead)))
  
  lane_count = 0 # Let's count how many lanes get mapped
  if(serotype != 'NT') obsframe = all_obs else obsframe <- NT_obs_with_state
  
  codenum_out <- obsframe[start_points,'codenum']
  age <- obsframe[start_points,'age_d'] # Is this correct for NT?
  
  # print(paste(serotype, length(start_points)))
  
  if(length(start_points>0)){
    carried <- length <- lanes <- spec_dates <- rep(0, length(start_points))
    starts <- ends <- c()

    for(i in 1:length(start_points)){
      codenum = as.character(obsframe[start_points[i],'codenum'])
      # carried[i] = getCarried(start_points[i], codenum, pathcopy[[serotype]], obsframe)
      start <- getDateBoundary(start_points[i], codenum, pathcopy[[serotype]], obsframe, direction=-1)
      end <- getDateBoundary(start_points[i], codenum, pathcopy[[serotype]], obsframe, direction=1) 
      # keeping start_points here to find the end is not a mistake. -1 scans upwards, towards the start and +1 scans downwards towards the end
      length[i] <- as.numeric(end$end-start$start)
      carried[i] <- start$carried
      # Lane Mapping
      lanes_<- c()
      spec_dates_ <- c() # for output only
      # Try to add the acute phenotype - 20201012
      # Acute is only needed if there's a lane mapped to it, easier to edit here
    
      for(specs in 1:length(start$specdates)) {
        spec_dates_ = c(spec_dates_, start$specdates[specs])
        lanes_ = c(lanes_,seq_data$lane[which(seq_data$specdate == start$specdates[specs] & seq_data$serotype == serotype & seq_data$codenum == codenum)])
      }
      for(specs in 1:length(end$specdates)) {
        spec_dates_ = c(spec_dates_, end$specdates[specs])
        lanes_ = c(lanes_,seq_data$lane[which(seq_data$specdate == end$specdates[specs] & seq_data$serotype == serotype & seq_data$codenum == codenum)])
      }
      if(length(unique(lanes_)) == 1) {
          lanes[i] = unique(lanes_)
          lane_count = lane_count + length(unique(lanes_))
          spec_dates[i] = paste(unique(as.character(zoo::as.Date(unique(spec_dates_)))), collapse = ",")
        } else if(length(unique(lanes_)) == 0) {
          lane_count = lane_count + length(unique(lanes_))
          lanes[i] = NA
          spec_dates[i] = paste(unique(as.character(zoo::as.Date(unique(spec_dates_)))), collapse = ",")
        } else {
          lane_count = lane_count + length(unique(lanes_))
          lanes[i] = paste(unique(lanes_), collapse = ",")
          spec_dates[i] = paste(unique(as.character(zoo::as.Date(unique(spec_dates_)))), collapse = ",")
        }
      # For output only
      starts <- c(starts, as.Date(start$start, "%y-%m-%d"))
      ends <- c(ends, as.Date(end$end, "%y-%m-%d"))
    }
    class(starts) = "Date"
    class(ends) = "Date"
    temp_df = data.frame(codenum = codenum_out, age = age, carried = carried, serotype = rep(serotype, length(start_points)), start = starts, end = ends, length = length, lanes = lanes, specdate = spec_dates)
    print(serotype)
    carry_output = rbind(carry_output, temp_df)
  }
}

#### Add the acute phenotype ###
# All carry_outputs must match to an observation in all_obs
acute_mapped = acute_all = rep(NA, dim(carry_output)[1])
for(i in 1:dim(carry_output)[1]){
  dates = unlist(strsplit(carry_output$specdate[i], ","))
  # a = which(carry_output$age_d[i] == all_obs$age_d) # Age will vary with disease duration, cannot match
  b = which(carry_output$codenum[i]  == all_obs$codenum)
  # c = which(carry_output$serotype[i] == all_obs$serotype) # A serotype won't be present if not genotyped, acute will still be written
  acute = c() # all acute episodes in order
  for(j in 1:length(dates)){
    d = which(dates[j] == all_obs$specdate)
    x = intersect(b,d)
    acute = c(acute, all_obs$acute[x])
  }
  acute_all[i] = paste(acute, collapse = ",")
  
  lanes = unlist(strsplit(carry_output$lanes[i], ",")) 
  if(!any(is.na(lanes))){ # There are matching lanes
    acute = c()
    for(k in 1:length(lanes)){ # Let's go through each lane
      d = which(seq_data$specdate[which(seq_data$lane == lanes[k])] == all_obs$specdate)
      x = intersect(b,d)
      acute = c(acute, all_obs$acute[x])
    }
    print(acute)
    acute_mapped[i] = paste(acute, collapse = ",")
  }
  
}
carry_output%<>%
  mutate(accute_all=acute_all) %>%
  mutate(accute_mapped=acute_mapped)

#############################
write.csv(carry_output, "TempData/carryOp_datefix.csv")

# Let's try to reformat carry_output
output = data.frame()

for(i in 1:dim(carry_output)[1]){
  if(!is.na(carry_output$lanes[i])){ # has a matched lane
    lanes = carry_output$lanes[i]
    A = as.numeric(unlist(strsplit(carry_output$accute_mapped[i], ',')))
    X = unlist(strsplit(toString(lanes), ','))
    y = carry_output$length[i]*rep(1,length(X))
    carried = carry_output$carried[i]*rep(1,length(X))
    Cat = rep(c('INFANT'),length(X))
    seros = rep(carry_output$serotype[i], length(X))
    output = rbind(output, data.frame(sampleID = X, category = Cat, serotype = seros, carried = carried, carriageDuration = y, acute = A))
  }
  
}

write.csv(output, "Out/carriageDuration_datefix.csv", row.names = FALSE) # output saved to Out
  