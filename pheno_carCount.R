if(rstudioapi::isAvailable()) setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # WORKING DIRECTORY
# Checked 20210621

# Identify all carriage durations 
# requires paths - Run code to generate paths.Rdata first (pheno_pathGen.R)
paths = readRDS('Out/paths.Rdata')
full_sero = c(obsSeroList, 'NT')
stend = data.frame()
for (i in 1:length(paths)){
  path = paths[[full_sero[i]]]
  dis = path$observed+path$fitted
  dis_s = rep(FALSE,length(dis))
  dis_s[which(dis>2)] = TRUE
  paths[[full_sero[i]]]$disease = dis_s
  path = paths[[full_sero[i]]]
  temp = carriage_st_end(paths[[full_sero[i]]])
  temp$sero = full_sero[i]
  stend = rbind(stend,temp)
}

# stend_f = filter(stend,type==0) # Remove unsure entry/exit carriage episodes
stend_f = stend
# Find carriage duration
Carry_d = 0 
for (i in 1:dim(stend_f)[1]){
  s = stend_f$start_id[i]
  e = stend_f$end_id[i]
  if(stend_f$sero[i] == 'NT'){
    dte1 = as.Date(NT_obs_with_state$specdate[e+1])
    dts1 = as.Date(NT_obs_with_state$specdate[s+1])
    dte = as.Date(NT_obs_with_state$specdate[e])
    dts = as.Date(NT_obs_with_state$specdate[s])
  } else{
    dte1 = as.Date(all_obs$specdate[e+1])
    dts1 = as.Date(all_obs$specdate[s+1])
    dte = as.Date(all_obs$specdate[e])
    dts = as.Date(all_obs$specdate[s])
  }
  if(length(dts) > 0){
    start_date = 0.5*(dts-dts1) + dts
    end_date =  0.5*(dte1-dte) + dte
    Carry_d[i] = end_date - start_date
  } 
}
stend_f$carry_d = Carry_d


carriage_st_end <- function(path){
  sub = ""
  start <- start_id <- entry <- 0
  end <- end_id <- exit <- 0
  k <- l <- r <- s <- 1
  sub_id = ""
  for (i in 1:dim(path)[1]){
    if(!(sub == path$subject[i])) { # New subject entry point
      if(start>0) { # if there's a start, previous subject was carrying
        end = i-1
        end_id[k] = end
        exit[s] = end # Lower edge
        k = k + 1
        s = s + 1
      }
      start = 0 # Reset start point for new subject
      if(path$disease[i]){
        entry[r] = i
        r = r + 1
      }
    }
    if(path$disease[i] == TRUE){
      if(start==0) { # if start != 0, start already identified
        start = i
        start_id[l] = start
        sub_id[l] = as.character(path$subject[i])
        l = l + 1
      }
    }
    if(path$disease[i] == FALSE){
      if(start > 0){ # if start == 0, no start identified
        end = i-1 # ended in the previous step
        end_id[k] = end
        k = k + 1
        start = 0 # reset start after finding the end
      }
    }
    sub = path$subject[i] # keep track of previous subject
  }
  type = rep(0,length(start_id))
  type[which(start_id%in%entry)] = 1 # entry unsure
  if(path$disease[i] == TRUE){
    end_id[k] = i # Last entry is TRUE
    type[k] = 2 # last exit is unsure
  }
  type[which(end_id%in%exit)] = 2 # exit unsure (can be entry unsure too!)
  start_end_list = data.frame(sub_id,start_id,end_id,type)
  
  return(start_end_list)
}

