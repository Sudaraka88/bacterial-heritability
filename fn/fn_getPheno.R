# iterate through path to get length
getDateBoundary <- function(obs_row, codenum, path, obs,direction=-1)
{ # -1 - towards start, +1 - towards end
  boundary <- obs_row # This is >0
  spec_dates <- c()
  # keep looping till beginning/end of data (NA) or until we go over the the next subject
  spec_dates <- c(spec_dates, obs$specdate[boundary])
  while(!is.na(path$subject[boundary+direction]) & path$subject[boundary+direction] == codenum)
  {    
    # Always presume a positive swab is carriage (sometimes Viterbi will output no carriage
    # due to low init and/or transition probs)
    # Reduce or increase the obs_row until we find the edge of carriage episode
    if(path$fitted[boundary+direction] == 2 || path$observed[boundary+direction] == 2) {
      boundary <- boundary + direction
      spec_dates <- c(spec_dates, obs$specdate[boundary])
    }
    else break # We are at the episode edge
    if(boundary+direction==0) { # Beginning of data -1
      break # path$subject[0] does not return NA
    }
  }
  # If we are not at beginning/end of data/subject, correct by +/- 0.5 specdate to boundary
  if(boundary+direction>0) { # Not at the beginning of the dataset
    if(!is.na(path$subject[boundary+direction]) & path$subject[boundary+direction] == codenum) { # Not at the end or the edge of a subject
      date_boundary <- obs$specdate[boundary+direction] +(obs$specdate[boundary]-obs$specdate[boundary+direction])/2
    } else {
      # date_boundary <- obs$specdate[boundary] # We are at the beginning of the dataset, have to use the boundary as is - Original method
      ### Modified 01042020 ###
      if(direction == -1) date_boundary <- obs$specdate[boundary] - 15 # We are at the end of data or beginning/end subject, have to use the boundary as is
      else date_boundary <- obs$specdate[boundary] + 15 # We are at the end of data or beginning/end subject, have to use the boundary as is
      #########################
    }
  } else {
    # date_boundary <- obs$specdate[boundary] # We are at the beginning of the dataset, have to use the boundary as is
    ### Modified 01042020 ###
    if(direction == -1) date_boundary <- obs$specdate[boundary] - 15 # We are at the end of data or beginning/end subject, have to use the boundary as is
    else date_boundary <- obs$specdate[boundary] + 15 # We are at the end of data or beginning/end subject, have to use the boundary as is
    #########################
  }
  if(direction==-1){# Going towards start
    if(boundary + direction == 0) carried = 0 # Beginning of dataset, there can be no previous carry
    else if(!is.na(path$subject[boundary+direction]) & path$subject[boundary+direction] == codenum) carried = obs$carried[boundary+direction]
    else carried = 0
    class(spec_dates) = "Date"
    return(list("start"=date_boundary,"carried"=carried, "specdates"=spec_dates))
  }
  else {
    class(spec_dates) = "Date"
    return(list("end"=date_boundary, "specdates"=spec_dates))
  }
}

getPheno <- function(obs,alt_obs,code,specdate,path,i,sero='NT'){
  obs_row = which(obs$codenum==code & obs$specdate==specdate)
  if(length(obs_row)>1) print(paste('Multiple observations idx:',i, 'obs:', obs_row, 'codenum:',code,'specdate:',specdate)) 
  # What is the best approach here?
  obs_sero = unlist(strsplit(obs$serotype[obs_row], split = ',')) # sero list in observations 
  
  if(!(sero %in% obs_sero)) { # sero in metadata mismatches with observation, what is the best approach here?
    print(paste('Sero mismatch idx:', i, 'obs:', obs_row, 'codenum:', code, 'specdate:', specdate, 'sero:', sero, 'obs_sero:', obs$serotype[obs_row]))
    # We can assume that serotype_metadata is wrong, and if possible, change serotype to match observed
    # Match if only 1 observed serotype is present, if multiple are present, cannot pick one randomly!
    if(length(obs_sero == 1)){ # There is exactly 1 mismatching sero
      if (sero=='NT') { # If the original metadata sero was NT, we have to change to all_obs list
        obs_row = which(alt_obs$codenum==code & alt_obs$specdate==specdate)
        sero = obs_sero
        obs = alt_obs
        print(paste('Substituted with sero:', sero, 'obs:', obs_row))
      }
      if('NT'%in%obs_sero){ # If the obs_sero is NT and there's a mismatch, original sero was from all_obs, change to NT_obs 
        obs_row = which(alt_obs$codenum==code & alt_obs$specdate==specdate)
        sero = obs_sero
        obs = alt_obs
        print(paste('Substituted with sero:', sero, 'obs:', obs_row))
      }
      # If obs_sero and sero are both not 'NT', we are already looking in the right list
    }
    else { # Multiple mismatching obs_seros or none
      print(paste('No or multiple alternatives, returning NAs.'))
      return(list("age"=NA,"pheno"=NA,"carried"=NA)) # If nothing works, return NAs/0s
    }
  }
  age = as.numeric(obs$age_d[obs_row])
  startANDcarried = getDateBoundary(obs_row=obs_row,codenum=code,obs=obs,path=path)
  start = startANDcarried$start
  carried = startANDcarried$carried
  end_ = getDateBoundary(obs_row=obs_row,codenum=code,obs=obs,path=path,direction = 1)
  end = end_$end
  pheno = as.numeric(end-start)
  return(list("age"=age,"pheno"=pheno,"carried"=carried))
}