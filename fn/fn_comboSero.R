# # Function to combine 15B/C, 6A/C and remove ""
# comboSero <- function(sero_list, latex = F){
#   sero_list <- replace(sero_list, sero_list=="15B", "15B/C")
#   sero_list <- replace(sero_list, sero_list=="15C", "15B/C")
#   sero_list <- replace(sero_list, sero_list=="6A", "6A/C")
#   sero_list <- replace(sero_list, sero_list=="6C", "6A/C")
#   if (latex == T) vec_ret <- na.omit(unique(sero_list[sero_list != "NT" & sero_list !=""])) # Latex NT issue
#   else vec_ret <- na.omit(unique(sero_list[sero_list !=""]))
#   return(vec_ret)
# }

# Function to combine 15B/C, 6A/C and remove ""
comboSero <- function(who_list, latex_list=''){
  if('NT' %in% latex_list) latex_list = latex_list[-which(latex_list=='NT')]
  sero_list = c(who_list,latex_list)
  sero_list <- replace(sero_list, sero_list=="15B", "15B/C")
  sero_list <- replace(sero_list, sero_list=="15C", "15B/C")
  sero_list <- replace(sero_list, sero_list=="6A", "6A/C")
  sero_list <- replace(sero_list, sero_list=="6C", "6A/C")
  vec_ret <- na.omit(unique(sero_list[sero_list !=""]))
  return(vec_ret)
}

