# Observe serotype distro
viewDistro = function(T, y = 0, col = NULL) {
  x = barplot(T, xaxt = "n", col = col)
  labs <- paste(names(T))
  text(cex=1, x=x, y=y, labs, xpd=TRUE, srt=90)
}

# e21 parameter
e21_Out <- function(modOut,view_vec){
  val = 0
  for (i in 1:length(view_vec)) val[i] = ematrix.msm(modOut[[paste(view_vec[i],'M')]])$estimate[2,1]
  names(val) = view_vec
  return(val)
}

# intensity 21 parameter
t21_Out <- function(modOut,view_vec){
  val = 0
  for (i in 1:length(view_vec)) {
    temp = qmatrix.msm(modOut[[paste(view_vec[i],'M')]])
    if (temp$estimate[2,1] < 10) val[i] = temp$estimate[2,1] # Remove high values
    else val[i] = 0
  }
  names(val) = view_vec
  return(val)
}

# table 21 parameter
counts21_Out <- function(modOut,view_vec){
  val = 0
  for (i in 1:length(view_vec)) {
    temp = modOut[[paste(view_vec[i],'T')]]
    if(dim(temp)[2] == 2) val[i] = temp[2,1]
    else val[i] = 0
  }
  names(val) = view_vec
  return(val)
}

#  hessian
h_Out <- function(modOut, view_vec){
  val = FALSE
  for (i in 1:length(view_vec)) val[i] = modOut[[paste(view_vec[i],'M')]]$foundse
  return(val)
}