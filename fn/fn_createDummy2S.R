# 2 state HMM
createDummy2S <- function(serotype, carrying_serotypes) {
  dummy <- 1
  carrying_serotypes <- unlist(strsplit(carrying_serotypes, split=","))
  if (serotype %in% carrying_serotypes) dummy <- 2
  return(dummy)
}