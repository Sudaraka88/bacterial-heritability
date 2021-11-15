# this is the analysis of true values
if(rstudioapi::isAvailable()) setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # WORKING DIRECTORY
library(ggplot2)

gap_freq = readRDS("cd_core_uqepi_gap_freq.rds")

dirs_ = c("RUN03", "RUN04", "RUN05")
tt = data.frame()
for(dirs in dirs_){
  ph_ = dir(dirs)[grep(pattern = "^phe", dir(dirs))]
  if(length(grep(ph_, pattern = "pkl")) > 0) ph_ = ph_[!(ph_ %in%  ph_[grep(ph_, pattern = "pkl")])] # remove pkl
  
  for (ph in ph_){
    tt_ = readRDS(paste(dirs, "/OUT/max/", ph, "_tt.rds", sep = ""))
    tt = rbind(tt, cbind(tt_, tts = nrow(tt_) , h2 = rep(unlist(strsplit(ph, "_"))[3], nrow(tt_))))
  }
}

tt$tts = as.factor(tt$tts)
tt$h2 = as.factor(tt$h2)

plt_cmp = ggplot(tt, aes(x = -log10(lrt.pvalue), y = -log10(p))) + geom_point(aes(col = gf, shape = Test), size = 3) +
  geom_abline(slope = 1, intercept = 0) + geom_abline(slope = 0, intercept = -log10(0.05/100000)) + geom_vline(xintercept = -log10(0.05/100000)) + 
  xlab("-log10(p) MA test") + ylab("-log10(p) gap/snp test") +
  theme(text = element_text(size=24)) 

plt_cmp




TP_gapsnp = length(which(-log10(tt$p) > -log10(0.05/100000)))/nrow(tt)
TP_ma = length(which(-log10(tt$lrt.pvalue) > -log10(0.05/100000)))/nrow(tt)

