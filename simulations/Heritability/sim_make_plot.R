if(rstudioapi::isAvailable()) setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # WORKING DIRECTORY
library(ggplot2)
library(dplyr)

mfolder_ = c("RUN04","RUN07") #getwd()
geno_type = factor(c("Low-LD (lgtRate=0.2)", "High-LD (lgtRate=0.1)"), levels = c("Low-LD (lgtRate=0.2)", "High-LD (lgtRate=0.1)"))
plt_df = c()
cnt = 1
for (mfolder in mfolder_){
  ldsc_out = readRDS(file.path(mfolder, "ldsc_out.rds"))
  
  h2diff = c()
  for(i in 1:nrow(ldsc_out)){
    h2diff[i] = abs(ldsc_out$h2[i] - as.numeric(unlist(strsplit(ldsc_out$folder[i], split = "h"))[2]))
  }
  ldsc_out = cbind(ldsc_out, h2dif = h2diff)
  ldsc_out_ = ldsc_out[which(ldsc_out$PCAs==1), ]
  
  enet_out = readRDS(file.path(mfolder, "enet_out.rds"))
  
  lmm_out = readRDS(file.path(mfolder, "lmm_out.rds"))
  
  plt_df = rbind(plt_df,
                 data.frame(h2 = ldsc_out_$h2, phen = ldsc_out_$phen, folder = ldsc_out_$folder, Method = "LDSC", RUN = geno_type[cnt]), 
                 data.frame(enet_out, Method = "wg-enet", RUN = geno_type[cnt]),
                 data.frame(lmm_out, Method = "LMM", RUN = geno_type[cnt]))
  cnt = cnt + 1
}

plt_df$h2[which(plt_df$h2>1)] = 1
plt_df$h2[which(plt_df$h<0)] = 0


plt_df$h2tru = sapply(plt_df$folder, function(x) as.numeric(unlist(strsplit(x, "h"))[2]))
cor(plt_df$h2, plt_df$h2tru)

plt = ggplot(plt_df, aes(x = factor(h2tru), y = h2, col = Method))  +  geom_abline(slope = 0.1, intercept = 0, col = "grey", lwd = 1.1) +
  stat_summary(fun = mean, geom="point", size=2, shape = 19) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.5) +
  geom_hline(yintercept = seq(0,10)/10, linetype = "dotted", color = "gray80", lwd = 1.2) +
  theme_linedraw(base_size = 18) +
  scale_y_continuous(breaks= seq(0,10)/10, labels= as.character(seq(0,10)/10)) + 
  ylab(expression("Estimated "*h^2*"")) + xlab(expression("Simulated "*h^2*"")) +
  theme(legend.position = "bottom") +
  facet_wrap(~RUN, ncol = 1)

plt

ggsave("h2sim_c.png", plot = plt, width = 8, height = 12)

# Overall error in estimate
mean_se(plt_df$h2[which(plt_df$Method=="LDSC")] - plt_df$h2tru[which(plt_df$Method=="LDSC")])
mean_se(plt_df$h2[which(plt_df$Method=="wg-enet")] - plt_df$h2tru[which(plt_df$Method=="wg-enet")])
mean_se(plt_df$h2[which(plt_df$Method=="LMM")] - plt_df$h2tru[which(plt_df$Method=="LMM")])

# Error in low-LD
mean_se(plt_df$h2[which(plt_df$Method=="LDSC" & plt_df$RUN == "Low-LD (lgtRate=0.2)")] - plt_df$h2tru[which(plt_df$Method=="LDSC" & plt_df$RUN == "Low-LD (lgtRate=0.2)")])
mean_se(plt_df$h2[which(plt_df$Method=="wg-enet"& plt_df$RUN == "Low-LD (lgtRate=0.2)")] - plt_df$h2tru[which(plt_df$Method=="wg-enet"& plt_df$RUN == "Low-LD (lgtRate=0.2)")])
mean_se(plt_df$h2[which(plt_df$Method=="LMM"& plt_df$RUN == "Low-LD (lgtRate=0.2)")] - plt_df$h2tru[which(plt_df$Method=="LMM"& plt_df$RUN == "Low-LD (lgtRate=0.2)")])

# Error in high-LD
mean_se(plt_df$h2[which(plt_df$Method=="LDSC" & plt_df$RUN == "High-LD (lgtRate=0.1)")] - plt_df$h2tru[which(plt_df$Method=="LDSC" & plt_df$RUN == "High-LD (lgtRate=0.1)")])
mean_se(plt_df$h2[which(plt_df$Method=="wg-enet"& plt_df$RUN == "High-LD (lgtRate=0.1)")] - plt_df$h2tru[which(plt_df$Method=="wg-enet"& plt_df$RUN == "High-LD (lgtRate=0.1)")])
mean_se(plt_df$h2[which(plt_df$Method=="LMM"& plt_df$RUN == "High-LD (lgtRate=0.1)")] - plt_df$h2tru[which(plt_df$Method=="LMM"& plt_df$RUN == "High-LD (lgtRate=0.1)")])


error_est = c()
error_ciL = c()
error_ciH = c()
for(i in unique(plt_df$RUN)){
  for(j in unique(plt_df$Method)){
    for(k in unique(plt_df$h2tru)){
      error_est = c(error_est, mean(plt_df$h2[which(plt_df$h2tru==k & plt_df$Method == j & plt_df$RUN == i)] - plt_df$h2tru[which(plt_df$h2tru==k & plt_df$Method == j & plt_df$RUN == i)]))
      error_ciL = c(error_ciL, unlist(unname(mean_se(plt_df$h2[which(plt_df$h2tru==k & plt_df$Method == j & plt_df$RUN == i)] - plt_df$h2tru[which(plt_df$h2tru==k & plt_df$Method == j & plt_df$RUN == i)])[2])))
      error_ciH = c(error_ciH, unlist(unname(mean_se(plt_df$h2[which(plt_df$h2tru==k & plt_df$Method == j & plt_df$RUN == i)] - plt_df$h2tru[which(plt_df$h2tru==k & plt_df$Method == j & plt_df$RUN == i)])[3])))
    }
  }
}
error_est = matrix(error_est, ncol = 9)
error_ciL = matrix(error_ciL, ncol = 9)
error_ciH = matrix(error_ciH, ncol = 9)

colnames(error_est) = rep(c("LDSC","wg-enet","LMM"), 3)
rownames(error_est) = paste("h", seq(1:9)/10, sep = "")

colnames(error_ciL) = rep(c("LDSC","wg-enet","LMM"), 3)
rownames(error_ciL) = paste("h", seq(1:9)/10, sep = "")

colnames(error_ciH) = rep(c("LDSC","wg-enet","LMM"), 3)
rownames(error_ciH) = paste("h", seq(1:9)/10, sep = "")


corrplot::corrplot(corr = error_est, is.corr = F, method = "shade", addCoef.col = T) 

