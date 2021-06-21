# Majority of this code is in clusterSeqs.R
if(rstudioapi::isAvailable()) setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # WORKING DIRECTORY
library(ggtree)
library(phytools)
library(DECIPHER)
library(dplyr)
library(ggplot2)
library(colorspace)
library(ggpubr)
library(gridExtra)


ridx = function(clstnm, sampleID, splt){
  # Function to map cluster labels with sample IDs
   temp = unlist(strsplit(clstnm, split = splt))
  nm = paste(temp[1],"_",temp[2],"#",temp[3], sep = "")
  return(which(sampleID %in% nm))
}
 
# PRE PROCESSING #
# tr = read.tree("cfrm_MIC_new.labelled_tree.newick")
tr = read.tree("MICtree_new.treefile")
# tr = read.tree("MICtree_outlierRmd.treefile")
# tr = read.tree("cfrmOut_MIC_outlierRmd.labelled_tree.newick")
# tr$tip.label = sapply(tr$tip.label, function(x) unlist(strsplit(x, "_S_"))[1])

plot(tr, cex = 0.4)

# These labels correspond to: 6731_1#13, 6755_6#9, 6731_1#21, 6731_2#11 - 20210114
# which(tr$tip.label == "6731_1_13")
# which(tr$tip.label == "6755_6_9")
# which(tr$tip.label == "6731_1_21")
# which(tr$tip.label == "6731_2_11")
###########################################################################
# tr_dropped = drop.tip(tr, tr$tip.label[c(484:487)]) # These are extreme outliers!
# plot(tr_dropped, cex = 0.4)
tr_dropped = tr # Hack to avoid dropping

# Try to map with serotype
# pheno = readRDS("../Results/treecluster/k_treecluster/mapped_pheno_Wclust.rds")
pheno = readRDS("../Results/fbcluster/k_fbcluster/mapped_pheno_Wclust_MIC.rds")
# pheno = readRDS("../Results/fbcluster/k_fbcluster/mapped_pheno_Wclust_MIC.rds") # For MIC phenotypes - need a different tree file
# pheno = readRDS("../../Genotype/Pensar2019/mapped_pheno_MIC.rds")
# Let's drop to match tree


# dropped sequences


getSero = function(lab, pheno){
  # lab in x_y_z format, convert to x_y#z
  t = unname(unlist(strsplit(lab, "_")))
  t = paste(t[1],"_", t[2], "#", t[3], sep = "")
  return(pheno$serotype[which(pheno$sampleID %in% t)])
}
getCD = function(lab, pheno){
  # lab in x_y_z format, convert to x_y#z
  t = unname(unlist(strsplit(lab, "_")))
  t = paste(t[1],"_", t[2], "#", t[3], sep = "")
  return(pheno$carriage_duration[which(pheno$sampleID %in% t)])
}

getpen.mic = function(lab, pheno){
  # lab in x_y_z format, convert to x_y#z
  t = unname(unlist(strsplit(lab, "_")))
  t = paste(t[1],"_", t[2], "#", t[3], sep = "")
  return(pheno$Penicillin.MIC[which(pheno$sampleID %in% t)])
}
getcef.mic = function(lab, pheno){
  # lab in x_y_z format, convert to x_y#z
  t = unname(unlist(strsplit(lab, "_")))
  t = paste(t[1],"_", t[2], "#", t[3], sep = "")
  return(pheno$Ceftriaxone.MIC[which(pheno$sampleID %in% t)])
}

getAcute = function(lab, pheno){
  # lab in x_y_z format, convert to x_y#z
  t = unname(unlist(strsplit(lab, "_")))
  t = paste(t[1],"_", t[2], "#", t[3], sep = "")
  return(pheno$acute[which(pheno$sampleID %in% t)])
}
getPen = function(lab, pheno){
  # lab in x_y_z format, convert to x_y#z
  t = unname(unlist(strsplit(lab, "_")))
  t = paste(t[1],"_", t[2], "#", t[3], sep = "")
  return(pheno$Penicillin[which(pheno$sampleID %in% t)])
}
getCef = function(lab, pheno){
  # lab in x_y_z format, convert to x_y#z
  t = unname(unlist(strsplit(lab, "_")))
  t = paste(t[1],"_", t[2], "#", t[3], sep = "")
  return(pheno$Ceftriaxone[which(pheno$sampleID %in% t)])
}
getChlora = function(lab, pheno){
  # lab in x_y_z format, convert to x_y#z
  t = unname(unlist(strsplit(lab, "_")))
  t = paste(t[1],"_", t[2], "#", t[3], sep = "")
  return(pheno$Chloramphenicol[which(pheno$sampleID %in% t)])
}
getClinda = function(lab, pheno){
  # lab in x_y_z format, convert to x_y#z
  t = unname(unlist(strsplit(lab, "_")))
  t = paste(t[1],"_", t[2], "#", t[3], sep = "")
  return(pheno$Clindamycin[which(pheno$sampleID %in% t)])
}
getEryth = function(lab, pheno){
  # lab in x_y_z format, convert to x_y#z
  t = unname(unlist(strsplit(lab, "_")))
  t = paste(t[1],"_", t[2], "#", t[3], sep = "")
  return(pheno$Erythromycin[which(pheno$sampleID %in% t)])
}
getStrim = function(lab, pheno){
  # lab in x_y_z format, convert to x_y#z
  t = unname(unlist(strsplit(lab, "_")))
  t = paste(t[1],"_", t[2], "#", t[3], sep = "")
  return(pheno$Sulpha.trimethoprim[which(pheno$sampleID %in% t)])
}
getTetra = function(lab, pheno){
  # lab in x_y_z format, convert to x_y#z
  t = unname(unlist(strsplit(lab, "_")))
  t = paste(t[1],"_", t[2], "#", t[3], sep = "")
  return(pheno$Tetracycline[which(pheno$sampleID %in% t)])
}
getCluster = function(lab, pheno){
  # lab in x_y_z format, convert to x_y#z
  t = unname(unlist(strsplit(lab, "_")))
  t = paste(t[1],"_", t[2], "#", t[3], sep = "")
  return(pheno$cluster[which(pheno$sampleID %in% t)])
}

seros = sapply(tr_dropped$tip.label, function(x) getSero(x, pheno))
cd = sapply(tr_dropped$tip.label, function(x) getCD(x, pheno))
pen.mic = sapply(tr_dropped$tip.label, function(x) getpen.mic(x, pheno))
cef.mic = sapply(tr_dropped$tip.label, function(x) getcef.mic(x, pheno))
acute = sapply(tr_dropped$tip.label, function(x) getAcute(x, pheno))
pen = sapply(tr_dropped$tip.label, function(x) getPen(x, pheno))
cef = sapply(tr_dropped$tip.label, function(x) getCef(x, pheno))
chlora = sapply(tr_dropped$tip.label, function(x) getChlora(x, pheno))
clinda = sapply(tr_dropped$tip.label, function(x) getClinda(x, pheno))
eryth = sapply(tr_dropped$tip.label, function(x) getEryth(x, pheno))
strim = sapply(tr_dropped$tip.label, function(x) getStrim(x, pheno))
tetra = sapply(tr_dropped$tip.label, function(x) getTetra(x, pheno))
cluster = sapply(tr_dropped$tip.label, function(x) getCluster(x, pheno))
# fbcluster = readRDS("../FastBaps/best.partition.rds")
# Map serotypes

maketreeplt_nl = function(tr_dropped, d, title){
  p = ggtree(tr_dropped, layout = "circular")
  if(is.factor(d$ph)){
    p_ = p %<+% d +
      geom_tippoint(aes(col = factor(ph)), cex = 0.7) +
      theme(legend.position = c("right"),
            legend.title=element_blank()) +
      ggtitle(title)
  } else {
    p_ = p %<+% d +
      geom_tippoint(aes(colour = ph), cex = 0.7)+
      scale_color_gradientn(colours = colorspace::qualitative_hcl(10)) +
      theme(legend.position=c("right"),
            legend.title=element_blank()) +
      ggtitle(title)
  }
  return(p_)
}

dir.create("mappedPlots")

d1 = data.frame(seq = tr_dropped$tip.label, ph = as.factor(seros))
p1 = maketreeplt_nl(tr_dropped, d1, "Serotype mapped to phylogeny")
p1
ggsave("mappedPlots/sero.pdf")

d2_1 = data.frame(seq = tr_dropped$tip.label, ph = pen.mic)
p2_1 = maketreeplt_nl(tr_dropped, d2_1, "Penicillin MIC mapped to phylogeny")
p2_1
ggsave("mappedPlots/pen.mic.png")

d2_2 = data.frame(seq = tr_dropped$tip.label, ph = cef.mic)
p2_2 = maketreeplt_nl(tr_dropped, d2_2, "Ceftriaxone MIC mapped to phylogeny")
p2_2
ggsave("mappedPlots/cef.mic.pdf")


d2 = data.frame(seq = tr_dropped$tip.label, ph = cd)
p2 = maketreeplt_nl(tr_dropped, d2, "Carriage Duration mapped to phylogeny")
p2
ggsave("mappedPlots/cd.pdf")

d3 = data.frame(seq = tr_dropped$tip.label, ph = as.factor(acute))
p3 = maketreeplt_nl(tr_dropped, d3, "Clinical Symptoms mapped to phylogeny")
p3
ggsave("mappedPlots/acute.pdf")

d4 = data.frame(seq = tr_dropped$tip.label, ph = as.factor(pen))
p4 = maketreeplt_nl(tr_dropped, d4, "Penicillin AMR mapped to phylogeny")
p4
ggsave("mappedPlots/pen.pdf")

d5 = data.frame(seq = tr_dropped$tip.label, ph = as.factor(cef))
p5 = maketreeplt_nl(tr_dropped, d5, "Ceftriaxone AMR mapped to phylogeny")
p5
ggsave("mappedPlots/cef.pdf")

d6 = data.frame(seq = tr_dropped$tip.label, ph = as.factor(chlora))
p6 = maketreeplt_nl(tr_dropped, d6, "Chloramphenicol AMR mapped to phylogeny")
p6
ggsave("mappedPlots/chlora.pdf")

d7 = data.frame(seq = tr_dropped$tip.label, ph = as.factor(clinda))
p7 = maketreeplt_nl(tr_dropped, d7, "Clindamycin AMR mapped to phylogeny")
p7
ggsave("mappedPlots/clinda.pdf")


d8 = data.frame(seq = tr_dropped$tip.label, ph = as.factor(eryth))
p8 = maketreeplt_nl(tr_dropped, d8, "Erythromycin AMR mapped to phylogeny")
p8
ggsave("mappedPlots/eryth.pdf")

d9 = data.frame(seq = tr_dropped$tip.label, ph = as.factor(strim))
p9 = maketreeplt_nl(tr_dropped, d9, "Sulpha trimethoprim AMR mapped to phylogeny")
p9
ggsave("mappedPlots/strim.pdf")

d10 = data.frame(seq = tr_dropped$tip.label, ph = as.factor(tetra))
p10 = maketreeplt_nl(tr_dropped, d10, "Tetracycline AMR mapped to phylogeny")
p10
ggsave("mappedPlots/tetra.pdf")

d11 = data.frame(seq = tr_dropped$tip.label, ph = as.factor(cluster))
p11 = maketreeplt_nl(tr_dropped, d11, "Cluster AMR mapped to phylogeny")
p11
ggsave("mappedPlots/cluster.pdf")
# p1 = ggtree(tr_dropped, layout = "circular")
# p1 %<+% d1 +
#   geom_tippoint(aes(col = factor(serotype))) +
#   theme(legend.title = element_blank(), # no title
#         legend.key = element_blank()) # no keys
# 
# # Map carriage duration
# p2 = ggtree(tr_dropped, layout = "circular")
# 
# # Map acute
# p3 = ggtree(tr_dropped, layout = "circular")
# d3 = data.frame(seq = tr_dropped$tip.label, 
#                 Acute = acute)
# p3 %<+% d3 + geom_tippoint(aes(colour = factor(Acute)))
# 
# 
# 
# 
# #############################3
# # chr = chronos(tr_dropped)
# # saveRDS(chr, "maelatree.chr")
# # chr = readRDS("maelatree.chr") # this order is compltely different from pheno order, must match!
# # chr_rooted = root(chr, 1)
# 
# 
# sampleID = pheno$sampleID
# k = 30
# clust = cutree(chr_rooted, k)
# ord2 = sapply(names(clust), function(x) ridx(x, sampleID, "_"))
# 
# cluster = rep(NA, length(sampleID))
# cluster[ord2] = clust
# 
# pheno$cluster = cluster # Cuts cluster order, not data order
# 
# pheno = pheno[-which(is.na(pheno$cluster)), ]
# 
# 
# nclsts = length(table(pheno$cluster)) # can be different from k at times!
# if(nclsts <= 20) glabstat = TRUE else glabstat = FALSE
# if(nclsts <= 50) {
#   lpos = "right"
#   bvis = TRUE
# } else {
#   lpos = "none"
#   bvis = FALSE
# }
# plot_ = ggplot(pheno, aes(x=factor(serotype), y=carriage_duration,  col=factor(cluster))) + geom_point() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = lpos)
# print(plot_)
# plot_ = ggplot(pheno, aes(x=factor(cluster), y=carriage_duration,  col=factor(serotype))) + geom_point() + theme(axis.text.x = element_text(vjust = 0.5, hjust=1))
# print(plot_)
# cols = rainbow_hcl(n=nclsts,c=90,l=50)
# hist(pheno$cluster, breaks = seq(1,length(table(pheno$cluster))), col = cols, border = bvis)
# 
# branches_color(as.dendrogram(chr_rooted), k = nclsts, col = cols, groupLabels = glabstat) %>%
#   set("labels_cex", c(.1)) %>% plot
# colored_bars(colors = cols[cuts], dend = as.dendrogram(chr_rooted), sort_by_labels_order = FALSE, rowLabels = NULL)
# # DENDOBJ approach
# 
# 
# dendObj = ReadDendrogram("maelatree.treefile")
# cuts = cutree(dend, 100, order_clusters_as_data = F)
# 
# 
# sampleID = pheno$sampleID
# ord = sapply(names(cuts), function(x) ridx(x, sampleID))
# cluster = rep(NA, length(sampleID))
# cluster[ord] = cuts
# pheno = mutate(pheno, cluster = cluster) # Cuts cluster order, not data order
# 
# 
# tr_um = force.ultrametric(tr_dropped)
# 
# 
# ggtree(tr_dropped, layout = "circular")
# 
# 
# 
# 
# 
# tr_test = groupClade(tr_um, c(5000))
# ggtree(tr_test, aes(color = group)) + geom_tiplab()
# 
