if(rstudioapi::isAvailable()) setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # WORKING DIRECTORY
# Checked 20210621

ph = "cd"

print(paste("Begin at:", Sys.time(), "Loading libraries..."))
library(treeWAS)
library(ape)
library(Rcpp)
library(foreach)
library(doParallel)

registerDoParallel(cores = 8L) # query of n_cores buggy? - change n_cores according to system

# Function to fix maela_tree$tip.label
refmt = function(nm){
  temp = unlist(strsplit(nm, split = "_"))
  return(paste(temp[1], "_", temp[2], "#", temp[3], sep = ""))
}

print("Loading files...")

############ Deal with the genotype first ###########
# snp_pos = readRDS("pos.rds") # not provided 
# genotype = readRDS("mapped_mx_CC.rds") # not provided
# snps.unique = get.unique.matrix(genotype, MARGIN = 2)$unique.data # These are now unique SNPs
# MA = apply(genotype, 2, function(x) names(sort(table(x), decreasing = T)[1]))
# genotype_bc = matrix(nrow = nrow(genotype), ncol = ncol(genotype), data = 0)
# colnames(genotype_bc) = colnames(genotype)
# rownames(genotype_bc) = rownames(genotype)
# for(i in 1:ncol(genotype)) genotype_bc[which(genotype[,i] == MA[i]),i] = 1
# snps.unique = get.unique.matrix(genotype_bc, MARGIN = 2)$unique.data # These are now unique SNPs
# saveRDS(snps.unique, "snps.unique_cd.rds") # unique snps

snps.unique = readRDS("snps.unique_cd.rds") 
print("Loaded SNPs: snps.unique_cd.rd")
nSEQs = dim(snps.unique)[1]
nSNPs = dim(snps.unique)[2]

pheno = readRDS("mapped_pheno_Wclust.rds")
if(ph == "cd") phenotype = log10(pheno$carriage_duration) # Load CD pheno

names(phenotype) = pheno$sampleID
print("Loaded Phenotype")
str(phenotype)


maela_tree = read.tree("Trees/cfrm_cd.labelled_tree.newick")
print("Loaded tree")

maela_tree$tip.label = unname(sapply(maela_tree$tip.label, function(x) refmt(x)))

# Check if tree lengths < 0  # None for maela_tree
toChange <- which(maela_tree$edge.length < 0) 
if(length(toChange) > 0){
  maela_tree$edge.length[toChange] <- 0
  cat("Setting", length(toChange), "negative branch lengths to zero.\n", sep=" ")
}

str(maela_tree)

#####################################################################
##################### RECONSTRUCTION ################################
st = nSEQs + 1
en = length(maela_tree$edge.length) + 1
if(!file.exists("snps.rec_stpd.rds")){
  t = Sys.time()
  snps.rec_stpd = foreach(i = 1:nSNPs, .combine = "cbind", .verbose = F) %dopar% {
    asr(var = snps.unique[,i], tree = maela_tree)[st:en]
  }
  print(paste("Done, ET:", Sys.time() - t))
  # saveRDS(snps.rec_stpd, "snps.rec_stpd.rds")
  colnames(snps.rec_stpd) = colnames(snps.unique)
  saveRDS(snps.rec_stpd, "snps.rec_stpd.rds")
  print("snps recon complete...")
} else {
  snps.rec_stpd = readRDS("snps.rec_stpd.rds")
  print("Loaded SNP RECON: snps.rec_stpd.rds")
}

####################### SIMULATION ########################
if(!file.exists("simsnps.unique.rds")){
  print("Begin snp sim")
  t = Sys.time()
  n.snps.sim <- ncol(snps.unique)*10
  print(paste("Simulating", n.snps.sim, "SNPs"))
  n.subs = get.fitch.n.mts(snps.unique,maela_tree)
  n.subs = table(n.subs)
  noms <- as.numeric(names(n.subs))
  temp <- rep(0, max(noms))
  for(i in 1:max(noms)){
    if(i %in% noms) temp[i] <- n.subs[which(noms==i)]
  }
  n.subs <- temp
  names(n.subs) <- 1:length(n.subs) # homoplasy distribution
  dist.dna.model = "JC69"
  seed = 1
  simsnps = snp.sim(n.snps = n.snps.sim,
                    n.subs = n.subs,
                    n.snps.assoc = 0,
                    assoc.prob = 100,
                    tree = maela_tree,
                    phen.loci = NULL,
                    heatmap = FALSE,
                    reconstruct = FALSE,
                    dist.dna.model = dist.dna.model,
                    row.names = rownames(snps.unique),
                    seed = seed)
  
  print(paste("Done, ET:", Sys.time() - t))
  print("Saving all simulated SNPs: simsnps.rds")
  saveRDS(simsnps, "simsnps.rds")
  # There could be non-unique SNPs here, remove them
  # simsnps = readRDS("simsnps.rds")
  print(str(simsnps$snps))
  simsnps.unique = get.unique.matrix(simsnps$snps, MARGIN = 2)$unique.data
  print(paste("After duplicate removal,", ncol(simsnps.unique), "SNPs were retained"))
  saveRDS(simsnps.unique, "simsnps.unique.rds")
  print(str(simsnps.unique))
} else {
  simsnps.unique = readRDS("simsnps.unique.rds")
  print("Loaded SIM SNPS: simsnps.unique.rds")
}

##################### RECONSTRUCTION ################################
# This is very time and memory consuming, avoid repeating at all costs!
if(!file.exists("simsnps.rec_stpd.rds")){
  print("Begin recon for sim snps")
  nSIMSNPs = ncol(simsnps.unique) # This is the real number of simulated SNPs after removing duplicates
  print(paste("Identified", ncol(nSIMSNPs), "simulated SNPs"))
  registerDoParallel(cores = 8L)
  t = Sys.time()
  simsnps.rec_stpd = foreach(i = 1:nSIMSNPs, .combine = "cbind", .verbose = T) %dopar% {
    asr(var = simsnps.unique[,i], tree = maela_tree)[st:en]
  }
  print(paste("Done, ET:", Sys.time() - t))
  # saveRDS(simsnps.rec_stpd, "simsnps.rec_strpd.rds")
  colnames(simsnps.rec_stpd) = colnames(simsnps.unique)
  saveRDS(simsnps.rec_stpd, "simsnps.rec_stpd.rds")
  print("simulated SNPs recon complete...")
} else {
  simsnps.rec_stpd = readRDS("simsnps.rec_stpd.rds")
  print("Loaded RECON SIM SNPS: simsnps.rec_stpd.rds")
}

#####################################################################
#####################################################################
#####################################################################

print("Starting sanity checks and reordering...")

# Set SNPs and pheno in the same order as maela_tree
if(!identical(as.character(rownames(snps.unique)), as.character(maela_tree$tip.label))){
  ord <- match(maela_tree$tip.label, rownames(snps.unique))
  snps.unique <- snps.unique[ord,]
  ## check:
  if(!identical(as.character(rownames(snps.unique)), as.character(maela_tree$tip.label))){
    stop("Unable to rearrange snps such that rownames(snps)
           match content and order of tree$tip.label.
           Please check that these match.")
  }
}

if(!identical(as.character(names(phenotype)), as.character(maela_tree$tip.label))){
  ord <- match(maela_tree$tip.label, names(phenotype))
  phenotype <- phenotype[ord]
  ## check:
  if(!identical(as.character(names(phenotype)), as.character(maela_tree$tip.label))){
    stop("Unable to rearrange phen such that names(phen)
           match content and order of tree$tip.label.
           Please check that these match.")
  }
}

#############################################################
########## RECON PART - LOAD AFTER RUNNING ##################
#############################################################
# Conduct recon only if file doesn't exist
if(!file.exists(paste("phen.rec_stpd_", ph, ".rds", sep = ""))){
  print("Begin phenotype recon")
  # Let's make the ancestral reconstructions
  if(ph == "cd") phen.rec = asr(var = phenotype, tree = maela_tree, type = "ML", method = "continuous")
  phen.rec_strpd = phen.rec[st:en]
  saveRDS(phen.rec_strpd, paste("phen.rec_stpd_", ph, ".rds", sep = ""))
  print("Complete...")
} else {
  phen.rec_strpd = readRDS(paste("phen.rec_stpd_", ph, ".rds", sep = ""))
  print(paste("Loaded phen.rec_stpd_", ph, ".rds", sep = ""))
}


##############################################################
#################### Testing Phase ###########################
##############################################################
print("Begin testing, previous test results will be overriden...")
test = c("terminal", "simultaneous", "subsequent") 
# Real SNPs
print("snps.unique")
print(str(snps.unique))

print("snps.rec_stpd")
print(str(snps.rec_stpd))

# Simulated SNPs
print("simsnps.unique")
print(str(simsnps.unique))

print("simsnps.rec_stpd")
print(str(simsnps.rec_stpd))

# phenotype
print("phen.rec_strpd")
print(str(phen.rec_strpd))

sig.list = list()
for(i in 1:length(test)){
  assoc.scores <- get.assoc.scores(snps = snps.unique, 
                                   snps.sim = simsnps.unique, 
                                   phen = phenotype,
                                   tree = maela_tree, 
                                   test = test[i],
                                   snps.reconstruction = rbind(snps.unique, snps.rec_stpd), 
                                   snps.sim.reconstruction = rbind(simsnps.unique, simsnps.rec_stpd), 
                                   phen.reconstruction = c(phenotype, phen.rec_strpd),
                                   unique.cols = TRUE)
  print(str(assoc.scores))
  saveRDS(assoc.scores, paste(ph, "_t_", test[i], "_sc.rds", sep = ""))
  
  sig.list[[i]] = get.sig.snps(corr.dat = assoc.scores$corr.dat,
                               corr.sim = assoc.scores$corr.sim,
                               snps.names = colnames(snps.unique),
                               test = test[i],
                               n.tests = length(test),
                               p.value = 0.01,
                               p.value.correct = "bonf",
                               p.value.by = "count")
  saveRDS(sig.list[[i]], paste(ph, "_t_", test[i], "_sigsnps.rds", sep = ""))
  ######################################################################
  ############################# Plot Results ###########################
  ######################################################################
  
  
  if(test[i] == "fisher"){
    ylab <- NULL # --> (uncorrected or) -log10 p-value
    log10 <- TRUE
    ttl <- paste("\n \n(", test[i], "test)")
  }else{
    ylab <- paste(test[i], "score", sep=" ")
    log10 <- FALSE
    ttl <- paste("\n \n(", test[i], "score)")
  }
  ## plot:
  png(file = paste("p_", ph, "_", test[i], "mhp.png", sep = ""), width = 1920, height = 1080, type = "cairo")
  manhattan.plot(p.vals = abs(sig.list[[i]]$corr.dat),
                 col = "funky",
                 transp = 0.25,
                 sig.thresh =  sig.list[[i]]$sig.thresh,
                 thresh.col="red",
                 snps.assoc = NULL,
                 snps.assoc.col = "red",
                 jitter.amount = 0.00001,
                 min.p = NULL,
                 log10=log10,
                 ylab=ylab)
  ## Add subtitle:
  title(ttl, cex.main=0.9)
  dev.off()
  
  ####################################
  ## 4) (B,C) Plot the distribution ##
  ####################################
  png(file = paste("p_", ph, "_", test[i], "null_dist.png", sep = ""), width = 1920, height = 1080, type = "cairo")
  if(test[i] == "fisher"){
    ## Generate one histogram per test:
    plot_sig_snps(corr.dat = -log10(abs(sig.list[[i]]$corr.dat)),
                  corr.sim = -log10(abs(sig.list[[i]]$corr.sim)),
                  corr.sim.subset = NULL,
                  sig.corrs = -log10(abs(sig.list[[i]]$corr.dat[sig.list[[i]]$sig.snps])),
                  sig.snps = sig.list[[i]]$sig.snps.names,
                  sig.thresh = -log10(sig.list[[i]]$sig.thresh),
                  test = test[i],
                  sig.snps.col = "black",
                  hist.col = rgb(0,0,1,0.5),
                  hist.subset.col = rgb(1,0,0,0.5),
                  thresh.col = "red",
                  snps.assoc = NULL,
                  snps.assoc.col = "blue",
                  bg = "lightgrey",
                  grid = TRUE,
                  freq = FALSE,
                  plot.null.dist = TRUE,
                  plot.dist = FALSE)
  }else{
    ## Generate one histogram per test:
    plot_sig_snps(corr.dat = abs(sig.list[[i]]$corr.dat),
                  corr.sim = abs(sig.list[[i]]$corr.sim),
                  corr.sim.subset = NULL,
                  sig.corrs = abs(sig.list[[i]]$corr.dat[sig.list[[i]]$sig.snps]),
                  sig.snps = sig.list[[i]]$sig.snps.names,
                  sig.thresh = sig.list[[i]]$sig.thresh,
                  test = test[i],
                  sig.snps.col = "black",
                  hist.col = rgb(0,0,1,0.5),
                  hist.subset.col = rgb(1,0,0,0.5),
                  thresh.col = "red",
                  snps.assoc = NULL,
                  snps.assoc.col = "blue",
                  bg = "lightgrey",
                  grid = TRUE,
                  freq = FALSE,
                  plot.null.dist = TRUE,
                  plot.dist = FALSE)
  }
  
  dev.off()
  png(file = paste("p_", ph, "_", test[i], "emp_dist.png", sep = ""), width = 1920, height = 1080, type = "cairo")
  if(test[i] == "fisher"){
    ## Generate one histogram per test:
    plot_sig_snps(corr.dat = -log10(abs(sig.list[[i]]$corr.dat)),
                  corr.sim = -log10(abs(sig.list[[i]]$corr.sim)),
                  corr.sim.subset = NULL,
                  sig.corrs = -log10(abs(sig.list[[i]]$corr.dat[sig.list[[i]]$sig.snps])),
                  sig.snps = sig.list[[i]]$sig.snps.names,
                  sig.thresh = -log10(sig.list[[i]]$sig.thresh),
                  test = test[i],
                  sig.snps.col = "black",
                  hist.col = rgb(0,0,1,0.5),
                  hist.subset.col = rgb(1,0,0,0.5),
                  thresh.col = "red",
                  snps.assoc = NULL,
                  snps.assoc.col = "blue",
                  bg = "lightgrey",
                  grid = TRUE,
                  freq = FALSE,
                  plot.null.dist = FALSE,
                  plot.dist = TRUE)
  }else{
    ## Generate one histogram per test:
    plot_sig_snps(corr.dat = abs(sig.list[[i]]$corr.dat),
                  corr.sim = abs(sig.list[[i]]$corr.sim),
                  corr.sim.subset = NULL,
                  sig.corrs = abs(sig.list[[i]]$corr.dat[sig.list[[i]]$sig.snps]),
                  sig.snps = sig.list[[i]]$sig.snps.names,
                  sig.thresh = sig.list[[i]]$sig.thresh,
                  test = test[i],
                  sig.snps.col = "black",
                  hist.col = rgb(0,0,1,0.5),
                  hist.subset.col = rgb(1,0,0,0.5),
                  thresh.col = "red",
                  snps.assoc = NULL,
                  snps.assoc.col = "blue",
                  bg = "lightgrey",
                  grid = TRUE,
                  freq = FALSE,
                  plot.null.dist = FALSE,
                  plot.dist = TRUE)
  }
  dev.off() 
  
  
}
names(sig.list) = test
saveRDS(sig.list, paste("sig.list_", ph, ".rds", sep = ""))

print(paste("Complete at", Sys.time()))



