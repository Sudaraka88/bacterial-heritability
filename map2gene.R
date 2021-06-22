library(data.table)
# Checked 20210621

# reference_genome.rds is a reduced version of the SPN23F ref. genome published here: https://journals.asm.org/doi/full/10.1128/JB.01343-08
# The gff file is available at: https://dx.doi.org/10.6084/m9.figshare.7588832 (via pyseer: https://pyseer.readthedocs.io/en/master/tutorial.html)

# This code quickly matches a snp position to a nearby gene in the reference genome
gff_data = readRDS("reference_genome.rds") # provided
gene_data = gff_data[!is.na(gff_data$gene),] # Only mapped genes in this df
dst = data.table(w = gene_data$start, val = gene_data$start)
setattr(dst, "sorted", "w")
dend = data.table(w = gene_data$end, val = gene_data$end)
setattr(dend, "sorted", "w")


map2gene = function(x, gene_data){
  TMP = which(gene_data$start <= x & gene_data$end >= x)
  if(length(TMP) > 0) {
    gd = gene_data[TMP,]
    overlap = 0 # Overlaps with a gene
  } else {
    sid = dst[J(x), .I, roll = "nearest", by = .EACHI]$I
    eid = dend[J(x), .I, roll = "nearest", by = .EACHI]$I
    ids = c(sid, eid)
    ds = c(abs(gene_data$start[sid] - x), abs(gene_data$end[eid] - x))
    minds = min(ds)
    gd = gene_data[ids[which(ds == minds)],]
    overlap = minds # No overlap with gene
  }
  return(list(x = x, gd = gd, dist_from_gene = overlap))
}