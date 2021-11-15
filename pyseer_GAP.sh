#!/bin/bash
# conda must be active to run this

pheno="cd"
SNPS="gap_cd.vcf"
PHYLOSIM="phylosim_cd.tsv"
COVARIATES="cov_cd"


echo "Running LMM with phylosims"
pyseer --phenotypes pheno_${pheno} --vcf ${SNPS} --similarity ${PHYLOSIM} --covariates ${COVARIATES} --use-covariates 2 --lmm --save-lmm _LMMCACHE_PENMIC_PSIM --cpu 4 --continuous > res_${pheno}_lmm_ps.txt
cat <(echo "#CHR SNP BP minLOG10(P) log10(p) r^2") <(paste <(sed '1d' res_${pheno}_lmm_ps.txt | cut -d "_" -f 2) <(sed '1d' res_${pheno}_lmm_ps.txt | cut -f 4) |awk '{p = -log($2)/log(10); print "26",".",$1,p,p,"0"}' ) |tr ' ' '\t' > plt_${pheno}_lmm_ps.plot


pheno="cef.mic"
SNPS="gap_mic.vcf"
PHYLOSIM="phylosim_mic.tsv"

echo "Running LMM with phylosims"
pyseer --phenotypes pheno_${pheno} --vcf ${SNPS} --similarity ${PHYLOSIM} --lmm --save-lmm _LMMCACHE_PENMIC_PSIM --cpu 4 --continuous > res_${pheno}_lmm_ps.txt
cat <(echo "#CHR SNP BP minLOG10(P) log10(p) r^2") <(paste <(sed '1d' res_${pheno}_lmm_ps.txt | cut -d "_" -f 2) <(sed '1d' res_${pheno}_lmm_ps.txt | cut -f 4) |awk '{p = -log($2)/log(10); print "26",".",$1,p,p,"0"}' ) |tr ' ' '\t' > plt_${pheno}_lmm_ps.plot


pheno="cef.mic"
printf "Now processing %s\n" ${pheno}

echo "Running LMM with phylosims"
pyseer --phenotypes pheno_${pheno} --vcf ${SNPS} --similarity ${PHYLOSIM} --lmm --save-lmm _LMMCACHE_CEFMIC_PSIM --cpu 4 --continuous > res_${pheno}_lmm_ps.txt
cat <(echo "#CHR SNP BP minLOG10(P) log10(p) r^2") <(paste <(sed '1d' res_${pheno}_lmm_ps.txt | cut -d "_" -f 2) <(sed '1d' res_${pheno}_lmm_ps.txt | cut -f 4) |awk '{p = -log($2)/log(10); print "26",".",$1,p,p,"0"}' ) |tr ' ' '\t' > plt_${pheno}_lmm_ps.plot
echo "###########################################################################################################################################################"
