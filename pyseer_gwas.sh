#!/bin/bash
# conda must be active to run this if installed via conda

#read pheno
pheno="cd"
SNPS="biallelicsnps_cdacute.vcf" # Filtered SNP file
PHYLODIST="phylodist_cdacute.tsv"
PHYLOSIM="phylosim_cdacute.tsv"
COVARIATES="cov.cdacute"
LINEAGECLUSTERS="fastbaps_cd.acute.txt"
printf "Now processing %s\n" ${pheno}

echo "Running FEM with phylodistances\n"
pyseer --phenotypes pheno_${pheno} --vcf ${SNPS} --distances ${PHYLODIST} --covariates ${COVARIATES} --use-covariates 3 --cpu 4 --continuous > res_${pheno}_fem_pd.txt
cat <(echo "#CHR SNP BP minLOG10(P) log10(p) r^2") <(paste <(sed '1d' res_${pheno}_fem_pd.txt | cut -d "_" -f 2) <(sed '1d' res_${pheno}_fem_pd.txt | cut -f 4) |awk '{p = -log($2)/log(10); print "26",".",$1,p,p,"0"}' ) |tr ' ' '\t' > plt_${pheno}_fem_pd.plot
echo "###########################################################################################################################################################"

echo "Running FEM with fastbaps clusters\n"
pyseer --phenotypes pheno_${pheno} --vcf ${SNPS} --no-distances --covariates ${COVARIATES} --use-covariates 2 3 --cpu 4 > res_${pheno}_fem_fb.txt
cat <(echo "#CHR SNP BP minLOG10(P) log10(p) r^2") <(paste <(sed '1d' res_${pheno}_fem_fb.txt | cut -d "_" -f 2) <(sed '1d' res_${pheno}_fem_fb.txt | cut -f 4) |awk '{p = -log($2)/log(10); print "26",".",$1,p,p,"0"}' ) |tr ' ' '\t' > plt_${pheno}_fem_fb.plot
echo "###########################################################################################################################################################"

echo "Running LMM with phylosims"
pyseer --phenotypes pheno_${pheno} --vcf ${SNPS} --similarity ${PHYLOSIM} --covariates ${COVARIATES} --lmm --save-lmm _LMMCACHE_CD_PSIM --use-covariates 3 --cpu 4 --continuous > res_${pheno}_lmm_ps.txt
cat <(echo "#CHR SNP BP minLOG10(P) log10(p) r^2") <(paste <(sed '1d' res_${pheno}_lmm_ps.txt | cut -d "_" -f 2) <(sed '1d' res_${pheno}_lmm_ps.txt | cut -f 4) |awk '{p = -log($2)/log(10); print "26",".",$1,p,p,"0"}' ) |tr ' ' '\t' > plt_${pheno}_lmm_ps.plot
echo "###########################################################################################################################################################"

pheno="pen.mic"
SNPS="biallelicsnps_mic.vcf"
PHYLODIST="phylodist_mic.tsv"
PHYLOSIM="phylosim_mic.tsv"
COVARIATES="cov.mic"
LINEAGECLUSTERS="fastbaps_mic.txt"

printf "Now processing %s\n" ${pheno}
echo "Running FEM with phylodistances"
pyseer --phenotypes pheno_${pheno} --vcf ${SNPS} --distances ${PHYLODIST} --covariates ${COVARIATES} --use-covariates 3 4 --cpu 4 --continuous > res_${pheno}_fem_pd.txt
cat <(echo "#CHR SNP BP minLOG10(P) log10(p) r^2") <(paste <(sed '1d' res_${pheno}_fem_pd.txt | cut -d "_" -f 2) <(sed '1d' res_${pheno}_fem_pd.txt | cut -f 3) |awk '{p = -log($2)/log(10); print "26",".",$1,p,p,"0"}' ) |tr ' ' '\t' > plt_${pheno}_fem_pd.plot
echo "###########################################################################################################################################################"

echo "Running FEM with fastbaps clusters"
pyseer --phenotypes pheno_${pheno} --vcf ${SNPS} --no-distances --covariates ${COVARIATES} --use-covariates 2 3 4 --cpu 4 > res_${pheno}_fem_fb.txt
cat <(echo "#CHR SNP BP minLOG10(P) log10(p) r^2") <(paste <(sed '1d' res_${pheno}_fem_fb.txt | cut -d "_" -f 2) <(sed '1d' res_${pheno}_fem_fb.txt | cut -f 4) |awk '{p = -log($2)/log(10); print "26",".",$1,p,p,"0"}' ) |tr ' ' '\t' > plt_${pheno}_fem_fb.plot
echo "###########################################################################################################################################################"

echo "Running LMM with phylosims"
pyseer --phenotypes pheno_${pheno} --vcf ${SNPS} --similarity ${PHYLOSIM} --lmm --save-lmm _LMMCACHE_PENMIC_PSIM --covariates ${COVARIATES} --use-covariates 3 4 --cpu 4 --continuous > res_${pheno}_lmm_ps.txt
cat <(echo "#CHR SNP BP minLOG10(P) log10(p) r^2") <(paste <(sed '1d' res_${pheno}_lmm_ps.txt | cut -d "_" -f 2) <(sed '1d' res_${pheno}_lmm_ps.txt | cut -f 4) |awk '{p = -log($2)/log(10); print "26",".",$1,p,p,"0"}' ) |tr ' ' '\t' > plt_${pheno}_lmm_ps.plot

pheno="cef.mic"
printf "Now processing %s\n" ${pheno}
echo "Running FEM with phylodistances"
pyseer --phenotypes pheno_${pheno} --vcf ${SNPS} --distances ${PHYLODIST} --covariates ${COVARIATES} --use-covariates 3 4 --cpu 4 --continuous > res_${pheno}_fem_pd.txt
cat <(echo "#CHR SNP BP minLOG10(P) log10(p) r^2") <(paste <(sed '1d' res_${pheno}_fem_pd.txt | cut -d "_" -f 2) <(sed '1d' res_${pheno}_fem_pd.txt | cut -f 4) |awk '{p = -log($2)/log(10); print "26",".",$1,p,p,"0"}' ) |tr ' ' '\t' > plt_${pheno}_fem_pd.plot
echo "###########################################################################################################################################################"

echo "Running FEM with fastbaps clusters"
pyseer --phenotypes pheno_${pheno} --vcf ${SNPS} --no-distances --covariates ${COVARIATES} --use-covariates 2 3 4 --cpu 4 > res_${pheno}_fem_fb.txt
cat <(echo "#CHR SNP BP minLOG10(P) log10(p) r^2") <(paste <(sed '1d' res_${pheno}_fem_fb.txt | cut -d "_" -f 2) <(sed '1d' res_${pheno}_fem_fb.txt | cut -f 4) |awk '{p = -log($2)/log(10); print "26",".",$1,p,p,"0"}' ) |tr ' ' '\t' > plt_${pheno}_fem_fb.plot
echo "###########################################################################################################################################################"

echo "Running LMM with phylosims"
pyseer --phenotypes pheno_${pheno} --vcf ${SNPS} --similarity ${PHYLOSIM} --lmm --save-lmm _LMMCACHE_CEFMIC_PSIM --covariates ${COVARIATES} --use-covariates 3 4 --cpu 4 --continuous > res_${pheno}_lmm_ps.txt
cat <(echo "#CHR SNP BP minLOG10(P) log10(p) r^2") <(paste <(sed '1d' res_${pheno}_lmm_ps.txt | cut -d "_" -f 2) <(sed '1d' res_${pheno}_lmm_ps.txt | cut -f 4) |awk '{p = -log($2)/log(10); print "26",".",$1,p,p,"0"}' ) |tr ' ' '\t' > plt_${pheno}_lmm_ps.plot
echo "###########################################################################################################################################################"
