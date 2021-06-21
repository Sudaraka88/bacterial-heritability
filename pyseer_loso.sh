#!/bin/bash
# conda must be active to run this if installed via conda

pheno="cd"
SNPS="biallelicsnps_cdacute.vcf" # Filtered SNP file
PHYLODIST="phylodist_cdacute.tsv"
PHYLOSIM="phylosim_cdacute.tsv"
COVARIATES="cov.cdacute"
#LINEAGECLUSTERS="fastbaps_cd.acute.txt"
printf "Now processing %s\n" ${pheno}

#echo "Running full-wg-enet with phylodistances"
#pyseer --phenotypes pheno_${pheno} --vcf ${SNPS} --distances ${PHYLODIST} --wg enet --cpu 4 --continuous --save-model ${pheno}_wgenet_pd --save-vars ${pheno}_wgenetsnps > res_${pheno}_wgenet_pd.txt
#echo "###########################################################################################################################################################"


for var in 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28
do
    pyseer --phenotypes LOSO/${pheno}_${var}_trn --vcf ${SNPS} --distances ${PHYLODIST} --wg enet --cpu 4 --continuous --save-model LOSO/${pheno}_${var} > LOSO/${pheno}_${var}.res
    cut -f 1 LOSO/${pheno}_${var}_tst | sed '1d' > LOSO/${pheno}_${var}_tst.samples
    enet_predict_pyseer --vcf ${SNPS} --true-values LOSO/${pheno}_${var}_tst LOSO/${pheno}_${var}.pkl LOSO/${pheno}_${var}_tst.samples > LOSO/${pheno}_${var}.preds
    echo "###########################################################################################################################################################"
done

pheno="pen.mic"
SNPS="biallelicsnps_mic.vcf"
PHYLODIST="phylodist_mic.tsv"
PHYLOSIM="phylosim_mic.tsv"
COVARIATES="cov.mic"
#LINEAGECLUSTERS="fastbaps_mic.txt"

printf "Now processing %s\n" ${pheno}

#echo "Running wg-enet with phylodistances"
#pyseer --phenotypes pheno_${pheno} --vcf ${SNPS} --distances ${PHYLODIST} --wg enet --cpu 4 --continuous --save-model ${pheno}_wgenet_pd --save-vars ${pheno}_wgenetsnps > res_${pheno}_wgenet_pd.txt
#echo "###########################################################################################################################################################"

for var in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
do
    pyseer --phenotypes LOSO/${pheno}_${var}_trn --vcf ${SNPS} --distances ${PHYLODIST} --wg enet --cpu 4 --continuous --save-model LOSO/${pheno}_${var} > LOSO/${pheno}_${var}.res
    cut -f 1 LOSO/${pheno}_${var}_tst | sed '1d' > LOSO/${pheno}_${var}_tst.samples
    enet_predict_pyseer --vcf ${SNPS} --true-values LOSO/${pheno}_${var}_tst LOSO/${pheno}_${var}.pkl LOSO/${pheno}_${var}_tst.samples > LOSO/${pheno}_${var}.preds
    echo "###########################################################################################################################################################"
done

pheno="cef.mic"
printf "Now processing %s\n" ${pheno}
#echo "Running wg-enet with phylodistances"
#pyseer --phenotypes pheno_${pheno} --vcf ${SNPS} --distances ${PHYLODIST} --wg enet --cpu 4 --continuous --save-model ${pheno}_wgenet_pd --save-vars ${pheno}_wgenetsnps > res_${pheno}_wgenet_pd.txt
#echo "###########################################################################################################################################################"

for var in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
do
    pyseer --phenotypes LOSO/${pheno}_${var}_trn --vcf ${SNPS} --distances ${PHYLODIST} --wg enet --cpu 4 --continuous --save-model LOSO/${pheno}_${var} > LOSO/${pheno}_${var}.res
    cut -f 1 LOSO/${pheno}_${var}_tst | sed '1d' > LOSO/${pheno}_${var}_tst.samples
    enet_predict_pyseer --vcf ${SNPS} --true-values LOSO/${pheno}_${var}_tst LOSO/${pheno}_${var}.pkl LOSO/${pheno}_${var}_tst.samples > LOSO/${pheno}_${var}.preds
    echo "###########################################################################################################################################################"
done
