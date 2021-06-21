#!/bin/bash
# conda must be active to run this if installed via conda

pheno="cd"
SNPS="biallelicsnps_cdacute.vcf" # Filtered SNP file
PHYLODIST="phylodist_cdacute.tsv"
PHYLOSIM="phylosim_cdacute.tsv"
COVARIATES="cov.cdacute"
#LINEAGECLUSTERS="fastbaps_cd.acute.txt"
printf "Now processing %s\n" ${pheno}

# Prediction only works with --wg enet for now, let's stick to that

echo "Running wg-enet with phylodistances"
pyseer --phenotypes pheno_${pheno} --vcf ${SNPS} --distances ${PHYLODIST} --wg enet --cpu 4 --continuous --save-model ${pheno}_wgenet_pd > res_${pheno}_wgenet_pd.txt
echo "###########################################################################################################################################################"

#echo "Running wg-enet with fastbaps clusters"
#pyseer --phenotypes pheno_${pheno} --vcf ${SNPS} --no-distances --lineage-clusters ${LINEAGECLUSTERS} --sequence-reweighting --wg enet --cpu 4 --continuous --save-model ${pheno}_wgenet_fb > res_${pheno}__wgenet_fb.txt
#echo "###########################################################################################################################################################"

pheno="pen.mic"
SNPS="biallelicsnps_mic.vcf"
PHYLODIST="phylodist_mic.tsv"
PHYLOSIM="phylosim_mic.tsv"
COVARIATES="cov.mic"
#LINEAGECLUSTERS="fastbaps_mic.txt"

printf "Now processing %s\n" ${pheno}

echo "Running wg-enet with phylodistances"
pyseer --phenotypes pheno_${pheno} --vcf ${SNPS} --distances ${PHYLODIST} --wg enet --cpu 4 --continuous --save-model ${pheno}_wgenet_pd > res_${pheno}_wgenet_pd.txt
echo "###########################################################################################################################################################"

#echo "Running wg-enet with fastbaps clusters"
#pyseer --phenotypes pheno_${pheno} --vcf ${SNPS} --no-distances --lineage-clusters ${LINEAGECLUSTERS} --sequence-reweighting --wg enet --cpu 4 --continuous --alpha 1 --save-model ${pheno}_wgenet_fb > res_${pheno}_wgenet_fb.txt
#echo "###########################################################################################################################################################"

pheno="cef.mic"
printf "Now processing %s\n" ${pheno}
echo "Running wg-enet with phylodistances"
pyseer --phenotypes pheno_${pheno} --vcf ${SNPS} --distances ${PHYLODIST} --wg enet --cpu 4 --continuous --save-model ${pheno}_wgenet_pd > res_${pheno}_wgenet_pd.txt
echo "###########################################################################################################################################################"

#echo "Running wg-enet with fastbaps clusters"
#pyseer --phenotypes pheno_${pheno} --vcf ${SNPS} --no-distances --lineage-clusters ${LINEAGECLUSTERS} --sequence-reweighting --wg enet --cpu 4 --continuous --alpha 1 --save-model ${pheno}_wgenet_fb > res_${pheno}_wgenet_fb.txt
#echo "###########################################################################################################################################################"

