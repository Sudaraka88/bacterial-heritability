#! /bin/bash

NAME=${1?Error: Filename required...}

# Check for snp-sites, install if missing
if [ $(dpkg-query -W -f='${Status}' snp-sites 2>/dev/null | grep -c "ok installed") -eq 0 ];
then 
	sudo apt-get install snp-sites;
fi
# Check for vcftools, install if missing
if [ $(dpkg-query -W -f='${Status}' vcftools 2>/dev/null | grep -c "ok installed") -eq 0 ];
then
	sudo apt-get install vcftools;
fi
# Check for samtools, install if missing
if [ $(dpkg=query -W -f='${Status}' samtools 2>/dev/null | grep -c "ok installed") -eq 0 ];
then
	sudo apt-get install samtools;
fi

# Now we can output frequency stats
echo "Input fasta file:" $NAME
echo "Running snp-sites..."
snp-sites -v $NAME -o $NAME".vcf"
vcftools --vcf $NAME".vcf" --out $NAME --freq
