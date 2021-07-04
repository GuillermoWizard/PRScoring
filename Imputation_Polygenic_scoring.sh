#!/usr/bin/env bash
# Author:       Luis Jaramillo and Willebaldo G.
# Email:        ljaramillo@codigo46.com.mx, wgarcia@codigo46.com.mx 
# Date:         2021/05/10
# Usage:        Scriptname.sh <folder>
# Copyright: 	Codigo46 S.A. de C.V.
#
# Este es el script front del desarrollo,
# activa el entorno de imp y ejecuta el script
#
#


################################################################################
Help()
{
   # Display Help
   echo 
   echo "The pipeline imputes Genotype Data, selects target variants and calculates Polygenic Score of a trait"
   echo 
   echo "Syntax: Imputation_Polygenic_scoring.sh [A] [B] [C]"
   echo "options:"
   echo "[A]     Directory with Genotype Data in Final report format (Illumina genotype format)"
   echo "[B]     Association tsv file (GWAS summary must be a .tsv) in PRSice format (With columns SNP, CHR, BP, A1, A2, OR/BETA, SE, P)"
   echo "[C]     low threshold of P"
   echo "[D]     high threshold of P"
   echo "[E]     step of p value to be used for PRS analysis from C to D"
   echo ""
}

################################################################################
################################################################################
# Main program                                                                 #
################################################################################
################################################################################

Help


SECONDS=0;
a="$(echo $CONDA_EXE | sed 's/bin\/conda/\t/')"   # CONDA_EXE es una variable predefinida en la instalación de Conda
CONDA_SH=$(echo "$a/etc/profile.d/conda.sh" | sed 's/\t//g')
echo "sourcing ... $CONDA_SH"
source $CONDA_SH


# Batch a partir de un folder con todos los reportes en txt
folder="$1"
#analysisvcf="$2"
assocfile="$2"
a="$3"
b="$4"
s="$5"


# Polygeic Risk score 
prciser=/media/datashare/data/PRSice_v1.25/PRSice_v1.25.R
plinkprs=/media/datashare/data/PRSice_v1.25/plink_1.9_linux_160914
randomphenotyper=/home/centos/Wily/bash_production_scripts/Random_Phenotype_table.R
assoc_cleaner=/home/centos/Wily/bash_production_scripts/Clean_Repeats_Keep_significant.R
snp_selector=/home/centos/Wily/PRScoring/Select_target_variants.R
normalizerPRS=/home/centos/Wily/PRScoring/Normalize_Recompute_Polygenic_Scores.r

#instrument_path="/media/datashare/data/gwas_catalog/covid19/"

runningfolder=$(pwd)

for file in $(ls $folder/*_all_variants.txt); do
  #praefix="$(echo ${file/.txt})"
  praefix=$(echo $file | sed 's/\.txt//g')
  echo "Praefix:  $praefix" 

  mkdir $praefix
  cp $file $praefix
  cd $praefix

  path=$(pwd)
  inputfile=$(ls $praefix)
  echo -e "Working directory: $path" 
  echo -e "File: $inputfile"
 
  #cd .. 
  #vcf_file= "$(echo $praefix.vcf)"
  #inputfile= $(echo "$praefix/") 

  finalreport-vcf-parser.R $inputfile

  #echo "Conversion a VCF: $vcf_file"

  echo "Output written at directory: $path"
  cd $folder

done



#### data imputation

conda activate c46-1.x
cd $folder
wholeVCFfiles=""
for dir in $(ls -d */); do 
cd $dir
	for vcf in $(ls *_all_variants.vcf); do
	vcfimputednames="";
		 for chr in $(seq 1 22); do \
		    praefix="$(echo ${vcf/.vcf})"
		    imp.py -i $vcf -c $chr  -e beagle -o $praefix.$chr 
    		    rm $praefix.${chr}.log
		    vcfcompressname=$(echo "$praefix.$chr.vcf.gz");
		    tabix -p vcf $vcfcompressname 
		    vcfimputednames=$(echo "$vcfimputednames $vcfcompressname");		
		    echo "$vcfimputednames"		
		  done
	outputimputed=$(echo "$praefix.complete.vcf");
	echo "Building $outputimputed";
	/media/datashare/custom-bins/bcftools/bcftools concat $vcfimputednames -o $outputimputed 
	bgzip -c $outputimputed >$outputimputed.gz 
	tabix -p vcf $outputimputed.gz
	mv $outputimputed.gz $folder   
	mv $outputimputed.gz.tbi $folder 
	wholeVCFfiles=$(echo "$wholeVCFfiles $outputimputed.gz")
	done
cd $folder
done
conda deactivate 

### End of imputation

cd $folder 
datetoday=$(date '+%Y%m%d')

### These are the names of the VCF files with imputation data 
echo "$wholeVCFfiles"

### This is the name of the whole big vcf file 
wholeimputedfile=$(echo "VCF_samples_file_$datetoday.vcf")
/media/datashare/custom-bins/bcftools/bcftools merge $wholeVCFfiles -o $wholeimputedfile
bgzip -c $wholeimputedfile >$wholeimputedfile.gz #zippedversion
tabix -p vcf $wholeimputedfile.gz #indexed version 

### selection variants in vcf for several samples
selectedvariantsfile=$(echo "VCF_samples_file_variants_selected_$datetoday.vcf")
Rscript --vanilla $snp_selector $assocfile

### Get the name of the SNP IDS of target trait
directoryassoc=$(dirname $assocfile)
listfile=$(echo "$directoryassoc/SNP_list_directory") # intermediate file 
snplistfile=$(cat $listfile) ## actual file with snps ids 


#### We use bcftools to perform the snp selection from the whole big VCF file 
/media/datashare/custom-bins/bcftools/bcftools view -i "ID=@$snplistfile" $wholeimputedfile -o $selectedvariantsfile 

## bgzip -c $analysisvcf >$analysisvcf.gz
## tabix -p vcf $analysisvcf.gz
## echo "$wholeVCFfiles"
## selectedvariantsfile=$(echo "VCF_samples_file_variants_selected_$datetoday.vcf")
## /media/datashare/custom-bins/bcftools/bcftools isec -c none -n=2 -w1 $wholeimputedfile.gz $analysisvcf.gz -o $selectedvariantsfile

### Create a new prefix for the plink files 
plinkfileprefix=$(echo "VCF_samples_file_variants_selected_$datetoday")
plink --vcf $selectedvariantsfile --recode tab --double-id --out $plinkfileprefix #recode vcf to plink 
plink --file $plinkfileprefix --make-bed --out $plinkfileprefix # recode plink to binary plink 

### For Pheno file Construction
#### Identified the ped file so we can build the pheno file for prcise
pedfile=$(echo "$plinkfileprefix.ped")
samplesIDS=$(cut -f2 $pedfile) ### select sample IDS
phenofile=$(echo "$plinkfileprefix.pheno") ### Pheno file name 
phenotinputs=$(echo "$samplesIDS $phenofile") ### sample IDs and pheno filename will be used as arguments 
echo "$pedfile"
echo "$phenotinputs" 

### Run script to build a random pheno assign to samples 
Rscript --vanilla $randomphenotyper $phenotinputs 


### Clean Assoc file to assoc clean file free of repeated variants 
#echo -e "$assoc_cleaner $selectedvariantsfile $assocfile"
#Rscript --vanilla $assoc_cleaner $selectedvariantsfile $assocfile # run asssoc cleaner
#newassocfile=$(echo "$assocfile" | sed 's/assoc/clean.assoc/g'); # name of new assoc file 

### Run PRSice analysis 
R --file=$prciser --args plink $plinkprs base $assocfile target $plinkfileprefix slower $a sinc $s supper $b no.regression T covary F allow.no.sex T pheno.file $phenofile debug.mode T 

resultspath=$(pwd); resultsfile=$(echo "$resultspath/PRSice_SCORES_AT_ALL_THRESHOLDS.txt");

Rscript --vanilla $normalizerPRS $resultsfile


### End message (running time)
echo "Running time: $SECONDS seconds";




