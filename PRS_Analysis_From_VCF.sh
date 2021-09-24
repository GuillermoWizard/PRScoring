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
   echo "[A]     File path with Genotype Data in vcf format"
   echo "[B]     Association tsv file (GWAS summary must be a .tsv) in PRSice format (With columns SNP, CHR, BP, A1, A2, OR/BETA, P)"
   echo "[C]     low threshold of P"
   echo "[D]     high threshold of P"
   echo "[E]     step of p value to be used for PRS analysis from C to D"
   echo ""
}


# Batch a partir de un folder con todos los reportes en txt
vcf="$1"
assocfile="$2"
a="$3"
b="$4"
s="$5"


if [ -z $5 ]
then
   Help
   echo "Please input all the arguments";
elif [ -n $5 ]
then
# Polygeic Risk score 
   prciser=/opt/bin/PRSice/bin/PRSice
   randomphenotyper=/opt/repo/Will/bash_production_scripts/Random_Phenotype_table.R
   #assoc_cleaner=/opt/repo/Will/bash_production_scripts/Clean_Repeats_Keep_significant.R
   #snp_selector=/opt/repo/Will/PRScoring/Select_target_variants.R
   normalizerPRS=/opt/repo/Will/PRScoring/Normalize_Recompute_Polygenic_Scores.r
   buildreferencesvcf=/opt/repo/Will/PRScoring/Build_references_on_vcf.r

   #date
   datetoday=$(date '+%Y%m%d')
   #dirname
   directoryassoc=$(dirname $vcf);
   cd $directoryassoc;

   vcf_with_references=$(Rscript --vanilla $buildreferencesvcf $assocfile $vcf);

   ### Create a new prefix for the plink files 
   plinkfileprefix=$(echo "VCF_file_with_variants_selected_on_$datetoday")
   listfile=$(echo "$directoryassoc/$plinkfileprefix") # intermediate file 
   plinkfileprefix=$(echo $listfile)


   plink --vcf $vcf_with_references --recode tab --double-id --out $plinkfileprefix #recode vcf to plink
   plink --file $plinkfileprefix --make-bed --out $plinkfileprefix # recode plink to binary plink 


   ### For Pheno file build
   #### build the pheno file for prcise
   pedfile=$(echo "$plinkfileprefix.ped")
   samplesIDS=$(cut -f2 $pedfile) ### select sample IDS
   phenofile=$(echo "$plinkfileprefix.pheno") ### Pheno file name 
   phenotinputs=$(echo "$samplesIDS $phenofile") ### sample IDs and pheno filename will be used as arguments

   echo "$pedfile"
   echo "$phenotinputs" 

   ### Run script to build a random pheno assign to samples 
   Rscript --vanilla $randomphenotyper $phenotinputs

   ### Run PRSice analysis (POLYGENIC RISK SCORE ANALYSIS)
   $prciser --base $assocfile --target $plinkfileprefix --lower $a --interval $s --upper $b --no-regress --pheno $phenofile --ignore-fid 

   resultspath=$(pwd); 
   resultsfile=$(echo "$resultspath/PRSice_SCORES_AT_ALL_THRESHOLDS.txt");

   Rscript --vanilla $normalizerPRS $resultsfile
   
   ### Run PRSice analysis 
   $prciser --base $assocfile --target $plinkfileprefix --lower $a --interval $s --upper $b --no-regress --pheno $phenofile --ignore-fid 
   resultspath=$(pwd); 
   resultsfile=$(echo "$resultspath/PRSice_SCORES_AT_ALL_THRESHOLDS.txt");

   Rscript --vanilla $normalizerPRS $resultsfile
 

   ### End message (running time)
   echo "Running time: $SECONDS seconds";


fi
