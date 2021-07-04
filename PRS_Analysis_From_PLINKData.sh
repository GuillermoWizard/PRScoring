#!/usr/bin/env bash
# Author:       Willebaldo G.
# Email:        wgarcia@codigo46.com.mx 
# Date:         2021/05/10
# Usage:        Scriptname.sh <folder>
# Copyright: 	Codigo46 S.A. de C.V.
#
# Este es el script front del desarrollo,
# activa el entorno de imp y ejecuta el script
#
#

#SECONDS=0;
#a="$(echo $CONDA_EXE | sed 's/bin\/conda/\t/')"   # CONDA_EXE es una variable predefinida en la instalación de Conda
#CONDA_SH=$(echo "$a/etc/profile.d/conda.sh" | sed 's/\t//g')
#echo "sourcing ... $CONDA_SH"
#source $CONDA_SH

Help()
{
   # Display Help
   echo "Runs a PRSice polygenic score analysis based on target and base data using a P threshold (range)"
   echo
   echo "Syntax: PRS_Analysis_From_PLINKData.sh [target genotypes] [base data] [lowest P] [highest P] [step in P units]"
   echo "options:"
   echo "[target genotypes]	Genotype data in plink binary format (PED,BIM,FAM). Use only the prefix."
   echo "[base data]   Base data / summary statistics in PRSice format."
   echo "[lowest P]    Lowest value of P value to compute PS"
   echo "[highest P]   Highest value of P value to compute PS"
   echo "[Step]        Step of the P thresholds"
   #echo "[phenotypes]  Phenotype file in plink format of cases and controls [0,1,2,-9]"
   echo ""
}


Help

# Batch a partir de un folder con todos los reportes en txt
plinkprefix="$1"
assocfile="$2"
#phenofile="$3"
plow="$3"
phigh="$4"
step="$5"


# Polygeic Risk score 
prciser=/media/datashare/data/PRSice_v1.25/PRSice_v1.25.R
plinkprs=/media/datashare/data/PRSice_v1.25/plink_1.9_linux_160914
# Pheno file 
randomphenotyper=/home/centos/Wily/bash_production_scripts/Random_Phenotype_table.R
normalizerPRS=/home/centos/Wily/PRScoring/Normalize_Recompute_Polygenic_Scores.r
#assoc_cleaner=/home/centos/Wily/bash_production_scripts/Clean_Repeats_Keep_significant.R

### For Pheno file Construction
#### Identified the ped file so we can build the pheno file for prcise
pedfile=$(echo "$plinkprefix.ped")
samplesIDS=$(cut -f2 $pedfile) ### select sample IDS
phenofile=$(echo "$plinkprefix.pheno") ### Pheno file name 
phenotinputs=$(echo "$samplesIDS $phenofile") ### sample IDs and pheno filename will be used as arguments 
#echo "$pedfile"
#echo "$phenotinputs" 

### Run script to build a random pheno assign to samples 
Rscript --vanilla $randomphenotyper $phenotinputs 

R --file=$prciser --args plink $plinkprs base $assocfile target $plinkprefix slower $plow sinc $step supper $phigh no.regression T covary F allow.no.sex T pheno.file $phenofile debug.mode T 

resultspath=$(pwd); resultsfile=$(echo "$resultspath/PRSice_SCORES_AT_ALL_THRESHOLDS.txt");

Rscript --vanilla $normalizerPRS $resultsfile



