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
   echo "The pipeline imputes Genotype Data from a VCF"
   echo 
   echo "Syntax: Imputation_Polygenic_scoring.sh [A] [B] [C]"
   echo "options:"
   echo "[A]     Directory with Genotype Data format VCF"
   echo ""
}

################################################################################
################################################################################
# Main program                                                                 #
################################################################################
################################################################################

Help

folder="$1"

SECONDS=0;
a="$(echo $CONDA_EXE | sed 's/bin\/conda/\t/')"   # CONDA_EXE es una variable predefinida en la instalación de Conda
CONDA_SH=$(echo "$a/etc/profile.d/conda.sh" | sed 's/\t//g')
echo "sourcing ... $CONDA_SH"
source $CONDA_SH

conda activate c46-1.x
cd $folder  ##Entra al folder del VCF 
#wholeVCFfiles=""; ### String for name storing  
#for dir in $(ls -d */); do ### dentro del directorio padre se van a leer todos los directorios, los cuales deben contener al menos un VCF
#cd $dir ## se accede a cada uno de esos subdirectorios 
	for vcf in $(ls *.vcf); do ### se lee cada uno de los vcf en el subdirectorio (debe ser uno por subdirectorio)
	vcfimputednames=""; ### los nombres de los archivos imputados se van a guardar en una variable   
		 for chr in $(seq 1 22); do \ 
		    praefix="$(echo ${vcf/.vcf})" ## se edita el prefijo del nombre 
		    imp.py -i $vcf -c $chr  -e beagle -o $praefix.$chr ##se van a imputar datos del cromosoma 1 al 22 para cada muestra
    		rm $praefix.${chr}.log ## elimina todos los logs generados por beagle 
		    vcfcompressname=$(echo "$praefix.$chr.vcf.gz"); ### identifica el nombre comprimido de los datos imputados 
		    tabix -p vcf $vcfcompressname #### genera los indices con tabix 
		    vcfimputednames=$(echo "$vcfimputednames $vcfcompressname"); ## Concatena todos los nombres de los comprimidos imputados para el vcf merge 		
		    #echo "$vcfimputednames"	### todos estos se van fusionar en un gran archivo gigante 	
		  done
	outputimputed=$(echo "$praefix.complete.vcf"); ## se genera el nombre del archivo gigante 
	#echo "Building $outputimputed"; ## se envia un mensaje al display
	/media/datashare/custom-bins/bcftools/bcftools concat $vcfimputednames -o $outputimputed ### aca se fusionan
	echo "Concatenating ... $vcfimputednames to $outputimputed"	### todos estos se van fusionar en un gran archivo gigante 	

	bgzip $outputimputed ### se comprime el vcf gigante
	tabix -p vcf $outputimputed.gz ## se genera su indice 
	mv $outputimputed.gz $folder   ## 
	mv $outputimputed.gz.tbi $folder 
	#wholeVCFfiles=$(echo "$wholeVCFfiles $outputimputed.gz")
	done
cd $folder
#done
conda deactivate 


