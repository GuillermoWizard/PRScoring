#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
stfu = suppressPackageStartupMessages
defaultW <- getOption("warn") 
options(warn = -1)
stfu(library(crayon))
stfu(library(data.table))
stfu(library(stringr))
stfu(library(tidyverse))
stfu(library(vcfR))
stfu(library(biomaRt))
stfu(library(dplyr))
stfu(library(fs))


display_help <- function(){
  message("
  ## Cambia VCF ID names from Illumina Names to DBSNP ID (RS)
	## Método de Uso:
        $ Rename_illuminanames_vcf.R [VCF] 
	Procesamiento de final reports individuales.
  ## Descripción
  ## Entradas:
      [VCF FILE]: Archivo VCF  
")}



if (length(args)==0) {
    display_help()
    stop("At least one arguments must be supplied (VCF file).", call.=FALSE)
    display_help()
} else if (length(args)==1) {
    # default output file
    file1 = args[1];
} 

main <- function(filevcf){
    filevcf="/Users/bioinfo/Documents/QUERIES_DB/V727007_all_variants.vcf"
    vcf <- read.vcfR( filevcf, verbose = FALSE )
    ########
    vcf_fix=as.data.frame(vcf@fix); 
    vcf_fix$ID=gsub('GSA-','',vcf_fix$ID); # 971 - 977 - 981
    vcf_fix$ID=gsub('seq-','',vcf_fix$ID); # 652, 318, 60 
    vcf_fix$ID=gsub('IlmnSeq_','',vcf_fix$ID); #980, 780, 781 
    vcf_fix$ID=gsub('Ilmnseq_','',vcf_fix$ID); # vcf_fix$ID[3000:4000] 9, 21, 33
    vcf_fix$ID=gsub('.1','',vcf_fix$ID); # 318 (vcf_fix$ID[0:1000]), 466, 474, 472 (vcf_fix$ID[4000:5000])
    vcf_fix$ID=gsub('.2','',vcf_fix$ID); 
    vcf@fix[,3]=vcf_fix$ID;
    
    outputfile<-rename_output(filevcf)
    print(outputfile)
    write.vcf( vcf, file = outputfile, mask = FALSE, APPEND = FALSE )
    
    #vcf_fix_char=as.character(vcf_fix)
    #print(vcf_fix$ID)
    #print(vcf_fix_char)
    
}


rename_output <- function(namef=FRfile){
  executiondate=format(Sys.time(), "%Y%m%d")
  sufix=gsub(' ','',paste("_",executiondate,".vcf.gz",collapser=""))
  Outputfilename=gsub(".vcf", sufix, namef)
  return(Outputfilename)
}


## MAIN
if (!interactive()){
  message("\n\n",yellow$underline$bold("Código 46 S.A. de C.V."),green("Rename Illumina Names"),"\n")
  arg=commandArgs(TRUE)
  if(length(arg)!=0){
    input_file=arg[1] # Finalreport de lectura
    tm = system.time(main(file1))
    message("\nTiempo de Procesamiento: ")
    print(tm)
    message("Proceso Completado\t",green("OK! :)"))
  } else {
    display_help()
  }
} else {
  message("Sourcing ...")
}
