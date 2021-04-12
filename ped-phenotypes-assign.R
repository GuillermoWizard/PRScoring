#!/usr/bin/env Rscript
# Author:       Guillermo G
# Email:        wgarcia@codigo46.com.mx
# Date:         2021/01/31
# Version:      1.0
# Usage:        
# Description:  Script

###
stfu = suppressPackageStartupMessages
stfu(library(crayon))
stfu(library(data.table))
stfu(library(stringr))
stfu(library(tidyverse))
stfu(library(ggplot2))
stfu(library(readr))
#library(reshape)
green = make_style(rgb(red=0,green=250,blue=0,maxColorValue=255))

args = commandArgs(trailingOnly=TRUE)


display_help <- function(){
  message("
  ## Assign Phenotypes to PED data rows
	## Usage:
        $ ped-phenotypes-assign.R [ped] [pheno]
  ## Description
  ## Inputs:
      [ped]:    PED file with individual IDS matching the individual IDS from Pheno file.
      [pheno]:  PHENO file containing two columns. First column has the individual ID from the PED and the second one contains the Phenotype (1,2,-9, etc)  
")}


if (length(args)==0) {
  #display_help()
  print(length(args))
  display_help()
  stop("At least two arguments must be supplied (ped and pheno files).", call.=FALSE)
  
  
  
} else if (length(args)==1) {
  # default output file
  #display_help()
  print(length(args))
  display_help()
  stop("At least two arguments must be supplied (ped and pheno files).", call.=FALSE)
  
  
} else if (length(args)==2) {
  file1 = args[1]; #ped
  file2 = args[2]; #pheno 
}







main<-function(file1,file2){
##############
#file1="/Users/Robot/Documents/Descargas_Codigo46/Investigación/PRS/Bipolar_disorder/ALL.SAMPLES.COHORT.HEAD.ped"
ped=as.data.frame(fread(file1,header=F))
##############
#file1="/Users/Robot/Documents/Descargas_Codigo46/Investigación/PRS/Bipolar_disorder/ALL.SAMPLES.COHORT.pheno"
pheno=as.data.frame(fread(file2,header=F))
ncol_ped=ncol(ped)
pedsamples=ped[,7:ncol_ped]
pedinfo=ped[,1:6]
#head(pheno)
ped_pheno=merge(pedinfo,pheno,by.x="V2",by.y="V1")
ped_pheno$V1=ped_pheno$V2
ped_pheno_modified  = ped_pheno %>% select('V1','V2','V3','V4','V5','V2.y')
ped_data = cbind(ped_pheno_modified,pedsamples)
#ncol(ped_data)
outname=rename_output(file1)
message("\n\n",yellow$underline$bold(paste("Writing output at:",outname)))
write.table(ped_data, file = outname,col.names = FALSE, row.names=F, sep="\t",quote = F)


}

rename_output <- function(PEDfile){
  #executiondate=format(Sys.time(), "%Y%m%d")
  #gsub("","ot",x)
  Outputfilename=gsub(".ped", ".Random.Pheno.ped", PEDfile)
  #Outputfilename=gsub(".PED", ".GENENAME.vcf.tsv", VCFFile)
  return(Outputfilename)
}

## MAIN
if (!interactive()){
  message("\n\n",yellow$underline$bold("Código 46 S.A. de C.V."),green("RANDOM ASSIGN PHENOTYPES TO PED DATA ROWS"),"\n")
  args=commandArgs(TRUE)
  if(length(args)==2){
    #PWD=arg[1] 
    #input_file_2=arg[2] 
    tm = system.time(main(file1,file2))
    message("\nRunning time: ")
    print(tm)
    message("Complete\t",green("OK! =✪ ᆺ ✪="))
  } else {
    print(length(args))
    display_help()
    stop("At least two arguments must be supplied (ped and pheno files).", call.=FALSE)
    
  }
} else {
  print(length(args))
  
  message("Sourcing ...")
}