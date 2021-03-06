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
green = make_style(rgb(red=0,green=250,blue=0,maxColorValue=255))
options(warn = defaultW)

helpmessage="
  Build_Reference_Genotypes.R [Assoc File] [VCF with samples genotypes]
    Inputs:
      [Assoc]   Tab separated file with (tsv extension) with summary statistics from GWAS in PRCSice format. 
      [VCF]     VCF file (vcf extension) with genotype data from sample individuals.      
  
  "

if (length(args)==0) {
  cat(helpmessage)
  stop("At least one argument must be supplied (assoc, vcf ).", call.=FALSE)
} else if (length(args)==1) {
  cat(helpmessage)
  stop("At least two arguments/files must be supplied.", call.=FALSE)
} else if (length(args)==2) {
  file1=args[1];
  file2=args[2];
}


main <- function(AssocFile, vcfFile){
  defaultW <- getOption("warn") 
  options(warn = -1)
  
  Assocdata = as.data.frame(fread(AssocFile,header=T))
  vcf <- read.vcfR(vcfFile, verbose = FALSE )
  
  ##### Metadata in header
  #metadata=(as.data.frame(vcf@meta))
  meta=vcf@meta
  
  ##### VCF genetic information 
  vcffix=(as.data.frame(vcf@fix));
  colnamesvcf=colnames(vcffix);
  ##### Samples Genotype data
  vcfGeno=(as.data.frame(vcf@gt))
  colnamesgenotypes=colnames(vcfGeno)
  
  #####
  vcfdata=cbind(vcffix,vcfGeno)
  n1=nrow(vcfdata)  
  #####
  fused_data=merge(vcfdata,Assocdata,by.x="ID",by.y="SNP")
  n2=nrow(fused_data)
  ###
  fused_data_neg=filter(fused_data,BETA<0)
  n3=nrow(fused_data_neg)
  ###
  fused_data_pos=filter(fused_data,BETA>0)
  n4=nrow(fused_data_pos)
  ### Generate simulated genotypes 
  fused_data_pos$AFFECTED1 = mapply(paste0,fused_data_pos$A1,rep("|",n4),fused_data_pos$A1)
  fused_data_neg$AFFECTED1 = mapply(paste0,fused_data_neg$A2,rep("|",n3),fused_data_neg$A2)
  fused_data_pos$AFFECTED2 = mapply(paste0,fused_data_pos$A1,rep("|",n4),fused_data_pos$A2)
  fused_data_neg$AFFECTED2 = mapply(paste0,fused_data_neg$A2,rep("|",n3),fused_data_neg$A2)
  fused_data_pos$NEUTRAL = mapply(paste0,fused_data_pos$A2,rep("|",n4),fused_data_pos$A2)
  fused_data_neg$NEUTRAL = mapply(paste0,fused_data_neg$A2,rep("|",n3),fused_data_neg$A2)
  fused_data_pos$PROTECTED1 = mapply(paste0,fused_data_pos$A2,rep("|",n4),fused_data_pos$A2)
  fused_data_neg$PROTECTED1 = mapply(paste0,fused_data_neg$A1,rep("|",n3),fused_data_neg$A1)
  fused_data_pos$PROTECTED2 = mapply(paste0,fused_data_pos$A2,rep("|",n4),fused_data_pos$A2)
  fused_data_neg$PROTECTED2 = mapply(paste0,fused_data_neg$A1,rep("|",n3),fused_data_neg$A2)
  ## bind negative and positive risk matrices of genotypes 
  vcf_simulated=rbind(fused_data_pos,fused_data_neg)
  n5=nrow(vcf_simulated)
  n6=ncol(vcf_simulated)
  #### Generate binary genoypes  
  vcf_simulated$AFECTADO1=matrix(mapply(gsub,vcf_simulated$REF,rep("0",n5),vcf_simulated$AFFECTED1),ncol = 1)### First 
  vcf_simulated$AFECTADO1=matrix(mapply(gsub,vcf_simulated$ALT,rep("1",n5),vcf_simulated$AFECTADO1),ncol = 1)
  vcf_simulated$AFECTADO1=matrix(mapply(gsub,"A",rep("1",n5),vcf_simulated$AFECTADO1),ncol = 1)
  vcf_simulated$AFECTADO1=matrix(mapply(gsub,"T",rep("1",n5),vcf_simulated$AFECTADO1),ncol = 1)
  vcf_simulated$AFECTADO1=matrix(mapply(gsub,"C",rep("1",n5),vcf_simulated$AFECTADO1),ncol = 1)
  vcf_simulated$AFECTADO1=matrix(mapply(gsub,"G",rep("1",n5),vcf_simulated$AFECTADO1),ncol = 1)
  ###
  vcf_simulated$AFECTADO2=matrix(mapply(gsub,vcf_simulated$REF,rep("0",n5),vcf_simulated$AFFECTED2),ncol = 1)### First 
  vcf_simulated$AFECTADO2=matrix(mapply(gsub,vcf_simulated$ALT,rep("1",n5),vcf_simulated$AFECTADO2),ncol = 1)
  vcf_simulated$AFECTADO2=matrix(mapply(gsub,"A",rep("1",n5),vcf_simulated$AFECTADO2),ncol = 1)
  vcf_simulated$AFECTADO2=matrix(mapply(gsub,"T",rep("1",n5),vcf_simulated$AFECTADO2),ncol = 1)
  vcf_simulated$AFECTADO2=matrix(mapply(gsub,"C",rep("1",n5),vcf_simulated$AFECTADO2),ncol = 1)
  vcf_simulated$AFECTADO2=matrix(mapply(gsub,"G",rep("1",n5),vcf_simulated$AFECTADO2),ncol = 1)
  ###
  vcf_simulated$NEUTRO=matrix(mapply(gsub,vcf_simulated$REF,rep("0",n5),vcf_simulated$NEUTRAL),ncol = 1)### First 
  vcf_simulated$NEUTRO=matrix(mapply(gsub,vcf_simulated$ALT,rep("1",n5),vcf_simulated$NEUTRO),ncol = 1)
  vcf_simulated$NEUTRO=matrix(mapply(gsub,"A",rep("1",n5),vcf_simulated$NEUTRO),ncol = 1)
  vcf_simulated$NEUTRO=matrix(mapply(gsub,"T",rep("1",n5),vcf_simulated$NEUTRO),ncol = 1)
  vcf_simulated$NEUTRO=matrix(mapply(gsub,"C",rep("1",n5),vcf_simulated$NEUTRO),ncol = 1)
  vcf_simulated$NEUTRO=matrix(mapply(gsub,"G",rep("1",n5),vcf_simulated$NEUTRO),ncol = 1)
  ###
  vcf_simulated$PROTEGIDO1=matrix(mapply(gsub,vcf_simulated$REF,rep("0",n5),vcf_simulated$PROTECTED1),ncol = 1)### First 
  vcf_simulated$PROTEGIDO1=matrix(mapply(gsub,vcf_simulated$ALT,rep("1",n5),vcf_simulated$PROTEGIDO1),ncol = 1)
  vcf_simulated$PROTEGIDO1=matrix(mapply(gsub,"A",rep("1",n5),vcf_simulated$PROTEGIDO1),ncol = 1)
  vcf_simulated$PROTEGIDO1=matrix(mapply(gsub,"T",rep("1",n5),vcf_simulated$PROTEGIDO1),ncol = 1)
  vcf_simulated$PROTEGIDO1=matrix(mapply(gsub,"C",rep("1",n5),vcf_simulated$PROTEGIDO1),ncol = 1)
  vcf_simulated$PROTEGIDO1=matrix(mapply(gsub,"G",rep("1",n5),vcf_simulated$PROTEGIDO1),ncol = 1)
  ###
  vcf_simulated$PROTEGIDO2=matrix(mapply(gsub,vcf_simulated$REF,rep("0",n5),vcf_simulated$PROTECTED2),ncol = 1)### First 
  vcf_simulated$PROTEGIDO2=matrix(mapply(gsub,vcf_simulated$ALT,rep("1",n5),vcf_simulated$PROTEGIDO2),ncol = 1)
  vcf_simulated$PROTEGIDO2=matrix(mapply(gsub,"A",rep("1",n5),vcf_simulated$PROTEGIDO2),ncol = 1)
  vcf_simulated$PROTEGIDO2=matrix(mapply(gsub,"T",rep("1",n5),vcf_simulated$PROTEGIDO2),ncol = 1)
  vcf_simulated$PROTEGIDO2=matrix(mapply(gsub,"C",rep("1",n5),vcf_simulated$PROTEGIDO2),ncol = 1)
  vcf_simulated$PROTEGIDO2=matrix(mapply(gsub,"G",rep("1",n5),vcf_simulated$PROTEGIDO2),ncol = 1)
  ####
  validcols=c(colnamesvcf,colnamesgenotypes,c("AFECTADO1","AFECTADO2","NEUTRO","PROTEGIDO1","PROTEGIDO2"))
  new_vcf = vcf_simulated %>% select(validcols)
  #### sort new vcf data 
  new_vcf=new_vcf[order(new_vcf[,1], new_vcf[,2] ),];
  colnames(new_vcf)[1] = paste0("#",validcols[1])
  outputfile<-rename_output(file2)
  writeLines(meta, outputfile)
  write.table(new_vcf, file = outputfile, sep = "\t",col.names = TRUE,row.names=FALSE,quote = FALSE,append = TRUE)
  options(warn = defaultW)
  cat(outputfile)
  
}


rename_output <- function(File){
  executiondate=format(Sys.time(), "%Y%m%d")
  sufix=gsub(' ','',paste(".Reference.vcf",collapser=""))
  Outputfilename=gsub(".vcf",sufix, File)
  return(Outputfilename)
}

## MAIN
if (!interactive()){
  #message("\n\n",yellow$underline$bold("C??digo 46 S.A. de C.V."),green("Generate SNP list from Assoc file (SNP column)"),"\n")
  arg=commandArgs(TRUE)
  if(length(arg)!=0){
    tm = system.time(main(file1,file2))
    #message("\n Running time: ")
    #print(tm)
    #message("Run complete\t",green("!"))
  } else {
    #display_help()
  }
} else {
  #message("Sourcing ...")
}


