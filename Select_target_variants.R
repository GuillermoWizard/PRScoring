#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
stfu = suppressPackageStartupMessages
stfu(library(crayon))
stfu(library(data.table))
stfu(library(stringr))
stfu(library(tidyverse))
green = make_style(rgb(red=0,green=250,blue=0,maxColorValue=255))
 

if (length(args)==0) {
  display_help()
  stop("At least one argument must be supplied (final report file, associacion table, variants table, barcode ).", call.=FALSE)
  display_help()
} else if (length(args)==1) {
  file1 = args[1]; # assoc 
}


main <- function(AssocFile=file1){
  Assoc = as.data.frame(fread(AssocFile,header=T))
  SNPlist=as.data.frame(Assoc$SNP)
  colnames(SNPlist)=c("SNP")
  #print(Assoc$SNP)
  outputfile<-rename_output(AssocFile)
  write.table(SNPlist, file = outputfile, sep = "\t",col.names = TRUE,row.names=FALSE,quote = FALSE)
  directory=dirname(outputfile)                                          
  SNP_list_filename <- file.path(directory, "SNP_list_directory");
  writeLines(outputfile, SNP_list_filename)
  message(paste("List of SNPS written at:",outputfile))
}



rename_output <- function(name){
  executiondate=format(Sys.time(), "%Y%m%d")
  sufix=gsub(' ','',paste(".SNP.IDs.",executiondate,".tsv",collapser=""))
  #Outputfilename=gsub(".txt", sufix, name)
  Outputfilename=gsub(".tsv", sufix, name)
  #Outputfilename=gsub(".assoc", sufix, name)
  return(Outputfilename)
}

## MAIN
if (!interactive()){
  #message("\n\n",yellow$underline$bold("Código 46 S.A. de C.V."),green("Generate SNP list from Assoc file (SNP column)"),"\n")
  arg=commandArgs(TRUE)
  if(length(arg)!=0){
    tm = system.time(main(file1))
    #message("\n Running time: ")
    #print(tm)
    #message("Run complete\t",green("OK!"))
  } else {
    display_help()
  }
} else {
  message("Sourcing ...")
}
