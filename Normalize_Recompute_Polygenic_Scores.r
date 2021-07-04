#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
stfu = suppressPackageStartupMessages
stfu(library(crayon))
stfu(library(data.table))
stfu(library(stringr))
stfu(library(tidyverse))


if (length(args)==0) {
  stop("At least one argument must be supplied (final report file, associacion table, variants table, barcode ).", call.=FALSE)
} else if (length(args)==1) {
  file1=args[1];
} 


main <- function(PrsiceScores){
  Scores = as.data.frame(fread(PrsiceScores,header=T,sep=" "))
  n1=nrow(Scores)
  n2=ncol(Scores)

  Percentages=Scores;
  for(i in 2:n2){
    min=min(Scores[,i])
    max=max(Scores[,i])
    #print(c(min,max))
    for(j in 1:n1){
    if(Scores[j,i]<0){
      Percentages[j,i]=((100)*(-1)*(Percentages[j,i]/min));
    }
    if(Scores[j,i]>0){
      Percentages[j,i]=(100)*(Percentages[j,i]/max);
    }    
    }
  }
  outputfile<-rename_output(PrsiceScores)
  write.table(Percentages, file = outputfile, sep = "\t",col.names = TRUE,row.names=FALSE,quote = FALSE)
  message("\n\n",yellow$underline$bold(paste("Writing output at:",outputfile)))
  
}

rename_output <- function(File){
  executiondate=format(Sys.time(),"%Y%m%d")
  sufix=gsub(' ','',paste(".PERCENTAGE.",executiondate,".txt",collapser=""))
  Outputfilename=gsub(".txt",sufix, File)
  return(Outputfilename)
}


## MAIN
if (!interactive()){
  #message("\n\n",yellow$underline$bold("Código 46 S.A. de C.V."),green("Generate SNP list from Assoc file (SNP column)"),"\n")
  arg=commandArgs(TRUE)
  if(length(arg)!=0){
    tm = system.time(main(file1))
    message("\n Running time: ")
    print(tm)
    message("Run complete\t",green("!"))
  } else {
    #display_help()
  }
} else {
  message("Sourcing ...")
}

