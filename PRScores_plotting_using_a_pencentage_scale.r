#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
stfu = suppressPackageStartupMessages
stfu(library(crayon))
stfu(library(data.table))
stfu(library(stringr))
stfu(library(tidyverse))
green = make_style(rgb(red=0,green=250,blue=0,maxColorValue=255))


if (length(args)==0) {
  stop("At least one argument must be supplied (final report file, associacion table, variants table, barcode ).", call.=FALSE)
} else if (length(args)==1) {
  file=args[1];
}

main <- function(file){
        ########
        #file="/Users/bioinfo/Documents/SCRIPTS/test_data/PRSice_SCORES_AT_ALL_THRESHOLDS.PERCENTAGE.20210704.txt"
        scores = as.data.frame(fread(file,header=T))
        ########
        Pscores = setDT(scores)
        mPscores = melt.data.table(Pscores,"IID")
        ########
        mPscores$color=mPscores$value>0
        mPscores$absolute=abs(mPscores$value)
        ########
        mPscores=filter(mPscores,variable!="pT_0")
        mPscores$IID=gsub("AFECTADO","POSITIVO",mPscores$IID)
        mPscores$IID=gsub("PROTEGIDO","NEGATIVO",mPscores$IID)
        ########
        PRS_percentage=ggplot(data=mPscores, aes(x=variable, y=value, fill=IID)) + labs(y="Polygenic Score (pencentage)", x = "Thresholds of p-values") +
          geom_bar(stat="identity", position=position_dodge()) + ggtitle("Normalized Polygenic Scores in a percentage scale using P-value thresholding")
        ########    
        outputfile<-rename_output(file)
        ggsave(filename=outputfile,plot=PRS_percentage,width = 28,height = 20,units = c("cm"))
}

rename_output <- function(File){
  executiondate=format(Sys.time(), "%Y%m%d")
  sufix=gsub(' ','',paste(".percentage.",executiondate,".png",collapser=""))
  Outputfilename=gsub(".txt",sufix, File)
  return(Outputfilename)
}


## MAIN
if (!interactive()){
  arg=commandArgs(TRUE)
  if(length(arg)!=0){
    tm = system.time(main(file))
    message("\n Running time: ")
    print(tm)
    message("Run complete\t",green("!"))
  } else {
    #display_help()
  }
} else {
  message("Sourcing ...")
}


