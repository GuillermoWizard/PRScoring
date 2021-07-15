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
  ######## Scores of a population
  plot2=ggplot(data=scores, aes(x=reorder(IID,-pT_0.000001), y=pT_0.000001, fill=IID)) + labs(y="Polygenic Score (pencentage)", x = "Threshold of p-value") + geom_bar(stat="identity", position=position_dodge()) + theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ######## Scores of a population in a Histogram 
  scores_hist=scores[grep(pattern="PROTEGIDO*|AFECTADO*|NEUTRO",scores$IID,perl=TRUE,invert = TRUE),]
  ni=nrow(scores_hist)
  mini=min(scores_hist$pT_0.000001)-10
  maxi=max(scores_hist$pT_0.000001)+10
  histogram1=ggplot(data=scores_hist, aes(x=pT_0.000001)) + labs(y="Prevalence (density function)", x = "Polygenic Score (pencentage)") + geom_histogram(aes(y=..density..), position="identity", alpha=0.5, bins=ni) + theme(legend.position = "none") + ggtitle("Histogram of Normalized Polygenic Scores in a percentage scale") + geom_density(alpha=.2, fill="#FF6666") + xlim(mini,maxi)

  ########    
  outputfile<-rename_output(file)
  outputfile2<-rename_output_hist(file)
  ggsave(filename=outputfile,plot=plot2,width = 28,height = 20,units = c("cm"))
  ggsave(filename=outputfile2,plot=histogram1,width = 28,height = 20,units = c("cm"))
  
}

rename_output <- function(File){
  executiondate=format(Sys.time(), "%Y%m%d")
  sufix=gsub(' ','',paste(".percentage.",executiondate,".png",collapser=""))
  Outputfilename=gsub(".txt",sufix, File)
  return(Outputfilename)
}

rename_output_hist <- function(File){
  executiondate=format(Sys.time(), "%Y%m%d")
  sufix=gsub(' ','',paste(".histogram.",executiondate,".png",collapser=""))
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


