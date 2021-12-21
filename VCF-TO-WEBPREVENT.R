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
  #filevcf="/Users/bioinfo/Documents/QUERIES_DB/V727007_all_variants.vcf"
  vcf <- read.vcfR( filevcf, verbose = FALSE )
  
  ######## FIXED GENETIC INFORMATION
  ######## NAME CORRECTION
  vcf_fix=as.data.frame(vcf@fix); 
  
  ######## GENOTYPE
  vcf_genotype=as.data.frame(vcf@gt);
  sample_names=colnames(vcf_genotype);
  print("sample_names")
  print(sample_names)
  
  ########
  n=nrow(vcf_genotype)
  m=ncol(vcf_genotype)
  
  ########
  genotype2=apply(vcf_genotype,2,str_split_fixed,rep(":",n),rep(2,n));
  n1=nrow(genotype2);
  genotype2=genotype2[1:n,];
  n2=nrow(genotype2) 
  
  print(c(n,n1,n2))
  
  ########
  for(i in 2:m){
    Genotype_VCF_table_selected=""
    SAMPLENAME=as.data.frame(rep(sample_names[i],n2))
    print(c("Processing:",i,SAMPLENAME[i,]))
    colnames(SAMPLENAME)=c("SAMPLEID")
    #print(genotype2[,i])  
    individual_genotype_splitted=t(as.data.frame(sapply(genotype2[,i],str_split_fixed,"\\|",2)));
    #print(t(as.data.frame(sapply(genotype2[,i],str_split_fixed,"\\|",2))))
    colnames(individual_genotype_splitted)=c("S1","S2")

    #print(colnames(SAMPLENAME))
    Genotype_VCF_table=cbind(vcf_fix,SAMPLENAME,individual_genotype_splitted);
    #print(colnames(Genotype_VCF_table))

    ###
    STRING1=as.data.frame(mapply(gsub,rep(0,n2),Genotype_VCF_table$REF,Genotype_VCF_table$S1))
    colnames(STRING1)=c("ST1_1")
    ALLELES_ALT=t(as.data.frame(sapply(Genotype_VCF_table$ALT,str_split_fixed,",",3)));
    Genotype_VCF_table=cbind(Genotype_VCF_table,STRING1,ALLELES_ALT);
    Genotype_VCF_table_selected=Genotype_VCF_table %>% dplyr::select("SAMPLEID","ID","CHROM","POS","REF","ALT","1","2","3","S1","S2","ST1_1")
    #####
    STRING1=as.data.frame(mapply(gsub,rep(1,n2),Genotype_VCF_table_selected$"1",Genotype_VCF_table_selected$ST1_1))
    colnames(STRING1)=c("ST1_2")
    Genotype_VCF_table_selected = cbind(Genotype_VCF_table_selected,STRING1);
    ###
    STRING1=as.data.frame(mapply(gsub,rep(2,n2),Genotype_VCF_table_selected$"2",Genotype_VCF_table_selected$ST1_2))
    colnames(STRING1)=c("ST1_3")
    Genotype_VCF_table_selected = cbind(Genotype_VCF_table_selected,STRING1);
    ###
    STRING1=as.data.frame(mapply(gsub,rep(3,n2),Genotype_VCF_table_selected$"3",Genotype_VCF_table_selected$ST1_3))
    colnames(STRING1)=c("ST1_4")
    Genotype_VCF_table_selected = cbind(Genotype_VCF_table_selected,STRING1);
  
    ###
    STRING1=as.data.frame(mapply(gsub,rep(0,n2),Genotype_VCF_table_selected$REF,Genotype_VCF_table$S2))
    colnames(STRING1)=c("ST2_1")
    #ALLELES_ALT=t(as.data.frame(sapply(Genotype_VCF_table$ALT,str_split_fixed,",",3)));
    Genotype_VCF_table_selected=cbind(Genotype_VCF_table_selected,STRING1);
    #Genotype_VCF_table_selected=Genotype_VCF_table %>% dplyr::select("CHROM","POS","ID","REF","ALT","1","2","3","ST1")
    #####
    STRING1=as.data.frame(mapply(gsub,rep(1,n2),Genotype_VCF_table_selected$"1",Genotype_VCF_table_selected$ST2_1))
    colnames(STRING1)=c("ST2_2")
    Genotype_VCF_table_selected = cbind(Genotype_VCF_table_selected,STRING1);
    ###
    STRING1=as.data.frame(mapply(gsub,rep(2,n2),Genotype_VCF_table_selected$"2",Genotype_VCF_table_selected$ST2_2))
    colnames(STRING1)=c("ST2_3")
    Genotype_VCF_table_selected = cbind(Genotype_VCF_table_selected,STRING1);
    ###
    STRING1=as.data.frame(mapply(gsub,rep(3,n2),Genotype_VCF_table_selected$"3",Genotype_VCF_table_selected$ST2_3))
    colnames(STRING1)=c("ST2_4")
    Genotype_VCF_table_selected = cbind(Genotype_VCF_table_selected,STRING1);
    ####
    Genotype_VCF_table_selected=Genotype_VCF_table_selected %>% dplyr::select("SAMPLEID","ID","CHROM","POS","ST1_4","ST2_4");
    #print(head(Genotype_VCF_table_selected))
    if(i==2){
      Genotype_VCF_table_selected_whole_population=Genotype_VCF_table_selected;
    }
    if(i>2){
      Genotype_VCF_table_selected_whole_population=rbind(Genotype_VCF_table_selected_whole_population,Genotype_VCF_table_selected);
    }
    
  }
  ####### WRITE OUTPUT
  outputfile<-rename_output(filevcf)
  write.table(Genotype_VCF_table_selected_whole_population, file = outputfile, sep = "\t",col.names = FALSE,row.names=FALSE,quote = FALSE)
  message("\n\n",yellow$underline$bold(paste("Writing output at:",outputfile)))
  
}


rename_output <- function(namef=FRfile){
  executiondate=format(Sys.time(), "%Y%m%d")
  sufix=gsub(' ','',paste("_CODIGOPREVENT_WEB_",executiondate,".tsv",collapser=""))
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