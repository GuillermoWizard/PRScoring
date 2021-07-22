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
  Classify_variants_of_a_trait.R [Assoc File] [VCF with samples genotypes]
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

AssocFile = "/Users/bioinfo/Documents/ALCHEMY_C46_UAQ_MDD_BPD/PRS_Broad_Depression_14_julio/2.41467_2018_3819_MOESM5_ESM.DATASET3.BETA.tsv"
vcfFile = "/Users/bioinfo/Documents/ALCHEMY_C46_UAQ_MDD_BPD/PRS_Broad_Depression_14_julio/Alchemy-C46-Samples-beagle.BROAD.DEPRESSION.SNV.Reference.vcf"
main <- function(AssocFile, vcfFile){
  defaultW <- getOption("warn") 
  options(warn = -1)
  ########
  Assocdata = as.data.frame(fread(AssocFile,header=T))
  vcf = read.vcfR(vcfFile, verbose = FALSE )
  ##### Metadata in header
  meta=vcf@meta
  ##### VCF info 
  vcffix=(as.data.frame(vcf@fix));
  colnamesvcf=colnames(vcffix);
  ##### Samples Genotype data
  vcfGeno=(as.data.frame(vcf@gt))
  colnamesgenotypes=colnames(vcfGeno)
  ##### Genotypes + vcf info 
  vcfdata=cbind(vcffix,vcfGeno)
  n1=nrow(vcfdata)  
  ##### 
  fused_data=merge(vcfdata,Assocdata,by.x="ID",by.y="SNP")
  n2=nrow(fused_data)
  ### Genotype data is split from other non useful data
  editedsamples=vcfGeno;
  n2=ncol(vcfGeno);
  
  for(i in 1:n2){
    edited_genotype=str_split_fixed(vcfGeno[,i],":",2)
    colnames(edited_genotype)=c("A","B")
    edited_genotype=as.data.frame(edited_genotype)
    editedsamples[,i]=edited_genotype$A;
  } ### bind info vcf and genotypes 
  vcfdata=cbind(vcffix,editedsamples)
  ###
  n1=ncol(vcfdata)
  n2=nrow(vcfdata)
  #print(c(n1,n2))
  numberofsamples=n1-9;
  ### Los genotipos en 0/1 se pasan a las versiones CGTA (ALT/REF)
  ### Samples data should start in column 10 
  genotype_samples=vcfdata
  for(i in 1:n2){
    #print(genotype_samples[i,10:n1])
    genotype_samples[i,10:n1] = gsub("0",genotype_samples$REF[i],genotype_samples[i,10:n1]);
    genotype_samples[i,10:n1]=gsub("1",genotype_samples$ALT[i],genotype_samples[i,10:n1])
    genotype_samples[i,10:n1]=gsub("2",str_split_fixed(genotype_samples$ALT,",",3)[i,2],genotype_samples[i,10:n1])
    genotype_samples[i,10:n1]=gsub("3",str_split_fixed(genotype_samples$ALT,",",3)[i,3],genotype_samples[i,10:n1])
    
  }
  ##### corre bien hasta aqui
  #head(genotype_samples)
  ### VCF and ASSOC merge. 
  merged_genotype_assoc=merge(Assocdata,genotype_samples,by.y="ID",by.x="SNP");
  #head(merged_genotype_assoc)
  #####
  n1=ncol(merged_genotype_assoc)
  n2=nrow(merged_genotype_assoc)
  #####
  #newadddedcols=merged_genotype_assoc[,(c(n1-(numberofsamples-1))):n1];
  newadddedcols=apply(merged_genotype_assoc[,(c(n1-(numberofsamples-1))):n1],2,str_count,merged_genotype_assoc$A1)
  merged_genotype_assoc=cbind(merged_genotype_assoc,newadddedcols)
  ##### -1, -2 indicate protection states on SNVs
  ##### 1, 2 indicate risk states on SNVS by person/sample
  merged_genotype_assoc[as.numeric(which(merged_genotype_assoc$BETA<0)),c((n1+1)):(c(n1+numberofsamples))]=merged_genotype_assoc[as.numeric(which(merged_genotype_assoc$BETA<0)),c((n1+1)):(c(n1+numberofsamples))]*(-1)
  table_of_variants=apply(merged_genotype_assoc[,c((n1+1)):(c(n1+numberofsamples))],2,table)
  names_table_of_variants=data.frame(rownames(table_of_variants))
  colnames(names_table_of_variants) = "State"
  zygosity=gsub("-2","Trait-Reducing Homozygous variants",rownames(table_of_variants))
  zygosity=gsub("-1","Trait-Reducing Heterozygous variants",zygosity)
  zygosity=gsub("0","Neutral Homozygous variants",zygosity)
  zygosity=gsub("1","Trait-Increasing Heterozygous variants",zygosity)
  zygosity=gsub("2","Trait-Increasing Homozygous variants",zygosity)
  zygosity=data.frame(zygosity)
  table_of_variants=cbind(zygosity,names_table_of_variants,table_of_variants)
  ##### rename output table 
  table_of_variants_file=rename_output(vcfFile);
  write.table(table_of_variants, file = table_of_variants_file, sep = "\t", col.names = TRUE, row.names=FALSE,quote = FALSE)
  #colnames(merged_genotype_assoc)
  n3=ncol(merged_genotype_assoc);
  merged_genotype_info =  merged_genotype_assoc %>% select("SNP","CHR","POS","A1","A2","BETA","P")
  genotypes_states = merged_genotype_assoc[,(n3+1-2*numberofsamples):n3]
  genotypes_info_states = cbind(merged_genotype_info,genotypes_states)
  states_file=rename_output_2(vcfFile);
  write.table(genotypes_info_states, file = states_file, sep = "\t", col.names = TRUE, row.names=FALSE,quote = FALSE)
  
}

rename_output <- function(File){
  executiondate=format(Sys.time(), "%Y%m%d")
  sufix=gsub(' ','',paste(".Variants.table.By.user.tsv",collapser=""))
  Outputfilename=gsub(".vcf",sufix, File)
  return(Outputfilename)
}

rename_output_2 <- function(File){
  executiondate=format(Sys.time(), "%Y%m%d")
  sufix=gsub(' ','',paste(".Genotype.States.By.User.tsv",collapser=""))
  Outputfilename=gsub(".vcf",sufix, File)
  return(Outputfilename)
}

## MAIN
if (!interactive()){
  #message("\n\n",yellow$underline$bold("Código 46 S.A. de C.V."),green("Generate SNP list from Assoc file (SNP column)"),"\n")
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
