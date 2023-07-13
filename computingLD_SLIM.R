#!/usr/bin/env Rscript
#############################
##
## This script takes the output of SLiM and builds a 
## genome image of the mutations: varPos on the columns, genomes on the rows
## it takes as input *.segMutations and *.genomes 
##
##
## by Isabel Alves - Fev 2023
##
#############################


args = commandArgs(trailingOnly=TRUE)
library(ggplot2, lib.loc = "/home/ialves/R/x86_64-pc-linux-gnu-library/4.1")
library(tidyverse, lib.loc = "/home/ialves/R/x86_64-pc-linux-gnu-library/4.1")
source(../sumStats_functions.R)

wrkDir <- args[1]
sceName <- args[2]
setwd(paste0(wrkDir, "/", sceName))

# examples : c("twoLoci_rec_cre100kb_s0.001","twoLoci_rec_cre100kb_s0.001_epist", "twoLoci_rec_cre100kb_s0.0425", "twoLoci_rec_cre100kb_s0.0425_epist","twoLoci_rec_cre100kb_s0.1", "twoLoci_rec_cre100kb_s0.1_epist", "twoLoci_rec_cre100kb_s0.3", "twoLoci_rec_cre100kb_s0.3_epist")
#scenario name

genNb <- 10000 # from 0 4999
nbOfSims <- 100
popSize <- 1000 

sumStats_LD_fName <- paste0("sumStats__", sceName, "_LD.txt")
cat(paste("CREPosition", "CSPosition", "D", "Rsquared", "Pvalue-ChiSqTest", sep = "\t"), 
    sep = "\n", file = sumStats_LD_fName)

for(simID in 1:nbOfSims) {
  
  #simID <- 10
  #open the mutation section of SLIM's output
  #twoLoci_rec5kb_cre1kb_s-0.001_epist_54_g10000.segMutations 
  #setwd(paste0(wrkDir, "/sims/", sceName))
  #setwd("/Users/ialves/Dropbox/SLIM/sims/twoLoci_rec_cre1kb_s0.3/")
  #simID <- 98
  segSitesTable <- openSegSites(sceName, simID, genNb)
  #dim(segSitesTable)
  #checking whether there're variants present more than once in the mut table
  # These sites are excluded from all the following analysis 
  cleanSegSitesTable <- detectMultiallelicSites(segSitesTable, TRUE)
  #dim(cleanSegSitesTable)
  nbOfMultiallelicSites <- nrow(segSitesTable)-nrow(cleanSegSitesTable)
  rm(segSitesTable)
  
  # reading the genomes 
  genMatrix <- getGenomesFromSLIMoutput(sceName, simID, genNb, popSize)
  #dim(genMatrix)
  df_genomes <- data.frame(genMatrix)
  
  nbofM3Mut <- nrow(cleanSegSitesTable %>% filter(V3 == "m3"))
  
  # computing linkage disequilibrium 
  if(nbofM3Mut > 0) {
    
    CREsites <- MutPositions("m3")
    #freqCREsites <- alleleCounts("m3")
    
    CS_delSites <- MutPositions("m2")
    #freqCS_delSites <- alleleCounts("m2")
    
    sapply(CREsites, computing_LD)
    
  } else {
    cat("No variable sites in the CRE. NO LD computation was done.")
  }
  
}