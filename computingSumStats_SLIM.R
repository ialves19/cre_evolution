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

# arguments to pass - TEST 
# sceName <- "twoLoci_rec_cre1kb_s0.3"
# genNb <- 10000 # from 0 4999
# nbOfSims <- 100
# wrkDir <- "/home/ialves/SLiM"
# popSize <- 1000 
# simID <- 98

#test pi ans S computation
# set.seed(189989)
# seqTest <- sample(c("A","T", "G", "C"), 50, replace = T)
# df_seq <- data.frame(rbind(seqTest,seqTest,seqTest,seqTest,seqTest,seqTest,seqTest,
#                            seqTest,seqTest,seqTest))
# df_seq[c(2,5,6,9),4] <- "G"
# df_seq[c(7),17] <- "T"
# df_seq[c(2,4,5,6,7,8,9,10),44] <- "A"
# df_seq[c(4,8),50] <- "C"
# 
# test_uniqHaps <- df_seq %>% group_by_all %>% count
# test_HapFreq <- test_uniqHaps$n/sum(test_uniqHaps$n)
# pitest <- computing_Pi(test_uniqHaps[,-ncol(test_uniqHaps)], test_HapFreq, 10)
# segTest <- 4/sum(1/rep(1:(10-1)))


#creating output files with headers
sumStats_perSite_fName <- paste0("sumStats_", sceName, "_perSite.txt")
sumStats_pi_fName <- paste0("pi_", sceName, ".txt")
sumStats_haps_fName <- paste0("sumStats_", sceName, "_haps.txt")

cat(paste("Nb_m1", "Nb_m2", "Nb_m3", "AveHet", "AveHet_m1", "AveHet_m2", "AveHet_m3", sep = "\t"), sep = "\n", file = sumStats_perSite_fName)
cat(paste("Overall_Nucleotide_diversity", "CS_Nucleotide_diversity", "CRE_Nucleotide_diversity", 
          "Overall_SegSites", "Overall_ThetaS", "CS_SegSites", "CS_ThetaS", 
          "CRE_SegSites", "CRE_ThetaS", "OverallTajimasD", "TajimasDCS", "TajimasDCRE", sep = "\t"), sep = "\n", file = sumStats_pi_fName)
cat(paste("avePerHap_nbMut2_ontoM3haps", "avePerHap_nbMut2_ontoNOM3haps", 
          "NbHaps_w_m3_w_m2", "NbHaps_w_m3_wo_m2", "NbHaps_wo_m3_w_m2", "NbHaps_wo_m3_wo_m2", "Odds_m3_m2","p-value_m3_m2",
          "avePerHap_nbMut1_ontoM3haps", "avePerHap_nbMut1_ontoNOM3haps",
          "NbHaps_w_m3_w_m1", "NbHaps_w_m3_wo_m1", "NbHaps_wo_m3_w_m1", "NbHaps_wo_m3_wo_m1", "Odds_m3_m1","p-value_m3_m1", sep = "\t"), 
    sep = "\n", file = sumStats_haps_fName)


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
  
  #heterozygosity 
  hetTable <- computingHeterozygosity(cleanSegSitesTable, popSize)
  mut_tblName <- paste0(sceName, "_mutProperties_", simID, ".txt")
  write.table(hetTable, file=mut_tblName, quote = F, row.names = F, col.names = T, sep = "\t")
  #cat(paste("Nb_m1", "Nb_m2", "Nb_m3", "AveHet", "AveHet_m1", "AveHet_m2", "AveHet_m3", sep = "\t"), sep = "\n", file = sumStats_perSite_fName)
  cat(paste(sum(hetTable$muType == "m1"), sum(hetTable$muType == "m2"), sum(hetTable$muType == "m3"), mean(hetTable$Heterozyg),
            mean(hetTable$Heterozyg[hetTable$muType == "m1"]),
            mean(hetTable$Heterozyg[hetTable$muType == "m2"]),
            mean(hetTable$Heterozyg[hetTable$muType == "m3"]), sep = "\t"), sep = "\n", file = sumStats_perSite_fName, append = T)
  # reading the genomes 
  genMatrix <- getGenomesFromSLIMoutput(sceName, simID, genNb, popSize)
  #dim(genMatrix)
  df_genomes <- data.frame(genMatrix)
  #retireving unique haplotypes 
  uniq_haps <- df_genomes %>% group_by_all %>% count
  #computing haplotype frequencies
  haps_freq <- uniq_haps$n/sum(uniq_haps$n)
  #computing pi
  nucleotideDiv <- computing_Pi(uniq_haps[,-ncol(uniq_haps)], haps_freq, popSize)
  # computing Theta S
  thetaS <- nrow(cleanSegSitesTable)/sum(1/rep(1:((popSize*2)-1)))
  
  ### Within CRE and CS 
  #computing pi 
  nbofM3Mut <- nrow(cleanSegSitesTable %>% filter(V3 == "m3"))
  if(nbofM3Mut > 0 ) {
    
    CREsites <- MutPositions("m3")
    CREgenoMatrix <- df_genomes %>% select(colnames(df_genomes)[colnames(df_genomes) %in% paste0("X",CREsites)])
    CSgenoMatrix <- df_genomes %>% select(colnames(df_genomes)[!colnames(df_genomes) %in% paste0("X",CREsites)])
    uniqHapsCRE <- CREgenoMatrix %>% group_by_all %>% count
    uniqHapsCS <- CSgenoMatrix %>% group_by_all %>% count
    hapsFreqCRE <- uniqHapsCRE$n/sum(uniqHapsCRE$n)
    hapsFreqCS <- uniqHapsCS$n/sum(uniqHapsCS$n)
    
    piCRE <- computing_Pi(uniqHapsCRE[,-ncol(uniqHapsCRE)], hapsFreqCRE, popSize)
    piCS <- computing_Pi(uniqHapsCS[,-ncol(uniqHapsCS)], hapsFreqCS, popSize)
    
  } else {
    CSgenoMatrix <- df_genomes 
    uniqHapsCS <- CSgenoMatrix %>% group_by_all %>% count
    hapsFreqCS <- uniqHapsCS$n/sum(uniqHapsCS$n)
    # popSize here is in nb of diploid individuals
    piCS <- computing_Pi(uniqHapsCS[,-ncol(uniqHapsCS)], hapsFreqCS, popSize)
    piCRE <- 0
    
  }

  #computing S 
  CRESegSites <- nrow(cleanSegSitesTable %>% filter(V3 == "m3"))
  CREThetaSegSites <- CRESegSites/sum(1/rep(1:((popSize*2)-1)))
  CSSegSites <- nrow(cleanSegSitesTable %>% filter(V3 == "m1" | V3 == "m2"))
  CSThetaSegSites <- CSSegSites/sum(1/rep(1:((popSize*2)-1)))
  
  #computing Tajimas'D 
  tajDCRE <- piCRE - CREThetaSegSites
  tajDCS <- piCS - CSThetaSegSites
  ### ------------
  overTajimasD <- nucleotideDiv - thetaS
  cat(paste(nucleotideDiv, piCS, piCRE, nrow(cleanSegSitesTable), thetaS, CSSegSites, CSThetaSegSites,
            CRESegSites, CREThetaSegSites, overTajimasD, tajDCS, tajDCRE, sep = "\t"), sep = "\n", file = sumStats_pi_fName, append = T)
  
  #computing the number of mut of type X per haplotype
  nb_m3_per_hap <- computing_nb_mutX_per_haplotype(df_genomes, "m3")
  nb_m2_per_hap <- computing_nb_mutX_per_haplotype(df_genomes, "m2")
  nb_m1_per_hap <- computing_nb_mutX_per_haplotype(df_genomes, "m1")
  
  #the vector below contains the ave nb of m2 (or m1) mut (>0) within haps containing or not m3 mut.
  m2_onm3 <- nb_mX_onto_m3_background(df_genomes, nb_m2_per_hap, nb_m3_per_hap, "m2")
  names(m2_onm3) <- c("mean_nb_mut2_ontoM3haps", "mean_nb_mut2_ontoNoM3haps", "nb_hap_w_m3_w_m2",
                      "nb_hap_w_m3_wo_m2","nb_hap_wo_m3_w_m2","nb_hap_wo_m3_wo_m2")
  m1_onm3 <- nb_mX_onto_m3_background(df_genomes, nb_m1_per_hap, nb_m3_per_hap, "m1")
  names(m1_onm3) <- c("mean_nb_mut1_ontoM3haps", "mean_nb_mut1_ontoNoM3haps", "nb_hap_w_m3_w_m1",
                      "nb_hap_w_m3_wo_m1","nb_hap_wo_m3_w_m1","nb_hap_wo_m3_wo_m1")
  
  #are there more haplotypes with m3 and m2 than m3 and no m2 ? 
  contingencyTable_m2_m3 <- matrix(m2_onm3[-c(1:2)], ncol = 2)
  # fisher.test(contingencyTable_m2_m3)$p.value
  # fisher.test(contingencyTable_m2_m3)$odd
  contingencyTable_m1_m3 <- matrix(m1_onm3[-c(1:2)], ncol = 2)
  cat(paste(paste(m2_onm3, collapse = "\t"), fisher.test(contingencyTable_m2_m3)$estimate,fisher.test(contingencyTable_m2_m3)$p.value,
            paste(m1_onm3, collapse = "\t"), fisher.test(contingencyTable_m1_m3)$estimate,fisher.test(contingencyTable_m1_m3)$p.value,
            sep = "\t"), sep = "\n", file = sumStats_haps_fName, append = T)

} 
  

  
  #####
  ##### to be finished - fev 7, 2023
  #####
  #m3 SFS entries - do a table of the variables below
  # hap_counts_per_m2 <- apply(df_genomes, 2, function(x) {sum(x=="m2")})
  # hap_counts_per_m2_condm3 <- apply(df_genomes[which(nb_haps_w_m3>0),], 2, function(x) {sum(x=="m2")})
  # plot(table(hap_counts_per_m2)/sum(table(hap_counts_per_m2)), type = "b")
  # lines(table(hap_counts_per_m2_condm3)/sum(table(hap_counts_per_m2_condm3)), type = "b")
  # 
  # hap_counts_m3 <- apply(df_genomes, 2, function(x) {sum(x=="m3")})
  # hap_counts_m1 <- apply(df_genomes, 2, function(x) {sum(x=="m1")})
  