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

#############################
##
## Functions
##
#############################

#open the mutation table 
openSegSites <- function(k, y, z) {
  # k = sceName
  # y = simID (replicate)
  # z = genNb
  x <- file(paste0(k, "_", y, "_g", z, ".segMutations"), "r")
  line <- readLines(x, 1)
  count <- 1
  while(!grepl("^Individuals:+", line)) {
    line <- readLines(x, 1)
    count <- count + 1
  }
  #! note that skip = 5 requires knowing that there are originally 5 lines of header in this file
  mut_tb <- read.table(paste0(k, "_", y, "_g", z, ".segMutations"),
                       skip = 5, nrows = count-6 , sep = " ")
  return(mut_tb)
}
# handle multiallelic sites
detectMultiallelicSites <- function(segM, clean) {
  # segM = segSitesTable
  #clean = T/F; if T then multiallelic sites will be removed
  #otherwise renamed by giving the physical position followed by X.1, X.2 etc
  reNameMut <- table(segM$V4)[table(segM$V4) > 1]
  if(length(reNameMut) > 0) {
    if(clean) {
      index_to_replace <- which(segM$V4 %in% as.numeric(names(reNameMut)))
      segM <- segM[-index_to_replace,]
    } else {
      for(i in 1:length(reNameMut)) {
        index_to_replace <- which(segM$V4 == as.numeric(names(reNameMut[i])))
        segM$V4[index_to_replace] <- paste0(names(reNameMut[i]), ".", 1:reNameMut[i])
      } #end of for along the multiallelic sites
    } #end of the if clean T/F
  } #end of the if there are multiallelic sites
  return(segM)
} 

# computing heterozygosity and adding it to hte clean table
#computing heterozygosity  
computingHeterozygosity <- function(segTbl, ps) {
  # segTbl = cleanSegSitesTable
  # ps = popSize
  
  p <- segTbl$V9/(ps*2) 
  q <- 1-p
  het <- 2*p*q
  newTbl <- segTbl %>% mutate(Heteroz = het)
  names(newTbl) <- c("SiteIDI", "SiteIDII", "muType", "Pos", "selCoeff", "domCoeff", "subPop",
                     "genOrigin", "absFreq", "Heterozyg")
  return(newTbl)
}
###------------

# read the genome outputed by slim
getGenomesFromSLIMoutput <- function(k, y, z, ps) {
  # k = sceName
  # y = simID (replicate)
  # z = genNb
  # ps = popSize
  
  gen_l <- list(ps*2)
  for(genome in 0:(ps*2-1)) {
    
    #genome <- 0
    tmp_v <- as.numeric(scan(paste0(k, "_", y, "_g", z, ".genomes"), skip = 1+genome, nlines = 1, sep = " ",                                      what = character())[-c(1,2)])
    cleanTmp_v <- tmp_v[tmp_v %in% hetTable$SiteIDI]
    
    genomeAnc <- rep("A", nrow(hetTable))
    genomeAnc[hetTable$SiteIDI %in% cleanTmp_v] <- hetTable$muType[hetTable$SiteIDI %in% cleanTmp_v]
    gen_l[[genome+1]] <- genomeAnc
    rm(tmp_v, cleanTmp_v)
  }
  
  # tranform list to matrix
  gen_m <- do.call(rbind,gen_l)
  gen_m <- gen_m[,order(hetTable$Pos, decreasing = F)]
  genNames <- paste0("G", 1:(ps*2))
  colnames(gen_m) <- as.character(hetTable$Pos[order(hetTable$Pos, decreasing = F)])
  rownames(gen_m) <- genNames
  return(gen_m)
}  
##--------------------

## Computing nucleotide diversity 
####
##
computing_Pi <- function(haps_m, haps_frq_v, ps) {
  
  # haps_m <- uniq_haps[,-ncol(uniq_haps)]
  # haps_frq_v <- haps_freq
  # ps <- popSize
  
  dissimilarity.matrix <- apply(t(haps_m),2,function(x)colSums(x!=t(haps_m)))
  diag(dissimilarity.matrix)<-0
  
  pi_tmp <- matrix(, ncol = ncol(dissimilarity.matrix), nrow = nrow(dissimilarity.matrix))
  for(h1 in 1:(nrow(haps_m)-1)) {
    for(h2 in (h1+1):nrow(haps_m)) {
      pi_tmp[h1, h2] <- haps_frq_v[h1]*haps_frq_v[h2]*dissimilarity.matrix[h1,h2]
    }
  }
  pi <- (ps/(ps-1))*2*sum(apply(pi_tmp, 1, sum, na.rm=T))
  return(pi)
}  
##--------------------
#computing the number of mut of type X per haplotype 
computing_nb_mutX_per_haplotype <- function(genomes,muType) {
  
  nb_haps_w_mutX <- apply(genomes, 1, function(x){ sum(x == muType)})
  return(nb_haps_w_mutX)
}

#test whether the average nb of m2 mut differs in haplotypes w or wo m3 mut
nb_mX_onto_m3_background <- function(genomes, haps_w_mX, haps_w_m3, muType) {
  
  index_m3 <- which(haps_w_m3 > 0)
  index_no_m3 <- which(haps_w_m3 == 0)
  
  index_m2 <- which(haps_w_mX > 0)
  index_no_m2 <- which(haps_w_mX == 0)
  
  #intersect haps with m3 and m2
  nb_haps_w_m3_w_m2 <- intersect(index_m3, index_m2)
  nb_haps_w_m3_wo_m2 <- intersect(index_m3, index_no_m2)
  nb_haps_wo_m3_w_m2 <- intersect(index_no_m3, index_m2)
  nb_haps_wo_m3_wo_m2 <- intersect(index_no_m3, index_no_m2)
  
  mean_nb_m2_within_m3_haps <- mean(computing_nb_mutX_per_haplotype(genomes[nb_haps_w_m3_w_m2,], muType))
  mean_nb_m2_within_nom3_haps <- mean(computing_nb_mutX_per_haplotype(genomes[nb_haps_wo_m3_w_m2,], muType))
  return(c(mean_nb_m2_within_m3_haps, mean_nb_m2_within_nom3_haps, 
           length(nb_haps_w_m3_w_m2), length(nb_haps_w_m3_wo_m2), length(nb_haps_wo_m3_w_m2), 
           length(nb_haps_wo_m3_wo_m2)))
}

#############################
################-- End of functions
######

args = commandArgs(trailingOnly=TRUE)
library(ggplot2, lib.loc = "/home/ialves/R/x86_64-pc-linux-gnu-library/4.1")
library(tidyverse, lib.loc = "/home/ialves/R/x86_64-pc-linux-gnu-library/4.1")


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

#creating output files with headers
sumStats_perSite_fName <- paste0("sumStats_", sceName, "_perSite.txt")
sumStats_pi_fName <- paste0("pi_", sceName, ".txt")
sumStats_haps_fName <- paste0("sumStats_", sceName, "_haps.txt")

cat(paste("Nb_m1", "Nb_m2", "Nb_m3", "AveHet", "AveHet_m1", "AveHet_m2", "AveHet_m3", sep = "\t"), sep = "\n", file = sumStats_perSite_fName)
cat("Nucleotide_diversity", sep = "\n", file = sumStats_pi_fName)
cat(paste("avePerHap_nbMut2_ontoM3haps", "avePerHap_nbMut2_ontoNOM3haps", 
          "NbHaps_w_m3_w_m2", "NbHaps_w_m3_wo_m2", "NbHaps_wo_m3_w_m2", "NbHaps_wo_m3_wo_m2", "Odds_m3_m2","p-value_m3_m2",
          "avePerHap_nbMut1_ontoM3haps", "avePerHap_nbMut1_ontoNOM3haps",
          "NbHaps_w_m3_w_m1", "NbHaps_w_m3_wo_m1", "NbHaps_wo_m3_w_m1", "NbHaps_wo_m3_wo_m1", "Odds_m3_m1","p-value_m3_m1", sep = "\t"), 
    sep = "\n", file = sumStats_haps_fName)

#mean_nb_mut2_ontoM3haps", "mean_nb_mut2_ontoNoM3haps", "nb_hap_w_m3_w_m2",
#"nb_hap_w_m3_wo_m2","nb_hap_wo_m3_w_m2","nb_hap_wo_m3_wo_m2"
for(simID in 1:nbOfSims) {

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
  cat(nucleotideDiv, sep = "\n", file = sumStats_pi_fName, append = T)
  
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
  