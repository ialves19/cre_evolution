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
  pi <- ((ps*2)/((ps*2)-1))*2*sum(apply(pi_tmp, 1, sum, na.rm=T))
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

# get allele counts from the mutation table 
alleleCounts <- function(mutType) {
  
  x <- as.numeric(unlist(cleanSegSitesTable %>% filter(V3 == mutType)  %>% select(V9)))
  
  return(x)
}
# get allele physical positions from the mutation table 
MutPositions <- function(mutType) {
  
  x <- as.numeric(unlist(cleanSegSitesTable %>% filter(V3 == mutType)  %>% select(V4)))
  
  return(x)
}

# Compute LD statistics for each pairwise combination of CREalleles and CS (mut2) alleles
computing_LD <- function(PhyPosCREallele) {
  
  #PhyPosCREallele <- 100552
  for(i in CS_delSites) {
    
    #i <- 47226
    df.tmp <- df_genomes %>% select(colnames(df_genomes)[colnames(df_genomes) %in% c(paste0("X",PhyPosCREallele), paste0("X",i))])
    countHapsLD <- df.tmp %>% group_by_all %>% count
    
    mObsFreq <- matrix(rep(0,4), ncol = 2)
    colnames(mObsFreq) <- c("A", "m3")
    rownames(mObsFreq) <- c("A", "m2")
    mObsFreq[1,1] <- as.numeric(countHapsLD[which(countHapsLD[,1] == "A" & countHapsLD[,2] == "A"),3]/(popSize*2))
    mObsFreq[1,2] <- as.numeric(countHapsLD[which(countHapsLD[,1] == "A" & countHapsLD[,2] == "m3"),3]/(popSize*2))
    mObsFreq[2,1] <- as.numeric(countHapsLD[which(countHapsLD[,1] == "m2" & countHapsLD[,2] == "A"),3]/(popSize*2))
    mObsFreq[2,2] <- as.numeric(countHapsLD[which(countHapsLD[,1] == "m2" & countHapsLD[,2] == "m3"),3]/(popSize*2))
    
    mExpFreqIndep <- matrix(rep(0,4), ncol = 2)
    colnames(mExpFreqIndep) <- c("A", "m3")
    rownames(mExpFreqIndep) <- c("A", "m2")
    freqM3 <- as.numeric(cleanSegSitesTable %>% filter(V4 == PhyPosCREallele) %>% select(V9))/(popSize*2)
    freqM2 <- as.numeric(cleanSegSitesTable %>% filter(V4 == i) %>% select(V9))/(popSize*2)
    
    mExpFreqIndep[1,1] <- (1-freqM2)*(1-freqM3)
    mExpFreqIndep[1,2] <- (1-freqM2)*(freqM3)
    mExpFreqIndep[2,1] <- (freqM2)*(1-freqM3)
    mExpFreqIndep[2,2] <- (freqM2)*(freqM3)
    
    #mExpFreqIndep
    mObsFreq[is.na(mObsFreq)] <- 0
    XiSqPvalue <- chisq.test(mObsFreq*(popSize*2), mExpFreqIndep*(popSize*2))$p.value
    D <- mObsFreq[1,1]*mObsFreq[2,2]-mObsFreq[1,2]*mObsFreq[2,1]
    rSq <- D^2/((1-freqM2)*freqM2*(1-freqM3)*freqM3)
    cat(paste(PhyPosCREallele, i,D, rSq, XiSqPvalue, sep = "\t"), sep = "\n", file = sumStats_LD_fName, append = T)
    
  }
}
#############################
################-- End of functions
######
