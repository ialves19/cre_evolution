#############################
##
## Functions
##
#############################
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
##===========================
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
##===========================
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
##===========================
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
##===========================
########################--------------------
#########
### End of functions


#############################
###
###
###       MAIN
###
###
#############################
library(ggplot2)
library(tidyverse)

wrkDir <- "/shared/home/ialves/CRE_evolution/SLiM/sims_Jun2025"
sceName <- "twoLoci_rec0.5kb_cre1kb_s-0.01"
setwd(paste0(wrkDir, "/", sceName))


genNb <- 10000 # from 0 4999
nbOfSims <- 100
popSize <- 1000 

dafListM1 <- list()
dafListM2 <- list()
dafListM3 <- list()
hapsFreqList <-  list()

for(simID in 1:nbOfSims) {
  
  #simID <- 1
  
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
  # reading the genomes 
  genMatrix <- getGenomesFromSLIMoutput(sceName, simID, genNb, popSize)
  #dim(genMatrix)
  df_genomes <- data.frame(genMatrix)
  #retireving unique haplotypes 
  uniq_haps <- df_genomes %>% group_by_all %>% count
  #computing haplotype frequencies
  haps_freq <- uniq_haps$n/sum(uniq_haps$n)
  
  hetTableM1 <- hetTable %>% filter(muType == "m1")
  hetTableM2 <- hetTable %>% filter(muType == "m2")
  hetTableM3 <- hetTable %>% filter(muType == "m3")
  
  # DAFs
  dafListM1[[simID]] <- hist(hetTableM1$absFreq/2000, breaks=seq(0,1, by=0.05))$counts
  dafListM2[[simID]] <- hist(hetTableM2$absFreq/2000, breaks=seq(0,1, by=0.05))$counts
  dafListM3[[simID]] <- hist(hetTableM3$absFreq/2000, breaks=seq(0,1, by=0.05))$counts

  midPoints <- hist(hetTableM1$absFreq/2000, breaks=seq(0,1, by=0.05))$mids
  
  # HAP freqs.
  hapsFreqList[[simID]] <- hist(haps_freq, breaks=seq(0,0.5, by=0.02))$counts
  names(hapsFreqList[[simID]]) <- hist(haps_freq, breaks=seq(0,0.5, by=0.02))$mids
}

# from list to df 
df_DAF_m1 <- data.frame(do.call(rbind, dafListM1))
df_DAF_m2 <- data.frame(do.call(rbind, dafListM2))
df_DAF_m3 <- data.frame(do.call(rbind, dafListM3))

# summing up over all sims to have enough numbers
# as if we would have simulated 100 genomic regions
DAFallSimM1 <- data.frame(apply(df_DAF_m1, 2, sum), row.names = midPoints)
DAFallSimM2 <- data.frame(apply(df_DAF_m2, 2, sum), row.names = midPoints)
DAFallSimM3 <- data.frame(apply(df_DAF_m3, 2, sum), row.names = midPoints)

# saving DAF
write.table(DAFallSimM1, file=paste0(wrkDir, "/", sceName, "/", sceName, "_DAF_m1.daf"), col.names = F, row.names = T, quote = F)
write.table(DAFallSimM2, file=paste0(wrkDir, "/", sceName, "/", sceName, "_DAF_m2.daf"), col.names = F, row.names = T, quote = F)
write.table(DAFallSimM3, file=paste0(wrkDir, "/", sceName, "/", sceName, "_DAF_m3.daf"), col.names = F, row.names = T, quote = F)

#---------------
#
# Haplotype freq. 
#
#---------------
hapFreq_df <- data.frame(do.call(cbind, hapsFreqList))
colnames(hapFreq_df) <- paste0("Sim",1:nbOfSims)

write.table(hapFreq_df, file=paste0(wrkDir, "/", sceName, "/", sceName, "_HAPS.freq"), col.names = F, row.names = T, quote = F)


