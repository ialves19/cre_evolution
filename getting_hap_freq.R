###############################
##
## Compute haplotype diversity
##
###############################
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
#.libPaths("C:/software/Rpackages") #need to change this
library(ggplot2)
library(tidyverse)

wrkDir <- args[1]
sceNames <- args[2]
nbReplicates <- args[3]

wrkDir <- "/Users/ialves/Dropbox/SLIM/sims"
sceNames <- "twoLoci_rec_cre100kb_s0.3_epist"
nbReplicates <- 25
setwd(wrkDir)

#the idea here is to create a table with the counts of haplotypes 
# per replicate

haploDist_v <- c()
genNb <- 4999 # from 0 4999
cat(paste("ReplicateNb", "HaplotypeFreq", sep = " "), file=paste0(sceNames, "/", sceNames, "_g", genNb, ".haps"), sep = "\n")

for(r in 1:nbReplicates) {
  mut_tb <- read.table(paste0(sceNames, "/", sceNames, "_", r, "_g", genNb, ".mutations"),
                       skip = 2, sep = " ")
  gen_l <- list(500)
  for(genome in 0:499) {
    
    tmp_v <- as.numeric(scan(paste0(sceNames, "/", sceNames, "_", r, "_g", 
                                    genNb, ".genomes"), skip = 1+genome, nlines = 1, sep = " ",                                      what = character())[-c(1,2)])
    genomeAnc <- rep("A", nrow(mut_tb))
    genomeAnc[mut_tb$V1 %in% tmp_v] <- mut_tb$V3[mut_tb$V1 %in% tmp_v]
    gen_l[[genome+1]] <- genomeAnc
  }
  
  gen_m <- as.data.frame(do.call(rbind,gen_l))
  
  genNames <- paste0("G", 1:500)
  
  #checking variants present more than once in the mut table
  # it happens that in some cases there are chromosomal positions with >1 mutation
  # in this cases I am creating X.1, X.2 and so on to represent multiple mutations
  reNameMut <- table(mut_tb$V4)[table(mut_tb$V4) > 1]
  if(length(reNameMut) > 0) {
    for(i in 1:length(reNameMut)) {
      index_to_replace <- which(mut_tb$V4 == as.numeric(names(reNameMut[i])))
      mut_tb$V4[index_to_replace] <- paste0(names(reNameMut[i]), ".", 1:reNameMut[i])
    }
  }
  
  # because of the previous cases it can happen that the vector of chromosomal positions
  # is character instead of numeric. 
  gen_m <- gen_m[,order(as.numeric(mut_tb$V4), decreasing = F)]
  
  
  #varPos <- factor(sort(mut_tb$V4), levels = sort(mut_tb$V4))
  varPos <- sort(as.numeric(mut_tb$V4))
  colnames(gen_m) <- varPos
  rownames(gen_m) <- genNames
  countGenomes <- gen_m %>% group_by_all%>%count
  haploDist_v <- countGenomes$n
  cat(paste(r, paste(haploDist_v, collapse = " "), sep = " "), file=paste0(sceNames, "/", sceNames, "_g", genNb, ".haps"), sep = "\n", append = T)
  
  #computing heterozygosity 
  mut_tb <- mut_tb %>% mutate(heteroz=2*mut_tb$V9/500*(1-mut_tb$V9/500))
  heterozSum <- mut_tb %>% group_by(V3) %>% summarise(mean(heteroz), median(heteroz), quantile(heteroz, probs=c(0.025)), quantile(heteroz, probs=c(0.25)), 
                                         quantile(heteroz, probs=c(0.75), quantile(heteroz, probs=c(0.975)))) 
  heterozSum <- left_join(heterozSum, ) mut_tb %>% group_by(V3) %>% count %>% select(n)
  names(heterozSum) <- c("Mean", "Median", "Q2.5%", "Q25%", "Q75%", "Q95%")
}

hapFreq <- list()
for(r in 1:nbReplicates-1) {
  tmp <- as.numeric(scan(paste0(sceNames, "/", sceNames, "_g", genNb, ".haps"), 
                skip = 1+r, nlines = 1, sep = " ", what = numeric())[-c(1)])
  hapFreq[[r+1]] <- tmp
}

#Computing histogram on the haps 
hapHist <- lapply(1:length(hapFreq), function(x) {hist(hapFreq[[x]], plot = F)$counts})
#Computing the number of haps per replicate and save in a file
nbHaps <- unlist(lapply(hapHist, sum))
cat(paste(nbHaps, collapse = " "), file=paste0(sceNames, "/", sceNames, "_g", genNb, ".hapNb"), sep = "\n")

#Creating a matrix with the counts (25rep x haplotype frequency)
hapFreq_sameLength <- list()
for(r in 1:nbReplicates) {
  emptyV <- rep(0, max(lengths(hapHist)))
  emptyV[1:length(hapHist[[r]])] <- hapHist[[r]]
  hapFreq_sameLength[[r]] <- emptyV
}

##transforming into a df
hapFreq_df <- do.call(rbind, hapFreq_sameLength)

#summarising the df (mean, median, 2.5%, 95%, 1Q, 3Q)
meanV <- apply(hapFreq_df, 2,mean)
medianV <- apply(hapFreq_df, 2,median)
CI_025V <- apply(hapFreq_df, 2, function(x) {quantile(x, probs=0.025)})
CI_975V <- apply(hapFreq_df, 2, function(x) {quantile(x, probs=0.975)})
CI_25V <- apply(hapFreq_df, 2, function(x) {quantile(x, probs=0.25)})
CI_75V <- apply(hapFreq_df, 2, function(x) {quantile(x, probs=0.75)})

cat(paste0("Mean: ", paste(meanV, collapse = " ")), file=paste0(sceNames, "/", sceNames, "_g", genNb, ".hapSumStats"), sep = "\n")
cat(paste0("Median: ", paste(medianV, collapse = " ")), file=paste0(sceNames, "/", sceNames, "_g", genNb, ".hapSumStats"), sep = "\n", append = T)
cat(paste0("CI_025V: ", paste(CI_025V, collapse = " ")), file=paste0(sceNames, "/", sceNames, "_g", genNb, ".hapSumStats"), sep = "\n", append = T)
cat(paste0("CI_975V: ", paste(CI_975V, collapse = " ")), file=paste0(sceNames, "/", sceNames, "_g", genNb, ".hapSumStats"), sep = "\n", append = T)
cat(paste0("CI_25V: ", paste(CI_25V, collapse = " ")), file=paste0(sceNames, "/", sceNames, "_g", genNb, ".hapSumStats"), sep = "\n", append = T)
cat(paste0("CI_75V: ", paste(CI_75V, collapse = " ")), file=paste0(sceNames, "/", sceNames, "_g", genNb, ".hapSumStats"), sep = "\n", append = T)


