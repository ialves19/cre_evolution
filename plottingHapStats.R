library(ggplot2)
library(tidyverse)
library(wesanderson)
#############################
###
###
###
### FUNCTIONS
###
###
#############################
plottingHaplotypeSumStats <- function(listINPUT, nameCol, pName, yAxisLabel) {
  # global variables: namesStats, selCoeff, scenarios, majorScenarios, aseScenario
  # listINPUT = pi
  # nameCol = sumStats
  # pName = plotTitle
  # yAxisLabel = yLabel
  
  selCoeffLabels <- paste0("s=", selCoeff)
  aveNbMutPerHap_m <- matrix(as.numeric(unlist(lapply(listINPUT, "[", which(namesStats == nameCol)))), 
                             ncol=length(selCoeff))
  colnames(aveNbMutPerHap_m) <- names(listINPUT)
  df_haps <- data.frame(aveNbMutPerHap_m)
  df_haps <- df_haps %>% select(names(df_haps)) %>%
    pivot_longer(., cols = names(df_haps), names_to = "Scenario", values_to = "sumStat") %>% 
    mutate(selectionCoeff=rep(selCoeff, replicates), aseModel=rep(aseScenario,replicates), model=rep(majorScenarios, replicates)) 
  df_haps$sumStat[is.na(df_haps$sumStat)] <- 0
  
  
  p <- ggplot(df_haps, aes(x=model, y=sumStat)) + 
    labs(title = plotTitle, x = "Selection model", y = yLabel) + 
    geom_violin(aes(color=Scenario), trim=FALSE, position = position_dodge(0.9)) + 
    geom_boxplot(aes(fill=Scenario), width=0.5, position = position_dodge(0.9)) + 
    scale_x_discrete(labels=unique(selCoeffLabels)) + 
    scale_fill_manual(values=colorsClasses) + scale_color_manual(values=colorsClasses) +
    guides(color = "none") + theme_classic()
  
  # p-values
  pValuesList <- list()
  for(s in unique(selCoeff)) { 
    #print(s)
    pValuesList[[s]] <- wilcox.test(unlist(df_haps %>% filter(selectionCoeff == s & aseModel == "noASE") %>% select("sumStat")), 
                                    unlist(df_haps %>% filter(selectionCoeff == s & aseModel == "ASE") %>% select("sumStat")), alternative = "two")$p.value
    
  }
  
  pValues_v <- as.matrix(unlist(pValuesList), byrow = F, rownames=unique(selCoeff))
  write.table(pValues_v, file = paste0(targetModel, "_", nameCol, "_pValues.txt"),
              quote = F, row.names = T, col.names = F, sep = "\t")
  print(p)
  
}
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

#setting working directory 
wrkDir <- "/Users/isabel/Dropbox/UnivSTRASBOURG/PROJECTS/CRE_evolution/05-Results/sims_Feb2025"
setwd(wrkDir)

# creating a pdf folder in case it doesn't exist to save pdf files with results
if(!dir.exists(paste0(wrkDir, "/pdfs"))) {
  dir.create(paste0(wrkDir, "/pdfs"))
  outputPDF <- paste0(wrkDir, "/pdfs")
} else {
  outputPDF <- paste0(wrkDir, "/pdfs")
}

#####-----------------
##
## VARIABLES
##
#####-----------------
# nb of replicates / model
replicates <- 100

# setting target model 
targetModel <- "twoLoci_rec0.5kb_cre1kb"

# selection coefficients
selCoeff <- rep(c("-0.001", 
              "-0.01", 
              "-0.0425",
              "-0.1", 
              "-0.3"), each=2)
# scenarios
scenarios <- paste0(targetModel, "_s", paste0(selCoeff, c("", "_epist")))
majorScenarios <- paste0(targetModel, "_s=", selCoeff)
aseScenario <- rep(c("noASE","ASE"), length(selCoeff)/2)
colorsClasses <- wes_palette("Zissou1", length(majorScenarios), type = "continuous")
#####-----------------
#####-------- End of variables
#####----

# opening the pdf 
pdf(file=paste0(outputPDF, "/", targetModel, "_HapsPlots_violins.pdf"), height = 6, width = 6)

# reading summary stats SLiM
hapStats <- list()
pi <- list()

for (s in scenarios) {
  print(s);
  
  hapStats[[s]] <- read.table(paste0(s, "/sumStats_", s, "_haps.txt"), header = T, sep = "\t", 
                              na.strings = "NaN")
  
}

# processing hap stats 
listHaps <- list()
namesStats <- names(hapStats[[1]])
listHaps <- lapply(hapStats, function(x) {x %>% mutate(NbofMut3Haps=apply(x[,c(3,4)], 1, sum))})
listHaps <- lapply(listHaps, function(x) {x %>% mutate(relNbHaps_m3_m2=NbHaps_w_m3_w_m2/NbofMut3Haps, 
                                                       relNbHaps_NOm3_m2=NbHaps_wo_m3_w_m2/(2000-NbofMut3Haps))})
#######################
##
### Plotting & computing p-Values
##
#######################
#Deleterious mutations onto overexpressed haps
sumStats <- "avePerHap_nbMut2_ontoM3haps"
plotTitle <- "Ave Nb of deleterious mutations onto overexpressed haplotypes"
yLabel <- "Number of mutations"

plottingHaplotypeSumStats(hapStats,sumStats,plotTitle,yLabel)

#Deleterious mutations onto underexpressed haps
sumStats <- "avePerHap_nbMut2_ontoNOM3haps"
plotTitle <- "Ave Nb of deleterious mutations onto underexpressed haplotypes"
yLabel <- "Number of mutations"

plottingHaplotypeSumStats(hapStats,sumStats,plotTitle,yLabel)

#nb overexpressed haps with deleterious mutations
sumStats <- "NbHaps_w_m3_w_m2"
plotTitle <- "Nb of overexpressed haps with deleterious muts"
yLabel <- "Number of haplotypes"

plottingHaplotypeSumStats(hapStats,sumStats,plotTitle,yLabel)

#nb underexpressed haps with deleterious mutations
sumStats <- "NbHaps_wo_m3_w_m2"
plotTitle <- "Nb of underexpressed haps with deleterious muts"
yLabel <- "Number of haplotypes"

plottingHaplotypeSumStats(hapStats,sumStats,plotTitle,yLabel)

dev.off()
