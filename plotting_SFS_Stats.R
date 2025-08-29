library(ggplot2)
library(tidyverse)
library(RColorBrewer)
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
  
  # this part performed a violin plot which is confusing to read
  # p <- ggplot(df_haps, aes(x=model, y=sumStat)) + 
  #   labs(title = plotTitle, x = "Selection model", y = yLabel) + 
  #   geom_violin(aes(color=Scenario), trim=FALSE, position = position_dodge(0.9)) + 
  #   geom_boxplot(aes(fill=Scenario), width=0.5, position = position_dodge(0.9)) + 
  #   scale_x_discrete(labels=unique(selCoeffLabels)) + 
  #   scale_fill_manual(values=colorsClasses) + scale_color_manual(values=colorsClasses) +
  #   guides(color = "none") + theme_classic()
  
  p <- ggplot(df_haps, aes(x=model, y=sumStat)) + 
    labs(title = plotTitle, x = "Selection model", y = yLabel) + 
    geom_boxplot(aes(color=Scenario)) + 
    scale_x_discrete(labels=unique(selCoeffLabels)) + 
    scale_color_manual(values=colorsClasses) + 
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
wrkDir <- "/Users/isabel/Dropbox/UnivSTRASBOURG/PROJECTS/CRE_evolution/05-Results/sims_antoine_all"
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
pdf(file=paste0(outputPDF, "/", targetModel, "_SFSstats_violins.pdf"), height = 6, width = 6)

# reading summary stats SLiM
perSiteStats <- list()
pi <- list()

for (s in scenarios) {
  print(s);
  
  perSiteStats[[s]] <- read.table(paste0(wrkDir, "/", s, "/sumStats_", s, "_perSite.txt"), header = T, sep = "\t", 
                                  na.strings = "NaN")
  pi[[s]] <- read.table(file=paste0(wrkDir, "/", s, "/pi_", s, ".txt"), header = T, sep = "\t", na.strings = NA)
}

# processing SFS stats 
namesStats <- names(perSiteStats[[1]])
namesPi <- names(pi[[1]])

#######################
##
### Plotting & computing p-Values
##
#######################
## plot Nb of deleterious mutations - coding sequence 
sumStats <- "Nb_m2"
plotTitle <- "Deleterious coding"
yLabel <- "Number of mutations"

plottingHaplotypeSumStats(perSiteStats, sumStats,plotTitle,yLabel)
#===================================
## plot nb of neutral coding
sumStats <- "Nb_m1"
plotTitle <- "Neutral coding"
yLabel <- "Number of mutations"
#===================================
plottingHaplotypeSumStats(perSiteStats, sumStats,plotTitle,yLabel)

## plot nb of neutral cis-regulatory
sumStats <- "Nb_m3"
plotTitle <- "Neutral regulatory"
yLabel <- "Number of mutations"

plottingHaplotypeSumStats(perSiteStats, sumStats,plotTitle,yLabel)
#===================================
## plot overall heterozygosity 
sumStats <- "AveHet"
plotTitle <- "Overall heterozygosity"
yLabel <- "Heterozygosity"

plottingHaplotypeSumStats(perSiteStats, sumStats,plotTitle,yLabel)
#===================================
## plot m1 heterozygosity 
sumStats <- "AveHet_m1"
plotTitle <- "Heterozygosity neutral coding"
yLabel <- "Heterozygosity"

plottingHaplotypeSumStats(perSiteStats, sumStats,plotTitle,yLabel)
#===================================
## plot m2 heterozygosity 
sumStats <- "AveHet_m2"
plotTitle <- "Heterozygosity deleterious coding"
yLabel <- "Heterozygosity"

plottingHaplotypeSumStats(perSiteStats, sumStats,plotTitle,yLabel)
#===================================
## plot m3 heterozygosity 
sumStats <- "AveHet_m3"
plotTitle <- "Heterozygosity neutral regulatory"
yLabel <- "Heterozygosity"

plottingHaplotypeSumStats(perSiteStats, sumStats,plotTitle,yLabel)
#===================================
###------------------------
###           Pi 
###------------------------
# plot overall pi
namesStats <- names(pi[[1]])
 
sumStats <- "Overall_Nucleotide_diversity"
plotTitle <- expression("Overall nucleotide diversity ("*pi*")")
yLabel <- expression("Nucleotide diversity ("*pi*")")

plottingHaplotypeSumStats(pi, sumStats, plotTitle, yLabel)
#===================================
## Coding region pi
sumStats <- "CS_Nucleotide_diversity"
plotTitle <- expression("Nucleotide diversity ("*pi*") - Coding seq.")
yLabel <- expression("Nucleotide diversity ("*pi*")")

plottingHaplotypeSumStats(pi, sumStats, plotTitle, yLabel)
#===================================
## Regulatory region pi
sumStats <- "CRE_Nucleotide_diversity"
plotTitle <- expression("Nucleotide diversity ("*pi*") - CRE")
yLabel <- expression("Nucleotide diversity ("*pi*")")

plottingHaplotypeSumStats(pi, sumStats, plotTitle, yLabel)
#===================================
## Overall_SegSites
sumStats <- "Overall_SegSites"
plotTitle <- "Overall segrating sites (S)"
yLabel <- "Number of mutations"

plottingHaplotypeSumStats(pi, sumStats, plotTitle, yLabel)
#===================================
## Overall_theta S
sumStats <- "Overall_ThetaS"
plotTitle <- expression("Overall "*theta*"(S)")
yLabel <- expression(theta*"(S)")

plottingHaplotypeSumStats(pi, sumStats, plotTitle, yLabel)
#===================================
## CS_SegSites 
sumStats <- "CS_SegSites"
plotTitle <- "Segrating sites (S) - Coding seq."
yLabel <- "Number of mutations"

plottingHaplotypeSumStats(pi, sumStats, plotTitle, yLabel)
#===================================
## CS_ThetaS
sumStats <- "CS_ThetaS"
plotTitle <-  expression(theta*"(S) - Coding Seq.")
yLabel <- expression(theta*"(S)")

plottingHaplotypeSumStats(pi, sumStats, plotTitle, yLabel)
#===================================
## CRE_SegSites
sumStats <- "CRE_SegSites"
plotTitle <-  "Segrating sites (S) - CRE"
yLabel <- "Number of mutations"

plottingHaplotypeSumStats(pi, sumStats, plotTitle, yLabel)
#===================================
## CRE_ThetaS
sumStats <- "CRE_ThetaS"
plotTitle <-  expression(theta*"(S) - Coding Seq.")
yLabel <- expression(theta*"(S)")

plottingHaplotypeSumStats(pi, sumStats, plotTitle, yLabel)
#===================================
## OverallTajimasD
sumStats <- "OverallTajimasD"
plotTitle <-  "Overall Tajima's D"
yLabel <- "Tajima's D"

plottingHaplotypeSumStats(pi, sumStats, plotTitle, yLabel)
#===================================
## TajimasD - CS
sumStats <- "TajimasDCS"
plotTitle <-  "Tajima's D - Coding seq."
yLabel <- "Tajima's D"

plottingHaplotypeSumStats(pi, sumStats, plotTitle, yLabel)
#===================================
##TajimasD - CRE
sumStats <- "TajimasDCRE"
plotTitle <-  "Tajima's D - CRE"
yLabel <- "Tajima's D"

plottingHaplotypeSumStats(pi, sumStats, plotTitle, yLabel)
#===================================
dev.off()