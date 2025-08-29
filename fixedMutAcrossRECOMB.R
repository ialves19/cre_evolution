library(ggplot2)
library(tidyverse)
library(ggpubr)
library(RColorBrewer)
library(wesanderson)
library(stringr)

wrkDir <- "/Users/isabel/Dropbox/UnivSTRASBOURG/PROJECTS/CRE_evolution/05-Results/sims_Jun2025_Antoine/antoineHoffner"
setwd(wrkDir)

#### added by july 2025
# setting target model 
suffixModel <- "twoLoci_rec"
creSize <- 50000/1000
creLabel <- paste0("kb_cre", creSize,"kb_s")
creTitle <- paste0("cre - ", creSize, "kb")


# recombination rates
phyDist <- as.numeric(c(0.5, 1, 5, 10, 500))
recRate <- 1.25e-7*phyDist*1000
names(recRate) <- paste0(as.character(phyDist), "kb")
phyDistLabel <- rep(phyDist, each=2)
# selection coefficients
selCoeff <- rep(c("-0.0425"), each=2)
# scenarios
scenarios <- paste0(suffixModel, phyDistLabel, creLabel, paste0(selCoeff, c("", "_epist")))
majorScenarios <- paste0(targetModel, "_s=", selCoeff)
aseScenario <- rep(c("noASE","ASE"), length(selCoeff)/2)
colorsClasses <- wes_palette("Zissou1", length(majorScenarios), type = "continuous")
####



# models <- c("twoLoci_rec0.5kb_cre100kb", "twoLoci_rec5kb_cre100kb", "twoLoci_rec10kb_cre100kb",
#             "twoLoci_rec50kb_cre100kb", "twoLoci_rec100kb_cre100kb", "twoLoci_rec1000kb_cre100kb")

# colorsRecVar <- brewer.pal(n = 8, name = "GnBu")

#scenarios <- paste0(targetModel, c("_s0.001","_s0.001_epist", "_s0.0425", "_s0.0425_epist","_s0.1", "_s0.1_epist", "_s0.3", "_s0.3_epist"))

# selcoeff <- -0.0425
# scenariosV1 <- c(paste0("_s", selcoeff*-1), paste0("_s", selcoeff*-1,"_epist"))
# scenariosV2 <-c(paste0("_s", selcoeff), paste0("_s", selcoeff,"_epist"))
# ASEscenario <- c("noASE", "ASE")

fullTableM3Het <- data.frame()
fullTableM2Het <- data.frame()
fullTableM2fixed <- data.frame()
fullTableM3fixed <- data.frame()

for(targetModel in scenarios) {
      
    hasASE <- unlist(gregexpr(pattern ='epist', targetModel))
    whatRecomb <- unlist(strsplit(unlist(strsplit(targetModel, split = "\\_"))[2], split = "rec"))[2]
    
    if(hasASE > 0) {
      ASEscenario <- "ASE"
    } else {
      ASEscenario <- "noASE"
    }
    # stargetModel <- paste0(targetModel, s)
    # modelWrkDir <- paste0(wrkDir, "/", stargetModel)
    # openTableM2Het <- read.table(paste0(modelWrkDir, "/sumStats_", stargetModel, "_perSite.txt"), sep = "\t", header = T, na.strings = NaN)
    # openTableM2Het <- openTableM2Het %>% mutate(Model=rep(stargetModel, nrow(openTableM2Het)), ASEModel = rep(ASEscenario[aseModelIndex], nrow(openTableM2Het)))
    # 
    # fullTableM2Het <- rbind(fullTableM2Het, openTableM2Het %>% mutate(Distance=rep(phyDist[which(models %in% targetModel)], nrow(openTableM2Het))))
      
  
    #fixed mutations
    openTableM2Fixed <- read.table(paste0(targetModel, "_g10000_m2_fMut.out"), sep = "\t", header = F)
    openTableM2Fixed <- openTableM2Fixed %>% mutate(Model=rep(targetModel, nrow(openTableM2Fixed)), ASEModel = rep(ASEscenario, nrow(openTableM2Fixed)))
    
    fullTableM2fixed <- rbind(fullTableM2fixed, openTableM2Fixed %>% mutate(Distance=rep(whatRecomb, nrow(openTableM2Fixed))))
    
    openTableM3Fixed <- read.table(paste0(targetModel, "_g10000_m3_fMut.out"), sep = "\t", header = F)
    openTableM3Fixed <- openTableM3Fixed %>% mutate(Model=rep(targetModel, nrow(openTableM3Fixed)), ASEModel = rep(ASEscenario, nrow(openTableM3Fixed)))
    
    fullTableM3fixed <- rbind(fullTableM3fixed, openTableM3Fixed %>% mutate(Distance=rep(whatRecomb, nrow(openTableM3Fixed))))
    
  }

dim(fullTableM2fixed)
dim(fullTableM3fixed)
#dim(fullTableM2Het)
fullTableM3fixed <- fullTableM3fixed %>% mutate(GenetDist=case_when(Distance == "0.5kb" ~ recRate[names(recRate) == "0.5kb"], 
                                                Distance == "1kb" ~ recRate[names(recRate) == "1kb"], 
                                                Distance == "5kb" ~ recRate[names(recRate) == "5kb"], 
                                                Distance == "10kb" ~ recRate[names(recRate) == "10kb"], 
                                                Distance == "500kb" ~ recRate[names(recRate) == "500kb"]))

###-----------------
## Regulatory mutations
###-----------------
coefsM3 <- coef(lm(V1 ~ GenetDist, data = fullTableM3fixed %>% filter(ASEModel == "ASE")))
fullTableM3fixed$Distance <- factor(fullTableM3fixed$Distance, levels=unique(fullTableM3fixed$Distance))
p <- fullTableM3fixed %>% ggplot(aes(x=Distance, y=V1, group = GenetDist)) + 
  geom_boxplot() + 
  facet_wrap(~ASEModel) + 
  labs(title=paste0("Regulatory mutations: ", creTitle)) +
  ylab("Number of fixed mutations") +
  theme_bw()
p + geom_smooth(method='lm', se=F,aes(group=1)) 
ggsave(file=paste0("CRE",creSize, "kb_m3_acrossGenetDist.pdf"), width = 20, height = 10, units = "cm")


