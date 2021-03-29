#BEE conference: Ancestral state reconstruction: Caste as a discrete trait
#Louis Bell-Roberts
#12/03/2021

setwd("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Trees/Nelsen_ultrametric_species/")

library(tidyverse)
library(ape)
library(phytools)
library(geiger)
library(phylolm)
library(corHMM)
library(caper)

#May be smart to treat caste as a categorical variable rather than as a numerical variable. Not sure yet
d <- read.csv("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Primary Dataset/Data_caste_mating_colonyS_WPM_QueenN_cleaned.csv", header = T)
d$Caste1 <- as.numeric(as.character(d$Caste1))
data <- d

anttree_species <- read.tree(file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Trees/Nelsen_ultrametric_species/ultrametric_Nelsen_sp.tre")

#Select only ant species
antdata<-filter(data, type=="ant")
#View(antdata)

#Remove missing data
antdata_caste <- filter(antdata, Caste1 >= 0)
#View(antdata_caste)

#Prune database and phylogeny

#Match the tip labels from the tree and the species from the database and prune the phylogeny accordingly
#It is the material that is in the first named set, that is not in the second named set
pruned.tree_sp<-drop.tip(anttree_species, setdiff(anttree_species$tip.label, antdata_caste$animal))
pruned.tree_sp
plotTree(pruned.tree_sp,ftype="i",fsize=0.4,lwd=1)

#Filter through my dataframe and select only the rows that match the tips of my tree
antdata_caste.1<-filter(antdata_caste, animal %in% pruned.tree_sp$tip.label)
View(antdata_caste.1)

#Select only the animal and caste columns
antdata_caste.2 <- antdata_caste.1 %>% dplyr::select(animal, Caste1)

#ASR
#ancestral_hrm1<-corHMM(phy = pruned.tree_sp, data = antdata_caste.2,rate.cat = 1)
#ancestral_hrm2<-corHMM(phy = pruned.tree_sp, data = antdata_caste.2,rate.cat = 2)
#ancestral_hrm3<-corHMM(phy = pruned.tree_sp, data = antdata_caste.2,rate.cat = 1, model = "ER")
#ancestral_hrm4<-corHMM(phy = pruned.tree_sp, data = antdata_caste.2,rate.cat = 1, model = "ARD")
ancestral_hrm5<-corHMM(phy = pruned.tree_sp, data = antdata_caste.2,rate.cat = 1, model = "SYM")


#Model selection - lower AIC is better
#ancestral_hrm1$AICc
#ancestral_hrm2$AICc
#ancestral_hrm3$AICc
#ancestral_hrm4$AICc
ancestral_hrm5$AICc #has the lowest AIC

#Plot ASR
plotvec<-as.factor(antdata_caste.2$Caste1[match(pruned.tree_sp$tip.label,table=antdata_caste.2$animal)])
plotRECON(phy=ancestral_hrm5$phy,likelihoods = ancestral_hrm5$states,pie.cex=0.3, show.tip.label = T,
          piecolors = c("red","black", "yellow", "purple"),
          tip.color = c("red","black", "yellow", "purple")[plotvec])
















