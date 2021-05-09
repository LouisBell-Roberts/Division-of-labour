#BEE conference: Ancestral state reconstruction: Caste as a discrete trait or MF as a discrete trait
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

#####Caste analysis

#May be smart to treat caste as a categorical variable rather than as a numerical variable. Not sure yet
d <- read.csv("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Primary Dataset/Data_caste_mating_colonyS_WPM_QueenN_cleaned.csv", header = T)
#d <- read.csv(file.choose(), header = T)
d$Caste3 <- as.numeric(as.character(d$Caste3))
data <- d

anttree_species <- read.tree(file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Trees/Nelsen_ultrametric_species/ultrametric_Nelsen_sp.tre")
anttree_species <- read.tree(file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Trees/Polytomy_tree/Genus_polytomy_tree.tre")

#anttree_species <- read.tree(file.choose())
  
#Select only ant species
antdata<-filter(data, type=="ant")
#View(antdata)

#Remove missing data
antdata_caste <- filter(antdata, Caste3 >= 0)
View(antdata_caste)

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
antdata_caste.2 <- antdata_caste.1 %>% dplyr::select(animal, Caste3)

#ASR
pruned.tree_sp.jkl<-multi2di(pruned.tree_sp)

ancestral_hrm1<-corHMM(phy = pruned.tree_sp, data = antdata_caste.2,rate.cat = 1)
ancestral_hrm2<-corHMM(phy = pruned.tree_sp, data = antdata_caste.2,rate.cat = 2)
ancestral_hrm3<-corHMM(phy = pruned.tree_sp, data = antdata_caste.2,rate.cat = 1, model = "ER")
ancestral_hrm4<-corHMM(phy = pruned.tree_sp, data = antdata_caste.2,rate.cat = 1, model = "ARD")
ancestral_hrm5<-corHMM(phy = pruned.tree_sp, data = antdata_caste.2,rate.cat = 1, model = "SYM")
ancestral_hrm6<-corHMM(phy = pruned.tree_sp.jkl, data = antdata_caste.2,rate.cat = 1, model = "SYM") #randomly resolved tree model
ancestral_hrm7<-corHMM(phy = pruned.tree_sp, data = antdata_caste.2,rate.cat = 2, model = "ER")
ancestral_hrm8<-corHMM(phy = pruned.tree_sp, data = antdata_caste.2,rate.cat = 2, model = "ARD")
ancestral_hrm9<-corHMM(phy = pruned.tree_sp, data = antdata_caste.2,rate.cat = 2, model = "SYM")
ancestral_hrm10<-corHMM(phy = pruned.tree_sp.jkl, data = antdata_caste.2,rate.cat = 2, model = "SYM") #randomly resolved tree model


#Model selection - lower AIC is better
ancestral_hrm1$AICc
ancestral_hrm2$AICc
ancestral_hrm3$AICc
ancestral_hrm4$AICc
ancestral_hrm5$AICc #lowest AICc for single rate category models
ancestral_hrm6$AICc #when randomly resolving the polytomies, the AIC result is hardly different
ancestral_hrm7$AICc
ancestral_hrm8$AICc
ancestral_hrm9$AICc #lowest AICc value by a considerable way
ancestral_hrm10$AICc #when randomly resolving the polytomies, the AIC result is hardly different

#Plot ASR
plotvec<-as.factor(antdata_caste.2$Caste3[match(pruned.tree_sp$tip.label,table=antdata_caste.2$animal)])
plotRECON(phy=ancestral_hrm9$phy,likelihoods = ancestral_hrm9$states,pie.cex=0.3, show.tip.label = T,
          piecolors = c("red","black", "yellow", "purple", "blue", "white", "orange", "grey"),
          tip.color = c("red","black", "yellow", "purple")[plotvec])

plotRECON(phy=ancestral_hrm5$phy,likelihoods = ancestral_hrm5$states,pie.cex=0.3, show.tip.label = T,
          piecolors = c("red","black", "yellow", "purple"),
          tip.color = c("red","black", "yellow", "purple")[plotvec])







###### MF analysis (binary)


#Read in data and tree
data <- read.csv("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Primary Dataset/Data_caste_mating_colonyS_WPM_QueenN_cleaned.csv", header = T)
anttree_species <- read.tree(file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Trees/Polytomy_tree/Genus_polytomy_tree.tre")

#Select only ant species
antdata<-filter(data, type=="ant")
View(antdata)

#Create MF categorical variable 'combined' based on polyandry.clean and eff.mating.freq.MEAN.harmonic
##Create bins for continuous MF variable
antdata<-antdata%>%mutate(MFbins = cut(eff.mating.freq.MEAN.harmonic, breaks = c(-Inf,1,2,Inf), c("0", "1", "2")))
View(select(antdata, eff.mating.freq.MEAN.harmonic, MFbins))
##Make polyandry.clean a factor
antdata$polyandry.clean<-as.factor(antdata$polyandry.clean)
##Combine the two columns. MFbins takes priority and if NA, then value taken from polyandry.clean
antdata$MFcombined<-coalesce(antdata$MFbins, antdata$polyandry.clean)
View(select(antdata, animal, polyandry.clean, MFbins, MFcombined))


#Remove missing data
##Set MFcombined to numeric variable
antdata$MFcombined<-as.numeric(as.character(antdata$MFcombined))

antdata_MF <- filter(antdata, MFcombined >= 0)
View(antdata_MF)
#Prune database and phylogeny

#Match the tip labels from the tree and the species from the database and prune the phylogeny accordingly
#It is the material that is in the first named set, that is not in the second named set
pruned.tree_sp<-drop.tip(anttree_species, setdiff(anttree_species$tip.label, antdata_MF$animal))
pruned.tree_sp
plotTree(pruned.tree_sp,ftype="i",fsize=0.4,lwd=1)

#Filter through my dataframe and select only the rows that match the tips of my tree
antdata_MF.1<-filter(antdata_MF, animal %in% pruned.tree_sp$tip.label)
View(antdata_MF.1)

#Select only the animal and caste columns
antdata_MF.2 <- antdata_MF.1 %>% dplyr::select(animal, MFcombined)

#ASR
pruned.tree_sp.jkl<-multi2di(pruned.tree_sp)

ancestral_hrm1<-corHMM(phy = pruned.tree_sp, data = antdata_MF.2,rate.cat = 1)
ancestral_hrm2<-corHMM(phy = pruned.tree_sp, data = antdata_MF.2,rate.cat = 2)
ancestral_hrm3<-corHMM(phy = pruned.tree_sp, data = antdata_MF.2,rate.cat = 1, model = "ER")
ancestral_hrm4<-corHMM(phy = pruned.tree_sp, data = antdata_MF.2,rate.cat = 1, model = "ARD")
ancestral_hrm5<-corHMM(phy = pruned.tree_sp, data = antdata_MF.2,rate.cat = 1, model = "SYM")
ancestral_hrm6<-corHMM(phy = pruned.tree_sp.jkl, data = antdata_MF.2,rate.cat = 1, model = "SYM") #randomly resolved tree model
ancestral_hrm7<-corHMM(phy = pruned.tree_sp, data = antdata_MF.2,rate.cat = 2, model = "ER")
ancestral_hrm8<-corHMM(phy = pruned.tree_sp, data = antdata_MF.2,rate.cat = 2, model = "ARD")
ancestral_hrm9<-corHMM(phy = pruned.tree_sp, data = antdata_MF.2,rate.cat = 2, model = "SYM")
ancestral_hrm10<-corHMM(phy = pruned.tree_sp.jkl, data = antdata_MF.2,rate.cat = 2, model = "SYM") #randomly resolved tree model


#Model selection - lower AIC is better
ancestral_hrm1$AICc
ancestral_hrm2$AICc
ancestral_hrm3$AICc
ancestral_hrm4$AICc # has the lowest AICc
ancestral_hrm5$AICc 
ancestral_hrm6$AICc #when randomly resolving the polytomies, the AIC result is hardly different
ancestral_hrm7$AICc
ancestral_hrm8$AICc
ancestral_hrm9$AICc 
ancestral_hrm10$AICc #when randomly resolving the polytomies, the AIC result is hardly different

#Plot ASR
plotvec<-as.factor(antdata_MF.2$MFcombined[match(pruned.tree_sp$tip.label,table=antdata_MF.2$animal)])
plotRECON(phy=ancestral_hrm4$phy,likelihoods = ancestral_hrm4$states,pie.cex=0.3, show.tip.label = T,
          piecolors = c("red","black", "yellow"),
          tip.color = c("red","black", "yellow")[plotvec])







###### Queen number analysis (binary)

#Read in data and tree
data <- read.csv("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Primary Dataset/Data_caste_mating_colonyS_WPM_QueenN_cleaned.csv", header = T)
anttree_species <- read.tree(file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Trees/Polytomy_tree/Genus_polytomy_tree.tre")

#Select only ant species
antdata<-filter(data, type=="ant")
View(antdata)

#Remove missing data
#antdata_PG <- filter(antdata, number.queens.MEAN >= 0, Reference.2 != "")
antdata_PG <- filter(antdata, polygyny.clean >= 0, Reference.2 != "")

View(antdata_PG)

#Prune database and phylogeny

#Match the tip labels from the tree and the species from the database and prune the phylogeny accordingly
#It is the material that is in the first named set, that is not in the second named set
pruned.tree_sp<-drop.tip(anttree_species, setdiff(anttree_species$tip.label, antdata_PG$animal))
pruned.tree_sp
plotTree(pruned.tree_sp,ftype="i",fsize=0.4,lwd=1)

#Filter through my dataframe and select only the rows that match the tips of my tree
antdata_PG.1<-filter(antdata_PG, animal %in% pruned.tree_sp$tip.label)
View(antdata_PG.1)

#Select only the animal and caste columns
antdata_PG.2 <- antdata_PG.1 %>% dplyr::select(animal, polygyny.clean)

#ASR
pruned.tree_sp.jkl<-multi2di(pruned.tree_sp)

ancestral_hrm1<-corHMM(phy = pruned.tree_sp, data = antdata_PG.2,rate.cat = 1)
ancestral_hrm2<-corHMM(phy = pruned.tree_sp, data = antdata_PG.2,rate.cat = 2)
ancestral_hrm3<-corHMM(phy = pruned.tree_sp, data = antdata_PG.2,rate.cat = 1, model = "ER")
ancestral_hrm4<-corHMM(phy = pruned.tree_sp, data = antdata_PG.2,rate.cat = 1, model = "ARD")
ancestral_hrm5<-corHMM(phy = pruned.tree_sp, data = antdata_PG.2,rate.cat = 1, model = "SYM")
ancestral_hrm6<-corHMM(phy = pruned.tree_sp.jkl, data = antdata_PG.2,rate.cat = 1, model = "SYM") #randomly resolved tree model
ancestral_hrm7<-corHMM(phy = pruned.tree_sp, data = antdata_PG.2,rate.cat = 2, model = "ER")
ancestral_hrm8<-corHMM(phy = pruned.tree_sp, data = antdata_PG.2,rate.cat = 2, model = "ARD")
ancestral_hrm9<-corHMM(phy = pruned.tree_sp, data = antdata_PG.2,rate.cat = 2, model = "SYM")
ancestral_hrm10<-corHMM(phy = pruned.tree_sp.jkl, data = antdata_PG.2,rate.cat = 2, model = "SYM") #randomly resolved tree model


#Model selection - lower AIC is better
ancestral_hrm1$AICc
ancestral_hrm2$AICc
ancestral_hrm3$AICc
ancestral_hrm4$AICc #lowest AICc for single rate category models
ancestral_hrm5$AICc 
ancestral_hrm6$AICc #when randomly resolving the polytomies, the AIC result is hardly different
ancestral_hrm7$AICc #lowest AICc value overall
ancestral_hrm8$AICc
ancestral_hrm9$AICc 
ancestral_hrm10$AICc #when randomly resolving the polytomies, the AIC result is hardly different



#Plot ASR
plotvec<-as.factor(antdata_PG.2$polygyny.clean[match(pruned.tree_sp$tip.label,table=antdata_PG.2$animal)])
plotRECON(phy=ancestral_hrm4$phy,likelihoods = ancestral_hrm4$states,pie.cex=0.3, show.tip.label = T,
          piecolors = c("red","black"),
          tip.color = c("red","black")[plotvec])



###### Worker policing analysis (binary)


#Read in data and tree
data <- read.csv("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Primary Dataset/Data_caste_mating_colonyS_WPM_QueenN_cleaned.csv", header = T)
anttree_species <- read.tree(file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Trees/Polytomy_tree/Genus_polytomy_tree.tre")

#Select only ant species
antdata<-filter(data, type=="ant")
#select only species with data on worker sons
antdata_WS<-filter(antdata, Worker.sons.clean >= 0)

#Select only worker policing and species columns
View(select(antdata_WS, animal, W.policing.clean, Worker.sons.clean))
View(antdata$Worker.sons.clean)

#Prune database and phylogeny

#Match the tip labels from the tree and the species from the database and prune the phylogeny accordingly
#It is the material that is in the first named set, that is not in the second named set
pruned.tree_sp<-drop.tip(anttree_species, setdiff(anttree_species$tip.label, antdata_WS$animal))
pruned.tree_sp
plotTree(pruned.tree_sp,ftype="i",fsize=0.4,lwd=1)













