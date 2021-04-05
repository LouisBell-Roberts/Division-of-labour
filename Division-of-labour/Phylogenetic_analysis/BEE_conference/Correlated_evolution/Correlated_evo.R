#BEE conference: Correlated evolution models
#Louis Bell-Roberts
#13/03/2021

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
#d <- read.csv(file.choose(), header = T)
d$Caste1 <- as.numeric(as.character(d$Caste1))
data <- d

anttree_species <- read.tree(file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Trees/Nelsen_ultrametric_species/ultrametric_Nelsen_sp.tre")
#anttree_species <- read.tree(file.choose())

#Select only ant species
antdata<-filter(data, type=="ant")
#View(antdata)

########

#Caste vs. MF filtering
antdata_MF <- filter(antdata, Caste1 >=1, eff.mating.freq.MEAN.harmonic >=1)

#Make data binary for MF
antdata_MF_bin <- antdata_MF

antdata_MF_bin$eff.mating.freq.MEAN.harmonic<- cut(antdata_MF$eff.mating.freq.MEAN.harmonic,
                    breaks=c(-Inf, 2, Inf),
                    labels=c("low","high"))

#Make data binary for Caste
antdata_MF_Caste_bin <- antdata_MF_bin

antdata_MF_Caste_bin$Caste1<- cut(antdata_MF_Caste_bin$Caste1,
                                                   breaks=c(-Inf, 1.99, Inf),
                                                   labels=c("low","high"))

#Prune tree
pruned.tree_sp<-drop.tip(anttree_species, setdiff(anttree_species$tip.label, antdata_MF_Caste_bin$animal))
pruned.tree_sp
plotTree(pruned.tree_sp,ftype="i",fsize=0.4,lwd=1)

#Filter through my dataframe and select only the rows that match the tips of my tree
antdata_MF_Caste_bin.1<-filter(antdata_MF_Caste_bin, animal %in% pruned.tree_sp$tip.label)
#View(antdata_MF_Caste_bin.1)

#Select columns of interest: animal, MF, caste
antdata_MF_Caste_bin.2 <- antdata_MF_Caste_bin.1 %>% dplyr::select(animal, eff.mating.freq.MEAN.harmonic, Caste1)

#Set dataframe as a matrix
antdata_MF_Caste_bin.3 <- antdata_MF_Caste_bin.2
as.matrix(antdata_MF_Caste_bin.3)

#Needs row names removed
antdata_MF_Caste_bin.3 <-antdata_MF_Caste_bin.3 %>% remove_rownames() %>% column_to_rownames(var="animal")

#Black is high, Caste is on the right
dotTree(pruned.tree_sp,x = antdata_MF_Caste_bin.3,data.type="discrete")

#Format data for fit.pagel
##SPECIES ORDER IN THE TREE DOES NOT NEED TO MATCH THE SPECIES ORDER OF MY MATRICES - I've checked this (hashed out below)

Tcaste <- antdata_MF_Caste_bin.2$Caste1
names(Tcaste)<-antdata_MF_Caste_bin.2$animal
TMF <- antdata_MF_Caste_bin.2$eff.mating.freq.MEAN.harmonic
names(TMF)<-antdata_MF_Caste_bin.2$animal

#Fit the model: Trait 1 = MF, Trait 2 = Caste
ant_correlated_MF_Caste<-fitPagel(pruned.tree_sp, x = TMF, y = Tcaste)
ant_correlated_MF_Caste
plot.fitPagel(ant_correlated_MF_Caste)


############
#Matching the order of the phylo vector with the order of species and trait
#antdata_MF_Caste_bin.4 <- antdata_MF_Caste_bin.2[match(pruned.tree_sp$tip.label, antdata_MF_Caste_bin.2$animal),]

#Format data for fit.pagel ---- (dataframe rows reordered to match phylo tip labels)
#Tcaste.1 <- antdata_MF_Caste_bin.4$Caste1
#names(Tcaste.1)<-antdata_MF_Caste_bin.4$animal
#TMF.1 <- antdata_MF_Caste_bin.4$eff.mating.freq.MEAN.harmonic
#names(TMF.1)<-antdata_MF_Caste_bin.4$animal

#Fit the model: Trait 1 = MF, Trait 2 = Caste ---- (dataframe rows reordered to match phylo tip labels)
#ant_correlated_MF_Caste.1<-fitPagel(pruned.tree_sp, x = TMF.1, y = Tcaste.1)
#ant_correlated_MF_Caste.1
#plot.fitPagel(ant_correlated_MF_Caste.1)
############





#Caste vs. CS filtering ------ STILL WORKING ON THIS
#### WHAT TO DO IF I HAVE CONTINUOUS VARIABLE? (e.g. colony size) - are there correlated evolution models for continuous traits or do I need to make arbitrary cutoff?
#Could use a cut of 1000 - values start to rise quite rapidly after this in my data
antdata_CS <- filter(antdata, Caste1 >=1, colony.size >=1)

#plot(ln(antdata_CS$colony.size))
#View(antdata_CS)
#mode(antdata_CS$colony.size)
#axis(2, at = seq(0, 20, by = 1))
#axis(side=2,at=c(1:20))

#Make data binary for CS
antdata_CS_bin <- antdata_CS
#low<=1000  high>1000
antdata_CS_bin$colony.size<- cut(antdata_CS$colony.size,
                                                   breaks=c(-Inf, 2000, Inf),
                                                   labels=c("low","high"))

#Make data binary for Caste
antdata_CS_Caste_bin <- antdata_CS_bin

antdata_CS_Caste_bin$Caste1<- cut(antdata_CS_Caste_bin$Caste1,
                                  breaks=c(-Inf, 1.99, Inf),
                                  labels=c("low","high"))

#Prune tree
pruned.tree_sp<-drop.tip(anttree_species, setdiff(anttree_species$tip.label, antdata_CS_Caste_bin$animal))
pruned.tree_sp
plotTree(pruned.tree_sp,ftype="i",fsize=0.4,lwd=1)

#Filter through my dataframe and select only the rows that match the tips of my tree
antdata_CS_Caste_bin.1<-filter(antdata_CS_Caste_bin, animal %in% pruned.tree_sp$tip.label)
#View(antdata_CS_Caste_bin.1)

#Select columns of interest: animal, CS, caste
antdata_CS_Caste_bin.2 <- antdata_CS_Caste_bin.1 %>% dplyr::select(animal, colony.size, Caste1)

#Set dataframe as a matrix
antdata_CS_Caste_bin.3 <- antdata_CS_Caste_bin.2
as.matrix(antdata_CS_Caste_bin.3)

#Needs row names removed
antdata_CS_Caste_bin.3 <-antdata_CS_Caste_bin.3 %>% remove_rownames() %>% column_to_rownames(var="animal")

#Red is high, Caste is on the right
dotTree(pruned.tree_sp,x = antdata_CS_Caste_bin.3,data.type="discrete")

#Format data for fit.pagel
Tcaste <- antdata_CS_Caste_bin.2$Caste1
names(Tcaste)<-antdata_CS_Caste_bin.2$animal
TCS <- antdata_CS_Caste_bin.2$colony.size
names(TCS)<-antdata_CS_Caste_bin.2$animal

#Fit the model: Trait 1 = MF, Trait 2 = Caste
ant_correlated_CS_Caste<-fitPagel(pruned.tree_sp, x = TCS, y = Tcaste)
ant_correlated_CS_Caste
plot.fitPagel(ant_correlated_CS_Caste)




