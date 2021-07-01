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
d$Caste3 <- as.numeric(as.character(d$Caste3))
data <- d

#anttree_species <- read.tree(file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Trees/Nelsen_ultrametric_species/ultrametric_Nelsen_sp.tre")
anttree_species <- read.tree(file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Trees/Polytomy_tree/Genus_polytomy_tree.tre")
#anttree_species <- read.tree(file.choose())

#Select only ant species
antdata<-filter(data, type=="ant")
#View(antdata)

########

#Caste vs. MF filtering
antdata_MF <- filter(antdata, Caste3 >=1, eff.mating.freq.MEAN.harmonic >=1) #Using Caste3 instead

#Make data binary for MF
antdata_MF_bin <- antdata_MF

antdata_MF_bin$eff.mating.freq.MEAN.harmonic<- cut(antdata_MF$eff.mating.freq.MEAN.harmonic,
                    breaks=c(-Inf, 2, Inf),
                    labels=c("low","high"))

#Make data binary for Caste
antdata_MF_Caste_bin <- antdata_MF_bin

#Alternatively, using Caste3
antdata_MF_Caste_bin$Caste3<- cut(antdata_MF_Caste_bin$Caste3,
                                  breaks=c(-Inf, 1.99, Inf),
                                  labels=c("low","high"))

#Prune tree
pruned.tree_sp<-drop.tip(anttree_species, setdiff(anttree_species$tip.label, antdata_MF_Caste_bin$animal))
pruned.tree_sp
plotTree(pruned.tree_sp,ftype="i",fsize=0.4,lwd=1)
dev.off()

#Filter through my dataframe and select only the rows that match the tips of my tree
antdata_MF_Caste_bin.1<-filter(antdata_MF_Caste_bin, animal %in% pruned.tree_sp$tip.label)
#View(antdata_MF_Caste_bin.1)

#Select columns of interest: animal, MF, caste
antdata_MF_Caste_bin.2 <- antdata_MF_Caste_bin.1 %>% dplyr::select(animal, eff.mating.freq.MEAN.harmonic, Caste3)

#Set dataframe as a matrix
antdata_MF_Caste_bin.3 <- antdata_MF_Caste_bin.2
as.matrix(antdata_MF_Caste_bin.3)

#Needs row names removed
antdata_MF_Caste_bin.3 <-antdata_MF_Caste_bin.3 %>% remove_rownames() %>% column_to_rownames(var="animal")

#Black is high, Caste is on the right
dotTree(pruned.tree_sp,x = antdata_MF_Caste_bin.3,data.type="discrete")

#Format data for fit.pagel
##SPECIES ORDER IN THE TREE DOES NOT NEED TO MATCH THE SPECIES ORDER OF MY MATRICES - I've checked this (hashed out below)

Tcaste <- antdata_MF_Caste_bin.2$Caste3
names(Tcaste)<-antdata_MF_Caste_bin.2$animal
TMF <- antdata_MF_Caste_bin.2$eff.mating.freq.MEAN.harmonic
names(TMF)<-antdata_MF_Caste_bin.2$animal

#Fit the model: Trait 1 = MF, Trait 2 = Caste
##Smaller AIC = better model
ant_correlated_MF_Caste<-fitPagel(pruned.tree_sp, x = TMF, y = Tcaste, model = "ARD")
ant_correlated_MF_Caste
plot.fitPagel(ant_correlated_MF_Caste)

#Fitting uni-directional models of binary trait evolution i.e. x depends on y, but not the converse, or y depends on x, but not the converse
ant_correlated_MF_Caste_CASdep_on_MF<-fitPagel(pruned.tree_sp, x = TMF, y = Tcaste, model = "ARD", dep.var = "y") #Substitution rate in y depends on x 
ant_correlated_MF_Caste_CASdep_on_MF
plot.fitPagel(ant_correlated_MF_Caste_CASdep_on_MF)
ant_correlated_MF_Caste_MFdep_on_CAS<-fitPagel(pruned.tree_sp, x = TMF, y = Tcaste, model = "ARD", dep.var = "x") #Substitution rate in x depends on y
ant_correlated_MF_Caste_MFdep_on_CAS
plot.fitPagel(ant_correlated_MF_Caste_MFdep_on_CAS)

aic<-setNames(c(ant_correlated_MF_Caste$independent.AIC,
                ant_correlated_MF_Caste_CASdep_on_MF$dependent.AIC,
                ant_correlated_MF_Caste_MFdep_on_CAS$dependent.AIC,
                ant_correlated_MF_Caste$dependent.AIC),
              c("independent","dependent y",
                "dependent x","dependent x&y"))
aic





############
#Matching the order of the phylo vector with the order of species and trait
#antdata_MF_Caste_bin.4 <- antdata_MF_Caste_bin.2[match(pruned.tree_sp$tip.label, antdata_MF_Caste_bin.2$animal),]

#Format data for fit.pagel ---- (dataframe rows reordered to match phylo tip labels)
#Tcaste.1 <- antdata_MF_Caste_bin.4$Caste3
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
antdata_CS <- filter(antdata, Caste3 >=1, colony.size >=1) #To calculate quartiles from polymorphic species only

#plot(ln(antdata_CS$colony.size))
#View(antdata_CS)
#mode(antdata_CS$colony.size)
#axis(2, at = seq(0, 20, by = 1))
#axis(side=2,at=c(1:20))

#Make data binary for CS
##Median value is 218.5
##25% quartile (59.35)
##Smallest CS at which caste can evolve (22)
##25% quartile for polymorphic species (161.4375)
##75% quartile (1151.50)
antdata_CS_bin <- antdata_CS
#low<=1000  high>1000
antdata_CS_bin$colony.size<- cut(antdata_CS$colony.size,
                                                   breaks=c(-Inf, 218.5, Inf),
                                                   labels=c("low","high"))

#Make data binary for Caste
antdata_CS_Caste_bin <- antdata_CS_bin

antdata_CS_Caste_bin$Caste3<- cut(antdata_CS_Caste_bin$Caste3,
                                  breaks=c(-Inf, 1.99, Inf),
                                  labels=c("low","high"))

#Prune tree
pruned.tree_sp<-drop.tip(anttree_species, setdiff(anttree_species$tip.label, antdata_CS_Caste_bin$animal))
pruned.tree_sp
plotTree(pruned.tree_sp,ftype="i",fsize=0.4,lwd=1)
dev.off()
#Filter through my dataframe and select only the rows that match the tips of my tree
antdata_CS_Caste_bin.1<-filter(antdata_CS_Caste_bin, animal %in% pruned.tree_sp$tip.label)
#View(antdata_CS_Caste_bin.1)

#Select columns of interest: animal, CS, caste
antdata_CS_Caste_bin.2 <- antdata_CS_Caste_bin.1 %>% dplyr::select(animal, colony.size, Caste3)

#Set dataframe as a matrix
antdata_CS_Caste_bin.3 <- antdata_CS_Caste_bin.2
as.matrix(antdata_CS_Caste_bin.3)

#Needs row names removed
antdata_CS_Caste_bin.3 <-antdata_CS_Caste_bin.3 %>% remove_rownames() %>% column_to_rownames(var="animal")

#Red is high, Caste is on the right
#dotTree(pruned.tree_sp,x = antdata_CS_Caste_bin.3,data.type="discrete")

#Format data for fit.pagel
Tcaste <- antdata_CS_Caste_bin.2$Caste3
names(Tcaste)<-antdata_CS_Caste_bin.2$animal
TCS <- antdata_CS_Caste_bin.2$colony.size
names(TCS)<-antdata_CS_Caste_bin.2$animal

#Fit the model: Trait 1 = CS, Trait 2 = Caste
ant_correlated_CS_Caste<-fitPagel(pruned.tree_sp, x = TCS, y = Tcaste)
ant_correlated_CS_Caste
plot.fitPagel(ant_correlated_CS_Caste)


#Fitting uni-directional models of binary trait evolution i.e. x depends on y, but not the converse, or y depends on x, but not the converse
ant_correlated_CS_Caste_CASdep_on_CS<-fitPagel(pruned.tree_sp, x = TCS, y = Tcaste, model = "ARD", dep.var = "y") #Substitution rate in y depends on x 
ant_correlated_CS_Caste_CASdep_on_CS
plot.fitPagel(ant_correlated_CS_Caste_CASdep_on_CS)

ant_correlated_CS_Caste_CSdep_on_CAS<-fitPagel(pruned.tree_sp, x = TCS, y = Tcaste, model = "ARD", dep.var = "x") #Substitution rate in x depends on y 
ant_correlated_CS_Caste_CSdep_on_CAS
plot.fitPagel(ant_correlated_CS_Caste_CSdep_on_CAS)


aic2<-setNames(c(ant_correlated_CS_Caste$independent.AIC,
                 ant_correlated_CS_Caste_CASdep_on_CS$dependent.AIC,
                 ant_correlated_CS_Caste_CSdep_on_CAS$dependent.AIC,
                ant_correlated_CS_Caste$dependent.AIC),
              c("independent","dependent y",
                "dependent x","dependent x&y"))
aic2
#####################






#MF vs. CS filtering ------ STILL WORKING ON THIS
#### WHAT TO DO IF I HAVE CONTINUOUS VARIABLE? (e.g. colony size) - are there correlated evolution models for continuous traits or do I need to make arbitrary cutoff?
#Could use a cut of 1000 - values start to rise quite rapidly after this in my data
antdata_MF_CS <- filter(antdata, eff.mating.freq.MEAN.harmonic >=1, colony.size >=1) #To calculate quartiles from polymorphic species only

#plot(ln(antdata_CS$colony.size))
#View(antdata_CS)
#mode(antdata_CS$colony.size)
#axis(2, at = seq(0, 20, by = 1))
#axis(side=2,at=c(1:20))

#Make data binary for CS
##Median value is 218.5
##25% quartile (59.35)
##Smallest CS at which caste can evolve (22)
##25% quartile for polymorphic species (161.4375)
##75% quartile (1151.50)
antdata_CS_bin <- antdata_MF_CS
#low<=1000  high>1000
antdata_CS_bin$colony.size<- cut(antdata_MF_CS$colony.size,
                                 breaks=c(-Inf, 161.4375, Inf),
                                 labels=c("low","high"))

#Make data binary for MF
antdata_MF_CS_bin <- antdata_CS_bin

antdata_MF_CS_bin$eff.mating.freq.MEAN.harmonic<- cut(antdata_MF_CS_bin$eff.mating.freq.MEAN.harmonic,
                                                   breaks=c(-Inf, 2, Inf),
                                                   labels=c("low","high"))

antdata_CS_Caste_bin

#Prune tree
pruned.tree_sp<-drop.tip(anttree_species, setdiff(anttree_species$tip.label, antdata_MF_CS_bin$animal))
pruned.tree_sp
plotTree(pruned.tree_sp,ftype="i",fsize=0.4,lwd=1)
dev.off()
#Filter through my dataframe and select only the rows that match the tips of my tree
antdata_MF_CS_bin.1<-filter(antdata_MF_CS_bin, animal %in% pruned.tree_sp$tip.label)
#View(antdata_MF_CS_bin.1)

#Select columns of interest: animal, CS, MF
antdata_MF_CS_bin.2 <- antdata_MF_CS_bin.1 %>% dplyr::select(animal, colony.size, eff.mating.freq.MEAN.harmonic)

#Set dataframe as a matrix
antdata_MF_CS_bin.3 <- antdata_MF_CS_bin.2
as.matrix(antdata_MF_CS_bin.3)

#Needs row names removed
antdata_MF_CS_bin.3 <-antdata_MF_CS_bin.3 %>% remove_rownames() %>% column_to_rownames(var="animal")

#Red is high, MF is on the right
#dotTree(pruned.tree_sp,x = antdata_MF_CS_bin.3,data.type="discrete")

#Format data for fit.pagel
TMF <- antdata_MF_CS_bin.2$eff.mating.freq.MEAN.harmonic
names(TMF)<-antdata_MF_CS_bin.2$animal
TCS <- antdata_MF_CS_bin.2$colony.size
names(TCS)<-antdata_MF_CS_bin.2$animal

#Fit the model: Trait 1 = CS, Trait 2 = Caste
ant_correlated_MF_CS<-fitPagel(pruned.tree_sp, x = TCS, y = TMF)
ant_correlated_MF_CS
plot.fitPagel(ant_correlated_MF_CS)


#Fitting uni-directional models of binary trait evolution i.e. x depends on y, but not the converse, or y depends on x, but not the converse
ant_correlated_CS_MF_MFdep_on_CS<-fitPagel(pruned.tree_sp, x = TCS, y = TMF, model = "ARD", dep.var = "y") #Substitution rate in y depends on x 
ant_correlated_CS_MF_MFdep_on_CS
plot.fitPagel(ant_correlated_CS_MF_MFdep_on_CS)

ant_correlated_CS_MF_CSdep_on_MF<-fitPagel(pruned.tree_sp, x = TCS, y = TMF, model = "ARD", dep.var = "x") #Substitution rate in x depends on y 
ant_correlated_CS_MF_CSdep_on_MF
plot.fitPagel(ant_correlated_CS_MF_CSdep_on_MF)


aic3<-setNames(c(ant_correlated_MF_CS$independent.AIC,
                 ant_correlated_CS_MF_MFdep_on_CS$dependent.AIC,
                 ant_correlated_CS_MF_CSdep_on_MF$dependent.AIC,
                 ant_correlated_MF_CS$dependent.AIC),
               c("independent","dependent y",
                 "dependent x","dependent x&y"))
aic3






