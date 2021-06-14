## Create polytomies on species-level ant phylogeny
#Louis Bell-Roberts
#02/04/2021

setwd("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Trees/Polytomy_tree/")


library(tidyverse)
library(phytools)
library(ape)
library(phylolm)
library(geiger)

#Data file
d <- read.csv("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Primary Dataset/Data_caste_mating_colonyS_WPM_QueenN_cleaned.csv", header = T)
d$Caste1 <- as.numeric(as.character(d$Caste1))
data <- d

#Tree file - species
anttree_species <- read.tree(file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Trees/Nelsen_ultrametric_species/ultrametric_Nelsen_sp.tre")

data$Worker.sons.clean

#Filtering for all variables of interest
antdata <- filter(data, type == 'ant', (Caste1 >= 1 | eff.mating.freq.MEAN.harmonic >= 1 | colony.size >= 0 | 
                                           polygyny.clean >= 0 | W.policing.clean >= 0 | Worker.sons.clean >= 0))

#Find the species of interest from the database that are not represented in the phylogeny
no_overlap <- setdiff(antdata$animal, anttree_species$tip.label)

#Select the species rows from the database that are not in the phylogeny using the "no_overlap" vector
sp_add_phylo <- antdata[antdata$animal %in% no_overlap,]
#View(sp_add_phylo)
sp_add_phylo$animal


####### Create polytomies on the phylogeny
#Simultaneously check if a tree is rooted and has no polytomies
is.binary(anttree_species)

## add a set of species in a vector
species <- sp_add_phylo$animal
anttree_species1 <- anttree_species

for(i in 1:length(species)) anttree_species1<-add.species.to.genus(anttree_species1,species[i],where="root")
plotTree(anttree_species1,ftype="i",fsize=0.4,lwd=1)
#######

#Prune tree
pruned.tree_sp<-drop.tip(anttree_species1, setdiff(anttree_species1$tip.label, antdata$animal))
pruned.tree_sp
plotTree(pruned.tree_sp,ftype="i",fsize=0.4,lwd=1)

#Write the tree to file - Newick tree
#write.tree(pruned.tree_sp,file = "Polytomy_tree.tre")

###### - TREE MAY HAVE DUPLICATE TAXA - ERROR FROM FIGTREE - This problem should be fixed now
duplicated(pruned.tree_sp$tip.label)












###############################

##Create polytomy tree where every genus is a polytomy i.e. remove add species and remove all species-level phylogenetic relationships



#Data file
d <- read.csv("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Primary Dataset/Data_caste_mating_colonyS_WPM_QueenN_cleaned.csv", header = T)
d$Caste1 <- as.numeric(as.character(d$Caste1))
data <- d

#Tree file - this reads in my genus-level tree
anttree_species <- read.tree(file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Trees/Polytomy_tree/Trees used to create polytomy tree/ultrametric_Single_Taxon_Representative.tre")
plotTree(anttree_species,ftype="i",fsize=0.4,lwd=1)

#Filtering for all variables of interest
#antdata <- filter(data, type == 'ant', (Caste1 >= 1 | eff.mating.freq.MEAN.harmonic >= 1 | colony.size >= 0 | 
#                                          polygyny.clean >= 0 | W.policing.clean >= 0 | Worker.sons.clean >= 0))
antdata <- filter(data, type == 'ant')

#Find the species of interest from the database that are not represented in the phylogeny
no_overlap <- setdiff(antdata$animal, anttree_species$tip.label)

#Select the species rows from the database that are not in the phylogeny using the "no_overlap" vector
sp_add_phylo <- antdata[antdata$animal %in% no_overlap,]
#View(sp_add_phylo)
sp_add_phylo$animal


####### Create polytomies on the phylogeny
#Simultaneously check if a tree is rooted and has no polytomies
is.binary(anttree_species)

## add a set of species in a vector
species <- sp_add_phylo$animal
anttree_species1 <- anttree_species

#Create tree where every genus is a polytomy
##Had to adapt Liam's function by adding 'force.ultrametric' to the function - otherwise throws error in relation to branch lengths
##This is discussed in more detail on Liam's Phytools blog
###"root" and "random" are used to change the position of where new species are added to the tree

#anttree_species1<-force.ultrametric(anttree_species1)

for(i in 1:length(species)) anttree_species1<-force.ultrametric(add.species.to.genus(anttree_species1,species[i],where="root"))
plotTree(anttree_species1,ftype="i",fsize=0.4,lwd=1)

#If I want to replicate this function many times and create a multi-phylo object, try: *****SECTION NOT WORKING YET*****
#trees<-replicate(10,multi2di(pruned.tree_sp),simplify=FALSE)
trees <- replicate(2, for(i in 1:length(species)) anttree_species1<-force.ultrametric(add.species.to.genus(anttree_species1,species[i],where="random")), 
                   simplify = FALSE)
class(trees)<-"multiPhylo"
#***************************


#Write tree to file
#write.tree(anttree_species1, file='Genus_polytomy_tree.tre')

###### - TREE MAY HAVE DUPLICATE TAXA - ERROR FROM FIGTREE - should be fixed now
duplicated(anttree_species1$tip.label)


















