#Multimodel inference - dealing with phylogenetic uncertainty
#Louis Bell-Roberts
#09/04/2021

library(tidyverse)
library(ape)
library(phylolm)
library(phytools)
library(phangorn)
library(nlme)
library(caper)
library(AICcmodavg)

#Data file
d <- read.csv("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Primary Dataset/Data_caste_mating_colonyS_WPM_QueenN_cleaned.csv", header = T)
#d <- read.csv(file.choose(), header = T)
d$Caste1 <- as.numeric(as.character(d$Caste1))
data <- d

#Tree file - species
tree<-read.tree(text="(((A1,A2),(B1,B2,B3),C,D),E,F);")
anttree_species <- read.tree(file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Trees/Polytomy_tree/Genus_polytomy_tree.tre")

#########
#Filter through data and create my 'pruned.tree'
#Caste vs. MF filtering
antdata_MF <- filter(data, type == 'ant', Caste1 >=1, eff.mating.freq.MEAN.harmonic >=1)

#Prune tree - 63 species, 33 genera OR 121 species, 52 genera
pruned.tree_sp<-drop.tip(anttree_species, setdiff(anttree_species$tip.label, antdata_MF$animal))
pruned.tree_sp
plotTree(pruned.tree_sp,ftype="i",fsize=0.4,lwd=1)

######
#Liam's own function to create a tree with randomly resolved polytomies - doesn't seem to work on larger phylogenies
resolveRandom<-function(tree){
  while(!is.binary(tree)){
    nodes<-1:tree$Nnode+Ntip(tree)
    Nchildren<-function(node,tree) length(Children(tree,node))
    nchilds<-sapply(nodes,Nchildren,tree=tree)
    node<-nodes[which(nchilds>2)[1]]
    tree<-sample(resolveNode(tree,node),1)[[1]]
  }
  tree
}
#Use the function a single time
plotTree(resolveRandom(tree),type="cladogram",nodes="centered")
######

#####
#Using ape function instead of Revell's own function
trees<-replicate(10,multi2di(pruned.tree_sp),simplify=FALSE)
class(trees)<-"multiPhylo"





















