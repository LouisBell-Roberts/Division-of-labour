#Shows your working directory
getwd()
#Sets the working directory

#######CHANGE WORKING DIRECTORY#######
setwd("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Trees/Nelsen")

library(phytools)
library(ape)
library(phangorn)
library(geiger)
library(dplyr)
##Read in the tree and database##
anttree <- read.tree(file = "Dryad_Supplementary_File_4_RAxML_bipartitions.result.ladderized.dropped.tre")
antdata <- read.csv("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Cleaned/Data.csv", header = T)

#Mucking about#
plot(anttree,no.margin=TRUE,edge.width=2)
anttree
str(anttree)
anttree$tip.label
plotTree(anttree,ftype="i",fsize=0.6,lwd=1)


#### ADVICE FROM AN OLD SCRIPT - NA values must be excluded from database by leaving the cell blank####
#Match the tip labels from the tree and the species from the database and prune the phylogeny accordingly
#It is the material that is in the first named set, that is not in the second named set
pruned.tree<-drop.tip(anttree, setdiff(anttree$tip.label, antdata$animal))
pruned.tree
plot.phylo(pruned.tree)
plotTree(pruned.tree,ftype="i",fsize=0.2,lwd=1)

#The intersection of two sets is the material that they have in common: intersect(setA, setB)
#There is also a built-in function setequal for testing if two sets are equal
#You can use %in% for comparing sets. The result is a logical vector whose length matches the vector on the left


####################################
####Creates high resolution TIFF file of the phylogeny so that I can read the tips - be careful not to create files that are too large####
#tiff("Plot2.tiff", width = 30, height = 37, units = 'in', res = 300)
#plotTree(pruned.tree,ftype="i",fsize=0.5,lwd=1)
#dev.off()

#Uses less memory when png:

#png("Plot_pruned.png", width = 30, height = 37, units = 'in', res = 300)
#plotTree(pruned.tree,ftype="i",fsize=0.5,lwd=1)
#dev.off()

#png("Nelsen.png", width = 30, height = 37, units = 'in', res = 300)
#plotTree(anttree,ftype="i",fsize=0.19,lwd=1)
#dev.off()
####################################


#"is.element" Tests if two vectors (lists) contain the same items (species in this case).
is.element(antdata$animal, pruned.tree$tip.label)

Harmonic <- antdata$eff.mating.freq.MEAN.harmonic

datafilter_All <- filter(antdata, Caste1 == 0 | Caste1 == 1 | Caste1 == 2 | Caste1 == 3 | Caste1 == 4, colony.size >= 0, Harmonic >=0, 
                     W.policing.clean == 0 | W.policing.clean == 1)
datafilter_Caste_CS_MF <- filter(antdata, Caste1 == 0 | Caste1 == 1 | Caste1 == 2 | Caste1 == 3 | Caste1 == 4, colony.size >= 0, Harmonic >=0)
datafilter_Caste_MF <- filter(antdata, Caste1 == 0 | Caste1 == 1 | Caste1 == 2 | Caste1 == 3 | Caste1 == 4, Harmonic >=0)
datafilter_Caste_CS <- filter(antdata, Caste1 == 0 | Caste1 == 1 | Caste1 == 2 | Caste1 == 3 | Caste1 == 4, colony.size >= 0)
datafilter_Caste_CatMF <- filter(antdata, Caste1 == 0 | Caste1 == 1 | Caste1 == 2 | Caste1 == 3 | Caste1 == 4, polyandry.clean >= 0)
datafilter_CS_MF <- filter(antdata, colony.size >= 0, Harmonic >=0)
datafilter_Caste_CatMF

##Plotting trees##
#Caste_CS_MF
#Database only overlaps with the Formicoids
pruned.tree1<-drop.tip(anttree, setdiff(anttree$tip.label, datafilter_Caste_CS_MF$animal))
pruned.tree1
plotTree(pruned.tree1,ftype="i",fsize=0.7,lwd=1)

#Caste_MF
#Database only overlaps with the Formicoids
pruned.tree2<-drop.tip(anttree, setdiff(anttree$tip.label, datafilter_Caste_MF$animal))
pruned.tree2
plotTree(pruned.tree2,ftype="i",fsize=0.7,lwd=1)

#Caste_CS
#Database overlaps with the Formicoids and the Poneroids
pruned.tree3<-drop.tip(anttree, setdiff(anttree$tip.label, datafilter_Caste_CS$animal))
pruned.tree3
plotTree(pruned.tree3,ftype="i",fsize=0.4,lwd=1)

#CS_MF
#Database only overlaps with the Formicoids
pruned.tree4<-drop.tip(anttree, setdiff(anttree$tip.label, datafilter_CS_MF$animal))
pruned.tree4
plotTree(pruned.tree4,ftype="i",fsize=0.7,lwd=1)

#Caste_CategoricalMF
#
pruned.tree5<-drop.tip(anttree, setdiff(anttree$tip.label, datafilter_Caste_CatMF$animal))
pruned.tree5
plotTree(pruned.tree5,ftype="i",fsize=0.7,lwd=1)




a <- filter(antdata, Caste1 == 0 | Caste1 == 1 | Caste1 == 2 | Caste1 == 3 | Caste1 == 4 | colony.size >= 0| Harmonic >=0| W.policing.clean == 0 | W.policing.clean == 1, type == "ant")
View(a %>% select(Genus, animal, colony.size, eff.mating.freq.MEAN.harmonic, Caste1))



###############################Ancestral state reconstruction_Mating Frequency########################################
#continuous trait: eff.mating.freq.MEAN.harmonic

#######CHANGE WORKING DIRECTORY#######
setwd("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Cleaned")
#Required packages
library(ape)
library(phytools)
library(geiger)
#doesn't seem to be working library(phylolm)
library(corHMM)
library(caper)
library(stringr)
library(dplyr)
library(tidyverse)

#Read in data file
antdata <- read.csv("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Cleaned/Data.csv", header = T)
#Read in genus-level tree file
anttree <- read.tree(file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Trees/Genus_tree/Genus_Level_ML_TREE_treepl_185_outsdroppeds.ladderized.increasing.tre")

#Creates a data frame which splits the "animal" column into a species and a genus name
#Dataframe with numberous useful columns - however this may screw up the analysis further on
#data <- as.data.frame(
#  cbind(
#    species = word(antdata$animal,2,sep = "_"),
#    genus = word(antdata$animal,sep = "_"),
#    mating_frequency = antdata$eff.mating.freq.MEAN.harmonic,
#    family = antdata$Family
#  )
#) 

#Dataframe with just genus and mating frequency
data <- as.data.frame(
  cbind(
    genus = word(antdata$animal,sep = "_"),
    mating_frequency = antdata$eff.mating.freq.MEAN.harmonic  
    )
) 

View(data)

#Mating frequency is currently a factor. Convert to a numeric
data$mating_frequency <- as.numeric(as.character(data$mating_frequency))

data_filtered <- filter(data, mating_frequency>= 0)
View(data_filtered)

###Create genus-level averages
library(data.table)
keys <- colnames(data_filtered)[!grepl('mating_frequency',colnames(data_filtered))]
X <- as.data.table(data_filtered)
data1<-X[,list(MF= mean(mating_frequency)),keys]
View(X[,list(MF= mean(mating_frequency)),keys])

###Prune tree
#Match the tip labels from the tree and the species from the database and prune the phylogeny accordingly
#It is the material that is in the first named set, that is not in the second named set
pruned.tree<-drop.tip(anttree, setdiff(anttree$tip.label, data1$genus))
pruned.tree
plot.phylo(pruned.tree)
plotTree(pruned.tree,ftype="i",fsize=0.4,lwd=1)


#Filter through my dataframe which has genus-level averages and select only the rows that match the tips of my tree
data2<-filter(data1, genus %in% pruned.tree$tip.label)

#Log transform the mating frequency data
data2$MF<-log(data2$MF)
View(data2)
#Convert the values in a column into row names in an existing data frame
#So that my data is formatted in the same way as the tutorial
data3<-data2 %>% remove_rownames() %>% column_to_rownames(var="genus")
View(data3)

###Ensure that phylogeny tip labels are exactly matched to the row names of the database
a<-name.check(pruned.tree,data3)

#Visualise as a dot tree
dotTree(pruned.tree,data3,length=20,ftype="i")

#Set the dataframe to a vector
data4<-as.matrix(data3)[,1]

#####Estimate ancestral states#####
fit<-fastAnc(pruned.tree,data4,vars=TRUE,CI=TRUE)
fit

fit$CI[1,]
range(data4)

##Visualisation: projection of the reconstruction onto the edges of the tree
obj<-contMap(pruned.tree,data4,plot=FALSE)
plot(obj,legend=0.7*max(nodeHeights(pruned.tree)),
     fsize=c(0.7,0.9))

phenogram(pruned.tree,data4,fsize=0.6,spread.costs=c(1,0))


########END################################################################################

#The intersection of two sets is the material that they have in common: intersect(setA, setB)
#There is also a built-in function setequal for testing if two sets are equal
#You can use %in% for comparing sets. The result is a logical vector whose length matches the vector on the left

#"is.element" Tests if two vectors (lists) contain the same items (species in this case).
is.element(data1$genus, pruned.tree$tip.label)



