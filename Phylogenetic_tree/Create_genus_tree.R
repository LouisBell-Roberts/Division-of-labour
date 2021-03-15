#Gijsbert Werner, Dep Zoology, Oxford
#Louis modifications - changed the code to produce a genus level tree for the 'outsdropped' species tree
#January 2021

library(ape)
library(phytools)
library(stringr)
library(dplyr)

setwd("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Trees/Nelsen")

###Load trees.
RAxML_bipartitions.result.ladderized.dropped<- read.tree(file = "Dryad_Supplementary_File_4_RAxML_bipartitions.result.ladderized.dropped.tre")
ML_TREE_treepl_185<- read.tree(file = "Dryad_Supplementary_File_7_ML_TREE_treepl_185.tre")
ML_TREE_treepl_185_outsdroppeds.ladderized.increasing.tre<- read.tree(file = "ML_TREE_treepl_185_outsdroppeds.ladderized.increasing.tre")

#RAxML_bipartitions.result.ladderized.dropped <-
#  read.tree(file = "../../Data/Dryad_Supplementary_File_4_RAxML_bipartitions.result.ladderized.dropped.tre")
#ML_TREE_treepl_185 <-
#  read.tree(file = "../../Data/Dryad_Supplementary_File_7_ML_TREE_treepl_185.tre")
#ML_TREE_treepl_185_outsdroppeds.ladderized.increasing.tre <-
#  read.tree(file = "../../Data/ML_TREE_treepl_185_outsdroppeds.ladderized.increasing.tre")

#Quick look at the trees
RAxML_bipartitions.result.ladderized.dropped
ML_TREE_treepl_185
ML_TREE_treepl_185_outsdroppeds.ladderized.increasing.tre
#About same numner of species, outsdropped will presumably have outgroups dropped (15 species fewer)

length(
  which(
    ML_TREE_treepl_185_outsdroppeds.ladderized.increasing.tre$tip.label %in% ML_TREE_treepl_185$tip.label
  )
)
length(
  which(
    ML_TREE_treepl_185$tip.label %in% ML_TREE_treepl_185_outsdroppeds.ladderized.increasing.tre$tip.label
  )
)
length(
  which(
    RAxML_bipartitions.result.ladderized.dropped$tip.label %in% ML_TREE_treepl_185_outsdroppeds.ladderized.increasing.tre$tip.label
  )
)
#Yes overlap for all of the smaller numbers.

#Species list Louis Sep 2020
spec_list <- read.csv("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Cleaned/Data.csv", header = T)
#spec_list <- read.csv("../../Data/Species.csv")
head(spec_list)
nrow(spec_list) #1536 rows/species

##Overlap with the tree
#Gijsbert creates a column that is functionally equivalent to my pre-made column 'animal'
#Add the 'genus' column to the 'species' column
spec_list$gs<-
  paste0(spec_list$Genus,"_",spec_list$species)

head(spec_list) #Remarks like 'could be synonym for xxx' can throw errors. Better to put them in a separate 'remarks' colum or similar

length(which(spec_list$gs %in% ML_TREE_treepl_185_outsdroppeds.ladderized.increasing.tre$tip.label))
#335 of the species from the data file are in the tree. Could potentially be increased a bit with some checking (synonyms, spelling etc)
#you could do thay manually, but there's also ways to automate parts of the process. Look at the package taxize (I can dig up some code too. )


###########Reduce the three trees to genus level

#First the RAXML tree

#Create a data frame with the full species and genus name for each species. 
spec.RAxML_bipartitions.result.ladderized.dropped <- as.data.frame(
  cbind(
    species = ML_TREE_treepl_185_outsdroppeds.ladderized.increasing.tre$tip.label,
    genus = word(
      ML_TREE_treepl_185_outsdroppeds.ladderized.increasing.tre$tip.label,
      sep = "_"
    )
  )
) 
head(spec.RAxML_bipartitions.result.ladderized.dropped) #yes, worked.

#Number of species per genus genus
spec_per_genus_RAxML_bipartitions.result.ladderized.dropped <-
  spec.RAxML_bipartitions.result.ladderized.dropped %>% group_by(genus) %>% summarise(num_genus =
                                                                                        n())

#Randomly select the first species for each genus. If the tree is properly formed (monophyletic genera), this should not matter.
spec.RAxML_bipartitions.result.ladderized.dropped.genus.level <-
  spec.RAxML_bipartitions.result.ladderized.dropped[which(!duplicated(spec.RAxML_bipartitions.result.ladderized.dropped$genus)), ]
spec.RAxML_bipartitions.result.ladderized.dropped.genus.level #So this is 299 rows, i.e. genera.

#Now actually drop the tree to genus level
genus.tree.RAxML_bipartitions.result.ladderized.dropped <-
  drop.tip(ML_TREE_treepl_185_outsdroppeds.ladderized.increasing.tre,
           tip = ML_TREE_treepl_185_outsdroppeds.ladderized.increasing.tre$tip.label[!ML_TREE_treepl_185_outsdroppeds.ladderized.increasing.tre$tip.label %in% 
                                                                          spec.RAxML_bipartitions.result.ladderized.dropped.genus.level$species])
#This drops all the species names that are not in the genus level data file we have created before. 

genus.tree.RAxML_bipartitions.result.ladderized.dropped #As expected, tree with 299 species remaining. 285 if we are using the outsdropped tree

#Assign to 'rename.genus.tree.RAxML_bipartitions.result.ladderized.dropped' so that I don't lose the original genus tree I created
rename.genus.tree.RAxML_bipartitions.result.ladderized.dropped <- genus.tree.RAxML_bipartitions.result.ladderized.dropped

#Optional: rename to genus name - since a tip now properly represents the genus as a whole, no longer the particular species we happened to select from among it. 
rename.genus.tree.RAxML_bipartitions.result.ladderized.dropped$tip.label<-
  word(genus.tree.RAxML_bipartitions.result.ladderized.dropped$tip.label,sep = "_")
rename.genus.tree.RAxML_bipartitions.result.ladderized.dropped

write.tree(rename.genus.tree.RAxML_bipartitions.result.ladderized.dropped,file = "Genus_Level_ML_TREE_treepl_185_outsdroppeds.ladderized.increasing.tre")

##Same procedure can be used for the other two trees... 
