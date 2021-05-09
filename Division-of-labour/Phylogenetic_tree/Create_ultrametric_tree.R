#Function to create ultrametric trees
#Louis Bell-Roberts
#05/04/2021

library(ggplot2)
library(ape)

#Tree file - species. This tree is in NEWICK format
anttree_species <- read.tree(file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Trees/Nelsen/ML_TREE_treepl_185_outsdroppeds.ladderized.increasing.tre")

#Make tree ultrametric
u.anttree_species <- chronoMPL(anttree_species)
is.ultrametric(u.anttree_species)

##Ensure that the branch lengths between the orginial and final trees are the same
plot(anttree_species$edge.length, u.anttree_species$edge.length, xlab="Original tree", ylab="Post-chronos tree", pch=20, bty="n")
abline(a=0, b=1, col="gray", lwd=0.5) #This works

#Write the ultrametric tree to file as NEWICK tree
write.tree(u.anttree_species, file='ultrametric_Single_Taxon_Representative.tre')