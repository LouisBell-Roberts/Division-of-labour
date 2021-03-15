# BEE conference: phylogenetically controlled regressions for Caste vs. MF or CS
#Louis Bell-Roberts
#11/03/2021

setwd("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Trees/Nelsen_ultrametric_species/")

library(tidyverse)
library(ape)
library(phylolm)

#Data file
d <- read.csv("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Cleaned/Data.csv", header = T)
d$Caste1 <- as.numeric(as.character(d$Caste1))
data <- d

#Tree file - species. This tree is in NEWICK format
anttree_species <- read.tree(file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Trees/Nelsen/ML_TREE_treepl_185_outsdroppeds.ladderized.increasing.tre")

#Make tree ultrametric
u.anttree_species <- chronoMPL(anttree_species)
is.ultrametric(u.anttree_species)

##Ensure that the branch lengths between the orginial and final trees are the same
plot(anttree_species$edge.length, u.anttree_species$edge.length, xlab="Original tree", ylab="Post-chronos tree", pch=20, bty="n")
abline(a=0, b=1, col="gray", lwd=0.5) #This works

#Write the ultrametric tree to file as NEWICK tree
#write.tree(u.anttree_species, file='ultrametric_Nelsen_sp.tre')

#########


#Caste vs. MF filtering
antdata_MF <- filter(data, type == 'ant', Caste1 >=1, eff.mating.freq.MEAN.harmonic >=1)

#Prune tree
pruned.tree_sp<-drop.tip(anttree_species, setdiff(anttree_species$tip.label, antdata_MF$animal))
pruned.tree_sp
plotTree(pruned.tree_sp,ftype="i",fsize=0.4,lwd=1)

#Filter through my dataframe and select only the rows that match the tips of my tree
antdata_MF.1<-filter(antdata_MF, animal %in% pruned.tree_sp$tip.label)
View(antdata_MF.1)

#Configure dataframe into correct formatting
antdata_MF.2 <- cbind(antdata_MF.1$animal, antdata_MF.1)
##Remove column name for 'animal'
antdata_MF.3 <-antdata_MF.2 %>% remove_rownames() %>% column_to_rownames(var="animal")
###Rename column
antdata_MF.4 <- antdata_MF.3 %>% 
  rename(
    animal = `antdata_MF.1$animal`
  )

plot(antdata_MF.4$eff.mating.freq.MEAN.harmonic, antdata_MF.4$Caste1)
abline(CasteVs.MF_Phy_lm)

CasteVs.MF_Phy_lm <- phylolm(Caste1~eff.mating.freq.MEAN.harmonic,data = antdata_MF.4,phy = pruned.tree_sp)
CasteVs.MF_Phy_glm <- phyloglm(Caste1~eff.mating.freq.MEAN.harmonic,data = antdata_MF.4, phy = pruned.tree_sp, method = "poisson_GEE")
CasteVs.MF_lm <- lm(Caste1 ~ eff.mating.freq.MEAN.harmonic, data = antdata_MF.4)

summary(CasteVs.MF_Phy_lm)
summary(CasteVs.MF_Phy_glm)
summary(CasteVs.MF_lm)

#Assumption diagnostics
par(mfrow=c(2,2))
plot(density(CasteVs.MF_Phy_lm$residuals))
qqnorm(CasteVs.MF_Phy_lm$residuals); qqline(CasteVs.MF_Phy_lm$residuals)
plot(CasteVs.MF_Phy_lm$fitted.values, CasteVs.MF_Phy_lm$residuals)
plot(CasteVs.MF_Phy_lm)
#########




#Caste vs. CS filtering
antdata_CS <- filter(data, type == 'ant', Caste1 >=1, colony.size >=1)

#Prune tree
pruned.tree_sp_CS<-drop.tip(anttree_species, setdiff(anttree_species$tip.label, antdata_CS$animal))
pruned.tree_sp_CS
plotTree(pruned.tree_sp_CS,ftype="i",fsize=0.4,lwd=1)

#Filter through my dataframe and select only the rows that match the tips of my tree
antdata_CS.1<-filter(antdata_CS, animal %in% pruned.tree_sp_CS$tip.label)
View(antdata_CS.1)

#Configure dataframe into correct formatting
antdata_CS.2 <- cbind(antdata_CS.1$animal, antdata_CS.1)
##Remove column name for 'animal'
antdata_CS.3 <-antdata_CS.2 %>% remove_rownames() %>% column_to_rownames(var="animal")
###Rename column
antdata_CS.4 <- antdata_CS.3 %>% 
  rename(
    animal = `antdata_CS.1$animal`
  )

plot(log(antdata_CS.4$colony.size), antdata_CS.4$Caste1)
abline(CasteVs.CS_lm)

CasteVs.CS_Phy_lm <- phylolm(formula = Caste1~colony.size,data = antdata_CS.4,phy = pruned.tree_sp_CS)
CasteVs.CS_Phy_glm <- phyloglm(Caste1~colony.size,data = antdata_CS.4, phy = pruned.tree_sp_CS, method = "poisson_GEE")
CasteVs.CS_lm <- lm(Caste1 ~ colony.size, data = antdata_CS.4)

summary(CasteVs.CS_Phy_lm)
summary(CasteVs.CS_Phy_glm)
summary(CasteVs.CS_lm)


#########





#Multiple regression
#CS and MF filtering
antdata_CS_MF <- filter(data, type == 'ant', Caste1 >=1, colony.size >=1, eff.mating.freq.MEAN.harmonic >=1)

#Prune tree
pruned.tree_sp_CS_MF<-drop.tip(anttree_species, setdiff(anttree_species$tip.label, antdata_CS_MF$animal))
pruned.tree_sp_CS_MF
plotTree(pruned.tree_sp_CS_MF,ftype="i",fsize=0.4,lwd=1)

#Filter through my dataframe and select only the rows that match the tips of my tree
antdata_CS_MF.1<-filter(antdata_CS_MF, animal %in% pruned.tree_sp_CS_MF$tip.label)
View(antdata_CS_MF.1)

#Configure dataframe into correct formatting
antdata_CS_MF.2 <- cbind(antdata_CS_MF.1$animal, antdata_CS_MF.1)
##Remove column name for 'animal'
antdata_CS_MF.3 <-antdata_CS_MF.2 %>% remove_rownames() %>% column_to_rownames(var="animal")
###Rename column
antdata_CS_MF.4 <- antdata_CS_MF.3 %>% 
  rename(
    animal = `antdata_CS_MF.1$animal`
  )


CasteVs.CS_MF_Phy_lm <- phylolm(formula = Caste1~colony.size+eff.mating.freq.MEAN.harmonic,data = antdata_CS_MF.4,phy = pruned.tree_sp_CS_MF)
CasteVs.CS_MF_Phy_glm <- phyloglm(Caste1~colony.size*eff.mating.freq.MEAN.harmonic,data = antdata_CS_MF.4, phy = pruned.tree_sp_CS_MF, method = "poisson_GEE")
CasteVs.CS_MF_lm <- lm(Caste1 ~ colony.size+eff.mating.freq.MEAN.harmonic, data = antdata_CS_MF.4)

summary(CasteVs.CS_MF_Phy_lm)
summary(CasteVs.CS_MF_Phy_glm)
summary(CasteVs.CS_MF_lm)











