# BEE conference: phylogenetically controlled regressions for Caste vs. MF or CS
## Using gls function rather than phylolm
#Louis Bell-Roberts
#17/03/2021

setwd("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Trees/Nelsen_ultrametric_species/")

library(ape)
library(geiger)
library(nlme)
library(phytools)
library(tidyverse)

#Data file
#otherd<-read.csv("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Cleaned/Barbetdata.csv", header = T)

d <- read.csv("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Primary Dataset/Data_caste_mating_colonyS_WPM_QueenN_cleaned.csv", header = T)
d$Caste1 <- as.numeric(as.character(d$Caste1))
data <- d

#Tree file - species
anttree_species <- read.tree(file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Trees/Nelsen_ultrametric_species/ultrametric_Nelsen_sp.tre")

#########
#Caste vs. MF filtering
antdata_MF <- filter(data, type == 'ant', Caste1 >=1, eff.mating.freq.MEAN.harmonic >=1)

#Prune tree
pruned.tree_sp<-drop.tip(anttree_species, setdiff(anttree_species$tip.label, antdata_MF$animal))
pruned.tree_sp
plotTree(pruned.tree_sp,ftype="i",fsize=0.4,lwd=1)

#Filter through my dataframe and select only the rows that match the tips of my tree
antdata_MF.1<-filter(antdata_MF, animal %in% pruned.tree_sp$tip.label)
#View(antdata_MF.1)

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

#PGLS model
pglsModel <- gls(Caste1 ~ eff.mating.freq.MEAN.harmonic, correlation = corBrownian(1, phy = pruned.tree_sp),
                 data = antdata_MF.4, method = "ML")
summary(pglsModel)
coef(pglsModel)

#Plot regression
plot(antdata_MF.4$Caste1 ~ antdata_MF.4$eff.mating.freq.MEAN.harmonic)
#abline(a = coef(pglsModel)[1], b = coef(pglsModel)[2])
abline(pglsModel)

##With caper
antdata_MF.5 <- dplyr::select(antdata_MF.1, animal, eff.mating.freq.MEAN.harmonic, Caste1)
antdata_MF.6 <- antdata_MF.5 %>% 
  rename(
    Species = animal
  )
pruned.tree_sp

#Model
comp.data<-comparative.data(pruned.tree_sp, antdata_MF.6, names.col="Species", vcv.dim=2, warn.dropped=TRUE)
modelo4<-pgls(Caste1~eff.mating.freq.MEAN.harmonic, data=comp.data)
summary(modelo4)
coef(modelo4)

#Plot regression
plot(antdata_MF.6$Caste1 ~ antdata_MF.6$eff.mating.freq.MEAN.harmonic)
abline(modelo4)

coef(pglsModel)
coef(modelo4)
coef(CasteVs.MF_Phy_lm)
coef(CasteVs.MF_Phy_glm)
coef(CasteVs.MF_lm)
###############



#########

#Caste vs. CS filtering
antdata_CS <- filter(data, type == 'ant', Caste1 >=1, colony.size >=1)

#Prune tree
pruned.tree_sp_CS<-drop.tip(anttree_species, setdiff(anttree_species$tip.label, antdata_CS$animal))
pruned.tree_sp_CS
plotTree(pruned.tree_sp_CS,ftype="i",fsize=0.4,lwd=1)

#Filter through my dataframe and select only the rows that match the tips of my tree
antdata_CS.1<-filter(antdata_CS, animal %in% pruned.tree_sp_CS$tip.label)
#View(antdata_CS.1)

#Configure dataframe into correct formatting
antdata_CS.2 <- cbind(antdata_CS.1$animal, antdata_CS.1)
##Remove column name for 'animal'
antdata_CS.3 <-antdata_CS.2 %>% remove_rownames() %>% column_to_rownames(var="animal")
###Rename column
antdata_CS.4 <- antdata_CS.3 %>% 
  rename(
    animal = `antdata_CS.1$animal`
  )


#PGLS model
pglsModel2 <- gls(Caste1 ~ log(colony.size), correlation = corBrownian(1, phy = pruned.tree_sp_CS),
                 data = antdata_CS.4, method = "ML")
summary(pglsModel2)
coef(pglsModel2)

#Plot regression
plot(antdata_CS.4$Caste1 ~ log(antdata_CS.4$colony.size))
abline(a = coef(pglsModel2)[1], b = coef(pglsModel2)[2])


##With caper
antdata_CS.5 <- dplyr::select(antdata_CS.1, animal, colony.size, Caste1)
antdata_CS.6 <- antdata_CS.5 %>% 
  rename(
    Species = animal
  )
pruned.tree_sp_CS

#Model
comp.data2<-comparative.data(pruned.tree_sp_CS, antdata_CS.6, names.col="Species", vcv.dim=2, warn.dropped=TRUE)
modelo5<-pgls(Caste1~log(colony.size), data=comp.data2)
summary(modelo5)
coef(modelo5)

#Plot regression
plot(antdata_CS.6$Caste1 ~ log(antdata_CS.6$colony.size))
abline(modelo5)












