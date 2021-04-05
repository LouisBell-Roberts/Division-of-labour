# BEE conference: phylogenetically controlled regressions for Caste vs. MF or CS
#Louis Bell-Roberts
#11/03/2021

setwd("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Trees/Nelsen_ultrametric_species/")

library(tidyverse)
library(ape)
library(phylolm)

#Data file
d <- read.csv("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Primary Dataset/Data_caste_mating_colonyS_WPM_QueenN_cleaned.csv", header = T)
#d <- read.csv(file.choose(), header = T)
d$Caste1 <- as.numeric(as.character(d$Caste1))
data <- d

#Tree file - species
anttree_species <- read.tree(file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Trees/Nelsen_ultrametric_species/ultrametric_Nelsen_sp.tre")
#anttree_species <- read.tree(file.choose())

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

#Models
CasteVs.MF_Phy_lm <- phylolm(Caste1~eff.mating.freq.MEAN.harmonic,data = antdata_MF.4,phy = pruned.tree_sp)
CasteVs.MF_Phy_glm <- phyloglm(Caste1~eff.mating.freq.MEAN.harmonic,data = antdata_MF.4, phy = pruned.tree_sp, method = "poisson_GEE")
CasteVs.MF_lm <- lm(Caste1 ~ eff.mating.freq.MEAN.harmonic, data = antdata_MF.4)

summary(CasteVs.MF_Phy_lm)
summary(CasteVs.MF_Phy_glm)
summary(CasteVs.MF_lm)

#Plot regressions
plot(antdata_MF.4$eff.mating.freq.MEAN.harmonic, antdata_MF.4$Caste1)
abline(CasteVs.MF_Phy_lm)
abline(CasteVs.MF_lm)


#Assumption diagnostics
par(mfrow=c(2,2))
plot(density(CasteVs.MF_Phy_lm$residuals))
qqnorm(CasteVs.MF_Phy_lm$residuals); qqline(CasteVs.MF_Phy_lm$residuals)
plot(CasteVs.MF_Phy_lm$fitted.values, CasteVs.MF_Phy_lm$residuals)
plot(CasteVs.MF_Phy_lm)




#ggplot2 - plotting raw data for MF and caste
ggplot(antdata_MF.4, aes(x = eff.mating.freq.MEAN.harmonic, y = Caste1)) + 
  geom_jitter(width = 0.05, height = 0.05, size = 2) +
  theme_classic() +
  ggtitle("Caste Number vs. Mating Frequency") +
  xlab("Mating Frequency") +
  ylab("Caste Number") +
  theme(plot.title = element_text(hjust = 0.5, size = 22, face = "bold"), axis.line = element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = "black", size = 1), 
        axis.text.x = element_text(face="bold", color="black", 
                                   size=22),
        axis.text.y = element_text(face="bold", color = "black", 
                                   size=22),
        axis.text=element_text(size=22),
        axis.title=element_text(size=22,face="bold"))




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

#Models
CasteVs.CS_Phy_lm <- phylolm(formula = Caste1~log(colony.size),data = antdata_CS.4,phy = pruned.tree_sp_CS)
CasteVs.CS_Phy_glm <- phyloglm(Caste1~log(colony.size),data = antdata_CS.4, phy = pruned.tree_sp_CS, method = "poisson_GEE")
CasteVs.CS_lm <- lm(Caste1 ~ log(colony.size), data = antdata_CS.4)

summary(CasteVs.CS_Phy_lm)
summary(CasteVs.CS_Phy_glm)
summary(CasteVs.CS_lm)
coef(CasteVs.CS_Phy_lm)
coef(CasteVs.CS_Phy_glm)
coef(CasteVs.CS_lm)

#Plot regressions
plot(log(antdata_CS.4$colony.size), antdata_CS.4$Caste1)
abline(CasteVs.CS_Phy_lm)
abline(CasteVs.CS_lm)



#Diagnostics
par(mfrow=c(2,2))
plot(density(CasteVs.CS_Phy_lm$residuals))
qqnorm(CasteVs.CS_Phy_lm$residuals); qqline(CasteVs.CS_Phy_lm$residuals)
plot(CasteVs.CS_Phy_lm$fitted.values, CasteVs.CS_Phy_lm$residuals)
plot(CasteVs.CS_Phy_lm)



#ggplot2 - plotting raw data for CS and caste
ggplot(antdata_CS.4, aes(x = log(colony.size), y = Caste1)) + 
  geom_point(size = 3) +
  theme_classic() +
  ggtitle("Caste Number vs. log Colony Size") +
  xlab("log Colony Size") +
  ylab("Caste Number") +
  theme(plot.title = element_text(hjust = 0.5, size = 22, face = "bold"), axis.line = element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = "black", size = 1), 
        axis.text.x = element_text(face="bold", color="black", 
                                   size=22),
        axis.text.y = element_text(face="bold", color = "black", 
                                   size=22),
        axis.text=element_text(size=22),
        axis.title=element_text(size=22,face="bold"))

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
#View(antdata_CS_MF.1)

#Configure dataframe into correct formatting
antdata_CS_MF.2 <- cbind(antdata_CS_MF.1$animal, antdata_CS_MF.1)
##Remove column name for 'animal'
antdata_CS_MF.3 <-antdata_CS_MF.2 %>% remove_rownames() %>% column_to_rownames(var="animal")
###Rename column
antdata_CS_MF.4 <- antdata_CS_MF.3 %>% 
  rename(
    animal = `antdata_CS_MF.1$animal`
  )

#KICKS UP AN ERROR IF I DON'T LOG TRANSFORM COLONY SIZE?
CasteVs.CS_MF_Phy_lm <- phylolm(formula = Caste1~colony.size*eff.mating.freq.MEAN.harmonic,data = antdata_CS_MF.4,phy = pruned.tree_sp_CS_MF)
CasteVs.CS_MF_Phy_glm <- phyloglm(Caste1~log(colony.size)*eff.mating.freq.MEAN.harmonic,data = antdata_CS_MF.4, phy = pruned.tree_sp_CS_MF, method = "poisson_GEE")
CasteVs.CS_MF_lm <- lm(Caste1 ~ colony.size*eff.mating.freq.MEAN.harmonic, data = antdata_CS_MF.4)

summary(CasteVs.CS_MF_Phy_lm)
summary(CasteVs.CS_MF_Phy_glm)
summary(CasteVs.CS_MF_lm)











