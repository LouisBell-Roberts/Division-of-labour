# BEE conference: phylogenetically controlled regressions for Caste vs. MF or CS
#Louis Bell-Roberts
#11/03/2021


#*****************************************************************

setwd("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Trees/Nelsen_ultrametric_species/") 
#I think this folder is not embedded within the Github, is it? Might be useful to also include it? 

library(tidyverse)
library(ape)
library(phylolm)
library(phytools)
library(ggplot2)

#Read in data file - ensure that caste variable is set a numeric variable
d <- read.csv("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Primary Dataset/Data_caste_mating_colonyS_WPM_QueenN_cleaned.csv", header = T)
#d <- read.csv(file.choose(), header = T)

d <- read.csv("./Data/Cleaned/Data_caste_mating_colonyS_WPM_QueenN_cleaned.csv", header = T) #This one works on Gijsbert Computer

d$Caste3 <- as.numeric(as.character(d$Caste3))
data <- d

#Tree files - 'Genus_polytomy_tree.tre' is the most up-to-date tree
#anttree_species <- read.tree(file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Trees/Nelsen_ultrametric_species/ultrametric_Nelsen_sp.tre")
anttree_species <- read.tree(file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Trees/Polytomy_tree/Genus_polytomy_tree.tre")
#anttree_species <- read.tree(file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Trees/Polytomy_tree/Polytomy_tree.tre")
#anttree_species <- read.tree(file.choose())

#########
#PGLS analysis - Caste ~ Mating frequency

#Caste vs. MF filtering - selecting only the ant species from the database along with species that have data on caste number and mating frequency
antdata_MF <- filter(data, type == 'ant', Caste3 >=1, eff.mating.freq.MEAN.harmonic >=1)

#Prune tree to match the species in my database - 130 species in pruned tree
pruned.tree_sp<-drop.tip(anttree_species, setdiff(anttree_species$tip.label, antdata_MF$animal))
pruned.tree_sp
plotTree(pruned.tree_sp,ftype="i",fsize=0.4,lwd=1)

#Prune the database to select only the rows that match the tips of my tree
antdata_MF.1<-filter(antdata_MF, animal %in% pruned.tree_sp$tip.label)
#View(antdata_MF.1)

#Configure dataframe into correct formatting for phylolm - row names must be equal to the species names 
##Duplicate the 'animal' column which contains the combined genus and species names
antdata_MF.2 <- cbind(antdata_MF.1$animal, antdata_MF.1)
##Remove column name for 'animal' - this labels each row as the species name
antdata_MF.3 <-antdata_MF.2 %>% remove_rownames() %>% column_to_rownames(var="animal")
###Rename one of the duplicate columns to 'animal'
antdata_MF.4 <- antdata_MF.3 %>% 
  rename(
    animal = `antdata_MF.1$animal`
  )


#Models: 1st phylogenetic general linear model, 2nd phylogenetic poisson regression, 3rd non-phylogenetic general linear model
CasteVs.MF_Phy_lm <- phylolm(Caste3~eff.mating.freq.MEAN.harmonic,data = antdata_MF.4,phy = pruned.tree_sp)
CasteVs.MF_Phy_glm <- phyloglm(Caste3~eff.mating.freq.MEAN.harmonic,data = antdata_MF.4, phy = pruned.tree_sp, method = "poisson_GEE")
CasteVs.MF_lm <- lm(Caste3 ~ eff.mating.freq.MEAN.harmonic, data = antdata_MF.4)
CasteVs.MF_lm.1 <- lm((Caste3)^2 ~ eff.mating.freq.MEAN.harmonic, data = antdata_MF.4)

summary(CasteVs.MF_Phy_lm)
summary(CasteVs.MF_Phy_glm)
summary(CasteVs.MF_lm)

#Plot the PGLS regression lines onto the data and the non-phylogenetic regression
plot(antdata_MF.4$eff.mating.freq.MEAN.harmonic, antdata_MF.4$Caste3, ylim=c(-2,10))
abline(CasteVs.MF_Phy_lm)
abline(CasteVs.MF_Phy_glm)
abline(CasteVs.MF_lm)

#*****************************************************************

############ Assumption diagnostics for the Caste ~ Mating frequency PGLS model

#par(mfrow=c(2,2))
## Residuals vs. fitted values plot - check for linearity
plot(CasteVs.MF_Phy_lm$fitted.values, CasteVs.MF_Phy_lm$residuals)

## Normal Q-Q plot - check for assumption of normality
mod_resids <- CasteVs.MF_Phy_lm$residuals
mod_resids <- mod_resids/sd(mod_resids)
resid_order <- order(mod_resids)
all_resids <- qqnorm(mod_resids, plot.it = FALSE)
all_resids <- as.data.frame(all_resids)
head(all_resids, 10)
ggplot(all_resids, aes(x = x, y = y)) + 
  geom_point() + geom_abline(intercept = 0, slope = 1) +
  xlab("Theoretical Value") + ylab("Standardised Residual")
####

## Scale location plot - Check for constant variance
# extract the residuals
sqrt_abs_resids <- CasteVs.MF_Phy_lm$residuals
# step 1. standardise them
sqrt_abs_resids <- sqrt_abs_resids / sd(sqrt_abs_resids)
# step 2. find their absolute value
sqrt_abs_resids <- abs(sqrt_abs_resids)
# step 3. square root these
sqrt_abs_resids <- sqrt(sqrt_abs_resids)

plt_data <- 
  data.frame(Fitted = CasteVs.MF_Phy_lm$fitted.values, Resids = sqrt_abs_resids)

ggplot(plt_data, aes(x = Fitted, y = Resids)) + 
  geom_point() + 
  xlab("Fitted values") + ylab("Square root of absolute residuals")
#####

# Excess assumption diagnostics code
#plot(density(CasteVs.MF_Phy_lm$residuals))
#qqnorm(CasteVs.MF_Phy_lm$residuals); qqline(CasteVs.MF_Phy_lm$residuals) # check for normality (code below does a better job)
#plot(CasteVs.MF_Phy_lm)
#

#ggplot2 - plotting raw data for MF and caste
ggplot(antdata_MF.4, aes(x = eff.mating.freq.MEAN.harmonic, y = Caste3)) + 
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
#PGLS analysis - Caste ~ Colony size

#Caste vs. CS filtering - selecting only the ant species from the database along with species that have data on caste number and colony size
antdata_CS <- filter(data, type == 'ant', Caste3 >=1, colony.size >=1)

#Prune tree - 520 species
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
CasteVs.CS_Phy_lm <- phylolm(formula = Caste3~log(colony.size),data = antdata_CS.4,phy = pruned.tree_sp_CS)
CasteVs.CS_Phy_glm <- phyloglm(Caste3~log(colony.size),data = antdata_CS.4, phy = pruned.tree_sp_CS, method = "poisson_GEE")
CasteVs.CS_lm <- lm(Caste3 ~ log(colony.size), data = antdata_CS.4)

summary(CasteVs.CS_Phy_lm)
summary(CasteVs.CS_Phy_glm)
summary(CasteVs.CS_lm)
coef(CasteVs.CS_Phy_lm)
coef(CasteVs.CS_Phy_glm)
coef(CasteVs.CS_lm)

#Plot regressions
plot(log(antdata_CS.4$colony.size), antdata_CS.4$Caste3, ylim = c(-2, 5))
abline(CasteVs.CS_Phy_lm)
abline(CasteVs.CS_lm)
abline(CasteVs.CS_Phy_glm)



#Diagnostics
par(mfrow=c(2,2))
plot(density(CasteVs.CS_Phy_lm$residuals))
qqnorm(CasteVs.CS_Phy_lm$residuals); qqline(CasteVs.CS_Phy_lm$residuals)
plot(CasteVs.CS_Phy_lm$fitted.values, CasteVs.CS_Phy_lm$residuals)
plot(CasteVs.CS_Phy_lm)



#ggplot2 - plotting raw data for CS and caste
ggplot(antdata_CS.4, aes(x = log(colony.size), y = Caste3)) + 
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
antdata_CS_MF <- filter(data, type == 'ant', Caste3 >=1, colony.size >=1, eff.mating.freq.MEAN.harmonic >=1, polygyny.clean >= 0, Reference.2 != "")
, polygyny.clean >= 0, Reference.2 != ""
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
CasteVs.CS_MF_Phy_lm <- phylolm(formula = Caste3~log(colony.size)*eff.mating.freq.MEAN.harmonic,data = antdata_CS_MF.4,phy = pruned.tree_sp_CS_MF)
CasteVs.CS_MF_Phy_glm <- phyloglm(Caste3~log(colony.size)*eff.mating.freq.MEAN.harmonic,data = antdata_CS_MF.4, phy = pruned.tree_sp_CS_MF, method = "poisson_GEE")
CasteVs.CS_MF_PG_Phy_lm <- phyloglm(Caste3~log(colony.size)+eff.mating.freq.MEAN.harmonic+polygyny.clean,data = antdata_CS_MF.4, 
                                     phy = pruned.tree_sp_CS_MF)
CasteVs.CS_MF_PG_Phy_glm <- phyloglm(Caste3~log(colony.size)+eff.mating.freq.MEAN.harmonic+polygyny.clean+eff.mating.freq.MEAN.harmonic:log(colony.size)+
                                       eff.mating.freq.MEAN.harmonic:polygyny.clean,
                                     data = antdata_CS_MF.4, 
                                     phy = pruned.tree_sp_CS_MF, method = "poisson_GEE")
CasteVs.CS_MF_PG_Phy_glm <- phyloglm(Caste3~log(colony.size)*eff.mating.freq.MEAN.harmonic+polygyny.clean,
                                     data = antdata_CS_MF.4, 
                                     phy = pruned.tree_sp_CS_MF, method = "poisson_GEE")

CasteVs.CS_MF_PG_Phy_glm <- phyloglm(Caste3~polygyny.clean*log(colony.size)+eff.mating.freq.MEAN.harmonic*log(colony.size),
                                     data = antdata_CS_MF.4, 
                                     phy = pruned.tree_sp_CS_MF, method = "poisson_GEE")
CasteVs.CS_MF_PG_Phy_glm <- phyloglm(Caste3~polygyny.clean,
                                     data = antdata_CS_MF.4, 
                                     phy = pruned.tree_sp_CS_MF, method = "poisson_GEE")

summary(CasteVs.CS_MF_Phy_lm)
summary(CasteVs.CS_MF_Phy_glm)
summary(CasteVs.CS_MF_lm)
summary(CasteVs.CS_MF_PG_Phy_glm)











