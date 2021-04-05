# BEE conference: non-phylogenetically controlled regressions for Caste vs. MF or CS
#Louis Bell-Roberts
#10/03/2021

library(tidyverse)

d <- read.csv("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Primary Dataset/Data_caste_mating_colonyS_WPM_QueenN_cleaned.csv", header = T)
#d <- read.csv(file.choose(), header = T)
d$Caste1 <- as.numeric(as.character(d$Caste1))
data <- d

#Caste vs. MF filtering
antdata_MF <- filter(data, type == 'ant', Caste1 >=1, eff.mating.freq.MEAN.harmonic >=1)

#Plot MF against Caste
plot(antdata_MF$eff.mating.freq.MEAN.harmonic, antdata_MF$Caste1)

#Create linear model for Caste vs. MF
CasteVs.MF_lm <- lm(Caste1 ~ eff.mating.freq.MEAN.harmonic, data = antdata_MF)
summary(CasteVs.MF_lm) #R-squared is 0.12 and gradient is significant

#Base R plot
plot(Caste1 ~ eff.mating.freq.MEAN.harmonic, antdata_MF, xlim=c(0,30), ylim=c(0,5))
abline(CasteVs.MF_lm)

#ggplot2
ggplot(antdata_MF, aes(x = eff.mating.freq.MEAN.harmonic, y = Caste1)) + 
  geom_point() +
  stat_smooth(method = "lm", col = "black", se = F) + 
  theme_classic() +
  ggtitle("Caste Number vs. Mating Frequency") +
  xlab("Mating Frequency") +
  ylab("Caste Number") +
  theme(plot.title = element_text(hjust = 0.5))

#Create generalised linear model with poisson distribution family 
CasteVs.MF_glm <- glm(Caste1 ~ eff.mating.freq.MEAN.harmonic, family = "poisson", data = antdata_MF)
summary(CasteVs.MF_glm)

#Base R plot
plot(Caste1 ~ eff.mating.freq.MEAN.harmonic, antdata_MF, xlim=c(0,30), ylim=c(0,5))
abline(CasteVs.MF_glm)

#ggplot2
#linear model
ggplot(antdata_MF, aes(x = eff.mating.freq.MEAN.harmonic, y = Caste1)) + 
  geom_point() +
  stat_smooth(method = "lm", col = "red")

#Draws line through all points
ggplot(data = antdata_MF, aes(x = eff.mating.freq.MEAN.harmonic, y = Caste1)) +
  geom_jitter(width = 0.05, height = 0.05) +
geom_smooth()

#Plotting poisson regression
ggplot(data = antdata_MF, aes(x = eff.mating.freq.MEAN.harmonic, y = Caste1)) +
  geom_jitter(width = 0.05, height = 0.05) +
  geom_smooth(method = 'glm', method.args = list(family = 'poisson')) +
  theme_classic() +
  ggtitle("Caste Number vs. Mating Frequency") +
  xlab("Mating Frequency") +
ylab("Caste Number") +
theme(plot.title = element_text(hjust = 0.5))





####################





#Caste vs. CS filtering
antdata_CS <- filter(data, type == 'ant', Caste1 >=1, colony.size >=1)
#View(antdata_CS)

#Plot CS against Caste
plot(antdata_CS$colony.size, antdata_CS$Caste1)
plot(log(antdata_CS$colony.size), antdata_CS$Caste1)

#Create linear model for Caste vs. Cs
CasteVs.CS_lm <- lm(Caste1 ~ log(colony.size), data = antdata_CS)
CasteVs.CS_lm <- lm(Caste1 ~ colony.size, data = antdata_CS)
summary(CasteVs.CS_lm) #R-squared is 0.21 when log(CS)

#Base R plot
plot(Caste1 ~ log(colony.size), antdata_CS, xlim=c(0,20), ylim=c(0,5))
abline(CasteVs.CS_lm)

#ggplot2
ggplot(antdata_CS, aes(x = log(colony.size), y = Caste1)) + 
  geom_point() +
  stat_smooth(method = "lm", col = "black", se = F) + 
  theme_classic() +
  ggtitle("Caste Number vs. log(Colony Size)") +
  xlab("Log number of ants in a colony") +
  ylab("Number of physical castes in a colony") +
  theme(plot.title = element_text(hjust = 0.5))

#Create generalised linear model with poisson distribution family 
CasteVs.CS_glm <- glm(Caste1 ~ log(colony.size), family = "poisson", data = antdata_CS)
CasteVs.CS_glm <- glm(Caste1 ~ colony.size, family = "poisson", data = antdata_CS)
summary(CasteVs.CS_glm)

#Base R plot
plot(Caste1 ~ log(colony.size), antdata_CS, xlim=c(0,30), ylim=c(0,5))
abline(CasteVs.CS_glm)

#ggplot2
#linear model
ggplot(antdata_CS, aes(x = log(colony.size), y = Caste1)) + 
  geom_point() +
  stat_smooth(method = "lm", col = "red")

#Draws line through all points
ggplot(data = antdata_CS, aes(x = log(colony.size), y = Caste1)) +
  geom_jitter(width = 0.05, height = 0.05) +
  geom_smooth()

#Plotting poisson regression
ggplot(data = antdata_CS, aes(x = log(colony.size), y = Caste1)) +
  geom_jitter(width = 0.05, height = 0.05) +
  geom_smooth(method = 'glm', method.args = list(family = 'poisson')) +
  theme_classic() +
ggtitle("Caste Number vs. log(Colony Size)") +
  xlab("log(Colony Size)") +
  ylab("Caste Number") +
  theme(plot.title = element_text(hjust = 0.5))














####Move all of this stuff over to another script for phylogenetically controlled analyses
setwd("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Trees/Nelsen_ultrametric_species/")

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
write.tree(u.anttree_species, file='ultrametric_Nelsen_sp.tre')

