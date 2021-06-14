#Advanced modelling course - start with cleaned up data
#Louis Bell-Roberts
#10/04/2021
#Variables of interest: CS, MF (continuous/categorical), QN (continuous/binary), Caste

setwd("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Trees/Nelsen_ultrametric_species/")

library(tidyverse)
library(ape)
library(phytools)
library(geiger)
library(phylolm)

#Read in data
d <- read.csv("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Primary Dataset/Data_caste_mating_colonyS_WPM_QueenN_cleaned.csv", header = T)
#d <- read.csv(file.choose(), header = T)
d$Caste3 <- as.numeric(as.character(d$Caste3))
d$number.queens.MEAN <- as.numeric(as.character(d$number.queens.MEAN))
data <- d

#Read in the tree file
anttree_species <- read.tree(file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Trees/Polytomy_tree/Genus_polytomy_tree.tre")


#Select only ants
antdata<-filter(data, type=="ant")


#Create MF categorical variable 'combined' based on polyandry.clean and eff.mating.freq.MEAN.harmonic
##Create bins for continuous MF variable
antdata<-antdata%>%mutate(MFbins = cut(eff.mating.freq.MEAN.harmonic, breaks = c(-Inf,1,2,Inf), c("0", "1", "2")))
#View(select(antdata, eff.mating.freq.MEAN.harmonic, MFbins))
##Make polyandry.clean a factor
antdata$polyandry.clean<-as.factor(antdata$polyandry.clean)
##Combine the two columns. MFbins takes priority and if NA, then value taken from polyandry.clean
antdata$MFcombined<-coalesce(antdata$MFbins, antdata$polyandry.clean)
View(select(antdata, animal, polyandry.clean, MFbins, MFcombined))

#Cleaning the PG columns
#Using ifelse - It takes three arguments - the condition, the if output, and the else output.
antdata$PG <- ifelse(antdata$polygyny.clean >=0 & antdata$Reference.2 != "", antdata$polygyny.clean, NA)
#View(select(antdata, animal, polygyny.clean, Reference.2, PG))
antdata$PG_cont <- ifelse(antdata$number.queens.MEAN >=0 & antdata$Reference.2 != "", antdata$number.queens.MEAN, NA)
#View(select(antdata, animal, number.queens.MEAN, Reference.2, PG_cont))

#Remove outlier species: Formica_yessensis, Lasius_neglectus, Pseudomyrmex_veneficus
antdata.1 <- filter(antdata, animal != "Formica_yessensis", animal != "Lasius_neglectus", animal != "Pseudomyrmex_veneficus")

#Select only the columns of interest
##antclean is the cleaned dataset
antclean <- dplyr::select(antdata.1, animal, colony.size, eff.mating.freq.MEAN.harmonic, MFcombined, PG, PG_cont, Caste3)
View(antclean)
str(antclean$Caste3)













