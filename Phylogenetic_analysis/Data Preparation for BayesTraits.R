#Preparing my data for analysis in Bayes Traits#
#25/02/2021 Louis Bell-Roberts

#Set the working directory
setwd("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/")

#Installing BayesTraits R-wrapper
library("devtools")
#install_github("rgriff23/btw")
library("btw")

library(phangorn)
library(nlme)
library(ape)
library(MCMCglmm)
library(phytools)
library(geiger)
library(phylolm)
library(corHMM)
library(caper)
library(stringr)
library(dplyr)

#To run analyses in BayesTraits we need two types of file: 1) a tree file containing a phylogenetic tree, or sample of trees, 
#2) a data file containing the trait data to be analyzed. 
#BayesTraits requires input phylogenetic trees to be in NEXUS format (i.e. it should have a #NEXUS tag at the top of the file)
#There are no column headings and missing data needs to be represented by a hyphen ( - ) - should I even include missing data in my analysis practice?

sf<-read.table(file.choose(),header= FALSE)
View(sf)
phylosf<-read.nexus(file.choose())
dfjk<-read.nexus(file.choose())
plot(dfjk)


#Make tree ultrametric
u.anttree <- chronoMPL(anttree)
is.ultrametric(u.anttree)
##Ensure that the branch lengths between the orginial and final trees are the same
plot(anttree$edge.length, u.anttree$edge.length, xlab="Original tree", ylab="Post-chronos tree", pch=20, bty="n")
abline(a=0, b=1, col="gray", lwd=0.5) #This works

##Read in the database##
antdata <- read.csv("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Cleaned/Data.csv", header = T, stringsAsFactors = FALSE)
head(antdata)

databayes<-filter(antdata, type=="ant")
View(databayes)

#Select traits: Caste and MF
databayes.1<- databayes %>% select(animal, eff.mating.freq.MEAN.harmonic, Caste1)
View(databayes.1)

#Make data binary
#databayes.4 <- databayes.3 %>% mutate(Harmonic_cat=cut(eff.mating.freq.MEAN.harmonic, breaks=c(-Inf, 1, 2, Inf), labels=c("Monogamy","low","high")))
databayes.1$Caste1 <- as.numeric((as.character(databayes.1$Caste1)))

databayes.1.1 <- databayes.1 %>% mutate(Harmonic_cat=cut(eff.mating.freq.MEAN.harmonic, breaks=c(-Inf, 2, Inf), labels=c("0","1")))
databayes.1.2 <- databayes.1.1 %>% mutate(Caste_cat=cut(Caste1, breaks=c(-Inf, 1.9, Inf), labels=c("0","1")))
View(databayes.1.2)
is.factor(databayes.1.2$Caste_cat)

#Remove all species where we don't have data for at least one of the traits
databayes.2 <- filter(databayes.1.2, Caste1>=0 | eff.mating.freq.MEAN.harmonic>=0)
View(databayes.2)

is.numeric(databayes.2$Caste_cat)
databayes.2$Caste_cat <- as.numeric((as.character(databayes.2$Caste_cat)))
databayes.2$Harmonic_cat <- as.numeric((as.character(databayes.2$Harmonic_cat)))
is.numeric(databayes.2$eff.mating.freq.MEAN.harmonic)

#Replace empty cells with NA's
#databayes.2<-databayes.1.2 %>%
#  mutate(Caste1 = na_if(Caste1, ""))
#View(databayes.2)

#Replace NA's with '-'
databayes.3 <- databayes.2 %>% replace_na(list(Harmonic_cat = "-", Caste_cat = "-"))
View(databayes.3)


#Select traits: Caste_cat and Harmonic_cat
databayes.4<- databayes.3 %>% select(animal, Harmonic_cat, Caste_cat)
View(databayes.4)

#Remove column names
databayes.5 <- databayes.4
colnames(databayes.5) <- NULL
View(databayes.5)

#Write database to .txt file
##Ensure that parameters in 'write.table' are set to remove quotation marks, row names and column names
write.table(databayes.5, file = "ants.txt", sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)









#Preparing tree for BayesTraits
##Read in species-level tree -> prune it to the database -> (make it ultrametric - optional) -> write out as nexus file
anttree_sp <- read.tree(file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Trees/Nelsen/ML_TREE_treepl_185_outsdroppeds.ladderized.increasing.tre")

#Prune tree
pruned.tree<-drop.tip(anttree_sp, setdiff(anttree_sp$tip.label, databayes.4$animal))
pruned.tree

#Match database to the tree
databayes.4.1<-filter(databayes.4, animal %in% pruned.tree$tip.label)
View(databayes.4.1)

#Remove column names
databayes.4.2 <- databayes.4.1
colnames(databayes.4.2) <- NULL
View(databayes.4.2)

#Write database to .txt file
##Ensure that parameters in 'write.table' are set to remove quotation marks, row names and column names
write.table(databayes.4.2, file = "ants.txt", sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

#Writing my tree into a nexus format - do this at the end
write.nexus(pruned.tree, file="antsMCMC.tre", translate = T)



#Using BayesTraits in R
setwd("~/Documents/DTP_1st_project_rotation/BayesTraits/BayesTraitsV3.0.2-OSX/") 
data(primates)

#Primate data
primates$tree
View(primates$trait)

#My data
databayes.5
anttree <- read.tree(file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Trees/Nelsen/Dryad_Supplementary_File_7_ML_TREE_treepl_185.tre")
write.nexus(anttree, file="antsMCMC2.tre", translate = T)

is.integer(primates$trait$T1)
is.character(primates$trait$T1[1])

command_vec1 <- c("1", "1") #for Multistate model in ML mode

results_1 <- bayestraits(primate.discrete1, primate.tree1, command_vec1, silent = F)
results_1$Log

results_2 <- bayestraits(databayes.5, )





####Testing convergence in my models####














