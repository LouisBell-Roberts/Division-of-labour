#Data imputation - Rphylopars
#Louis Bell-Roberts
#05/04/2021

setwd("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Trees/")

library(Rphylopars)
library(phytools)
library(tidyverse)
library(ape)
library(geiger)

#Data file
d <- read.csv("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Primary Dataset/Data_caste_mating_colonyS_WPM_QueenN_cleaned.csv", header = T)
d$Caste1 <- as.numeric(as.character(d$Caste1))
data <- d

#Tree file - species
#anttree_species <- read.tree(file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Trees/Nelsen_ultrametric_species/ultrametric_Nelsen_sp.tre")
#anttree_species <- read.tree(file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Trees/Polytomy_tree/Genus_polytomy_tree.tre")
anttree_species <- read.tree(file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Trees/Polytomy_tree/Polytomy_tree.tre")

#Caste, MF and CS filtering
antdata <- filter(data, type == 'ant', (Caste1 >=1 | eff.mating.freq.MEAN.harmonic >=1 | colony.size >= 1))

#Create row names
antdata1 <-antdata %>% remove_rownames() %>% column_to_rownames(var="animal")

#Select caste, MF and CS columns
antdata2 <- antdata1 %>% dplyr::select(eff.mating.freq.MEAN.harmonic, colony.size, Caste1)

#Formatting: to columns for species name (one is a rownames column). Name one of the columns 'species'
antdata3 <- cbind(antdata2, antdata$animal)
antdata4 <- antdata3 %>% 
                    rename(
                      species = `antdata$animal`)
antdata4 <- antdata4 %>% dplyr::select(species, eff.mating.freq.MEAN.harmonic, colony.size, Caste1)

#Match tree to database
#Prune tree
pruned.tree_sp<-drop.tip(anttree_species, setdiff(anttree_species$tip.label, antdata4$species))
pruned.tree_sp
plotTree(pruned.tree_sp,ftype="i",fsize=0.4,lwd=1)

#Filter through my dataframe and select only the rows that match the tips of my tree
antdata5<-filter(antdata4, species %in% pruned.tree_sp$tip.label)

#Check if tree matches database
name.check(pruned.tree_sp, antdata5)

#Run Rphylopars function
p_BM <- phylopars(trait_data = antdata5, tree = pruned.tree_sp)

#Imputed values - SEEMS TO HAVE PRODUCED NEGATIVE COLONY SIZES...
imputed_data_MF <- p_BM$anc_recon[1:983, 1]
imputed_data_CS <- p_BM$anc_recon[1:983, 2]
imputed_data_Caste <- p_BM$anc_recon[1:983, 3]

#Linear models and plots
MFVs.CS_lm <- lm((imputed_data_CS)^0.1 ~ imputed_data_MF)
summary(MFVs.CS_lm)
plot(imputed_data_MF, (imputed_data_CS)^0.1)
abline(MFVs.CS_lm)

MFVs.Caste_lm <- lm(imputed_data_Caste ~ imputed_data_MF)
summary(MFVs.Caste_lm)
plot(imputed_data_MF, imputed_data_Caste)
abline(MFVs.Caste_lm)

CSVs.Caste_lm <- lm(imputed_data_Caste ~ ((imputed_data_CS)^2))
summary(CSVs.Caste_lm)
plot(imputed_data_Caste, imputed_data_CS)
abline(CSVs.Caste_lm)

p_BM$anc_recon


p_BM$anc_recon[1:983,] # Data with imputed species means
p_BM$anc_var[1:20,] # Variances for each estimate
p_BM$anc_recon[1:20,] - sqrt(p_BM$anc_var[1:20,])*1.96 # Lower 95% CI
p_BM$anc_recon[1:20,] + sqrt(p_BM$anc_var[1:20,])*1.96 # Upper 95% CI

plot.phylo(reorder(tree,"postorder"))
nodelabels()

p_BM$anc_recon[21:39,] # Reconstructed ancestral states for each trait
p_BM$anc_var[21:39,] # Variances for each estimate
p_BM$anc_recon[21:39,] - sqrt(p_BM$anc_var[21:39,])*1.96 # Lower 95% CI
p_BM$anc_recon[21:39,] + sqrt(p_BM$anc_var[21:39,])*1.96 # Upper 95% CI




library(Rphylopars)
library(phytools) # for simulating pure-birth phylogenies
set.seed(21) # Set the seed for reproducible results
trait_cov <- matrix(c(4,2,2.2,2,3,1.5,2.2,1.5,2.5),nrow = 3,ncol = 3)
trait_cov # Phylogenetic trait covariance
tree <- pbtree(n = 20)
sim_data <- simtraits(v = trait_cov,tree = tree,nmissing = 10)

sim_data$trait_data # Note that 10 observations are missing completely at random

p_BM <- phylopars(trait_data = sim_data$trait_data,tree = sim_data$tree)
p_BM # Estimated trait covariance
trait_cov # Simulated trait covariance

p_BM$anc_recon[1:20,] # Data with imputed species means
p_BM$anc_var[1:20,] # Variances for each estimate
p_BM$anc_recon[1:20,] - sqrt(p_BM$anc_var[1:20,])*1.96 # Lower 95% CI
p_BM$anc_recon[1:20,] + sqrt(p_BM$anc_var[1:20,])*1.96 # Upper 95% CI

plot.phylo(reorder(tree,"postorder"))
nodelabels()

p_BM$anc_recon[21:39,] # Reconstructed ancestral states for each trait
p_BM$anc_var[21:39,] # Variances for each estimate
p_BM$anc_recon[21:39,] - sqrt(p_BM$anc_var[21:39,])*1.96 # Lower 95% CI
p_BM$anc_recon[21:39,] + sqrt(p_BM$anc_var[21:39,])*1.96 # Upper 95% CI
