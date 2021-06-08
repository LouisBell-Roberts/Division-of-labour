#Phylogenetic path analysis
#Louis Bell-Roberts
#25/05/2021


library(tidyverse)
library(ape)
library(phylolm)
library(phytools)
library(ggplot2)
library(phylopath)

#Read in data file - ensure that caste variable is set a numeric variable
d <- read.csv("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Primary Dataset/Data_caste_mating_colonyS_WPM_QueenN_cleaned.csv", header = T)

#d <- read.csv("./Data/Cleaned/Data_caste_mating_colonyS_WPM_QueenN_cleaned.csv", header = T) #This one works on Gijsbert Computer

d$Caste3 <- as.numeric(as.character(d$Caste3))
d$number.queens.MEAN <- as.numeric(as.character(d$number.queens.MEAN))
data <- d

#Tree files - 'Genus_polytomy_tree.tre' is the most up-to-date tree
#anttree_species <- read.tree(file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Trees/Nelsen_ultrametric_species/ultrametric_Nelsen_sp.tre")
anttree_species <- read.tree(file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Trees/Polytomy_tree/Genus_polytomy_tree.tre")

#Filter data
antdata_multiple_regression <- filter(data, type == 'ant', Caste3 >=1, eff.mating.freq.MEAN.harmonic >=1, colony.size >=1, polygyny.clean >= 0, Reference.2 != "")

#Prune tree
pruned.tree_sp_multiple_regression<-drop.tip(anttree_species, setdiff(anttree_species$tip.label, antdata_multiple_regression$animal))
pruned.tree_sp_multiple_regression
plotTree(pruned.tree_sp_multiple_regression,ftype="i",fsize=0.4,lwd=1)
dev.off()

#Filter through my dataframe and select only the rows that match the tips of my tree
antdata_multiple_regression.1<-filter(antdata_multiple_regression, animal %in% pruned.tree_sp_multiple_regression$tip.label)

#Configure dataframe into correct formatting
antdata_multiple_regression.2 <- cbind(antdata_multiple_regression.1$animal, antdata_multiple_regression.1)
##Remove column name for 'animal'
antdata_multiple_regression.3 <-antdata_multiple_regression.2 %>% remove_rownames() %>% column_to_rownames(var="animal")
###Rename column
antdata_multiple_regression.4 <- antdata_multiple_regression.3 %>% 
  rename(
    animal = `antdata_multiple_regression.1$animal`
  )
#log transform all values in the colony.size column - transforming MF has significant effect
antdata_multiple_regression.5 <- antdata_multiple_regression.4
antdata_multiple_regression.5$colony.size <- log(antdata_multiple_regression.5$colony.size)
antdata_multiple_regression.5$eff.mating.freq.MEAN.harmonic <- log(antdata_multiple_regression.5$eff.mating.freq.MEAN.harmonic)



#Define the models for the analysis: response before the ~, effector after the ~
models <- define_model_set(
  one   = c(Caste3 ~ eff.mating.freq.MEAN.harmonic + colony.size),
  two   = c(Caste3 ~ eff.mating.freq.MEAN.harmonic),
  three = c(Caste3 ~ colony.size),
  four  = c(Caste3 ~ eff.mating.freq.MEAN.harmonic, colony.size ~ eff.mating.freq.MEAN.harmonic),
  five  = c(Caste3 ~ eff.mating.freq.MEAN.harmonic, colony.size ~ eff.mating.freq.MEAN.harmonic, Caste3 ~ colony.size),
  six   = c(Caste3 ~ eff.mating.freq.MEAN.harmonic + colony.size + polygyny.clean),
  seven = c(Caste3 ~ eff.mating.freq.MEAN.harmonic + colony.size, polygyny.clean ~ eff.mating.freq.MEAN.harmonic),
  eight = c(Caste3 ~ eff.mating.freq.MEAN.harmonic + colony.size, eff.mating.freq.MEAN.harmonic ~ polygyny.clean),
  nine  = c(Caste3 ~ eff.mating.freq.MEAN.harmonic, colony.size ~ eff.mating.freq.MEAN.harmonic),
  ten  = c(Caste3 ~ eff.mating.freq.MEAN.harmonic, eff.mating.freq.MEAN.harmonic ~ colony.size, Caste3 ~ colony.size),
  eleven = c(eff.mating.freq.MEAN.harmonic ~ Caste3),
  twelve = c(eff.mating.freq.MEAN.harmonic ~ Caste3, Caste3 ~ colony.size),
  thirteen = c(eff.mating.freq.MEAN.harmonic ~ Caste3, eff.mating.freq.MEAN.harmonic ~ colony.size),
  fourteen = c(eff.mating.freq.MEAN.harmonic ~ Caste3, eff.mating.freq.MEAN.harmonic ~ colony.size, Caste3 ~ colony.size),
  fifteen = c(colony.size ~ polygyny.clean, Caste3 ~ colony.size),
  sixteen = c(polygyny.clean ~ colony.size, Caste3 ~ colony.size),
  seventeen = c(polygyny.clean ~ colony.size, colony.size ~ Caste3),
  eighteen = c(eff.mating.freq.MEAN.harmonic ~ colony.size, Caste3 ~ eff.mating.freq.MEAN.harmonic),
  nineteen = c(colony.size ~ Caste3, eff.mating.freq.MEAN.harmonic ~ colony.size),
  twenty = c(colony.size ~ Caste3, Caste3 ~ eff.mating.freq.MEAN.harmonic)
  )  
plot_model_set(models)

result <- phylo_path(models, data = antdata_multiple_regression.5, tree = pruned.tree_sp_multiple_regression, model = 'lambda')
result
plot(summary(result))
plot(models$ten)
plot(models$fourteen)
plot(average(result))


best_model <- best(result)
plot(best_model)
average_model <- average(result)
plot(average_model, algorithm = 'mds', curvature = 0.1)

average_model_full <- average(result, avg_method = "full")
plot(average_model_full, algorithm = 'mds', curvature = 0.1)
coef_plot(average_model)

result$d_sep$ten

########
,
six   = c(Caste3 ~ eff.mating.freq.MEAN.harmonic + colony.size + polygyny.clean)
nine  = c(Caste3 ~ RS, RS ~ LS),
  .common = c(LS ~ BM, NL ~ BM, DD ~ NL)
)












