#Figures for TEEB talk 2021
#Louis Bell-Roberts
#22/06/2021


library(tidyverse)
library(ape)
library(phylolm)
library(phytools)
library(ggplot2)
library(brms)
library(MCMCglmm)



#Read in data file - ensure that caste variable is set a numeric variable
d <- read.csv("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Primary Dataset/Data_caste_mating_colonyS_WPM_QueenN_cleaned.csv", header = T)
#d <- read.csv(file.choose(), header = T)

#d <- read.csv("./Data/Cleaned/Data_caste_mating_colonyS_WPM_QueenN_cleaned.csv", header = T) #This one works on Gijsbert Computer

d$Caste3 <- as.numeric(as.character(d$Caste3))
d$number.queens.MEAN <- as.numeric(as.character(d$number.queens.MEAN))
data <- d

#Tree files - 'Genus_polytomy_tree.tre' is the most up-to-date tree and should be used
#anttree_species <- read.tree(file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Trees/Nelsen_ultrametric_species/ultrametric_Nelsen_sp.tre")
anttree_species <- read.tree(file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Trees/Polytomy_tree/Genus_polytomy_tree.tre")

##CS, MF and PG (binary)
antdata_multiple_regression <- filter(data, type == 'ant', Caste3 >=1, eff.mating.freq.MEAN.harmonic >=1, colony.size >=1, polygyny.clean >= 0, Reference.2 != "")


#Prune tree
pruned.tree_sp_multiple_regression<-drop.tip(anttree_species, setdiff(anttree_species$tip.label, antdata_multiple_regression$animal))
pruned.tree_sp_multiple_regression
plotTree(pruned.tree_sp_multiple_regression,ftype="i",fsize=0.4,lwd=1)
dev.off()

#Filter through my dataframe and select only the rows that match the tips of my tree
antdata_multiple_regression.1<-filter(antdata_multiple_regression, animal %in% pruned.tree_sp_multiple_regression$tip.label)
#View(antdata_CS_MF.1)


#Configure dataframe into correct formatting
antdata_multiple_regression.2 <- cbind(antdata_multiple_regression.1$animal, antdata_multiple_regression.1)
##Remove column name for 'animal'
antdata_multiple_regression.3 <-antdata_multiple_regression.2 %>% remove_rownames() %>% column_to_rownames(var="animal")
###Rename column
antdata_multiple_regression.4 <- antdata_multiple_regression.3 %>% 
  rename(
    animal = `antdata_multiple_regression.1$animal`
  )


#Data formatting for brms
antdata_multiple_regression.4$obs <- 1:nrow(antdata_multiple_regression.4)
antdata_multiple_regression.5 <- dplyr::select(antdata_multiple_regression.4, Caste3, colony.size, eff.mating.freq.MEAN.harmonic, polygyny.clean, animal, obs)
antdata_multiple_regression.5$eff.mating.freq.MEAN.harmonic <- log(antdata_multiple_regression.5$eff.mating.freq.MEAN.harmonic)
antdata_multiple_regression.5$colony.size <- log(antdata_multiple_regression.5$colony.size)
animal <- pruned.tree_sp_multiple_regression
A <- ape::vcv.phylo(animal)


#Negative binomial with log transformation on MF.
model_pois_MF <-brm(
  Caste3 ~ eff.mating.freq.MEAN.harmonic + (1|gr(animal, cov = A)),
  data = antdata_multiple_regression.5, family = negbinomial(), 
  data2 = list(A = A),
  chains = 2, cores = 2, iter = 4000,
  control = list(adapt_delta = 0.95))

print(model_pois_MF, digits = 10)
me <- conditional_effects(model_pois_MF)
plot(me, points = T)[[1]] + ggplot2::theme_classic() +
  ggtitle("Caste Number vs. Mating Frequency") +
  xlab("Log(Mating Frequency)") +
  ylab("Caste Number") +
  theme(plot.title = element_text(hjust = 0.5, size=30),
    axis.title.x = element_text(size=27),
    axis.title.y = element_text(size=27),
    axis.text.y = element_text(size=27),
    axis.text.x = element_text(size=27))


#To calculate pMCMC
Table=as.data.frame(summary(model_pois_MF)$fixed)
samps = as.matrix(as.mcmc(model_pois_MF))
pMCMC=data.frame()
for(i in 1:length(Table[,1])) {
  temp=mean(samps[,i] < 0)
  if(Table[i,1] < 0) {temp=mean(samps[,i] > 0)}
  pMCMC[i,1]=temp
  colnames(pMCMC)="pMCMC"
}
Table=cbind(Table, pMCMC)

#Negative binomial with log transformation on CS
model_pois_CS <-brm(
  Caste3 ~ colony.size + (1|gr(animal, cov = A)),
  data = antdata_multiple_regression.5, family = negbinomial(), 
  data2 = list(A = A),
  chains = 2, cores = 2, iter = 4000,
  control = list(adapt_delta = 0.95))

print(model_pois_CS, digits = 10)
me <- conditional_effects(model_pois_CS)
plot(me, points = T)[[1]] + ggplot2::theme_classic() +
  ggtitle("Caste Number vs. Colony Size") +
  xlab("Log(Colony Size)") +
  ylab("Caste Number") +
  theme(plot.title = element_text(hjust = 0.5, size=30),
        axis.title.x = element_text(size=27),
        axis.title.y = element_text(size=27),
        axis.text.y = element_text(size=27),
        axis.text.x = element_text(size=27))

#To calculate pMCMC
Table=as.data.frame(summary(model_pois_CS)$fixed)
samps = as.matrix(as.mcmc(model_pois_CS))
pMCMC=data.frame()
for(i in 1:length(Table[,1])) {
  temp=mean(samps[,i] < 0)
  if(Table[i,1] < 0) {temp=mean(samps[,i] > 0)}
  pMCMC[i,1]=temp
  colnames(pMCMC)="pMCMC"
}
Table=cbind(Table, pMCMC)

#Negative binomial with log transformation on CS and MF and includes polygyny
model_pois_MF_CS_PG <-brm(
  Caste3 ~ eff.mating.freq.MEAN.harmonic + colony.size + polygyny.clean + (1|gr(animal, cov = A)),
  data = antdata_multiple_regression.5, family = negbinomial(), 
  data2 = list(A = A),
  chains = 2, cores = 2, iter = 4000,
  control = list(adapt_delta = 0.95))

print(model_pois_MF_CS_PG, digits = 10)
me <- conditional_effects(model_pois_MF_CS_PG)
plot(me, points = T)[[1]] + ggplot2::theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))

#Plot poisson regression manually
xvalues <- sort(antdata_multiple_regression.5$eff.mating.freq.MEAN.harmonic)
log.means <- 0.1106692361+0.0428126557*xvalues

## Un-log the values to get to the lambdas:
mean.values <- exp(log.means)
plot(antdata_multiple_regression.5$eff.mating.freq.MEAN.harmonic, antdata_multiple_regression.5$Caste3)
lines(xvalues,mean.values)

#


#Plotting phylogeny with traits at tips

#Using phytools - quite basic and ugly plots
antdata_multiple_regression.6 <- dplyr::select(antdata_multiple_regression.4, Caste3, colony.size, eff.mating.freq.MEAN.harmonic, polygyny.clean)
antdata_multiple_regression.6$colony.size <- log(antdata_multiple_regression.6$colony.size)
antdata_multiple_regression.6$eff.mating.freq.MEAN.harmonic <- log(antdata_multiple_regression.6$eff.mating.freq.MEAN.harmonic)
#Rename columns
antdata_multiple_regression.6 <- antdata_multiple_regression.6 %>% 
  dplyr::rename(
    "Caste number" = Caste3,
    "Colony size" = colony.size,
    "Mating frequency" = eff.mating.freq.MEAN.harmonic,
    "Queen number" = polygyny.clean
  )
dotTree(pruned.tree_sp_multiple_regression,antdata_multiple_regression.6,standardize=TRUE,length=10)
phylo.heatmap(pruned.tree_sp_multiple_regression,antdata_multiple_regression.6,standardize=T,length=10, split=c(0.7,0.3),fsize=c(0.00001,0.8,0.8))

phylo.heatmap(anoletree,anole.resids,
              split=c(0.7,0.3),fsize=c(0.4,0.8,0.8),
              standardize=TRUE)












