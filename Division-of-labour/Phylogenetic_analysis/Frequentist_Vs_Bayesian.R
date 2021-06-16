#Frequentist vs Bayesian methods
#Louis Bell-Roberts
#16/06/2021


library(tidyverse)
library(ape)
library(phylolm)
library(phytools)
library(ggplot2)
library(MCMCglmm)

#Read in data file - ensure that caste variable is set a numeric variable
d <- read.csv("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Primary Dataset/Data_caste_mating_colonyS_WPM_QueenN_cleaned.csv", header = T)

#d <- read.csv("./Data/Cleaned/Data_caste_mating_colonyS_WPM_QueenN_cleaned.csv", header = T) #This one works on Gijsbert Computer

d$Caste3 <- as.numeric(as.character(d$Caste3))
d$number.queens.MEAN <- as.numeric(as.character(d$number.queens.MEAN))
data <- d

#Tree files - 'Genus_polytomy_tree.tre' is the most up-to-date tree and should be used
#anttree_species <- read.tree(file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Trees/Nelsen_ultrametric_species/ultrametric_Nelsen_sp.tre")
anttree_species <- read.tree(file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Trees/Polytomy_tree/Genus_polytomy_tree.tre")

#For Gijsbert Computer: anttree_species <- read.tree(file = "./Data/Trees//Genus_polytomy_tree.tre") #Am I using the right one here? - yes

#Caste vs. MF filtering - selecting only the ant species from the database along with species that have data on caste number and mating frequency
model_selection_data <- dplyr::filter(data, type == 'ant', Caste3 >=1, eff.mating.freq.MEAN.harmonic >=1, colony.size >=1, polygyny.clean >= 0, Reference.2 != "")

#Prune tree
pruned.tree<-drop.tip(anttree_species, setdiff(anttree_species$tip.label, model_selection_data$animal))
pruned.tree
plotTree(pruned.tree,ftype="i",fsize=0.4,lwd=1)
dev.off()

#Prune the database to select only the rows that match the tips of my tree
model_selection_data.1<-dplyr::filter(model_selection_data, animal %in% pruned.tree$tip.label)


#Configure dataframe into correct formatting for phylolm - row names must be equal to the species names 
##Duplicate the 'animal' column which contains the combined genus and species names
model_selection_data.2 <- cbind(model_selection_data.1$animal, model_selection_data.1)
##Remove column name for 'animal' - this labels each row as the species name
model_selection_data.3 <-model_selection_data.2 %>% remove_rownames() %>% column_to_rownames(var="animal")
###Rename one of the duplicate columns to 'animal'
model_selection_data.4 <- model_selection_data.3 %>% 
  rename(
    animal = `model_selection_data.1$animal`
  )

#Graphically check for collinearity and perform vif test
plot(log(model_selection_data.4$eff.mating.freq.MEAN.harmonic), log(model_selection_data.4$colony.size)) #Positive correlation is present
plot(model_selection_data.4$polygyny.clean, log(model_selection_data.4$eff.mating.freq.MEAN.harmonic)) #Negative correlation is present
abline(MF_PGlm)
plot(model_selection_data.4$polygyny.clean, log(model_selection_data.4$colony.size)) #No correlation
abline(CS_PGlm)
MF_PGlm <- lm(log(model_selection_data.4$eff.mating.freq.MEAN.harmonic) ~ model_selection_data.4$polygyny.clean)
CS_PGlm <- lm(log(model_selection_data.4$colony.size) ~ model_selection_data.4$polygyny.clean)
All_lm <- lm(model_selection_data.4$Caste3 ~ log(model_selection_data.4$eff.mating.freq.MEAN.harmonic) + log(model_selection_data.4$colony.size) + model_selection_data.4$polygyny.clean)
summary(MF_PGlm)
summary(CS_PGlm)
car::vif(All_lm)


#Frequentist models: 1st: normal linear model, 2nd: glm poisson, 3rd: phylogenetic linear model, 4th: phylogenetic poisson regression
linear_model <- lm(Caste3 ~ log(colony.size) + log(eff.mating.freq.MEAN.harmonic) + polygyny.clean, data = model_selection_data.4)
glm_poisson <- glm(Caste3 ~ log(colony.size) + log(eff.mating.freq.MEAN.harmonic) + polygyny.clean, data = model_selection_data.4,
                   family = "poisson")
Phy_lm <- phylolm(Caste3 ~ log(colony.size) + log(eff.mating.freq.MEAN.harmonic) + polygyny.clean, data = model_selection_data.4, phy = pruned.tree)
Phy_glm <- phyloglm(Caste3 ~ log(colony.size) + log(eff.mating.freq.MEAN.harmonic) + polygyny.clean, data = model_selection_data.4, phy = pruned.tree, method = "poisson_GEE")


summary(linear_model)
summary(glm_poisson)
summary(Phy_lm)
summary(Phy_glm)

#Plotting technique for a poisson glm
coeff <- coef(Phy_glm)
xvalues <- sort(log(model_selection_data.4$eff.mating.freq.MEAN.harmonic))
log.means <- coeff[1]+coeff[3]*xvalues #coeff[3] plots MF coefficient

## Un-log the values to get to the lambdas:
mean.values <- exp(log.means)

plot(log(model_selection_data.4$eff.mating.freq.MEAN.harmonic), model_selection_data.4$Caste3)
lines(xvalues,mean.values)

#Diagnostics
##Distribution of predictors
hist(log(model_selection_data.4$colony.size))
hist(log(model_selection_data.4$eff.mating.freq.MEAN.harmonic))
hist(sqrt(model_selection_data.4$eff.mating.freq.MEAN.harmonic))
hist(sqrt(model_selection_data.4$eff.mating.freq.MEAN.harmonic-min(model_selection_data.4$eff.mating.freq.MEAN.harmonic)))
hist(model_selection_data.4$Caste3)
table(model_selection_data.4$polygyny.clean)

##Check distribution of model residuals - result: not normally distributed. Poisson model is more appropriate. Log transforming Caste3 doesn't really help
hist(Phy_lm$residuals)
qqnorm(Phy_lm$residuals)
qqline(Phy_lm$residuals)
plot(x=fitted(Phy_lm), y=Phy_lm$residuals, pch=19)


######
#Bayesian models

#Prepare tree for the model
inv.pruned.tree_sp<-inverseA(pruned.tree,nodes="TIPS",scale=TRUE)

#Have created a number of different priors - I have read that at large sample sizes the priors specified don't really matter as they will not qualitatively change the results
##Using a couple of different priors here suggests that we probably have a good enough sample size as results don't change very much based on the priors used
prior<-list(G=list(G1=list(V=1,nu=0.002)),R=list(V=1,nu=0.002))
prior2<-list(G=list(G1=list(V=1,nu=0.02)),R=list(V=1,nu=0.02))
prior_ext<-list(R=list(V=1,nu=1), G=list(G1=list(V=1,nu=1,alpha.mu=0,alpha.V=1000)))
prior_exp<-list(R=list(V=1, fix=1), G=list(G1=list(V=1, nu=1, alpha.mu=0, alpha.V=1000))) #quoted from the course notes
prior_exp.1<-list(R=list(V=10, fix=1), G=list(G1=list(V=1, nu=1, alpha.mu=0, alpha.V=1000))) #quoted from the course notes but with higher V
prior.iw<-list(R=list(V=1, nu=1), G=list(G1=list(V=1, nu=1)))
prior.chi=list(R=list(V=1, fix=1), G=list(G1=list(V=1, nu=1000, alpha.mu=0,alpha.V=1)))
p.var<-var(model_selection_data.4$Caste3,na.rm=TRUE) 
prior1.1<-list(G=list( G1=list(V=matrix(p.var/2),n=1)), R=list(V=matrix(p.var/2),n=1))
p.var <- var(model_selection_data.4$Caste3, na.rm = TRUE)
prior.proper <- list(G = list(G1 = list(V = matrix(p.var * 0.05), nu = 1)), R = list(V = matrix(p.var * 0.95), nu = 1))

#Model
CasteVs.MF_CS_PG_glm_MCMCglmm<-MCMCglmm(Caste3~log(colony.size)+log(eff.mating.freq.MEAN.harmonic)+polygyny.clean,random=~animal,
                                     family="poisson",ginverse=list(animal=inv.pruned.tree_sp$Ainv),
                                     prior=prior,data=model_selection_data.4,nitt=250000,burnin=10000,thin=500)

CasteVs.MF_CS_PG_glm_MCMCglmm.1<-MCMCglmm(Caste3~log(colony.size)+log(eff.mating.freq.MEAN.harmonic)+polygyny.clean,random=~animal,
                                          family="poisson",ginverse=list(animal=inv.pruned.tree_sp$Ainv),
                                          prior=prior_exp,data=model_selection_data.4,nitt=250000,burnin=10000,thin=500)

#Comparing summary outputs from Bayesian and frequentist models:
summary(CasteVs.MF_CS_PG_glm_MCMCglmm) #Despite using quite different priors, coefficients remain roughly the same (sample size must be large enough that prior specification doesn't matter too much)
summary(Phy_glm) #phylogenetic poisson regression
summary(glm_poisson) #non-phylogenetic
#Coefficients differ by approximately a factor of nearly 2 for both colony size and MF when comparing between frequentist and bayesian approaches
#Despite coefficients being lower in PGLS model, the coefficients are more significant
#Coefficients from Bayesian model and non-phylogenetic poisson glm are more similar. Levels of significance are also more similar




#Bayesian diagnostics
plot(CasteVs.MF_CS_PG_glm_MCMCglmm.1)
autocorr.diag(CasteVs.MF_CS_PG_glm_MCMCglmm$Sol)
autocorr.diag(CasteVs.MF_CS_PG_glm_MCMCglmm$VCV)
posterior.mode(CasteVs.MF_CS_PG_glm_MCMCglmm$VCV)
HPDinterval(CasteVs.MF_CS_PG_glm_MCMCglmm$VCV)
effectiveSize(CasteVs.MF_CS_PG_glm_MCMCglmm$Sol)
effectiveSize(CasteVs.MF_CS_PG_glm_MCMCglmm$VCV)

heidel.diag(CasteVs.MF_CS_PG_glm_MCMCglmm$VCV)

#Phylogenetic signal
herit <- CasteVs.MF_CS_PG_glm_MCMCglmm$VCV[, "animal"]/(CasteVs.MF_CS_PG_glm_MCMCglmm$VCV[, "animal"] + CasteVs.MF_CS_PG_glm_MCMCglmm$VCV[, "units"])
effectiveSize(herit)
mean(herit)
HPDinterval(herit) #Display 95% credible interval
plot(herit)
















