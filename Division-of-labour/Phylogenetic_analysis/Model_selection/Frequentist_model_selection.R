#Model selection - PGLS glm's
#Louis Bell-Roberts
#07/06/2021

library(tidyverse)
library(ape)
library(phylolm)
library(phytools)
library(ggplot2)
library(geiger)



#Read in data file - ensure that caste variable is set a numeric variable
d <- read.csv("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Primary Dataset/Data_caste_mating_colonyS_WPM_QueenN_cleaned.csv", header = T)
#d <- read.csv(file.choose(), header = T)

#d <- read.csv("./Data/Cleaned/Data_caste_mating_colonyS_WPM_QueenN_cleaned.csv", header = T) #This one works on Gijsbert Computer

d$Caste3 <- as.numeric(as.character(d$Caste3))
d$number.queens.MEAN <- as.numeric(as.character(d$number.queens.MEAN))
data <- d

#Tree files - 'Genus_polytomy_tree.tre' is the most up-to-date tree
#anttree_species <- read.tree(file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Trees/Nelsen_ultrametric_species/ultrametric_Nelsen_sp.tre")
anttree_species <- read.tree(file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Trees/Polytomy_tree/Genus_polytomy_tree.tre")
#anttree_species <- read.tree(file = "./Data/Trees//Genus_polytomy_tree.tre") #For Gijsbert computer. I think you are loading it from somewhere else, but I presume it's the same file? - yes
#anttree_species <- read.tree(file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Trees/Polytomy_tree/Polytomy_tree.tre")
#anttree_species <- read.tree(file.choose())


#PGLS analysis - Caste ~ Mating frequency

#Caste vs. MF filtering - selecting only the ant species from the database along with species that have data on caste number and mating frequency
model_selection_data <- dplyr::filter(data, type == 'ant', Caste3 >=1, eff.mating.freq.MEAN.harmonic >=1, colony.size >=1, polygyny.clean >= 0, Reference.2 != "")


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

#Graphically check for collinearity
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


#Models: 1st phylogenetic general linear model, 2nd phylogenetic poisson regression, 3rd non-phylogenetic general linear model
linear_model <- lm(Caste3 ~ log(colony.size) + log(eff.mating.freq.MEAN.harmonic) + polygyny.clean, data = model_selection_data.4)
glm_poisson <- glm(Caste3 ~ log(colony.size) + log(eff.mating.freq.MEAN.harmonic) + polygyny.clean, data = model_selection_data.4,
                  family = "poisson")
Phy_lm <- phylolm(Caste3 ~ log(colony.size) + log(eff.mating.freq.MEAN.harmonic) + polygyny.clean, data = model_selection_data.4, phy = pruned.tree)
Phy_glm <- phyloglm(Caste3 ~ log(colony.size) + log(eff.mating.freq.MEAN.harmonic) + polygyny.clean, data = model_selection_data.4, phy = pruned.tree, method = "poisson_GEE")


summary(linear_model)
summary(glm_poisson)
summary(Phy_lm)
summary(Phy_glm)

#Correct plotting technique
coeff <- coef(Phy_glm)
xvalues <- sort(log(model_selection_data.4$eff.mating.freq.MEAN.harmonic))
log.means <- coeff[1]+coeff[3]*xvalues

## Un-log the values to get to the lambdas:
mean.values <- exp(log.means)

plot(log(model_selection_data.4$eff.mating.freq.MEAN.harmonic), model_selection_data.4$Caste3)
lines(xvalues,mean.values)

#Diagnostics
hist(log(model_selection_data.4$colony.size))
hist(log(model_selection_data.4$eff.mating.freq.MEAN.harmonic))
hist(sqrt(model_selection_data.4$eff.mating.freq.MEAN.harmonic))
hist(sqrt(model_selection_data.4$eff.mating.freq.MEAN.harmonic-min(model_selection_data.4$eff.mating.freq.MEAN.harmonic)))
hist(model_selection_data.4$Caste3)
table(model_selection_data.4$polygyny.clean)

##Check distribution of model residuals - result: not normally distributed. Poisson model is more appropriate
hist(Phy_lm$residuals)
qqnorm(Phy_lm$residuals)
qqline(Phy_lm$residuals)
plot(x=fitted(Phy_lm), y=Phy_lm$residuals, pch=19)#looks pretty good
#

#Models of interest - could only do phylolm rather than phyloglm becaues phyloglm objects do not contain AIC scores
mf <- phylolm(Caste3 ~ log(eff.mating.freq.MEAN.harmonic), data = model_selection_data.4, 
                    phy = pruned.tree)
cs <- phylolm(Caste3 ~ log(colony.size), data = model_selection_data.4, 
               phy = pruned.tree)
pg <- phylolm(Caste3 ~ polygyny.clean, data = model_selection_data.4, 
               phy = pruned.tree)
mf_cs <- phylolm(Caste3 ~ log(colony.size) + log(eff.mating.freq.MEAN.harmonic), data = model_selection_data.4, 
                  phy = pruned.tree)
mf_cs_int <- phylolm(Caste3 ~ log(colony.size) * log(eff.mating.freq.MEAN.harmonic), data = model_selection_data.4, 
                  phy = pruned.tree)
mf_pg <- phylolm(Caste3 ~ log(eff.mating.freq.MEAN.harmonic) + polygyny.clean, data = model_selection_data.4, 
                  phy = pruned.tree)
mf_pg_int <- phylolm(Caste3 ~ log(eff.mating.freq.MEAN.harmonic) * polygyny.clean, data = model_selection_data.4, 
                  phy = pruned.tree)
cs_pg <- phylolm(Caste3 ~ log(colony.size) + polygyny.clean, data = model_selection_data.4, 
                  phy = pruned.tree)
cs_pg_int <- phylolm(Caste3 ~ log(colony.size) * polygyny.clean, data = model_selection_data.4, 
                  phy = pruned.tree)
mf_cs_pg <- phylolm(Caste3 ~ log(colony.size) + log(eff.mating.freq.MEAN.harmonic) + polygyny.clean, data = model_selection_data.4, 
                     phy = pruned.tree)
mf_cs_int_pg <- phylolm(Caste3 ~ log(colony.size) * log(eff.mating.freq.MEAN.harmonic) + polygyny.clean, data = model_selection_data.4, 
                     phy = pruned.tree)
mf_cs_pg_int <- phylolm(Caste3 ~ log(colony.size) * polygyny.clean + log(eff.mating.freq.MEAN.harmonic), data = model_selection_data.4, 
                        phy = pruned.tree)

summary(mf)
summary(cs) #2nd best
summary(pg)
summary(mf_cs) #3rd best
summary(mf_cs_int)
summary(mf_pg)
summary(mf_pg_int)
summary(cs_pg)
summary(cs_pg_int) #best
summary(mf_cs_pg)
summary(mf_cs_int_pg)
summary(mf_cs_pg_int)

aicw(c(mf$aic, cs$aic, pg$aic, mf_cs$aic, mf_cs_int$aic, mf_pg$aic, mf_pg_int$aic, cs_pg$aic, cs_pg_int$aic, mf_cs_pg$aic, mf_cs_int_pg$aic, mf_cs_pg_int$aic))
#

####GW Notes
#So looking at the above I guess these are the reasonable (>10%) candidate models from best to worst. With the top two very close to each other (~22%, and number 3 and 4 too ~13%)
summary(cs_pg_int)
summary(cs)
summary(mf_cs)
summary(mf_cs_pg_int)
#So I guess this tells us: CS is definitely important (no big surprise), because always included. Polygyny probably, because in best model, and fourth, and potentiallly MF. 
#If Mf inluced however, it's not really significant. So really only CS and potentially polygyny? 
#I guess this is not entirely in line with expectations, right? 
#So do we trust this result and why not/yes. 
#Have you done any model anaylysis on these / the best model? Larger colony size = more castest make sense, more polygyny more castes less so? 

#Other approach can we calculate aic for the poissoin regressions? 
#Look at model from above? 
summary(Phy_glm)
#In principle should be easy to calculate aic with basic equation, if we have loglikelihood? 
Phy_glm$logLik
#It doesn't give it, why? 
#Looking at original paper for the poisson_GEE method (https://pubmed.ncbi.nlm.nih.gov/12381290/), I think that's because it's actually using a quasi-likelihood estimator, rather than a likelihood estimator.
#Can you calculate 'quasi-aic' for these? Perhaps: useful Ben Bolker Vignette (do you know him?): https://cran.r-project.org/web/packages/bbmle/vignettes/quasi.pdf


###########Old stuff after here and some code for using 'loo'############






##CS, MF and PG (binary)
model_selection_data <- dplyr::filter(data, type == 'ant', Caste3 >=1, eff.mating.freq.MEAN.harmonic >=1, colony.size >=1, polygyny.clean >= 0, Reference.2 != "")
##CS, MF and PG (continuous)
#model_selection_data <- filter(data, type == 'ant', Caste3 >=1, eff.mating.freq.MEAN.harmonic >=1, colony.size >=1, number.queens.MEAN >= 0, Reference.2 != "")



#Prune tree - 521 species
pruned.tree<-drop.tip(anttree_species, setdiff(anttree_species$tip.label, model_selection_data$animal))
pruned.tree
plotTree(pruned.tree,ftype="i",fsize=0.4,lwd=1)

#Filter through my dataframe and select only the rows that match the tips of my tree
model_selection_data.1<-filter(model_selection_data, animal %in% pruned.tree$tip.label)
#View(model_selection_data.1)

#Configure dataframe into correct formatting
model_selection_data.2 <- cbind(model_selection_data.1$animal, model_selection_data.1)
##Remove column name for 'animal'
model_selection_data.3 <-model_selection_data.2 %>% remove_rownames() %>% column_to_rownames(var="animal")
###Rename column
model_selection_data.4 <- model_selection_data.3 %>% 
  rename(
    animal = `model_selection_data.1$animal`
  )

#Models
linear_model <- lm(Caste3 ~ log(colony.size), data = model_selection_data.4)
glm_poisson <- glm(Caste3 ~ log(colony.size) + log(eff.mating.freq.MEAN.harmonic) + polygyny.clean, family="poisson",data = model_selection_data.4)
Phy_lm <- phylolm(Caste3~log(colony.size) + log(eff.mating.freq.MEAN.harmonic) + polygyny.clean, data = model_selection_data.4, phy = pruned.tree)
Phy_glm <- phyloglm(Caste3~log(colony.size) + log(eff.mating.freq.MEAN.harmonic) + polygyny.clean, data = model_selection_data.4, 
                    phy = pruned.tree, method = "poisson_GEE")

CasteVs.MF_lm <- lm(Caste3 ~ eff.mating.freq.MEAN.harmonic, data = model_selection_data.4)


Summary(linear_model)
Summary(glm_poisson)
Summary(Phy_lm)
Summary(Phy_glm)
Summary(CasteVs.MF_lm)
















antdata_multiple_regression <- dplyr::filter(data, type == 'ant', Caste3 >=1, eff.mating.freq.MEAN.harmonic >=1, colony.size >=1, polygyny.clean >= 0, Reference.2 != "")

pred.vars = c("log.eff.mating.freq.MEAN.harmonic", "log.colony.size", 
              "polygyny.clean")
m.mat = permutations(n = 2, r = 3, v = c(F, T), repeats.allowed = T)
models = apply(cbind(T, m.mat), 1, function(xrow) {
  paste(c("1", pred.vars)[xrow], collapse = "+")
})
models = paste("(1|gr(animal, cov = C)", models, sep = "+")

models = paste("Caste3", models, sep = "~")


# AIC of models
all.waic = rep(NA, length(models))
# Estiamted lambdas
all.lambda = rep(NA, length(models))
# Which predictors are estimated in the models beside the intercept
m.mat = cbind(1, m.mat)
colnames(m.mat) = c("(Intercept)", pred.vars)
# number of parameters estimated in the models
n.par = 2 + apply(m.mat, 1, sum)

ant.mods <- comparative.data(phy = pruned.tree_sp_multiple_regression, data = antdata_multiple_regression.4, names.col = animal, 
                            vcv = TRUE, na.omit = FALSE, warn.dropped = TRUE)
coefs=m.mat# define an object to store the coefficients   
for (k in 1:length(models)) {
  res = try(brm(as.formula(models[k]), data = antdata_multiple_regression.4, family = negbinomial(),
                data2 = list(C = C),
                chains = 2, cores = 2, iter = 4000,
                control = list(adapt_delta = 0.95)
                ))
  }
}


brms_negbinom_MF_CS_PG_bin <-brm(
  Caste3 ~ log(eff.mating.freq.MEAN.harmonic) + log(colony.size) + polygyny.clean + (1|gr(animal, cov = C)),
  data = antdata_multiple_regression.4, family = negbinomial(), 
  data2 = list(C = C),
  chains = 2, cores = 2, iter = 4000,
  control = list(adapt_delta = 0.95))















