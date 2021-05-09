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
anttree_species <- read.tree(file = "./Data/Trees//Genus_polytomy_tree.tre") #For Gijsbert computer. I think you are loading it from somewhere else, but I presume it's the same file?
#anttree_species <- read.tree(file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Trees/Polytomy_tree/Polytomy_tree.tre")
#anttree_species <- read.tree(file.choose())

#########
#PGLS analysis - Caste ~ Mating frequency

#Caste vs. MF filtering - selecting only the ant species from the database along with species that have data on caste number and mating frequency
antdata_MF <- filter(data, type == 'ant', Caste3 >=1, eff.mating.freq.MEAN.harmonic >=1)
#GW: I presume you restrict to ants because you only have a tree for them? Otherwise, particularly when doing ASR, having an outgroup is actually very good, so could be worthwhile to include a few 'outgroup' species from wasp/bees?

#Prune tree to match the species in my database - 130 species in pruned tree
pruned.tree_sp<-drop.tip(anttree_species, setdiff(anttree_species$tip.label, antdata_MF$animal))
pruned.tree_sp
plotTree(pruned.tree_sp,ftype="i",fsize=0.4,lwd=1)

#Prune the database to select only the rows that match the tips of my tree
antdata_MF.1<-filter(antdata_MF, animal %in% pruned.tree_sp$tip.label)
#View(antdata_MF.1)
nrow(antdata_MF.1) #Sanity check, looks like the numbers match :-)

antdata_MF.1$animal[80]
antdata_MF.1$animal[80]<-"Cardiocondyla_argyrotricha" #Better to get rid of " within names. 
antdata_MF.1$animal[80]

pruned.tree_sp$tip.label[102]
pruned.tree_sp$tip.label[102]<-"Cardiocondyla_argyrotricha"
pruned.tree_sp$tip.label[102]
#Otherwise problems with mcmcglmm later down the line. 

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

#Gijsbert's code
#What's the range of the two variables? 
summary(antdata_MF.4$Caste3)
summary(antdata_MF.4$eff.mating.freq.MEAN.harmonic)
#Replotting the regression lines expanding the y-axis a bit more to better see the variation. 
plot(antdata_MF.4$eff.mating.freq.MEAN.harmonic, antdata_MF.4$Caste3, ylim=c(0.5,5))
abline(CasteVs.MF_Phy_lm)
abline(CasteVs.MF_Phy_glm)
abline(CasteVs.MF_lm)
##First reflection: the relationship doesn't seem to be very strong right? 
#There may be a weak effect where species with a higher mat freq are more likely to have more castes, but it would only be weak I guess. 
#After all, for most caste numbers most mating frequency value are pretty close to 1. For 3 castes they are a bit higher, but few species. 

#Let's first look at it another way. 
boxplot(eff.mating.freq.MEAN.harmonic~Caste3,data=antdata_MF.4)
table(antdata_MF.4$Caste3)
#This shows that mean/median is pretty much the same for 1 vs 2 castes, but it's higher for 3 and 4. 
#But, there only very few species, particulalry for 4 (only 1) - so won't weigh superheavy in regresssion. 

###Let's look back at the regression 
#Let's first have a look at the phylogenetic logistic regression (nr 2)
summary(CasteVs.MF_Phy_glm)
#This regression line drops out of the range of the data, which can't really be right. 

#Why not? Weird? Let's do a normal poisson regression to compare
CasteVs.MF_glm_poisson <- glm(Caste3 ~ eff.mating.freq.MEAN.harmonic, family="poisson",data = antdata_MF.4)
summary(CasteVs.MF_glm_poisson)
plot(antdata_MF.4$eff.mating.freq.MEAN.harmonic, antdata_MF.4$Caste3, ylim=c(0.5,5))
abline(CasteVs.MF_glm_poisson$coefficients)
#Still not great, but slightly better at least. 

#Only other thing I can think of is try MCMC and fit a poission model.
#Code from MPCM book: http://www.mpcm-evolution.com/OPM/Chapter11_OPM/4_nongaussian.html
library(MCMCglmm)
inv.pruned.tree_sp<-inverseA(pruned.tree_sp,nodes="TIPS",scale=TRUE)
prior<-list(G=list(G1=list(V=1,nu=0.02)),R=list(V=1,nu=0.02))
CasteVs.MF_glm_MCMCglmm<-MCMCglmm(Caste3~eff.mating.freq.MEAN.harmonic,random=~animal,
                     family="poisson",ginverse=list(animal=inv.pruned.tree_sp$Ainv),
                     prior=prior,data=antdata_MF.4,nitt=250000,burnin=10000,thin=500)
heidel.diag(CasteVs.MF_glm_MCMCglmm$Sol) #GW: Have to look closer at how these work. 
plot(CasteVs.MF_glm_MCMCglmm)
summary(CasteVs.MF_glm_MCMCglmm)
#Ok, this is qualitatively pretty similar to the previous non-phylogenetic regrssion, i think
CasteVs.MF_glm_poisson$coefficients

#Other option, could the problem be the class of the response variable? 
class(antdata_MF.4$Caste3)
antdata_MF.4$Caste3<-as.integer(antdata_MF.4$Caste3)
class(antdata_MF.4$Caste3)
#Nope, no effect. As it shouldn't have I guess. 

###Ok, so this result seem pretty robust now, but we don't really understand it yet. 
#Have to read more about poisson regressions, which I admit I know very little about to understand thisbetter. 

#Other option, with phyloglm at least, is to turn it into a logistic regression. 
#Ia logistic regression model should always take a binary response variable (yes/no, 1/0, success/fail), of course.
#Let's see what happens if we 'binarise' the response variable, will the model behave? 
table(antdata_MF.4$Caste3)
antdata_MF.4$Caste3_binarised<-
  ifelse(antdata_MF.4$Caste3>1,1,0)
table(antdata_MF.4$Caste3_binarised) #So, 1 castee is now represented by '0', more than one by '1'
CasteVs.MF_Phy_glm_binarised <- 
  phyloglm(Caste3_binarised~eff.mating.freq.MEAN.harmonic,data = antdata_MF.4, 
           phy = pruned.tree_sp, method = "logistic_MPLE",btol = 100)
summary(CasteVs.MF_Phy_glm_binarised)


plot(antdata_MF.4$eff.mating.freq.MEAN.harmonic,antdata_MF.4$Caste3_binarised)
newdat <- data.frame(eff.mating.freq.MEAN.harmonic=seq(min(antdata_MF.4$eff.mating.freq.MEAN.harmonic), 
                            max(antdata_MF.4$eff.mating.freq.MEAN.harmonic),len=100))
newdat$Caste3_binarised<-
  1/
  (1+exp(-(CasteVs.MF_Phy_glm_binarised$coefficients[[2]]*newdat$eff.mating.freq.MEAN.harmonic-CasteVs.MF_Phy_glm_binarised$coefficients[[1]])))
lines(Caste3_binarised ~ eff.mating.freq.MEAN.harmonic, newdat, col="green4", lwd=2)
#This makes some sense I guess, some effect, but quite weak of mating freq on case number, when treated as a binary variable (more than 1 yes/no)

####Ok let's now look back at the linear regression
plot(antdata_MF.4$eff.mating.freq.MEAN.harmonic, antdata_MF.4$Caste3, ylim=c(0.5,5))
abline(CasteVs.MF_Phy_lm,lty="dashed")
abline(CasteVs.MF_lm)
summary(CasteVs.MF_Phy_lm)
summary(CasteVs.MF_lm)
#Ok, so here we have a result which is very significant without phylogeny, but not sig with phylogeny. 
#Is that a bug, or shouldn't we surprised. Let's plot it on the phylogeny

head(antdata_MF.4)
antdata_MF.4_barplot<-antdata_MF.4 %>% select(eff.mating.freq.MEAN.harmonic,Caste3)
head(antdata_MF.4_barplot)
#Let's first have a look at the distribution of the caste numbers. 
plotTree.barplot(pruned.tree_sp,
                 antdata_MF.4_barplot %>% select(Caste3),
                 args.plotTree=list(fsize=0.25))
#Ok, so very strong clustering of larger number of castes than 1. 
#For instance, look at the camponotus polytomy in the centre And the big cluster at the bottom. 
#I guess this really reduces the effective sample size, because there is really only a few transitions from 1 to more castes. 
#Probably about 5-6? 

plotTree.barplot(pruned.tree_sp,
                 antdata_MF.4_barplot %>% select(eff.mating.freq.MEAN.harmonic),
                 args.plotTree=list(fsize=0.25))
#Bit similar story for mating frequency. 

#Plotting them both at the same time. 
plotTree.barplot(pruned.tree_sp,
                 antdata_MF.4_barplot,
                 args.plotTree=list(fsize=0.25),
                 args.barplot=list(beside=TRUE,xlim=c(0,26)))
#Difficult to see, but if you blow it up and look closely you also see quite a lot correlation.

#Overall thoughts by Gijsbert
#So there is a strong effect without phylogeny.
#But, with phylogeny it seems pretty weak
#Unfortunately, I would really be surprised if that's not really an artefact, but just a result from the clustering of the higher caste clades. 
#In other words, the effects seems really weak.
#I suppose this is disappointing, although it's what we do phylogenetic controls for of course. 
#Some things you could do/try now.

#1. Dig deeper into the poisson regression, why do we get these weird parameter values? 
#My hypothesis: abline is not suitable for plotting these coefficients, because it's essentially assuming a linear model (y = a + b*x), which the poission regression isn't
#So if that's true, in fact there's nothing wrong with the regression, but with our attempt to plot it. 
#Btw: Same story with logistic regression, you can't plot it with abline either. 
#But, of course even if the mcmcm parameter values turn out sensible, once you do this (which I would suspect), they are still not signficant.
#Can we get more data? Particulalry, more transitions to separate high caste number clades, to remedy this? 

#2. Do ASRs to get a better feel of where/when high case numbers and/or mating frequency evolve. 
#It's a different question, I guess, but may still be interesint? 


#GW: Btw, I did not look beyond here. Can do, but wanted to send this to you already :-)

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











