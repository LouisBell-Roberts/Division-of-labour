# BEE conference: phylogenetically controlled regressions for Caste vs. MF or CS
#Louis Bell-Roberts
#11/03/2021


#*****************************************************************

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
#anttree_species <- read.tree(file = "./Data/Trees//Genus_polytomy_tree.tre") #For Gijsbert computer. I think you are loading it from somewhere else, but I presume it's the same file? - yes
#anttree_species <- read.tree(file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Trees/Polytomy_tree/Polytomy_tree.tre")
#anttree_species <- read.tree(file.choose())

#########
#PGLS analysis - Caste ~ Mating frequency

#Caste vs. MF filtering - selecting only the ant species from the database along with species that have data on caste number and mating frequency
antdata_MF <- filter(data, type == 'ant', Caste3 >=1, eff.mating.freq.MEAN.harmonic >=1)
#GW: I presume you restrict to ants because you only have a tree for them? Otherwise, particularly when doing ASR, having an outgroup is actually very good, so could be worthwhile to include a few 'outgroup' species from wasp/bees?

a <- ggplot(antdata_MF, aes(x = sqrt(eff.mating.freq.MEAN.harmonic-min(eff.mating.freq.MEAN.harmonic))))
a <- ggplot(antdata_MF, aes(x = -((eff.mating.freq.MEAN.harmonic)^-1)))
reciprocal

a + geom_histogram(bins = 30, color = "black", fill = "gray")

dplot(sqrt(antdata_MF$eff.mating.freq.MEAN.harmonic-min(antdata_MF$eff.mating.freq.MEAN.harmonic)), antdata_MF$Caste3)

#Prune tree to match the species in my database - 131 species in pruned tree **Dorylus_nigricans_molestus?**
pruned.tree_sp<-drop.tip(anttree_species, setdiff(anttree_species$tip.label, antdata_MF$animal))
pruned.tree_sp
plotTree(pruned.tree_sp,ftype="i",fsize=0.4,lwd=1)
dev.off()
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
CasteVs.MF_Phy_lm <- phylolm(log(Caste3)~log(eff.mating.freq.MEAN.harmonic),data = antdata_MF.4, phy = pruned.tree_sp)
CasteVs.MF_Phy_glm <- phyloglm(Caste3~log(eff.mating.freq.MEAN.harmonic),data = antdata_MF.4, phy = pruned.tree_sp, method = "poisson_GEE")
CasteVs.MF_lm <- lm(Caste3 ~ eff.mating.freq.MEAN.harmonic, data = antdata_MF.4)

summary(CasteVs.MF_Phy_lm)
summary(CasteVs.MF_Phy_glm)
summary(CasteVs.MF_lm)

## Get and plot the fitted mean structure:
coeff <- coef(CasteVs.MF_Phy_glm)
xvalues <- sort(antdata_MF.4$eff.mating.freq.MEAN.harmonic)
log.means <- coeff[1]+coeff[2]*xvalues

## Un-log the values to get to the lambdas:
mean.values <- exp(log.means)

plot(antdata_MF.4$eff.mating.freq.MEAN.harmonic, antdata_MF.4$Caste3)
lines(xvalues,mean.values)


#Plot the PGLS regression lines onto the data and the non-phylogenetic regression
plot(antdata_MF.4$eff.mating.freq.MEAN.harmonic, antdata_MF.4$Caste3, ylim=c(1,4))
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
#But, there only very few species, particularly for 4 (only 1) - so won't weigh super-heavy in regression. 

###Let's look back at the regression 
#Let's first have a look at the phylogenetic logistic regression (nr 2)
summary(CasteVs.MF_Phy_glm)
#This regression line drops out of the range of the data, which can't really be right. 

########LOOK MORE HERE - this is much more similar to my poisson phylogenetic regressions with MCMC and brms - slight effect of phylogeny that reduces the effect size of MF but not so dramatic anymore!! Plus, in G's plot, he plots it incorrectly which also gives the wrong impression. Make the slope coefficient of 0.03 seem small - but really this is on a log scale! actually not such a small effect size after all! try plotting non-phylogenetic poisson regression with my new system ##########
#Why not? Weird? Let's do a normal poisson regression to compare
CasteVs.MF_glm_poisson <- glm(Caste3 ~ eff.mating.freq.MEAN.harmonic, family="poisson",data = antdata_MF.4)
summary(CasteVs.MF_glm_poisson)
# Incorrect plotting technique: plot(antdata_MF.4$eff.mating.freq.MEAN.harmonic, antdata_MF.4$Caste3, ylim=c(0.5,5))
#abline(CasteVs.MF_glm_poisson$coefficients)
#Still not great, but slightly better at least. 

#Correct plotting technique
coeff <- coef(CasteVs.MF_glm_poisson)
xvalues <- sort(antdata_MF.4$eff.mating.freq.MEAN.harmonic)
log.means <- coeff[1]+coeff[2]*xvalues

## Un-log the values to get to the lambdas:
mean.values <- exp(log.means)

plot(antdata_MF.4$eff.mating.freq.MEAN.harmonic, antdata_MF.4$Caste3)
lines(xvalues,mean.values)

#Only other thing I can think of is try MCMC and fit a poisson model.
#Code from MPCM book: http://www.mpcm-evolution.com/OPM/Chapter11_OPM/4_nongaussian.html
library(MCMCglmm)
inv.pruned.tree_sp<-inverseA(pruned.tree_sp,nodes="TIPS",scale=TRUE)
prior<-list(G=list(G1=list(V=1,nu=0.02)),R=list(V=1,nu=0.02))
CasteVs.MF_glm_MCMCglmm<-MCMCglmm(Caste3~eff.mating.freq.MEAN.harmonic,random=~animal,
                     family="poisson",ginverse=list(animal=inv.pruned.tree_sp$Ainv),
                     prior=prior,data=antdata_MF.4,nitt=250000,burnin=10000,thin=500)

CasteVs.MF_glm_MCMCglmm.1<-MCMCglmm(Caste3~log(eff.mating.freq.MEAN.harmonic),random=~animal,
                                  family="poisson",ginverse=list(animal=inv.pruned.tree_sp$Ainv),
                                  prior=prior,data=antdata_MF.4,nitt=250000,burnin=10000,thin=500)

heidel.diag(CasteVs.MF_glm_MCMCglmm$Sol) #GW: Have to look closer at how these work. 
plot(CasteVs.MF_glm_MCMCglmm)
summary(CasteVs.MF_glm_MCMCglmm)

## Get and plot the fitted mean structure:
xvalues <- sort(antdata_MF.4$eff.mating.freq.MEAN.harmonic)
log.means <- 0.155200+0.030549*xvalues

## Un-log the values to get to the lambdas:
mean.values <- exp(log.means)

plot(antdata_MF.4$eff.mating.freq.MEAN.harmonic, antdata_MF.4$Caste3)
lines(xvalues,mean.values)

#Ok, this is qualitatively pretty similar to the previous non-phylogenetic regression, I think
CasteVs.MF_glm_poisson$coefficients

#Other option, could the problem be the class of the response variable? 
class(antdata_MF.4$Caste3)
antdata_MF.4$Caste3<-as.integer(antdata_MF.4$Caste3)
class(antdata_MF.4$Caste3)
#Nope, no effect. As it shouldn't have I guess. 

###Ok, so this result seem pretty robust now, but we don't really understand it yet. 
#Have to read more about poisson regressions, which I admit I know very little about to understand this better. 

#Other option, with phyloglm at least, is to turn it into a logistic regression. 
#In a logistic regression model should always take a binary response variable (yes/no, 1/0, success/fail), of course.
#Let's see what happens if we 'binarise' the response variable, will the model behave? 
table(antdata_MF.4$Caste3)
antdata_MF.4$Caste3_binarised<-
  ifelse(antdata_MF.4$Caste3>1,1,0)
table(antdata_MF.4$Caste3_binarised) #So, 1 caste is now represented by '0', more than one by '1'
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
#This makes some sense I guess, some effect, but quite weak of mating freq on caste number, when treated as a binary variable (more than 1 yes/no)

####Ok let's now look back at the linear regression
plot(antdata_MF.4$eff.mating.freq.MEAN.harmonic, antdata_MF.4$Caste3, ylim=c(0.5,5))
abline(CasteVs.MF_Phy_lm,lty="dashed")
abline(CasteVs.MF_lm)
summary(CasteVs.MF_Phy_lm)
summary(CasteVs.MF_lm)
#Ok, so here we have a result which is very significant without phylogeny, but not sig with phylogeny. 
#Is that a bug, or shouldn't we be surprised? Let's plot it on the phylogeny

head(antdata_MF.4)
antdata_MF.4_barplot<-antdata_MF.4 %>% dplyr::select(eff.mating.freq.MEAN.harmonic,Caste3)
head(antdata_MF.4_barplot)
#Let's first have a look at the distribution of the caste numbers. 
plotTree.barplot(pruned.tree_sp,
                 antdata_MF.4_barplot %>% dplyr::select(Caste3),
                 args.plotTree=list(fsize=0.25))
#Ok, so very strong clustering of larger number of castes than 1. 
#For instance, look at the camponotus polytomy in the centre And the big cluster at the bottom. 
#I guess this really reduces the effective sample size, because there is really only a few transitions from 1 to more castes. 
#Probably about 5-6? 

plotTree.barplot(pruned.tree_sp,
                 antdata_MF.4_barplot %>% dplyr::select(eff.mating.freq.MEAN.harmonic),
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
#Unfortunately, I would really be surprised if that's not really an artifact, but just a result from the clustering of the higher caste clades. 
#In other words, the effects seems really weak.
#I suppose this is disappointing, although it's what we do phylogenetic controls for of course. 
#Some things you could do/try now.

#1. Dig deeper into the poisson regression, why do we get these weird parameter values? 
#My hypothesis: abline is not suitable for plotting these coefficients, because it's essentially assuming a linear model (y = a + b*x), which the poisson regression isn't
#So if that's true, in fact there's nothing wrong with the regression, but with our attempt to plot it. 
#Btw: Same story with logistic regression, you can't plot it with abline either. 
#But, of course even if the mcmcm parameter values turn out sensible, once you do this (which I would suspect), they are still not significant.
#Can we get more data? Particularly, more transitions to separate high caste number clades, to remedy this? 

#2. Do ASRs to get a better feel of where/when high caste numbers and/or mating frequency evolve. 
#It's a different question, I guess, but may still be interesting? 


#GW: Btw, I did not look beyond here. Can do, but wanted to send this to you already :-)

###Analyses following the comments from Gijsbert: brms
#Data formatting for brms
antdata_MF.4$obs <- 1:nrow(antdata_MF.4)
antdata_MF.5 <- dplyr::select(antdata_MF.4, Caste3, eff.mating.freq.MEAN.harmonic, animal, obs)
animal <- pruned.tree_sp
A <- ape::vcv.phylo(animal)

#The '(1|obs)' bit is required as the Poisson distribution does not have a natural overdispersion parameter. 
##We model the residual variance via the group-level effects of obs
#1st model: poisson
model_pois <-brm(
  Caste3 ~ eff.mating.freq.MEAN.harmonic + (1|gr(animal, cov = A)) + (1|obs),
  data = antdata_MF.5, family = poisson("log"), 
  data2 = list(A = A),
  chains = 2, cores = 2, iter = 4000,
  control = list(adapt_delta = 0.95)
  )

#Not sure yet how to change axis labels and plotting preferences - might use similar technique to ggplot?
plot(conditional_effects(model_pois), points = TRUE)
summary(model_pois)

#Plot manually
## Get and plot the fitted mean structure:
xvalues <- sort(antdata_MF.5$eff.mating.freq.MEAN.harmonic)
log.means <- 0.17+0.03*xvalues

## Un-log the values to get to the lambdas:
mean.values <- exp(log.means)

plot(antdata_MF.5$eff.mating.freq.MEAN.harmonic, antdata_MF.5$Caste3)
lines(xvalues,mean.values)

#2nd model: negative binomial with log transformation on MF.
model_pois.1 <-brm(
  Caste3 ~ log(eff.mating.freq.MEAN.harmonic) + (1|gr(animal, cov = A)),
  data = antdata_MF.5, family = negbinomial(), 
  data2 = list(A = A),
  chains = 2, cores = 2, iter = 4000,
  control = list(adapt_delta = 0.95))

summary(model_pois.1)
#plot(conditional_effects(model_pois.1), points = TRUE)  ---- this doesn't seem to work because there's probably an issue with MF not being on log-scale

#Plot manually
xvalues <- sort(log(antdata_MF.5$eff.mating.freq.MEAN.harmonic))
log.means <- 0.14+0.21*xvalues

## Un-log the values to get to the lambdas:
mean.values <- exp(log.means)
plot(log(antdata_MF.5$eff.mating.freq.MEAN.harmonic), antdata_MF.5$Caste3)
lines(xvalues,mean.values)

#
pp_check(model_pois.1)
loo(model_pois)
########
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
ggplot(antdata_MF.4, aes(x = log(eff.mating.freq.MEAN.harmonic), y = Caste3)) + 
  geom_point(size = 2, alpha = 0.4) +
  geom_smooth(method = "loess") +
  theme_classic() +
  ggtitle("Number of ant castes vs. Female Mating Frequency") +
  xlab("Mating Frequency") +
  ylab("Number of castes") +
  theme(plot.title = element_text(hjust = 0.5, size = 22, face = "bold"), axis.line = element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = "black", size = 1), 
        axis.text.x = element_text(face="bold", color="black", 
                                   size=22),
        axis.text.y = element_text(face="bold", color = "black", 
                                   size=22),
        axis.text=element_text(size=22),
        axis.title=element_text(size=22,face="bold"))
#  geom_jitter(width = 0.05, height = 0.05, size = 2) +







#########
#PGLS analysis - Caste ~ Colony size

#Caste vs. CS filtering - selecting only the ant species from the database along with species that have data on caste number and colony size
antdata_CS <- filter(data, type == 'ant', Caste3 >=1, colony.size >=1)

#Prune tree - 521 species
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
CasteVs.CS_glm_poisson <- glm(Caste3 ~ log(colony.size), family="poisson",data = antdata_CS.4)


summary(CasteVs.CS_Phy_lm)
summary(CasteVs.CS_Phy_glm)
summary(CasteVs.CS_lm)
summary(CasteVs.CS_glm_poisson)

#CS poisson glm plot
#Plot manually
xvalues <- sort(log(antdata_CS.4$colony.size))
log.means <- -0.10706+0.05065*xvalues

## Un-log the values to get to the lambdas:
mean.values <- exp(log.means)
plot(log(antdata_CS.4$colony.size), antdata_CS.4$Caste3)
lines(xvalues,mean.values)

#CS lm plot
plot(log(antdata_CS.4$colony.size), antdata_CS.4$Caste3, ylim = c(1, 4))
abline(CasteVs.CS_lm)

##2nd model: negative binomial with log transformation on MF.
model_pois_CS.1 <-brm(
  Caste3 ~ log(colony.size) + (1|gr(animal, cov = B)),
  data = antdata_CS.5, family = negbinomial(), 
  data2 = list(B = B),
  chains = 2, cores = 2, iter = 4000,
  control = list(adapt_delta = 0.95))
summary(model_pois_CS.1)
#MCMCglmm models - taking a very long time (only at 54000 after 2 hours)
inv.pruned.tree_sp_CS<-inverseA(pruned.tree_sp_CS,nodes="TIPS",scale=TRUE)
prior<-list(G=list(G1=list(V=1,nu=0.02)),R=list(V=1,nu=0.02))
CasteVs.CS_Phy_MCMCglm<-MCMCglmm(Caste3~log(colony.size),random=~animal,
                                    family="poisson",ginverse=list(animal=inv.pruned.tree_sp_CS$Ainv),
                                    prior=prior,data=antdata_CS.4,nitt=250000,burnin=10000,thin=500)
summary(CasteVs.CS_Phy_MCMCglm)

#Plot manually - very similar to brms result but took a lot longer to run
xvalues <- sort(log(antdata_CS.4$colony.size))
log.means <- -0.16364+0.05057*xvalues

## Un-log the values to get to the lambdas:
mean.values <- exp(log.means)
plot(log(antdata_CS.4$colony.size), antdata_CS.4$Caste3)
lines(xvalues,mean.values)

#brms models
antdata_CS.4$obs <- 1:nrow(antdata_CS.4)
antdata_CS.5 <- dplyr::select(antdata_CS.4, Caste3, colony.size, animal, obs)
animal <- pruned.tree_sp_CS
B <- ape::vcv.phylo(animal)

#The '(1|obs)' bit is required as the Poisson distribution does not have a natural overdispersion parameter. 
##We model the residual variance via the group-level effects of obs
##1st model: poisson
model_pois_CS <-brm(
  Caste3 ~ log(colony.size) + (1|gr(animal, cov = B)) + (1|obs),
  data = antdata_CS.5, family = poisson("log"), 
  data2 = list(B = B),
  chains = 2, cores = 2, iter = 4000,
  control = list(adapt_delta = 0.95)
)

#Summary
plot(conditional_effects(model_pois_CS), points = TRUE)
summary(model_pois_CS)

#Plot manually
xvalues <- sort(log(antdata_CS.5$colony.size))
log.means <- -0.15+0.05*xvalues

## Un-log the values to get to the lambdas:
mean.values <- exp(log.means)
plot(log(antdata_CS.5$colony.size), antdata_CS.5$Caste3)
lines(xvalues,mean.values)

##2nd model: negative binomial with log transformation on MF.
model_pois_CS.1 <-brm(
  Caste3 ~ log(colony.size) + (1|gr(animal, cov = B)),
  data = antdata_CS.5, family = negbinomial(), 
  data2 = list(B = B),
  chains = 2, cores = 2, iter = 4000,
  control = list(adapt_delta = 0.95))


#Summary
plot(conditional_effects(model_pois_CS.1), points = TRUE)
summary(model_pois_CS.1)

#Plot manually
xvalues <- sort(log(antdata_CS.5$colony.size))
log.means <- -0.15+0.05*xvalues

## Un-log the values to get to the lambdas:
mean.values <- exp(log.means)
plot(log(antdata_CS.5$colony.size), antdata_CS.5$Caste3)
lines(xvalues,mean.values)



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
  geom_point(size = 3, alpha = 0.1) +
  geom_smooth(method = "loess") +
  theme_classic() +
  ggtitle("Number of ant castes vs. log Colony Size") +
  xlab("log Colony Size") +
  ylab("Number of castes") +
  theme(plot.title = element_text(hjust = 0.5, size = 22, face = "bold"), axis.line = element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = "black", size = 1), 
        axis.text.x = element_text(face="bold", color="black", 
                                   size=22),
        axis.text.y = element_text(face="bold", color = "black", 
                                   size=22),
        axis.text=element_text(size=22),
        axis.title=element_text(size=22,face="bold"))




#########





#Multiple regression - testing various combinations of predictors to identify their effects on each other
#Predictor filtering - variables of interest: caste, MF, CS, PG (binary), PG (categorical)

##CS and MF
antdata_multiple_regression <- filter(data, type == 'ant', Caste3 >=1, eff.mating.freq.MEAN.harmonic >=1, colony.size >=1)
##MF and PG (binary)
antdata_multiple_regression <- filter(data, type == 'ant', Caste3 >=1, eff.mating.freq.MEAN.harmonic >=1, polygyny.clean >= 0, Reference.2 != "")
##MF and PG (continuous)
antdata_multiple_regression <- filter(data, type == 'ant', Caste3 >=1, eff.mating.freq.MEAN.harmonic >=1, number.queens.MEAN >= 0, Reference.2 != "")
##CS and PG (binary)
antdata_multiple_regression <- filter(data, type == 'ant', Caste3 >=1, colony.size >=1, polygyny.clean >= 0, Reference.2 != "")
##CS and PG (continuous)
antdata_multiple_regression <- filter(data, type == 'ant', Caste3 >=1, colony.size >=1, number.queens.MEAN >= 0, Reference.2 != "")
##CS, MF and PG (binary)
antdata_multiple_regression <- filter(data, type == 'ant', Caste3 >=1, eff.mating.freq.MEAN.harmonic >=1, colony.size >=1, polygyny.clean >= 0, Reference.2 != "")
##CS, MF and PG (continuous)
antdata_multiple_regression <- filter(data, type == 'ant', Caste3 >=1, eff.mating.freq.MEAN.harmonic >=1, colony.size >=1, number.queens.MEAN >= 0, Reference.2 != "")


#Remove some outliers (optional)
#antdata_CS_MF <- filter(antdata_CS_MF, animal != "Formica_yessensis", animal != "Lasius_neglectus", animal != "Pseudomyrmex_veneficus")

#Prune tree
pruned.tree_sp_multiple_regression<-drop.tip(anttree_species, setdiff(anttree_species$tip.label, antdata_multiple_regression$animal))
pruned.tree_sp_multiple_regression
plotTree(pruned.tree_sp_multiple_regression,ftype="i",fsize=0.4,lwd=1)
dev.off()

#Filter through my dataframe and select only the rows that match the tips of my tree
antdata_multiple_regression.1<-filter(antdata_multiple_regression, animal %in% pruned.tree_sp_multiple_regression$tip.label)
#View(antdata_CS_MF.1)

#Change these species names
antdata_CS_MF.1$animal[72]
antdata_CS_MF.1$animal[72]<-"Cardiocondyla_argyrotricha" #Better to get rid of " within names. 
antdata_CS_MF.1$animal[72]

pruned.tree_sp_CS_MF$tip.label[94]
pruned.tree_sp_CS_MF$tip.label[94]<-"Cardiocondyla_argyrotricha"
pruned.tree_sp_CS_MF$tip.label[94]
#Otherwise problems with mcmcglmm later down the line. 


#Configure dataframe into correct formatting
antdata_multiple_regression.2 <- cbind(antdata_multiple_regression.1$animal, antdata_multiple_regression.1)
##Remove column name for 'animal'
antdata_multiple_regression.3 <-antdata_multiple_regression.2 %>% remove_rownames() %>% column_to_rownames(var="animal")
###Rename column
antdata_multiple_regression.4 <- antdata_multiple_regression.3 %>% 
  rename(
    animal = `antdata_multiple_regression.1$animal`
  )

###Models
#MCMC test - see if concur with brms. Passes the test
inv.pruned.tree_sp<-inverseA(pruned.tree_sp_multiple_regression,nodes="TIPS",scale=TRUE)
#When fitting priors in THIS format, it is only refering to the random effects. Fitting priors for fixed effects is rare
prior<-list(G=list(G1=list(V=1,nu=0.002)),R=list(V=1,nu=0.002))
prior_ext<-list(R=list(V=1,nu=1), G=list(G1=list(V=1,nu=1,alpha.mu=0,alpha.V=1000)))
prior_exp<-list(R=list(V=1, fix=1), G=list(G1=list(V=1, nu=1, alpha.mu=0, alpha.V=1000))) #quoted from the course notes
prior_exp.1<-list(R=list(V=10, fix=1), G=list(G1=list(V=1, nu=1, alpha.mu=0, alpha.V=1000))) #quoted from the course notes but with higher V
prior.iw<-list(R=list(V=1, nu=1), G=list(G1=list(V=1, nu=1)))
prior.chi=list(R=list(V=1, fix=1), G=list(G1=list(V=1, nu=1000, alpha.mu=0,alpha.V=1)))
p.var<-var(antdata_multiple_regression.4$Caste3,na.rm=TRUE) 
prior1.1<-list(G=list( G1=list(V=matrix(p.var/2),n=1)), R=list(V=matrix(p.var/2),n=1))
#Defined as 'proper priors'
p.var <- var(antdata_multiple_regression.4$Caste3, na.rm = TRUE)
prior.proper <- list(G = list(G1 = list(V = matrix(p.var * 0.05), nu = 1)), R = list(V = matrix(p.var * 0.95), nu = 1))

CasteVs.MF_CS_glm_MCMCglmm<-MCMCglmm(Caste3~log(eff.mating.freq.MEAN.harmonic)+log(colony.size),random=~animal,
                                    family="gaussian",ginverse=list(animal=inv.pruned.tree_sp$Ainv),
                                    prior=prior,data=antdata_multiple_regression.4,nitt=250000,burnin=10000,thin=500)

CasteVs.MF_CS_glm_MCMCglmm.a<-MCMCglmm(Caste3~log(eff.mating.freq.MEAN.harmonic)+log(colony.size),random=~animal,
                                     family="poisson",ginverse=list(animal=inv.pruned.tree_sp$Ainv),
                                     prior=prior,data=antdata_multiple_regression.4,nitt=250000,burnin=10000,thin=500) #nu made smaller

CasteVs.MF_CS_glm_MCMCglmm.b<-MCMCglmm(Caste3~log(eff.mating.freq.MEAN.harmonic)+log(colony.size),random=~animal,
                                       family="poisson",ginverse=list(animal=inv.pruned.tree_sp$Ainv),
                                       prior=prior,data=antdata_multiple_regression.4,nitt=250000,burnin=10000,thin=500) #nu made bigger

CasteVs.MF_CS_glm_MCMCglmm.c<-MCMCglmm(Caste3~log(eff.mating.freq.MEAN.harmonic)+log(colony.size),random=~animal,
                                       family="poisson",ginverse=list(animal=inv.pruned.tree_sp$Ainv),
                                       prior=prior,data=antdata_multiple_regression.4,nitt=1500000,burnin=10000,thin=500) #very large number of iterations with small nu. nit 1000000 was not enough

CasteVs.MF_CS_glm_MCMCglmm.d<-MCMCglmm(Caste3~log(eff.mating.freq.MEAN.harmonic)+log(colony.size),random=~animal,
                                       family="poisson",ginverse=list(animal=inv.pruned.tree_sp$Ainv),
                                       prior=prior,data=antdata_multiple_regression.4,nitt=2000000,burnin=10000,thin=500)
CasteVs.MF_CS_glm_MCMCglmm.e<-MCMCglmm(Caste3~log(eff.mating.freq.MEAN.harmonic)+log(colony.size),random=~animal,
                                       family="poisson",ginverse=list(animal=inv.pruned.tree_sp$Ainv),
                                       prior=prior_ext,data=antdata_multiple_regression.4,nitt=250000,burnin=10000,thin=500) #using prior extension
CasteVs.MF_CS_glm_MCMCglmm.f<-MCMCglmm(Caste3~log(eff.mating.freq.MEAN.harmonic)+log(colony.size),random=~animal,
                                       family="poisson",ginverse=list(animal=inv.pruned.tree_sp$Ainv),
                                       prior=prior_ext,data=antdata_multiple_regression.4,nitt=10000000,burnin=10000,thin=700)#using prior extension and more nit
CasteVs.MF_CS_glm_MCMCglmm.g<-MCMCglmm(Caste3~log(eff.mating.freq.MEAN.harmonic)+log(colony.size),random=~animal,
                                       family="poisson",ginverse=list(animal=inv.pruned.tree_sp$Ainv),
                                       prior=prior.iw,data=antdata_multiple_regression.4,nitt=5000000,burnin=10000,thin=500) #inv.wishart prior - distributions look a lot better!
CasteVs.MF_CS_glm_MCMCglmm.h<-MCMCglmm(Caste3~log(eff.mating.freq.MEAN.harmonic)+log(colony.size),random=~animal,
                                       family="poisson",ginverse=list(animal=inv.pruned.tree_sp$Ainv),
                                       prior=prior.proper,data=antdata_multiple_regression.4,nitt=250000,burnin=10000,thin=500) #using 'proper' priors
CasteVs.MF_CS_glm_MCMCglmm.h<-MCMCglmm(Caste3~log(eff.mating.freq.MEAN.harmonic)+log(colony.size),random=~animal,
                                       family="poisson",ginverse=list(animal=inv.pruned.tree_sp$Ainv),
                                       prior=prior1.1,data=antdata_multiple_regression.4,nitt=250000,burnin=10000,thin=500) #prior1.1
CasteVs.MF_CS_glm_MCMCglmm.i<-MCMCglmm(Caste3~log(eff.mating.freq.MEAN.harmonic)+log(colony.size),random=~animal,
                                       family="poisson",ginverse=list(animal=inv.pruned.tree_sp$Ainv),
                                       prior=prior_exp,data=antdata_multiple_regression.4,nitt=1000000,burnin=10000,thin=500) #prior from course notes
CasteVs.MF_CS_glm_MCMCglmm.j<-MCMCglmm(Caste3~log(eff.mating.freq.MEAN.harmonic)+log(colony.size),random=~animal,
                                       family="poisson",ginverse=list(animal=inv.pruned.tree_sp$Ainv),
                                       prior=prior.chi,data=antdata_multiple_regression.4,nitt=250000,burnin=10000,thin=500) #prior chi
CasteVs.MF_CS_glm_MCMCglmm.k<-MCMCglmm(Caste3~log(eff.mating.freq.MEAN.harmonic)+log(colony.size),random=~animal,
                                       family="poisson",ginverse=list(animal=inv.pruned.tree_sp$Ainv),
                                       prior=prior_exp,data=antdata_multiple_regression.4,nitt=250000,burnin=10000,thin=500, slice = T) #prior from course notes with slice
CasteVs.MF_CS_glm_MCMCglmm.l<-MCMCglmm(Caste3~log(eff.mating.freq.MEAN.harmonic)+log(colony.size),random=~animal,
                                       family="poisson",ginverse=list(animal=inv.pruned.tree_sp$Ainv),
                                       prior=prior_exp.1,data=antdata_multiple_regression.4,nitt=250000,burnin=10000,thin=500, slice = T) #prior_exp with V=10

#Checking for convergence = we should see no trend in the trace (i.e. still increasing)
#Checking for autocorrelation = values in the trace are widely spread (goes up and down) + autocorr.diag() function gives values all below 0.1 from the 1st lag
##When autocorrelation is high, the chain needs to be run for longer and thin increased
###nitt/thin = 1000-2000 and Lag < 0.1
#Autocorrelation is not the problem per say, but it reduces the effective size parameter which is a problem
#Curves of the trace should be unimodal, symmetrical and not aberrant
#To get a larger effective sample size, run the chain for longer (>10,000 sample size)
#Posterior distribution refers to the plots that you get from the MCMCglmm model (lines and histogram - give the values in the summary output!)
#Posterior distributions of the variance components refers to the plots titled 'animal' and 'units' - these are more interesting than the 
#The variance components:
##Va = animal = additive genetic variance
##Vr = residual variance = units
#The influence of the prior distribution fade away with a sufficient sample size.
#Model selection with DIC- DIC are not well understood and that the DIC may be focused at the wrong level for most people’s intended level of inference - particularly with non-Gaussian responses
#This prior specification used to be used a lot because it was believed to be relatively uninformative, and is equivalent to an inverse-gamma prior with shape and scale equal to 0.001. In many cases it is relatively uninformative but when the posterior distribution for the variances has suport close to zero it can behave poorly.

summary(CasteVs.MF_CS_glm_MCMCglmm.i) #Despite using quite different priors, coefficients remain roughly the same (sample size must be large enough that prior specification doesn't matter too much)
autocorr.diag(CasteVs.MF_CS_glm_MCMCglmm.i$Sol)
autocorr.diag(CasteVs.MF_CS_glm_MCMCglmm.i$VCV)
posterior.mode(CasteVs.MF_CS_glm_MCMCglmm.i$VCV)
HPDinterval(CasteVs.MF_CS_glm_MCMCglmm.i$VCV)
effectiveSize(CasteVs.MF_CS_glm_MCMCglmm.i$Sol)
effectiveSize(CasteVs.MF_CS_glm_MCMCglmm.i$VCV)

heidel.diag(CasteVs.MF_CS_glm_MCMCglmm.g$VCV)

herit <- CasteVs.MF_CS_glm_MCMCglmm.g$VCV[, "animal"]/(CasteVs.MF_CS_glm_MCMCglmm.g$VCV[, "animal"] + CasteVs.MF_CS_glm_MCMCglmm.g$VCV[, "units"])
effectiveSize(herit)
mean(herit)
HPDinterval(herit) #Display 95% credible interval
plot(herit)

##brms models:
animal <- pruned.tree_sp_multiple_regression
C <- ape::vcv.phylo(animal)


#MF
brms_negbinom_MF <-brm(
  Caste3 ~ log(eff.mating.freq.MEAN.harmonic) + (1|gr(animal, cov = C)),
  data = antdata_multiple_regression.4, family = negbinomial(), 
  data2 = list(C = C),
  chains = 2, cores = 2, iter = 4000,
  control = list(adapt_delta = 0.95))
summary(brms_negbinom_MF)

#CS
brms_negbinom_CS <-brm(
  Caste3 ~ log(colony.size) + (1|gr(animal, cov = C)),
  data = antdata_multiple_regression.4, family = negbinomial(), 
  data2 = list(C = C),
  chains = 2, cores = 2, iter = 4000,
  control = list(adapt_delta = 0.95))
summary(brms_negbinom_CS)


#MF + CS
brms_negbinom_MF_CS <-brm(
  Caste3 ~ log(eff.mating.freq.MEAN.harmonic) + log(colony.size) + (1|gr(animal, cov = C)),
  data = antdata_multiple_regression.4, family = negbinomial(), 
  data2 = list(C = C),
  prior = c(
    prior(normal(0, 10), "b"),
    prior(normal(0, 50), "Intercept"),
    prior(student_t(3, 0, 20), "sd")
    ),
  chains = 2, cores = 2, iter = 10000, warmup = 1000,
  control = list(adapt_delta = 0.95))
summary(brms_negbinom_MF_CS, waic = TRUE)

prior <- get_prior(Caste3 ~ log(eff.mating.freq.MEAN.harmonic) + log(colony.size) + (1|gr(animal, cov = C)),
                    data = antdata_multiple_regression.4, family = gaussian(), 
                   data2 = list(C = C)
)

#MF + CS - normal distribution
brms_norm_MF_CS <-brm(
  Caste3 ~ log(eff.mating.freq.MEAN.harmonic) + log(colony.size) + (1|gr(animal, cov = C)),
  data = antdata_multiple_regression.4, family = gaussian(), 
  data2 = list(C = C),
  chains = 2, cores = 2, iter = 4000,
  control = list(adapt_delta = 0.95))
summary(brms_negbinom_MF_CS, waic = TRUE)

#MF + CS
brms_negbinom_MF_CS <-brm(
  Caste3 ~ log(eff.mating.freq.MEAN.harmonic) + log(colony.size) + (1|gr(animal, cov = C)),
  data = antdata_multiple_regression.4, family = negbinomial(), 
  data2 = list(C = C),
  chains = 2, cores = 2, iter = 4000,
  control = list(adapt_delta = 0.95))
summary(brms_negbinom_MF_CS)

#MF * CS
brms_negbinom_MF_CS_interaction <-brm(
  Caste3 ~ log(eff.mating.freq.MEAN.harmonic) * log(colony.size) + (1|gr(animal, cov = C)),
  data = antdata_multiple_regression.4, family = negbinomial(), 
  data2 = list(C = C),
  chains = 2, cores = 2, iter = 4000,
  control = list(adapt_delta = 0.95))
summary(brms_negbinom_MF_CS_interaction)

#MF + PG (binary)
brms_negbinom_MF_PG_bin <-brm(
  Caste3 ~ log(eff.mating.freq.MEAN.harmonic) + polygyny.clean + (1|gr(animal, cov = C)),
  data = antdata_multiple_regression.4, family = negbinomial(), 
  data2 = list(C = C),
  chains = 2, cores = 2, iter = 4000,
  control = list(adapt_delta = 0.95))
summary(brms_negbinom_MF_PG_bin)

#MF * PG (binary)
brms_negbinom_MF_PG_bin_interaction <-brm(
  Caste3 ~ log(eff.mating.freq.MEAN.harmonic) * polygyny.clean + (1|gr(animal, cov = C)),
  data = antdata_multiple_regression.4, family = negbinomial(), 
  data2 = list(C = C),
  chains = 2, cores = 2, iter = 4000,
  control = list(adapt_delta = 0.95))
summary(brms_negbinom_MF_PG_bin_interaction)

#MF + PG (continuous)
brms_negbinom_MF_PG_cont <-brm(
  Caste3 ~ log(eff.mating.freq.MEAN.harmonic) + number.queens.MEAN + (1|gr(animal, cov = C)),
  data = antdata_multiple_regression.4, family = negbinomial(), 
  data2 = list(C = C),
  chains = 2, cores = 2, iter = 4000,
  control = list(adapt_delta = 0.95))
summary(brms_negbinom_MF_PG_cont)

#MF * PG (continuous)
brms_negbinom_MF_PG_cont_interaction <-brm(
  Caste3 ~ log(eff.mating.freq.MEAN.harmonic) * number.queens.MEAN + (1|gr(animal, cov = C)),
  data = antdata_multiple_regression.4, family = negbinomial(), 
  data2 = list(C = C),
  chains = 2, cores = 2, iter = 4000,
  control = list(adapt_delta = 0.95))
summary(brms_negbinom_MF_PG_cont_interaction)

#CS + PG (binary)
brms_negbinom_CS_PG_bin <-brm(
  Caste3 ~ log(colony.size) + polygyny.clean + (1|gr(animal, cov = C)),
  data = antdata_multiple_regression.4, family = negbinomial(), 
  data2 = list(C = C),
  chains = 2, cores = 2, iter = 4000,
  control = list(adapt_delta = 0.95))
summary(brms_negbinom_CS_PG_bin)

#CS * PG (binary)
brms_negbinom_CS_PG_bin_interaction <-brm(
  Caste3 ~ log(colony.size) * polygyny.clean + (1|gr(animal, cov = C)),
  data = antdata_multiple_regression.4, family = negbinomial(), 
  data2 = list(C = C),
  chains = 2, cores = 2, iter = 4000,
  control = list(adapt_delta = 0.95))
summary(brms_negbinom_CS_PG_bin_interaction)

#CS + PG (continuous)
brms_negbinom_CS_PG_cont <-brm(
  Caste3 ~ log(colony.size) + number.queens.MEAN + (1|gr(animal, cov = C)),
  data = antdata_multiple_regression.4, family = negbinomial(), 
  data2 = list(C = C),
  chains = 2, cores = 2, iter = 4000,
  control = list(adapt_delta = 0.95))
summary(brms_negbinom_CS_PG_cont)

#CS * PG (continuous)
brms_negbinom_CS_PG_cont_interaction <-brm(
  Caste3 ~ log(colony.size) * number.queens.MEAN + (1|gr(animal, cov = C)),
  data = antdata_multiple_regression.4, family = negbinomial(), 
  data2 = list(C = C),
  chains = 2, cores = 2, iter = 4000,
  control = list(adapt_delta = 0.95))
summary(brms_negbinom_CS_PG_cont_interaction)

#MF + CS + PG (binary)
brms_negbinom_MF_CS_PG_bin <-brm(
  Caste3 ~ log(eff.mating.freq.MEAN.harmonic) + log(colony.size) + polygyny.clean + (1|gr(animal, cov = C)),
  data = antdata_multiple_regression.4, family = negbinomial(), 
  data2 = list(C = C),
  chains = 2, cores = 2, iter = 4000,
  control = list(adapt_delta = 0.95))
summary(brms_negbinom_MF_CS_PG_bin)

#MF + CS * PG (binary)
brms_negbinom_MF_CS_int_PG_bin <-brm(
  Caste3 ~ log(eff.mating.freq.MEAN.harmonic) + log(colony.size) * polygyny.clean + (1|gr(animal, cov = C)),
  data = antdata_multiple_regression.4, family = negbinomial(), 
  data2 = list(C = C),
  chains = 2, cores = 2, iter = 4000,
  control = list(adapt_delta = 0.95))
summary(brms_negbinom_MF_CS_int_PG_bin)

#MF * CS * PG (binary)
brms_negbinom_MF_CS_PG_bin_interaction <-brm(
  Caste3 ~ log(eff.mating.freq.MEAN.harmonic) * log(colony.size) + polygyny.clean + (1|gr(animal, cov = C)),
  data = antdata_multiple_regression.4, family = negbinomial(), 
  data2 = list(C = C),
  chains = 2, cores = 2, iter = 4000,
  control = list(adapt_delta = 0.95))
summary(brms_negbinom_MF_CS_PG_bin_interaction)

#MF + CS + PG (continuous)
brms_negbinom_MF_CS_PG_cont <-brm(
  Caste3 ~ log(eff.mating.freq.MEAN.harmonic) + log(colony.size) + number.queens.MEAN + (1|gr(animal, cov = C)),
  data = antdata_multiple_regression.4, family = negbinomial(), 
  data2 = list(C = C),
  chains = 2, cores = 2, iter = 4000,
  control = list(adapt_delta = 0.95))
summary(brms_negbinom_MF_CS_PG_cont)

#MF * CS * PG (continuous)
brms_negbinom_MF_CS_PG_cont_interaction <-brm(
  Caste3 ~ log(eff.mating.freq.MEAN.harmonic) * number.queens.MEAN + log(colony.size) * number.queens.MEAN + (1|gr(animal, cov = C)),
  data = antdata_multiple_regression.4, family = negbinomial(), 
  data2 = list(C = C),
  chains = 2, cores = 2, iter = 4000,
  control = list(adapt_delta = 0.95))
summary(brms_negbinom_MF_CS_PG_cont_interaction)


#Model selection using looAic - Bayesian version of Aic
##All very similar but Caste ~ MF + CS seems to be the best
loo(brms_negbinom_MF, brms_negbinom_CS, brms_negbinom_MF_CS, brms_negbinom_MF_PG_bin, brms_negbinom_CS_PG_bin, brms_negbinom_MF_CS_PG_bin)
loo_compare(loo(brms_negbinom_MF), loo(brms_negbinom_CS), loo(brms_negbinom_MF_CS), loo(brms_negbinom_MF_PG_bin), loo(brms_negbinom_CS_PG_bin), loo(brms_negbinom_MF_CS_PG_bin), loo(brms_negbinom_MF_CS_int_PG_bin))

loo(brms_negbinom_MF, brms_negbinom_CS, brms_negbinom_MF_CS, brms_negbinom_MF_CS_interaction, brms_negbinom_MF_PG_bin, brms_negbinom_MF_PG_bin_interaction, brms_negbinom_MF_PG_cont, 
    brms_negbinom_MF_PG_cont_interaction, brms_negbinom_CS_PG_bin, brms_negbinom_CS_PG_bin_interaction, brms_negbinom_CS_PG_cont, 
    brms_negbinom_CS_PG_cont_interaction, brms_negbinom_MF_CS_PG_bin, brms_negbinom_MF_CS_PG_bin_interaction, brms_negbinom_MF_CS_PG_cont, 
    brms_negbinom_MF_CS_PG_cont_interaction)



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
CasteVs.CS_MF_PG_Phy_glm <- phyloglm(Caste3~eff.mating.freq.MEAN.harmonic*log(colony.size),
                                     data = antdata_CS_MF.4, 
                                     phy = pruned.tree_sp_CS_MF, method = "poisson_GEE")
CasteVs.CS_MF_PG_Phy_lm <- phylolm(Caste3~log(colony.size) * eff.mating.freq.MEAN.harmonic + polygyny.clean,
                                     data = antdata_CS_MF.4, 
                                     phy = pruned.tree_sp_CS_MF)

summary(CasteVs.CS_MF_Phy_lm)
summary(CasteVs.CS_MF_Phy_glm)
summary(CasteVs.CS_MF_lm)
summary(CasteVs.CS_MF_PG_Phy_glm)
summary(CasteVs.CS_MF_PG_Phy_lm)

#brms: caste ~ MF + PG
#brms models
antdata_CS_MF.4$obs <- 1:nrow(antdata_CS_MF.4)
antdata_CS.5 <- dplyr::select(antdata_CS_MF.4, Caste3, eff.mating.freq.MEAN.harmonic, polygyny.clean, animal, obs)
animal <- pruned.tree_sp_CS_MF
C <- ape::vcv.phylo(animal)

#The '(1|obs)' bit is required as the Poisson distribution does not have a natural overdispersion parameter. 
##We model the residual variance via the group-level effects of obs
##1st model: poisson
model_pois_MF_PG <-brm(
  Caste3 ~ eff.mating.freq.MEAN.harmonic * polygyny.clean + (1|gr(animal, cov = C)) + (1|obs),
  data = antdata_CS.5, family = poisson("log"), 
  data2 = list(C = C),
  chains = 2, cores = 2, iter = 4000,
  control = list(adapt_delta = 0.95)
)

summary(model_pois_MF_PG)
plot(conditional_effects(model_pois_MF_PG), points = TRUE)


#VIF test - VIF scores seem to be low...
##MF ~ CS analysis
CasteVs.CS_MF_lm <- lm(log(colony.size) ~ eff.mating.freq.MEAN.harmonic, data = antdata_multiple_regression)
CasteVs.CS_MF_glm_poisson <- glm(log(colony.size) ~ eff.mating.freq.MEAN.harmonic, family="poisson",data = antdata_CS_MF.4)

plot(antdata_CS_MF.4$eff.mating.freq.MEAN.harmonic, log(antdata_CS_MF.4$colony.size))
abline(CasteVs.CS_MF_lm)
summary(CasteVs.CS_MF_lm)

## Caste ~ MF + CS
CasteVs.CS_MF_lm.1 <- lm(Caste3 ~ log(colony.size) + log(eff.mating.freq.MEAN.harmonic) + number.queens.MEAN, data = antdata_multiple_regression)
CasteVs.CS_MF_glm_poisson.1 <- glm(Caste3 ~ log(colony.size) + log(eff.mating.freq.MEAN.harmonic), family="poisson",data = antdata_multiple_regression)
summary(CasteVs.CS_MF_lm.1)
summary(CasteVs.CS_MF_glm_poisson.1)
car::vif(CasteVs.CS_MF_lm.1)
car::vif(CasteVs.CS_MF_glm_poisson.1)

## Caste ~ MF + CS + PG (binary)
CasteVs.CS_MF_lm.2 <- lm(Caste3 ~ log(colony.size) + log(eff.mating.freq.MEAN.harmonic) + polygyny.clean, data = antdata_multiple_regression)
CasteVs.CS_MF_glm_poisson.2 <- glm(Caste3 ~ log(colony.size) + log(eff.mating.freq.MEAN.harmonic) + polygyny.clean, family="poisson",data = antdata_multiple_regression)
plot(log(antdata_multiple_regression$eff.mating.freq.MEAN.harmonic), log(antdata_multiple_regression$colony.size))
CasteVs.PG_MF_lm <- lm(polygyny.clean ~ log(eff.mating.freq.MEAN.harmonic), data = antdata_multiple_regression)
plot(log(antdata_multiple_regression$eff.mating.freq.MEAN.harmonic), antdata_multiple_regression$polygyny.clean)
abline(CasteVs.PG_MF_lm)
CasteVs.PG_CS_lm <- lm(polygyny.clean ~ log(colony.size), data = antdata_multiple_regression)
plot(log(antdata_multiple_regression$colony.size), antdata_multiple_regression$polygyny.clean)
abline(CasteVs.PG_CS_lm)
summary(CasteVs.CS_MF_lm.2)
summary(CasteVs.CS_MF_glm_poisson.2)
car::vif(CasteVs.CS_MF_lm.2)
car::vif(mod = CasteVs.CS_MF_glm_poisson.2)


#MCMCglmm
inv.pruned.tree_sp.1<-inverseA(pruned.tree_sp_CS_MF,nodes="TIPS",scale=TRUE)
prior<-list(G=list(G1=list(V=1,nu=0.02)),R=list(V=1,nu=0.02))
CasteVs.MF_glm_MCMCglmm.1<-MCMCglmm(Caste3~eff.mating.freq.MEAN.harmonic+log(colony.size),random=~animal,
                                    family="poisson",ginverse=list(animal=inv.pruned.tree_sp.1$Ainv),
                                    prior=prior,data=antdata_CS_MF.4,nitt=250000,burnin=10000,thin=500)

CasteVs.MF_glm_MCMCglmm.2<-MCMCglmm(Caste3~eff.mating.freq.MEAN.harmonic*log(colony.size),random=~animal,
                                    family="poisson",ginverse=list(animal=inv.pruned.tree_sp.1$Ainv),
                                    prior=prior,data=antdata_CS_MF.4,nitt=250000,burnin=10000,thin=500)

CasteVs.MF_glm_MCMCglmm.3<-MCMCglmm(Caste3~eff.mating.freq.MEAN.harmonic*polygyny.clean,random=~animal,
                                    family="poisson",ginverse=list(animal=inv.pruned.tree_sp.1$Ainv),
                                    prior=prior,data=antdata_CS_MF.4,nitt=250000,burnin=10000,thin=500)

CasteVs.MF_glm_MCMCglmm.4<-MCMCglmm(Caste3~eff.mating.freq.MEAN.harmonic+log(colony.size)+polygyny.clean,random=~animal,
                                    family="poisson",ginverse=list(animal=inv.pruned.tree_sp.1$Ainv),
                                    prior=prior,data=antdata_CS_MF.4,nitt=250000,burnin=10000,thin=500)

CasteVs.MF_glm_MCMCglmm.5<-MCMCglmm(Caste3~eff.mating.freq.MEAN.harmonic+polygyny.clean,random=~animal,
                                    family="poisson",ginverse=list(animal=inv.pruned.tree_sp.1$Ainv),
                                    prior=prior,data=antdata_CS_MF.4,nitt=250000,burnin=10000,thin=500)

CasteVs.MF_glm_MCMCglmm.6<-MCMCglmm(Caste3~eff.mating.freq.MEAN.harmonic*polygyny.clean+log(colony.size),random=~animal,
                                    family="poisson",ginverse=list(animal=inv.pruned.tree_sp.1$Ainv),
                                    prior=prior,data=antdata_CS_MF.4,nitt=250000,burnin=10000,thin=500)

CasteVs.MF_glm_MCMCglmm.7<-MCMCglmm(Caste3~polygyny.clean+log(colony.size),random=~animal,
                                    family="poisson",ginverse=list(animal=inv.pruned.tree_sp.1$Ainv),
                                    prior=prior,data=antdata_CS_MF.4,nitt=250000,burnin=10000,thin=500)

CasteVs.MF_glm_MCMCglmm.8<-MCMCglmm(Caste3~polygyny.clean*log(colony.size),random=~animal,
                                    family="poisson",ginverse=list(animal=inv.pruned.tree_sp.1$Ainv),
                                    prior=prior,data=antdata_CS_MF.4,nitt=250000,burnin=10000,thin=500)

CasteVs.MF_glm_MCMCglmm.9<-MCMCglmm(Caste3~polygyny.clean*log(colony.size)+eff.mating.freq.MEAN.harmonic,random=~animal,
                                    family="poisson",ginverse=list(animal=inv.pruned.tree_sp.1$Ainv),
                                    prior=prior,data=antdata_CS_MF.4,nitt=250000,burnin=10000,thin=500)

CasteVs.MF_glm_MCMCglmm.10<-MCMCglmm(Caste3~log(eff.mating.freq.MEAN.harmonic)+polygyny.clean,random=~animal,
                                    family="poisson",ginverse=list(animal=inv.pruned.tree_sp.1$Ainv),
                                    prior=prior,data=antdata_CS_MF.4,nitt=250000,burnin=10000,thin=500)


heidel.diag(CasteVs.MF_glm_MCMCglmm$Sol) #GW: Have to look closer at how these work. 
plot(CasteVs.MF_glm_MCMCglmm)
summary(CasteVs.MF_glm_MCMCglmm.1)
summary(CasteVs.MF_glm_MCMCglmm.2)
summary(CasteVs.MF_glm_MCMCglmm.3)
summary(CasteVs.MF_glm_MCMCglmm.4)
summary(CasteVs.MF_glm_MCMCglmm.5)
summary(CasteVs.MF_glm_MCMCglmm.6)
summary(CasteVs.MF_glm_MCMCglmm.7)
summary(CasteVs.MF_glm_MCMCglmm.8)
summary(CasteVs.MF_glm_MCMCglmm.9)
summary(CasteVs.MF_glm_MCMCglmm.10)

#Significant results arise when controlling for polygyny but only when either CS or MF are excluded from the model as well.

#Plotting the interaction between MF and CS on Caste
A=25
B=0:15

y=-0.0985621+0.0990318*A+0.0434442*B+(-0.0095733*A*B)
plot(y, type="l", col="blue", lwd=1, xlab="CS", ylab="Caste", ylim = c(-5,5))
lines(y, col="pink", lwd=1)





