# MCMC Phylogenetic Generalised Linear Mixed Models (MCMCglmm)
#17/02/2021 Louis Bell-Roberts

setwd("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Trees/Genus_tree/")

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

##Read in the tree and database##
fulltree <- read.tree(file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Trees/Nelsen/ML_TREE_treepl_185_outsdroppeds.ladderized.increasing.tre")
anttree <- read.tree(file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Trees/Genus_tree/Genus_Level_ML_TREE_treepl_185_outsdroppeds.ladderized.increasing.tre")
antdata <- read.csv("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Cleaned/Data.csv", header = T)

#Check if the original tree is ultrametric - visually it looks ultrametric
plotTree(anttree,ftype="i",fsize=0.4,lwd=1)
#However, returns false
is.ultrametric(fulltree)

#Dataframe with genus added as additional column
data_with_genus <- as.data.frame(
  cbind(
    genus = word(antdata$animal,sep = "_"),
    antdata
  )
) 

View(data_with_genus)

#Select the columns that I want to use
alltaxa <- data_with_genus %>% select(genus, eff.mating.freq.MEAN.harmonic, colony.size, Caste1)
head(alltaxa)

#Don't need to select only the ants from the dataframe as all other taxa will be removed when matching with my ant phylogeny

#Check if my variables are set as numeric - change if not
is.numeric(alltaxa$eff.mating.freq.MEAN.harmonic)
is.numeric(alltaxa$colony.size)
is.numeric(alltaxa$Caste1)

#Set caste as numeric
alltaxa1 <- alltaxa
alltaxa1$Caste1 <- as.numeric(as.character(alltaxa1$Caste1))
View(alltaxa1)

#Match the tip labels from the tree and the species from the database and prune the phylogeny accordingly
#It is the material that is in the first named set, that is not in the second named set
pruned.tree<-drop.tip(anttree, setdiff(anttree$tip.label, alltaxa1$genus))
pruned.tree
plotTree(pruned.tree,ftype="i",fsize=0.4,lwd=1)

#Filter through my dataframe and select only the rows that match the tips of my tree
data<-filter(alltaxa1, genus %in% pruned.tree$tip.label)
View(data)



##########Creating the model#############
#Looking at the textbook data files:
a<-read.table(file.choose())
View(a)
b<-read.table(file.choose())
View(b)


#Ensure that tree is ultrametric
u.tree <- force.ultrametric(pruned.tree, method = c("extend"))

#Create inverse matrix of phylogeny
inv.phylo<-inverseA(u.tree,nodes="TIPS",scale=TRUE)

#Ensure all of my data is named correctly so that it matches
##the model would only run if I rename the column heading 'genus' to phylo and 'pruned.tree' to phylo - maybe both of these things need to have the same name
data1<-data %>% 
  rename(
    phylo = genus,
    MF = eff.mating.freq.MEAN.harmonic,
    CS = colony.size
  )
phylo <- pruned.tree

#select only the phylo, Caste1 and MF column
data2 <- data1 %>%
  select(phylo, MF, Caste1)

#Create my default prior
prior<-list(G=list(G1=list(V=1,nu=0.02)),R=list(V=1,nu=0.02)) #G - random, R - fixed

#Ensure that there is no missing data in the predictor variables
data3 <- filter(data2, MF>=0)

####Convert data to genus-level averages##### - For use in the model where I create my own genus-level averages and it has single random effect
#Example code from textbook: data4$spec_mean_cf<-sapply(split(data3$MF,data3$phylo),mean)[data3$phylo]
l<-sapply(split(data3$MF,data3$phylo), mean)
k<-sapply(split(data3$Caste1,data3$phylo),mean, na.rm=TRUE)
j<-cbind(l,k)
h<-as.data.frame(j)
g<-rownames_to_column(h, var = "phylo")
data11<-g %>% 
  rename(
    MF = l,
    Caste = k
  )
View(data11)
data12<-data11
data12$Caste <- na_if(data11$Caste, "NaN")
head(data12)
#############################################

####IMPORTANT - I don't need to create genus-level averages when running this type of model: MCMCglmm will do that for me if I want####
#1ST MODEL - Genus-level average; Gaussian; single fixed effect; single random effect
model_simple<-MCMCglmm(Caste~MF,random=~phylo,
                       family="gaussian",ginverse=list(phylo=inv.phylo$Ainv),prior=prior,
                       data=data12,nitt=1100000,burnin=10000,thin=1000)

summary(model_simple)
posterior.mode(model_simple$Sol)
plot(model_simple)
heidel.diag(model_simple$Sol)
############################################

#2ND MODEL - Genus-level average; Poisson; single fixed effect; single random effect - DID NOT WORK
##Using the same prior and phylo as in model 1

#model_simple_poiss<-MCMCglmm(Caste~MF,random=~phylo,
 #                      family="poisson",ginverse=list(phylo=inv.phylo$Ainv),prior=prior,
  #                     data=data12,nitt=1100000,burnin=10000,thin=1000)

#Run again with gaussian but with log transformation this time - THIS DOES WORK
model_simple_poiss<-MCMCglmm(log(Caste)~log(MF),random=~phylo,
                             family="gaussian",ginverse=list(phylo=inv.phylo$Ainv),prior=prior,
                             data=data12,nitt=1100000,burnin=10000,thin=1000)


summary(model_simple_poiss)
posterior.mode(model_simple_poiss$Sol)
plot(model_simple_poiss)
heidel.diag(model_simple_poiss$Sol)

############################################

#3RD MODEL - Multiple measurements per genus; Gaussian; single fixed effect; 2 random effects (phylo + species)
View(data3)

#Add additional column identical to the phylo column and name it species 
species<-data3$phylo
data4<-cbind(species, data3)
View(data4)

#Calculate the specific mean of the predictor variable (MF)
data4$spec_mean_MF<-sapply(split(data4$MF,data4$phylo),mean)[data4$phylo]

#Inverse phylo matrix created previously in the script. Set the new priors with two random effects parameters
prior3<-list(G=list(G1=list(V=1,nu=0.02),G2=list(V=1,nu=0.02)),
             R=list(V=1,nu=0.02)) #G - random effects; R - fixed effects

mcmcmodel_3<-MCMCglmm(Caste1~spec_mean_MF,random=~phylo+species,
                        family="gaussian",ginverse=list(phylo=inv.phylo$Ainv),
                        prior=prior3,data=data4,nitt=1100000,burnin=10000,thin=1000)

summary(mcmcmodel_3)
posterior.mode(mcmcmodel_3$Sol)
plot(mcmcmodel_3)
heidel.diag(mcmcmodel_3$Sol)

#In order to get an estimate for the ``between-species'' and ``within-species'', we need to use the within-group centering technique:
data5<-data4
data5$within_spec_MF<-data5$MF-data4$spec_mean_MF

#Elaborated model
mcmcmodel_3_elaborated<-MCMCglmm(Caste1~spec_mean_MF+within_spec_MF,random=~phylo+species,
                        family="gaussian",ginverse=list(phylo=inv.phylo$Ainv),
                        prior=prior3,data=data5,nitt=1100000,burnin=10000,thin=1000)

summary(mcmcmodel_3_elaborated)
posterior.mode(mcmcmodel_3_elaborated$Sol)
plot(mcmcmodel_3_elaborated)
heidel.diag(mcmcmodel_3_elaborated$Sol)

#Calculate Lambda using:
lambda <- mcmcmodel_3_elaborated$VCV[,'phylo']/
  (mcmcmodel_3_elaborated$VCV[,'phylo']+mcmcmodel_3_elaborated$VCV[,'species']+
     mcmcmodel_3_elaborated$VCV[,'units'])

summary(lambda)
mean(lambda)
posterior.mode(lambda)
HPDinterval(lambda)
############################################

#4TH MODEL - Multiple measurements per genus; Poisson; single fixed effect; 2 random effects (phylo + species)
##Run the model 3 again but change the family to poisson
mcmcmodel_4_poisson<-MCMCglmm(Caste1~spec_mean_MF,random=~phylo+species,
                      family="poisson",ginverse=list(phylo=inv.phylo$Ainv),
                      prior=prior3,data=data4,nitt=1100000,burnin=10000,thin=1000)

summary(mcmcmodel_4_poisson)
posterior.mode(mcmcmodel_4_poisson$Sol)
plot(mcmcmodel_4_poisson)
heidel.diag(mcmcmodel_4_poisson$Sol)

#Run model 3_elaborated again but change the family to poisson

#Elaborated model
mcmcmodel_4_poisson_elaborated<-MCMCglmm(Caste1~spec_mean_MF+within_spec_MF,
                                 random=~phylo+species,family="poisson",
                                 ginverse=list(phylo=inv.phylo$Ainv),prior=prior3,data=data5,
                                 nitt=1100000,burnin=10000,thin=1000)

summary(mcmcmodel_4_poisson_elaborated)
posterior.mode(mcmcmodel_4_poisson_elaborated$Sol)
plot(mcmcmodel_4_poisson_elaborated)
heidel.diag(mcmcmodel_4_poisson_elaborated$Sol)

#Calculate Lambda using:
lambda <- mcmcmodel_4_poisson_elaborated$VCV[,'phylo']/
  (mcmcmodel_4_poisson_elaborated$VCV[,'phylo']+mcmcmodel_4_poisson_elaborated$VCV[,'species']+
     mcmcmodel_4_poisson_elaborated$VCV[,'units'])

summary(lambda)
mean(lambda)
posterior.mode(lambda)
HPDinterval(lambda)
############################################

#5th MODEL - genus-level average; complete case analysis; Gaussian; two fixed effects (MF + CS); 1 random effect (phylo)
View(data12)

#Create dataframe with genus-level averages for MF CS and Caste: data14

##Filter the data so that I only have species with complete information
data13<-filter(data1, MF>=0, CS>=0, Caste1>=0)
View(data13)

#Create genus-level averages
data13.1<-sapply(split(data13$MF,data13$phylo), mean)
data13.2<-sapply(split(data13$CS,data13$phylo), mean)
data13.3<-sapply(split(data13$Caste1,data13$phylo), mean)

#Combine into new dataframe
data13.4<-as.data.frame(cbind(data13.1, data13.2, data13.3))
data13.5<-rownames_to_column(data13.4, var = "phylo")

#Rename the columns
data14<-data13.5 %>% 
  rename(
    MF = data13.1,
    CS = data13.2,
    Caste = data13.3
  )

prior<-list(G=list(G1=list(V=1,nu=0.02)),R=list(V=1,nu=0.02)) #G - random, R - fixed
mcmcmodel_5<-MCMCglmm(Caste~MF+CS,random=~phylo,
                                 family="gaussian",ginverse=list(phylo=inv.phylo$Ainv),
                                 prior=prior,data=data14,nitt=1100000,burnin=10000,thin=1000)

summary(mcmcmodel_5)
posterior.mode(mcmcmodel_5$Sol)
plot(mcmcmodel_5)
heidel.diag(mcmcmodel_5$Sol)
############################################

#6th MODEL - genus-level average; complete case analysis; Poisson; two fixed effect (MF + CS); 1 random effect (phylo) - DID NOT WORK
##Run the same model as model 5, but change the distribution family to poisson

plot(data14$MF, data14$Caste)
plot(log(data14$CS), data14$Caste)
hist(data14$Caste)

lm<-lm(Caste^5~MF, data = data14)
anova(lm)
summary(lm)

pred_data <- data.frame(MF = seq(0, 20, length.out = 20))
predict(lm, pred_data)
pred_data <- mutate(pred_data, Caste = predict(lm, pred_data))

ggplot(pred_data, aes(x = MF, y = Caste)) + 
  geom_line() + geom_point(data = data14) + 
  xlab("MF") + ylab("Caste") + 
  theme_bw(base_size = 22)

#Checking the linearity assumption
plt_data <- 
  data.frame(Fitted = fitted(lm), 
             Resids =  resid(lm))

ggplot(plt_data, aes(x = Fitted, y = Resids)) + 
  geom_point() + 
  xlab("Fitted values") + ylab("Residuals")

#Checking the assumption of normality
mod_resids <- resid(lm)
mod_resids <- mod_resids / sd(mod_resids)
resid_order <- order(mod_resids)
resid_order

all_resids <- qqnorm(mod_resids, plot.it = FALSE)
all_resids <- as.data.frame(all_resids)

#Create normal probability plot
ggplot(all_resids, aes(x = x, y = y)) + 
  geom_point() + geom_abline(intercept = 0, slope = 1) +
  xlab("Theoretical Value") + ylab("Standardised Residual")

#Check for Constant variance
# extract the residuals
sqrt_abs_resids <- resid(lm)
# step 1. standardise them
sqrt_abs_resids <- sqrt_abs_resids / sd(sqrt_abs_resids)
# step 2. find their absolute value
sqrt_abs_resids <- abs(sqrt_abs_resids)
# step 3. square root these
sqrt_abs_resids <- sqrt(sqrt_abs_resids)

plt_data <- 
  data.frame(Fitted = fitted(lm), Resids = sqrt_abs_resids)

ggplot(plt_data, aes(x = Fitted, y = Resids)) + 
  geom_point() + 
  xlab("Fitted values") + ylab("Square root of absolute residuals")

a.data<-mutate(data14, logCaste=Caste^3)
lm_log <- lm(logCaste ~ MF, data = a.data)

plot(lm_log, which = 3)

#linear model
lm<-lm(Caste^3.5~MF, data = data14)
lm<-lm(Caste~MF, data = data14)
lm<-lm(log(Caste)~MF, data = data14)

#Create histogram plot of my residuals
g3 <- qplot(lm$residuals,
            geom = "histogram",
            bins = 10) +
  labs(title = "Histogram of residuals",
       x = "residual")
g3


##Change from poisson distribution to gaussian but log transform the variables
prior<-list(G=list(G1=list(V=1,nu=0.02)),R=list(V=1,nu=0.02)) #G - random, R - fixed
mcmcmodel_6_poisson<-MCMCglmm(Caste~MF+CS,random=~phylo,
                      family="poisson",ginverse=list(phylo=inv.phylo$Ainv),
                      prior=prior,data=data14,nitt=1100000,burnin=10000,thin=1000)

summary(mcmcmodel_6_poisson)
posterior.mode(mcmcmodel_6_poisson$Sol)
plot(mcmcmodel_6_poisson)
heidel.diag(mcmcmodel_6_poisson$Sol)

############################################

#7th MODEL - multiple measurements per genus; Gaussian; two fixed effect (MF + CS); 2 random effects (phylo + species)
#Cannot have missing data in the predictor variables - remove these species

############################################

#8th MODEL - multiple measurements per genus; Poisson; two fixed effect (MF + CS); 2 random effects (phylo + species)

############################################



yu<-filter(data4, Caste1>=0)
yu
mean(yu$Caste1)
var(yu$Caste1)
plot(log(data13$CS), data13$Caste1)

ggplot(data=data13, aes(x=CS, y=Caste1))+
  geom_jitter(width = 0.05, height = 0.05)+
  geom_smooth(method = 'glm')



myglm<-glm(Caste1 ~ MF, family = "poisson", data = data13)
myglm.1<-glm(Caste1 ~ log(CS), family = "poisson", data = data13)
myglm1<-glm(Caste1 ~ MF + log(CS) + MF*log(CS), family = "poisson", data = data13)
summary(myglm)
summary(myglm.1)
summary(myglm1)

cor(data13$MF, data13$CS, method = 'pearson')
plot(myglm.1)
vif(myglm1)



