#Phylogenetic Generalised Least Squares Regression (PGLS)
#03/02/2021 Louis Bell-Roberts

setwd("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Trees/Genus_tree/")

library(phangorn)
library(ggpubr)
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
library(data.table)
library(tidyverse)
library(rcompanion)

##Read in the tree and database##
anttree <- read.tree(file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Trees/Genus_tree/Genus_Level_ML_TREE_treepl_185_outsdroppeds.ladderized.increasing.tre")
antdata <- read.csv("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Cleaned/Data.csv", header = T)

#The tree needs to be ultrametric
u.tree<-force.ultrametric(pruned.tree, method=c("extend"))
#Need to calculate the inverse of the matrix of phylogenetic correlation:
inv.phylo<-inverseA(u.tree,nodes="TIPS",scale=TRUE)

#For some reason, the model would only run if I rename the column heading 'genus' to phylo and 'pruned.tree' to phylo - maybe both of these things need to have the same name?
data7<-data2 %>% 
  rename(
    phylo = genus,
  )
phylo<-pruned.tree

#Create our prior - this is a default prior?
prior<-list(G=list(G1=list(V=1,nu=0.02)),R=list(V=1,nu=0.02))


####IMPORTANT - I don't need to create genus-level averages when running this type of model: MCMCglmm will do that for me####
#Create model
model_simple<-MCMCglmm(CS~MF,random=~phylo,
                       family="gaussian",ginverse=list(phylo=inv.phylo$Ainv),prior=prior,
                       data=data7,nitt=5000000,burnin=1000,thin=500)

summary(model_simple)
plot(model_simple)
heidel.diag(model_simple$Sol)

#Calculate the posterior probability of the phylogenetic signal 'lambda'
lambda <- model_simple$VCV[,'phylo']/
  (model_simple$VCV[,'phylo']+model_simple$VCV[,'units'])
mean(lambda)
posterior.mode(lambda)
HPDinterval(lambda)

#####################
xx <- 0:20
xx <- 0:271
counts <- c(49, 36, 42, 26, 22, 22, 8, 12, 2, 4, 7, 0, 1, 1, 1, 1, 2, 1, 0, 1, 0)
obs <- rep(xx,counts)
obs
hist(obs)

poisson.density <- length(kl)*dpois(xx,mean(kl))
nb <- fitdistr(kl,"negative binomial")
nb.density <- length(kl)*dnbinom(xx,size=nb$estimate["size"],mu=nb$estimate["mu"])

foo <- barplot(counts,names.arg=xx,ylim=range(c(counts,poisson.density)))
lines(foo[,1],poisson.density,lwd=2)
lines(foo[,1],nb.density,lwd=2,col="red")
legend("topright",lwd=2,col=c("black","red"),legend=c("Poisson","Negative Binomial"))

#I want to test if MF is normally distributed - it is very skewed and log transformations do not seem to help. Is this important?
ggdensity((antdata$eff.mating.freq.MEAN.harmonic), 
          main = "Density plot of MF",
          xlab = "MF")
ggqqplot((antdata$eff.mating.freq.MEAN.harmonic))
shapiro.test((antdata$eff.mating.freq.MEAN.harmonic))

l <- filter(antdata, eff.mating.freq.MEAN.harmonic>= 0)
kl<-l$eff.mating.freq.MEAN.harmonic
jkl<-sqrt(kl-min(kl))

fitdistr(kl, "Poisson")

library(fitdistrplus)
library(logspline)

descdist((T_box), discrete = FALSE)
fit.gamma <- fitdist(kl, "Poisson", )
plot(fit.gamma)

plotNormalHistogram(jkl)

Box = boxcox(kl ~ 1,              # Transform Turbidity as a single vector
             lambda = seq(-10,30,1)      # Try values -6 to 6 by 0.1
)
Cox = data.frame(Box$x, Box$y)            # Create a data frame with the results
Cox2 = Cox[with(Cox, order(-Cox$Box.y)),] # Order the new data frame by decreasing y
Cox2[1,]                                  # Display the lambda with the greatest
#    log likelihood
lambda = Cox2[1, "Box.x"]                 # Extract that lambda
T_box = (kl ^ lambda - 1)/lambda   # Transform the original data
plotNormalHistogram(sqrt(T_box))
library(rcompanion)
plotNormalHistogram(T_box)
########################

###Phylogenetic Least Squares Regression PGLS###

#Dataframe with just genus and mating frequency
data <- as.data.frame(
  cbind(
    genus = word(antdata$animal,sep = "_"),
    mating_frequency = antdata$eff.mating.freq.MEAN.harmonic
  )
) 

View(data)

#Mating frequency is currently a factor. Convert to a numeric
data$mating_frequency <- as.numeric(as.character(data$mating_frequency))

#Filter mating frequency data to remove all species without data
data_filtered <- filter(data, mating_frequency>= 0)

View(data_filtered)

#Create genus-level averages for MF
keys <- colnames(data_filtered)[!grepl('mating_frequency', colnames(data_filtered))]
X <- as.data.table(data_filtered)
X[,list(MF= mean(mating_frequency)),keys]
data1<-X[,list(MF= mean(mating_frequency)),keys]
View(data1)

###Now do the same for colony size
data_c <- as.data.frame(
  cbind(
    genus = word(antdata$animal,sep = "_"),
    colony_size = antdata$colony.size
  )
) 

#Colony size is currently a factor. Convert to a numeric
data_c$colony_size <- as.numeric(as.character(data_c$colony_size))

#Filter colony size to remove all species without data
data_filtered_c <- filter(data_c, colony_size>=0)

#Create genus-level averages for colony size
keys_c <- colnames(data_filtered_c)[!grepl('colony_size', colnames(data_filtered_c))]
X_c <- as.data.table(data_filtered_c)
X_c[,list(CS= mean(colony_size)),keys_c]
data1_c<-X_c[,list(CS= mean(colony_size)),keys_c]
View(data1_c)

###Combine and filter averaged genus-level dataframes for CS and MF: data1 and data1_c
matching<-intersect(data1$genus, data1_c$genus)
match_vector<-data1$genus %in% matching
match_vector_filtered<-filter(data1, match_vector)

match_vector2<-data1_c$genus %in% matching
match_vector2_filtered<-filter(data1_c, match_vector2)

#Bind the CS column to my MF dataframe
bound<-cbind(match_vector_filtered,
             CS = match_vector2_filtered$CS)


####Prune tree####

### bound$genus appeared to retain all of the previous 252 levels which had been lost. This may be a solution?
#However, this may not even be a problem at all
#bound1<-bound
#as.factor(as.character(bound1$genus))
#str(bound1$genus)

#Match the tip labels from the tree and the species from the database and prune the phylogeny accordingly
#It is the material that is in the first named set, that is not in the second named set
pruned.tree<-drop.tip(anttree, setdiff(anttree$tip.label, bound$genus))
pruned.tree
plot.phylo(pruned.tree)
plotTree(pruned.tree,ftype="i",fsize=0.4,lwd=1)

#Filter through my dataframe which has genus-level averages and select only the rows that match the tips of my tree
data2<-filter(bound, genus %in% pruned.tree$tip.label)

#Convert the values in a column into row names in an existing data frame
#So that my data is formatted in the same way as the tutorial
#data must be formatted in this way in order to match with the phylogeny tips
data3<-data2 %>% remove_rownames() %>% column_to_rownames(var="genus")
View(data3)

#Adds a 'Species' column which ape seems to want
data4<-cbind(data3,
             Species = data2$genus)
View(data4)
###Ensure that phylogeny tip labels are exactly matched to the row names of the database
a<-name.check(pruned.tree,data4)
a

#Phylogenetic GLS is basically a linear model in which the covariance (correlation) structure between species is permitted 
#to match that expected under a Brownian motion process* of evolution on the tree. (*Or other processes.) Consequently, 
#the first step is to define this covariance structure.
bm<-corBrownian(1, pruned.tree, form = ~Species)
bm

#A model investigating the relationship between altitude at which a species is found and the length of the note in its song.
modelo1<-gls(log(CS)~MF, data=data4, correlation=bm, method = "ML")
modelo1<-gls(log(CS)~MF, data4, correlation=corBrownian(1, pruned.tree, form = ~Species), method = "ML")
summary(modelo1)
plot(modelo1)


#####Using caper instead of ape####
data2<-filter(bound, genus %in% pruned.tree$tip.label)
data20<-data2
data20$genus<-as.factor(as.character(data20$genus))

data20$genus
pruned.tree$tip.label

tip.order <- match(pruned.tree$tip.label, data20$genus)

###In the code below, Gijsbert matches the order of the species in the tree to the order of the species in the database and plots the 
#However, this code does not re-order the database which could cause further problems downstream
comp_auch_plotvec<-as.factor(auchdata$comp.ident[match(pruned.tree$tip.label,table=auchdata$hemipteran.species)])
jkl<-data20[match(pruned.tree$tip.label, data20$genus)]


###To solve problem - does the order of the tip labels and the order of the genus names in the database matter?
comp.data<-comparative.data(pruned.tree, jkl, names.col="genus", vcv.dim=2, warn.dropped=TRUE)
modelo4<-pgls(log(CS)~MF, data=comp.data)
summary(modelo4)


datos<-read.csv(file.choose(),header=TRUE,row.names=1)
View(datos)

datos1<-read.csv(file.choose(),header=TRUE)
View(datos1)
