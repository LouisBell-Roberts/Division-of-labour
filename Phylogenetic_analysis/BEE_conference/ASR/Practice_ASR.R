#Ancestral State Reconstructions for BEE presentation using Phytools
##02/03/2021 Louis Bell-Roberts

library(ape)
library(phytools)
library(geiger)
library(phylolm)
library(corHMM)
library(caper)

anttree_genus <- read.tree(file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Trees/Genus_tree/Genus_Level_ML_TREE_treepl_185_outsdroppeds.ladderized.increasing.tre")
anttree_species <- read.tree(file = "/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Trees/Nelsen/ML_TREE_treepl_185_outsdroppeds.ladderized.increasing.tre")

data <- read.csv("/Users/louis.bell-roberts/Documents/DTP_1st_project_rotation/Data/Shared database_Juliet/Shared_Juliet.csv", header = T)

#Tree is rooted but not ultrametric
is.rooted(anttree_species)
is.ultrametric(anttree_species)

#Select only ant species
antdata<-filter(data, type=="ant")
View(antdata)

#Make the data binary
is.numeric(antdata$eff.mating.freq.MEAN.harmonic)
antdata$Caste1 <- as.numeric((as.character(antdata$Caste1)))

antdata1.1 <- antdata %>% mutate(Harmonic_cat=cut(eff.mating.freq.MEAN.harmonic, breaks=c(-Inf, 2, Inf), labels=c("1","2")))
antdata1.2 <- antdata1.1 %>% mutate(Caste_cat=cut(Caste1, breaks=c(-Inf, 1.9, Inf), labels=c("1","2")))

antdata2 <- antdata1.2 %>% select(animal, Harmonic_cat, Caste_cat)
View(antdata2)

#Set caste and MF as a numeric
antdata2$Caste_cat <- as.numeric(as.character(antdata2$Caste_cat))
antdata2$Harmonic_cat <- as.numeric(as.character(antdata2$Harmonic_cat))

#Remove missing data
antdata3 <- filter(antdata2, Caste_cat >= 0)
View(antdata3)

#Prune database and phylogeny

#Match the tip labels from the tree and the species from the database and prune the phylogeny accordingly
#It is the material that is in the first named set, that is not in the second named set
pruned.tree_sp<-drop.tip(anttree_species, setdiff(anttree_species$tip.label, antdata3$animal))
pruned.tree_sp
plotTree(pruned.tree_sp,ftype="i",fsize=0.4,lwd=1)

#Filter through my dataframe and select only the rows that match the tips of my tree
antdata4<-filter(antdata3, animal %in% pruned.tree_sp$tip.label)
View(antdata4)

antdata4.1 <- antdata4 %>% select(animal, Caste_cat)

######This one doesn't work#######
#ASR:
ancestral<-ace(phy=pruned.tree_sp,x = antdata4.1$Caste_cat,type="d")

#However, this code does not re-order the database which could cause further problems downstream
plotvec.1<-as.factor(antdata4.1$Caste_cat[match(pruned.tree_sp$tip.label,table=antdata4.1$animal)])
plotvec<-as.factor(antdata4$Caste_cat[match(pruned.tree_sp$tip.label,table=antdata4$animal)])

#Visualisation
plot.phylo(x = pruned.tree_sp,show.tip.label = T, cex = 0.4, tip.color = c("red","blue")[plotvec.1])
nodelabels(thermo = ancestral$lik.anc,piecol = c("red","blue"),cex=0.1)
####################################

plotvec<-as.factor(antdata4$Caste_cat[match(pruned.tree_sp$tip.label,table=antdata4$animal)])
##corHMM analysis with ASR ###
#To run analyses with corHMM, we need to give a two dimensional table with species names 
ant_Caste_hrm <- antdata4 %>% select(animal, Caste_cat)

#Let's try a hidden rate model with 1 and 2 potential speeds of evolution (phylogeny is too small for this 
#to make much sense)
ancestral_hrm1<-corHMM(phy = pruned.tree_sp, data = ant_Caste_hrm,rate.cat = 1)
ancestral_hrm2<-corHMM(phy = pruned.tree_sp, data = ant_Caste_hrm,rate.cat = 2)

#Let's consider the AICc-values  - lower value for model 2: model 2 is better
ancestral_hrm1$AICc
ancestral_hrm2$AICc

plotRECON(phy=ancestral_hrm2$phy,likelihoods = ancestral_hrm2$states,pie.cex=0.5,
          piecolors = c("red","blue"), tip.color = c("red","blue")[plotvec])


####Correlated evolution among two binary traits
#Let's explore correlated/uncorrelated models of evolution among two binary traits
antdata5 <- filter(antdata2, Caste_cat >= 0, Harmonic_cat >=0)

#Prune tree
#Prune database and phylogeny

#Match the tip labels from the tree and the species from the database and prune the phylogeny accordingly
#It is the material that is in the first named set, that is not in the second named set
pruned.tree_sp2<-drop.tip(anttree_species, setdiff(anttree_species$tip.label, antdata5$animal))
pruned.tree_sp2
plotTree(pruned.tree_sp2,ftype="i",fsize=0.4,lwd=1)

#Filter through my dataframe and select only the rows that match the tips of my tree
antdata6<-filter(antdata5, animal %in% pruned.tree_sp$tip.label)
View(antdata6)

#Set dataframe as a matrix
antdata7<-antdata6
as.matrix(antdata7)

#Needs column names removed
antdata8 <-antdata7 %>% remove_rownames() %>% column_to_rownames(var="animal")

#Black is high, Caste is on the right
dotTree(pruned.tree_sp2,x = antdata8,data.type="discrete")

T1<-antdata8$Caste_cat
names(T1)<-pruned.tree_sp2$tip.label
T2<-antdata8$Harmonic_cat
names(T2)<-pruned.tree_sp2$tip.label

ant_correlated_T1T2<-fitPagel(pruned.tree_sp2, x = T1, y = T2)
ant_correlated_T1T2
plot.fitPagel(ant_correlated_T1T2)
#Lacks support for a dependent model of evolution.
##However, transition rates for dependent model are very interesting
########

#Alternative ASR - will I get the same results: appears to match corHMM hidden rates model quite closely. Ace just doesn't seem to be working in comparison
antdata4

#Run the simulation
ant_dat_simul_ER<-rayDISC(phy = pruned.tree_sp,data = antdata4,charnum = 2,ntraits = 1,
                               model="ER",node.states = "marginal",root.p = "maddfitz")
ant_dat_simul_SYM<-rayDISC(phy = pruned.tree_sp,data = antdata4,charnum = 2,ntraits = 1,
                                model="SYM",node.states = "marginal",root.p = "maddfitz")
ant_dat_simul_ARD<-rayDISC(phy = pruned.tree_sp,data = antdata4,charnum = 2,ntraits = 1,
                                model="ARD",node.states = "marginal",root.p = "maddfitz")

ant_dat_simul_ER$AICc
ant_dat_simul_SYM$AICc
ant_dat_simul_ARD$AICc

#For colouring the figure correctly
plotvec<-as.factor(antdata4$Caste_cat[match(pruned.tree_sp$tip.label,table=antdata4$animal)])

plotRECON(phy=ant_dat_simul_ER$phy,likelihoods = ant_dat_simul_ER$states,pie.cex=0.5,
          piecolors = c("red","blue"), tip.color = c("red","blue")[plotvec])
ant_dat_simul_ER$solution



###########Genus level ASR and Correlated evolution models - Caste and MF##### ---------- NOT WORKING AT THE MOMENT
#Select only ant species
antdata<-filter(data, type=="ant")
View(antdata)

#Set caste and MF as a numeric
antdata$Caste1 <- as.numeric(as.character(antdata$Caste1))
antdata$eff.mating.freq.MEAN.harmonic <- as.numeric(as.character(antdata$eff.mating.freq.MEAN.harmonic))

#Select only the species, MF and Caste columns
antdata1 <- antdata %>% select(animal, eff.mating.freq.MEAN.harmonic, Caste1)
View(antdata1)

#Dataframe with genus added as additional column
data_with_genus <- as.data.frame(
  cbind(
    genus = word(antdata1$animal,sep = "_"),
    antdata1
  )
) 

#Create genus-level averages
l<-sapply(split(data_with_genus$eff.mating.freq.MEAN.harmonic,data_with_genus$genus), mean, na.rm=T)
k<-sapply(split(data_with_genus$Caste1,data_with_genus$genus),mean, na.rm=TRUE)
j<-cbind(l,k)
h<-as.data.frame(j)
g<-rownames_to_column(h, var = "phylo")
antdata2<-g %>% 
  rename(
    MF = l,
    Caste = k
  )
View(antdata2)
antdata2$Caste <- na_if(antdata2$Caste, "NaN")
antdata2$MF <- na_if(antdata2$MF, "NaN")

#Remove all of the missing data for Caste
antdata2.1<-filter(antdata2, Caste>=0)
View(antdata2.1)

#Make data binary - Anything <= 1.5 is turned into 1
antdata2.2 <- antdata2.1 %>% mutate(Caste_cat=cut(Caste, breaks=c(-Inf, 1.5, Inf), labels=c("1","2")))

#To run analyses with corHMM, we need to give a two dimensional table with species names 
genus_Caste_hrm <- antdata2.2 %>% select(phylo, Caste_cat)

#Prune database and phylogeny using genus-level tree
#Match the tip labels from the tree and the species from the database and prune the phylogeny accordingly
#It is the material that is in the first named set, that is not in the second named set
pruned.tree_gs<-drop.tip(anttree_genus, setdiff(anttree_genus$tip.label, genus_Caste_hrm$phylo))
pruned.tree_gs
plotTree(pruned.tree_gs,ftype="i",fsize=0.4,lwd=1)

#Filter through my dataframe and select only the rows that match the tips of my tree
genus_Caste_hrm_pruned<-filter(genus_Caste_hrm, phylo %in% pruned.tree_gs$tip.label)

#Let's try a hidden rate model with 1 and 2 potential speeds of evolution (phylogeny is too small for this 
#to make much sense)
ancestral_hrm1<-corHMM(phy = pruned.tree_gs, data = genus_Caste_hrm_pruned,rate.cat = 1)
ancestral_hrm2<-corHMM(phy = pruned.tree_gs, data = genus_Caste_hrm_pruned,rate.cat = 2)

#Let's consider the AICc-values  - lower value for model 1: model 1 is better
ancestral_hrm1$AICc
ancestral_hrm2$AICc

plotvec2<-as.factor(genus_Caste_hrm_pruned$Caste_cat[match(pruned.tree_gs$tip.label,table=genus_Caste_hrm_pruned$phylo)])

plotRECON(phy=ancestral_hrm2$phy,likelihoods = ancestral_hrm2$states,pie.cex=0.5,
          piecolors = c("red","blue"), tip.color = c("red","blue")[plotvec2])






