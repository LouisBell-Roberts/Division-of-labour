#Create polytomies function

#Simultaneously check if a tree is rooted and has no polytomies
is.binary()

## add a single species to a genus at the root of the genus
species<-c("Pheidole_dentata","Pheidole_sp")
t1<-add.species.to.genus(pruned.tree_sp,species)
plotTree(t1,ftype="i")


## add a set of species in a vector
species<-c("Pheidole_dentata","Pheidole_sp")
for(i in 1:length(species)) pruned.tree_sp<-add.species.to.genus(pruned.tree_sp,species[i],where="root")
plotTree(pruned.tree_sp,ftype="i",fsize=0.4,lwd=1)