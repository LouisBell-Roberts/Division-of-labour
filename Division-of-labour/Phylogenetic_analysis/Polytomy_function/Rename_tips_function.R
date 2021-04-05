library(ape)
library(phytools)


tree1 <- stree(4)
plot(a)
plotTree(tree1,node.numbers=T)
a$tip.label

b <- stree(3)
plotTree(b, node.numbers=T)
a$tip.label

#where=which(a$tip.label=="t3")
which(node==a$node.label)


tree2<-bind.tree(a,"t21",where="t3",position=tree$edge.length[which(a$edge[,2]==6)])
#Binding to tree
tree1<-pbtree(n=5)
plot(tree1,no.margin=T,label.offset=0.1)
nodelabels()
tiplabels()

#the bit in quotation is the name of the new branch
##where = the branch it comes off of
tree2<-bind.tip(tree1,"t19",where=1,position= 0.5*tree1$edge.length[which(tree1$edge[,2]==1)])
plot(tree2,no.margin=T,label.offset=0.1)
nodelabels()
tiplabels()


tree3<-pbtree(n=20)
plot(tree3,no.margin=T,label.offset=0.1)
nodelabels()
tiplabels()

tree4<-bind.tip(tree3,"t21",where=23,position= 0.5*tree3$edge.length[which(tree3$edge[,2]==23)])
plot(tree4,no.margin=T,label.offset=0.1)
nodelabels()
tiplabels()


tree<-pbtree(n=20)
plot(tree,no.margin=T,label.offset=0.1) # offset may vary
nodelabels()
tiplabels()


which(node==tree$node.label)
which(tree$tip.label=="t20")



#Function for binding tips
bind.tip<-function(tree,tip.label,edge.length=NULL,where=NULL){
  if(is.null(where)) where<-length(tree$tip)+1
  tip<-list(edge=matrix(c(2,1),1,2),
            tip.label=tip.label,
            edge.length=edge.length,
            Nnode=1)
  class(tip)<-"phylo"
  obj<-bind.tree(tree,tip,where=where)
  return(obj)
}






############################
#Function for renaming tips
a$tip.label[a$tip.label=="t4"] <- "something_else"

#Function renames tip labels
rename.tips <- function(phy, old_names, new_names) {
  mpos <- match(old_names,phy$tip.label)
  phy$tip.label[mpos] <- new_names
  return(phy)
}

b<-c("q","w","e","r")

rename.tips(a,a$tip.label,b)

old_names<-a$tip.label
match(old_names,a$tip.label)
