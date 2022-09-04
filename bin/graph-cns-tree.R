#!/usr/bin/env Rscript

library("RColorBrewer")
library("codetools")
library("phytools")
packageVersion("phytools")

#import command line arguments
args <- commandArgs(TRUE)
tree_fp = as.character(args[1])
cns_table_fp = as.character(args[2])
outgroup_name = as.character(args[3])
main_gene_name = as.character(args[4])
colorful = as.character(args[5])

#read in tree and cns data
tree <- read.tree(tree_fp)
x <- read.csv(cns_table_fp)

#root the tree on the outgroup
tree <- root(tree, outgroup_name)

#process dataframe for phytools
x_matrix <- as.matrix(x)
rownames(x_matrix) <- x$tree_genes
x_matrix<-x_matrix[,-2] # delete column 2
x_matrix<-x_matrix[,-1]
x <- as.data.frame(x_matrix)
par(xpd = NA)

#make grey and white graph
if (colorful=="false") {
save_name = paste0(main_gene_name,".cns.tree.pdf")
#convert from true false dataframe to 1s and 0s, this is for color mapping purposes
x[] <- lapply( x, factor)
x[] <- lapply( x, as.logical)
x[] <- lapply( x, as.integer)
pdf(file = save_name)
colors<-setNames(replicate(ncol(x),setNames(c("white","darkgrey"),0:1),simplify=FALSE),colnames(x))
plotTree.datamatrix(tree,x,colors=colors,length=4,ftype="i",fsize=.35,sep=0,xlim=10,ylim=-10)
nodelabels(tree$node.label, node=2:tree$Nnode+Ntip(tree), adj=c(1,-0.2),frame="none",cex=.5)
add.simmap.legend(colors=setNames(c("white","darkgrey"),c("absent","present")),shape="square",prompt=FALSE,x=0,y=3,fsize=0.9)
dev.off()
} else if (colorful=="true") {
save_name = paste0(main_gene_name,".colorful.cns.tree.pdf")
#convert to factor
x[] <- lapply( x, factor)
pdf(file = save_name)
plotTree.datamatrix(tree,x,length=4,ftype="i",fsize=.35,sep=0,xlim=10,ylim=-10)
nodelabels(tree$node.label, node=2:tree$Nnode+Ntip(tree), adj=c(1,-0.2),frame="none",cex=.5)
dev.off()
} else {
print("Invalid colorful parameter. Options are false or true (lowercase).")
}

