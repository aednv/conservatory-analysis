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

save_name = paste0(main_gene_name,".cns.tree.pdf")

tree <- read.tree(tree_fp)
x <- read.csv(cns_table_fp)
tree <- root(tree, outgroup_name)
x_matrix <- as.matrix(x)
rownames(x_matrix) <- x$tree_genes
x_matrix<-x_matrix[,-2] # delete column 2
x_matrix<-x_matrix[,-1]
x <- as.data.frame(x_matrix)
#convert from true false to 1s and 0s
x[] <- lapply( x, factor)
x[] <- lapply( x, as.logical)
x[] <- lapply( x, as.integer)
par(xpd = NA)

pdf(file = save_name)
colors<-setNames(replicate(ncol(x),setNames(c("white","darkgrey"),0:1),simplify=FALSE),colnames(x))
plotTree.datamatrix(tree,x,colors=colors,length=4,ftype="i",fsize=.35,sep=0,xlim=10,ylim=-10)
nodelabels(tree$node.label, node=2:tree$Nnode+Ntip(tree), adj=c(1,-0.2),frame="none",cex=.5)
add.simmap.legend(colors=setNames(c("white","darkgrey"),c("absent","present")),shape="square",prompt=FALSE,x=0,y=3,fsize=0.9)
dev.off()
