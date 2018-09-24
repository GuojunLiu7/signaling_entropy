### DoIntegPIN.R
### Author: Andrew Teschendorff
### Date: 31 Dec 2013.

### This function takes as input
### (i) adj.m: an adjacency matrix representing a PPI network, with rownames/colnames of the matrix annotated to Entrez gene IDs.
### (ii) exp.m: an expression data matrix with rows labeling the genes (also annotated to entrez gene IDs.

### The output of the function is a list with following entries:
### a: the maximally connected network upon integrating the input PPI with the gene expression data matrix.
### e: the corresponding gene expression data matrix (same number of rows as $a).
### gr.o: the igraph R object corresponding to $a.

DoIntegPIN <- function(adj.m,exp.m){

require(igraph);

### find genes in network with gene expression profiles
commonEID.v <- intersect(rownames(adj.m),rownames(exp.m));
match(commonEID.v,rownames(exp.m)) -> map.idx;
expPIN.m <- exp.m[map.idx,];
match(commonEID.v,rownames(adj.m)) -> map1.idx;
adjEXP.m <- adj.m[map1.idx,map1.idx];

gr.o <- graph.adjacency(adjEXP.m,mode="undirected");
comp.l <- clusters(gr.o)
cd.v <- summary(factor(comp.l$member));
mcID <- as.numeric(names(cd.v)[which.max(cd.v)]);
maxc.idx <- which(comp.l$member==mcID);
adjMC.m <- adjEXP.m[maxc.idx,maxc.idx];
expMC.m <- expPIN.m[maxc.idx,];
grMC.o <- graph.adjacency(adjMC.m,mode="undirected");

return(list(a=adjMC.m,e=expMC.m,gr=grMC.o));

}


