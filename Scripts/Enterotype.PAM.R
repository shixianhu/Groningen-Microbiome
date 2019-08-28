library(ggplot2)
library(cluster)
library(factoextra)
library(ade4)
library(clusterSim)

#============================================================================
dist.JSD <- function(inMatrix, pseudocount=0.000001, ...) {
  KLD <- function(x,y) sum(x *log(x/y))
  JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
  matrixColSize <- length(colnames(inMatrix))
  matrixRowSize <- length(rownames(inMatrix))
  colnames <- colnames(inMatrix)
  resultsMatrix <- matrix(0, matrixColSize, matrixColSize)

  inMatrix = apply(inMatrix,1:2,function(x) ifelse (x==0,pseudocount,x))

  for(i in 1:matrixColSize) {
    for(j in 1:matrixColSize) {
      resultsMatrix[i,j]=JSD(as.vector(inMatrix[,i]),
                             as.vector(inMatrix[,j]))
    }
  }
  colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
  as.dist(resultsMatrix)->resultsMatrix
  attr(resultsMatrix, "method") <- "dist"
  return(resultsMatrix)
 }

pam.clustering=function(x,k) { # x is a distance matrix and k the number of clusters
  cluster = as.vector(pam(as.dist(x), k, diss=TRUE)$clustering)
}

# import metaphlan data of dag3, filtered, re-scaled (recalculate relative abundance)

metaphlan=read.table("../../DAG3_metaphlan_merged.txt",header = T,
                     stringsAsFactors = F,sep = "\t",check.names = F,row.names = 1)
taxa=as.data.frame(t(metaphlan))
genus=taxa[,grep("g__",colnames(taxa))]
genus=genus[,grep("s__",colnames(genus),invert = T)]
colnames(genus)=lapply(colnames(genus),function(x){
  strsplit(x,"g__")[[1]][2]
})
genus=as.data.frame(t(genus))
genus.dist=dist.JSD(genus)
genus.cluster=pam.clustering(genus.dist, k=2)

nclusters = index.G1(t(genus), genus.cluster, d = genus.dist, centrotypes = "medoids")
nclusters=NULL

for (k in 1:10) {
  if (k==1) {
    nclusters[k]=NA
  } else {
    genus.cluster_temp=pam.clustering(genus.dist, k)
    nclusters[k]=index.G1(t(genus),genus.cluster_temp,  d = genus.dist,
                          centrotypes = "medoids")
  }
}

pdf("CHindex.pdf")
plot(nclusters, type="h", xlab="k clusters", ylab="CH index")
dev.off()

cluster=data.frame(row.names = colnames(genus),Cluster=genus.cluster)
write.table(cluster,file = "Clusters.txt",sep = "\t",row.names = T,quote = F)

pdf("PCoA.pdf")
beta_diversity=vegdist(as.data.frame(t(genus)),method = "bray")
pcoa_analysis=as.data.frame(cmdscale(beta_diversity,k=4))
table=merge(cluster,pcoa_analysis,by="row.names",all = F)
ggplot (table, aes(V1,V2)) + geom_point(aes(colour = as.factor(Cluster)),size=1) + theme_bw()
dev.off()

