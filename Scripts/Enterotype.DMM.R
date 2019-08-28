# the DMM ======================
# author: shixian hu

library(DirichletMultinomial)
library(reshape2)
library(magrittr)
library(dplyr)
library(vegan)
library(ggplot2)

uncorrect=read.table("../../DAG3_metaphlan_merged.txt",sep = "\t",check.names = F,
                     stringsAsFactors = F,header = T,row.names = 1)
uncorrect=as.data.frame(t(uncorrect))
species=uncorrect[,grep("s__",colnames(uncorrect))]
#species=species[,grep("s__",colnames(species),invert = T)]
colnames(species)=lapply(colnames(species),function(x){
  strsplit(x,"s__")[[1]][2]
})
species = species[rowSums(species > 0) > 0,]
beta_diversity=vegdist(species,method = "bray")
pcoa_analysis=as.data.frame(cmdscale(beta_diversity,k=4))
species=species*1000000
write.table(species,"Count.txt",quote=F,row.names=T,sep="\t")

fit <- mclapply(1:6, dmn, count = as.matrix(species), verbose=TRUE,mc.cores=1)

# lplc
lplc <- sapply(fit, laplace)
print(lplc)

pdf('lplc.model.pdf')
plot(lplc, type="b", xlab="Number of Dirichlet Components", ylab="Model Fit")
dev.off()

best <- fit[[which.min(lplc)]]
cluster <- data.frame(Cluster=apply(mixture(best), 1, which.max))

table=merge(cluster,pcoa_analysis,by="row.names",all = F)
write.table(table,"lplc.cluster.txt",quote=F,row.names=T,sep="\t")
ggplot (table, aes(V1,V2)) + geom_point(aes(colour = as.factor(Cluster)),size=1) + theme_bw()
ggsave("lplc.PCoA.pdf")

for (k in seq(ncol(fitted(best)))) {
  d <- melt(fitted(best))
  colnames(d) <- c("Species", "cluster", "value")
  d <- subset(d, cluster == k) %>%
    arrange(value) %>%
    mutate(Species = factor(Species, levels = unique(Species))) %>%
    filter(abs(value) > quantile(abs(value), 0.8))
  
  p <- ggplot(d, aes(x = Species, y = value)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(title = paste("Top drivers: community type", k))
  print(p)
  ggsave(paste(k,"lplc.core.species.pdf",sep="_"))
}

heatmapdmn(as.matrix(species), fit[[1]], best, 10)
ggsave("lplc.Heatmap.pdf")

# aic
aic  <- sapply(fit, AIC)
print(aic)

pdf('aic.model.pdf')
plot(aic, type="b", xlab="Number of Dirichlet Components", ylab="Model Fit")
dev.off()


best <- fit[[which.min(aic)]]
cluster <- data.frame(Cluster=apply(mixture(best), 1, which.max))

table=merge(cluster,pcoa_analysis,by="row.names",all = F)
write.table(table,"aic.cluster.txt",quote=F,row.names=T,sep="\t")
ggplot (table, aes(V1,V2)) + geom_point(aes(colour = as.factor(Cluster)),size=1) + theme_bw()
ggsave("aic.PCoA.pdf")

for (k in seq(ncol(fitted(best)))) {
  d <- melt(fitted(best))
  colnames(d) <- c("Species", "cluster", "value")
  d <- subset(d, cluster == k) %>%
    arrange(value) %>%
    mutate(Species = factor(Species, levels = unique(Species))) %>%
    filter(abs(value) > quantile(abs(value), 0.8))
  
  p <- ggplot(d, aes(x = Species, y = value)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(title = paste("Top drivers: community type", k))
  print(p)
  ggsave(paste(k,"aic.core.species.pdf",sep="_"))
}

heatmapdmn(as.matrix(species), fit[[1]], best, 10)
ggsave("aic.Heatmap.pdf")


# bic
bic  <- sapply(fit, BIC)
print(bic)

pdf('bic.model.pdf')
plot(bic, type="b", xlab="Number of Dirichlet Components", ylab="Model Fit")
dev.off()

best <- fit[[which.min(bic)]]
cluster <- data.frame(Cluster=apply(mixture(best), 1, which.max))

table=merge(cluster,pcoa_analysis,by="row.names",all = F)
write.table(table,"bic.cluster.txt",quote=F,row.names=T,sep="\t")
ggplot (table, aes(V1,V2)) + geom_point(aes(colour = as.factor(Cluster)),size=1) + theme_bw()
ggsave("bic.PCoA.pdf")

for (k in seq(ncol(fitted(best)))) {
  d <- melt(fitted(best))
  colnames(d) <- c("Species", "cluster", "value")
  d <- subset(d, cluster == k) %>%
    arrange(value) %>%
    mutate(Species = factor(Species, levels = unique(Species))) %>%
    filter(abs(value) > quantile(abs(value), 0.8))
  
  p <- ggplot(d, aes(x = Species, y = value)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(title = paste("Top drivers: community type", k))
  print(p)
  ggsave(paste(k,"bic.core.species.pdf",sep="_"))
}

heatmapdmn(as.matrix(species), fit[[1]], best, 10)
ggsave("bic.Heatmap.pdf")
