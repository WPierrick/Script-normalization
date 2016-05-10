# Clustering

# Using Ward's distance between samples and euclidean distance per row

setwd("/clusterdata/uqpwains/ibscratch/QC-IDATS")

load("mValsSw.RObject")
revmval  <- t(mValsSw)
cluster1 = hclust(dist(revmval, method = "euclidean"), method = "ward.D2")
png("cluster.png")
plot(cluster1,"Clustering")
dev.off()
