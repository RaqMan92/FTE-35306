

#############
# Libraries #

# Used for the plotting of the PCA.
library(plotrix)

library(randomForest)
library(tree)

#############
# Functions #

#' Cols
#' Creates a range of colors that can be used for plotting.
#' @param vec: is a vector that will be given colors ID's.
#'
#' @return Color ID's for the use in plotting
Cols <- function(vec) {
  cols <- rainbow(length(unique(vec)))
  
  return(cols[as.numeric(as.factor(vec))])
  
}


#############
# Main code #

# 1. Apply PCA to the expression data from the brain subregions. Perform PCA
# both with and without scaling to unit standard deviation. Address the
# following questions:

# Perform a PCA on the data, x must be numeric so remove the tissues from the PCA.
pr.out.scale <- prcomp(mydat[-152], scale = TRUE)
pr.out <- prcomp(mydat[-152], scale = FALSE)

# a. Analyze how much of the variation is explained by the first few PCs.

# Calculate the variance explained by each principal component, do this for the scale and non-scale.
pr.var.scale <- pr.out.scale$sdev ^ 2
pr.var <- pr.out$sdev ^ 2

# Calculate the proportion of the variance that is being explained by each principal component.
pve.scale <- pr.var.scale / sum(pr.var.scale)
pve <- pr.var / sum(pr.var)

# Plot the PVE for each component and the cumulative PVE as well.
par(mfrow = c(2,2))
# Plot the scaled data.
plot(
  pve.scale,
  xlab = "Principal Component",
  ylab = "Proportion of Variance Explained",
  main = "Proportion of Variance Explained Scaled PCA",
  ylim = c(0, 1),
  type = 'b'
)
plot(
  cumsum(pve.scale),
  xlab = "Principal Component",
  ylab = "Cumulative Proportion of Variance Explained",
  main = "Cumulative Proportion of Variance Explained Scaled PCA",
  ylim = c(0,1) ,
  type = 'b'
)
# Plot the non-scaled data.
plot(
  pve,
  xlab = "Principal Component",
  ylab = "Proportion of Variance Explained",
  main = "Proportion of Variance Explained Non-scaled PCA",
  ylim = c(0, 1),
  type = 'b'
)
plot(
  cumsum(pve),
  xlab = "Principal Component",
  ylab = "Cumulative Proportion of Variance Explained",
  main = "Cumulative Proportion of Variance Explained Non-scaled PCA",
  ylim = c(0,1) ,
  type = 'b'
)

# There is quite a difference between the scale and the non-scale, the first PC
# in the non-scale explains 81.18%. While the scale only explains 27.68%.

# b. How are the different tissues placed in the space defined by
# these PCs? Comment on how this relates to observations when you classified
# these tissues compared to each other (week 1,2).

# See the number of unique tissues
unique(expr4T.filtered$tissue)

# Each entry corrosponds to the same index as the tissue, so first tissue has the first color.
sapply(Cols(unique(expr4T.filtered$tissue)), color.id)

#plotting samples after transformation:

plot(pr.out.scale$x[,1:2],col = Cols(expr4T.filtered$tissue),pch = 19)

pr.out.scale$rotation

cbind(pr.out.scale$x[,1:2], expr4T.filtered$tissue)


plot(pr.out$x[,1:2],col = Cols(expr4T.filtered$tissue),pch = 19)

pr.out$rotation

cbind(pr.out$x[,1:2], expr4T.filtered$tissue)

# Which genes contribute substantially to first or second PC:

pr.out.scale$rotation[which(abs(pr.out.scale$rotation[,1]) > 0.2 |
                              abs(pr.out.scale$rotation[,2]) > 0.2),1:2]

# Which genes contribute substantially to first or second PC:

pr.out$rotation[which(abs(pr.out$rotation[,1]) > 0.2 |
                        abs(pr.out$rotation[,2]) > 0.2),1:2]


# @ TODO: old code remove it.
# # Set a fixed seed so that the results are reproducable
# set.seed(42)
# random.sample.data <- sample(1:nrow(expr4T.filtered), round(nrow(expr4T.filtered) * 0.15))
# random.column.data <- sample(1:ncol(expr4T.filtered) - 1, round(ncol(expr4T.filtered) * 0.15))
# # Decrease the number of rows
# random.data.expr4T.filtered <- expr4T.filtered[random.sample.data,random.column.data]
# pr.out.scale.expr4T.filtered <- prcomp(random.data.expr4T, scale = TRUE)
# pr.out.expr4T.filtered <- prcomp(random.data.expr4T, scale = FALSE)
#
#
# library(BiplotGUI)
# Biplots(random.data.expr4T.filtered)
#
#
# biplot(prcomp(mydat[which(mydat$tissue == tissue1),-152], scale = T),
#        scale = 0)
# biplot(prcomp(t(mydat[which(mydat$tissue == tissue1),-152]), scale = T),
#        scale = 0)
#
# par(mfrow=c(1,1))
# biplot(pr.out.scale, scale =0)
# biplot(pr.out, scale =0)

# c. Use principal component loading vectors to identify genes which might be relevant for separation of
# different tissues. Compare these with your observation on informative features in week 2.

# @ TODO: Compare with week 2

# Look at which genes contribute substantially to each principal component loading vectors.
pr.out$rotation[which(abs(pr.out$rotation[,1]) > 0.2),]
pr.out.scale$rotation[which(abs(pr.out.scale$rotation[,1]) > 0.2),]

# 2. Apply hierarchical clustering, testing two different distance methods and
# two different linkages. Also apply K-means clustering. Finally, use the PCA
# results from step 1 (above), pick the first few principal components, and
# apply clustering to the principal component score vectors. Address the
# following questions:

# a. Compare the different clustering results. Do different methods give very
# different or very similar results?


# Check if the kmeans is going well.
clustering.data <- expr4T.filtered[-152]
ncenters <- 13

km.out.scale.13 <-
  kmeans(scale(clustering.data), centers = ncenters)
km.out.13 <- kmeans(clustering.data, centers = ncenters)
km.out.nstart.50 <-
  kmeans(clustering.data, centers = ncenters, nstart = 50)
km.out.nstart.50.scale <-
  kmeans(scale(clustering.data), centers = ncenters, nstart = 50)
km.out.13.pca <- kmeans(pr.out$x[,1:2], centers = ncenters)
km.out.nstart.50.pca <-
  kmeans(pr.out$x[,1:2], centers = ncenters, nstart = 20)

table(km.out.scale.13$cluster)
table(km.out.13$cluster)
# They give very different results for each cluster.

table(km.out.nstart.50$cluster)
table(km.out.nstart.50.scale$cluster)
# They also give very different results for each cluster.

table(km.out.13.pca$cluster)
table(km.out.nstart.50.pca$cluster)
# They also give very different results for each cluster.

hc.complete.pca <- hclust(dist(pr.out$x[,1:2]), method = "complete")
hc.average.pca <- hclust(dist(pr.out$x[,1:2]), method = "average")
hc.single.pca <- hclust(dist(pr.out$x[,1:2]), method = "single")


# Do a hierarchical clustering on the data
hc.complete <- hclust(dist(clustering.data), method = "complete")
hc.average <- hclust(dist(clustering.data), method = "average")
hc.single <- hclust(dist(clustering.data), method = "single")


cutree(hc.complete , ncenters)
cutree(hc.average , ncenters)
cutree(hc.single , ncenters)

# b. Compare the clustering results with the results obtained when classifying
# tissues in weeks 1 and 2.


# Hint: You can use the table() function in R to compare the true
# class labels to the class labels obtained by clustering. Be careful
# how you interpret the results: K-means clustering will arbitrarily
# number the clusters, so you cannot simply check whether the true
# class labels and clustering labels are the same.


tissues.clustering.assigned <-
  unique(expr4T.filtered$tissue)[as.vector(km.out.13$cluster)]

dat <- data.frame(expr4T.filtered, K = tissues.clustering.assigned)
table(dat$K)
table(dat$tissue)
plot(
  as.vector(table(dat$tissue)) - as.vector(table(dat$K)),
  main = "Differences between the class labels",
  ylab = "Difference true class labels and clustering labels",
  xlab = "Class label",
  type = "b"
)

# @ TODO: compare the results with the week 1 and 2.

# c. Can you say something on how many clusters there are in this dataset?

#http://stackoverflow.com/questions/15376075/cluster-analysis-in-r-determine-the-optimal-number-of-clusters
wss <-
  (nrow(clustering.data) - 1) * sum(apply(clustering.data, 2, var))
for (i in 1:ncenters) {
  wss[i] <- sum(kmeans(clustering.data,
                       centers = i)$withinss)
}
plot(
  1:ncenters,
  wss,
  type = "b",
  main = "Compare SSE for a number of clusters",
  xlab = "Number of Clusters",
  ylab = "Within groups sum of squares"
)

# There seems to be clusters in the data, this can be seen in the decrease of the within groups sum of squares.
# The number of clusters seem to be 6, since there is no longer a sharp decrease of SSE. AFter ti the reduction slows down quite a bit.

# 3. Use one of the approaches from step 1 or 2 above to pre-process the gene
# expression data. After this preprocessing, re-train at least two of the
# classification approaches applied in week 1 or 2. Compare the prediction
# performance with what was obtained previously.

# Use clustering to find a intreseting group of tissues, you can train a model of that and use test data to validate it.


interesting.genes <- c()
# Loop trough all the centers in the tree.
for (i in 1:ncenters) {
  # If there are just 2 group members within a single cut, than save it.
  if (length(which(cutree(hc.average , ncenters) == i)) <= 2) {
    # Take the name of the  interesting gene.
    interesting.genes <-
      c(interesting.genes, names(which(cutree(
        hc.average , ncenters
      ) == i)))
  }
}

p <- ncol(expr4T.filtered) - 1
p.2 <- p / 2
p.3 <- sqrt(p)


rf.tissues.clustering.preprocess <- randomForest(
  tissue ~ ENSG00000225972.1_MTND1P23 + ENSG00000225630.1_MTND2P28 + ENSG00000237973.1_hsa.mir.6723 + ENSG00000229344.1_RP5.857K21.7 + ENSG00000154146.8_NRGN + ENSG00000131095.7_GFAP + ENSG00000266844.1_RP11.862L9.3 + ENSG00000197971.10_MBP + ENSG00000123560.9_PLP1 + ENSG00000226958.1_CTD.2328D6.1 + ENSG00000198695.2_MT.ND6 + ENSG00000101405.3_OXT + ENSG00000101200.5_AVP,
  data = expr4T.filtered,
  subset = train.expr4T.data,
  mtry = (length(interesting.genes) - 1) / 2,
  ntree = 500,
  importance = TRUE
)
# OOB estimate of  error rate: 43.17%

rf.clus.pred <-
  predict (rf.tissues.clustering.preprocess, newdata = expr4T.filtered[-train.expr4T.data,])
1 - mean(rf.clus.pred == expr4T.filtered[-train.expr4T.data,]$tissue)


rf.tissues.all.genes <- randomForest(
  tissue ~ .,
  data = expr4T.filtered,
  subset = train.expr4T.data,
  mtry = p.2,
  ntree = 500,
  importance = TRUE
)
# OOB estimate of  error rate: 16.67%

rf.pred <-
  predict (rf.tissues.all.genes, newdata = expr4T.filtered[-train.expr4T.data,])
1 - mean(rf.pred == expr4T.filtered[-train.expr4T.data,]$tissue)



tree.expr.clus <-
  tree(
    tissue ~ ENSG00000225972.1_MTND1P23 + ENSG00000225630.1_MTND2P28 + ENSG00000237973.1_hsa.mir.6723 + ENSG00000229344.1_RP5.857K21.7 + ENSG00000154146.8_NRGN + ENSG00000131095.7_GFAP + ENSG00000266844.1_RP11.862L9.3 + ENSG00000197971.10_MBP + ENSG00000123560.9_PLP1 + ENSG00000226958.1_CTD.2328D6.1 + ENSG00000198695.2_MT.ND6 + ENSG00000101405.3_OXT + ENSG00000101200.5_AVP,
    data = expr4T.filtered,
    subset = train.expr4T.data
  )
yhat.tree <-
  predict(tree.expr.clus, newdata = expr4T.filtered[-train.expr4T.data,], type = "class")
newdata.test <- expr4T.filtered[-train.expr4T.data,]$tissue
error.tree.clus <- 1 - mean(yhat.tree == newdata.test)

cv.tissue <- cv.tree(tree.expr.clus, FUN = prune.misclass)
best.set <- min.set.selection(cv.tissue)

prune.expr <- prune.misclass(tree.expr.clus, best = best.set)
tree.pred <-
  predict(prune.expr, expr4T.filtered[-train.expr4T.data,], type = "class")
error.prune.tree.clust <-
  1 - mean(tree.pred == expr4T.filtered[-train.expr4T.data,]$tissue) #Test error rate


tree.all.genes <- tree(tissue ~ .,
                       data = expr4T.filtered,
                       subset = train.expr4T.data)
yhat.tree <-
  predict(tree.all.genes, newdata = expr4T.filtered[-train.expr4T.data,], type = "class")
newdata.test <- expr4T.filtered[-train.expr4T.data,]$tissue
error.tree <- 1 - mean(yhat.tree == newdata.test)

cv.tissue <- cv.tree(tree.all.genes, FUN = prune.misclass)
best.set <- min.set.selection(cv.tissue)

prune.expr <- prune.misclass(tree.all.genes, best = best.set)
tree.pred <-
  predict(prune.expr, expr4T.filtered[-train.expr4T.data,], type = "class")
error.prune.tree <-
  1 - mean(tree.pred == expr4T.filtered[-train.expr4T.data,]$tissue) #Test error rate


error.tree
error.tree.clus
error.prune.tree
error.prune.tree.clust
