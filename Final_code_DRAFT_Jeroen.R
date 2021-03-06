#########################################
# Name: Raquel Manzo & Jeroen Lodewijk  #
# Student numbers:
#########################################
#test

#############
# Packages #
#############

# install.packages("car")
# 
# install.packages("class")
# 
# install.packages("leaps")
# 
# install.packages("MASS")

#############
# Libraries #
#############
# Library used for doing a variance inflation factors test.
library(car)
#
library(class)
#
library(leaps)
#
library(MASS)

# Used for the plotting of the PCA.
library(plotrix)

# Library used for performing a randdomForest.
library(randomForest)

# Library is used for the creation of trees.
library(tree)


####################
#### Functions #####
####################
#' GlmPredictionErrorRate
#' Perform a prediction on the test data and calculate the error rate of that prediction.
#' @param glm.fit: Fitted generalized linear model that has been trained by a trainings set. 
#' @param tissue.1: First tissue you want to use for prediction. 
#' @param tissue.2: Second tissue you want to use for prediction. 
#' @param test.data: A test dataset, that doesn't contain any trainings data. 
#' @param type.prediction: Type of prediction that is beinf performed by the predict(). Default on response. 
#'
#' @return Table of predictions and the error rate of the predictions against the test data.
GlmPredictionErrorRate <-
  function(glm.fit, tissue.1, tissue.2, test.data, type.prediction = "response") {
    glm.probs <-  predict(glm.fit, test.data, type = type.prediction)
    glm.pred <- rep(tissue.1, nrow(test.data))
    glm.pred[glm.probs>0.5] <- tissue.2
    print(table(glm.pred, test.data$tissue))
    return(1 - mean(glm.pred == test.data$tissue))
  }

### PredictivePerformanceLm ###
#' Check how well the predictive performance is of the linear model.
#' @param y: the y used in the creation of the linear model.
#' @param data.set: the dataset that has been used for the creation of the linear model.
#' @param training.data: the trainings data indexes.
#' @param lm.training: the linear model, where the trainings data has been apllied to.
#'
#' @return R-squeare and fractions of variability explained by the model.
PredictivePerformanceLm <-
  function(y, data.set, training.data, lm.training) {
    # Create a test dataset to check how much of the variability is explained by the current model.
    # Take all the entries that do not belong to the trainings dataset.
    test <- -training.data
    
    # Create a prediction using the lm of the trainings data on the test data.
    test.pred <- predict(lm.training, newdata = data.set[test, ])
    
    # See how well it predicts the y.
    test.y    <- data.set[test, y]
    
    # Calculate the Sum of Squares total, residual, regression and total. In order to know the fraction of variability explained by the model.
    SS.total      <- sum((test.y - mean(test.y)) ^ 2)
    SS.residual   <- sum((test.y - test.pred) ^ 2)
    SS.regression <- sum((test.pred - mean(test.y)) ^ 2)
    
    # Calculate the R-square (NOT the fraction of variability explained by the model)
    cat("R-square = ", 1 - SS.residual / SS.total)
    cat("\nFraction of variability explained by the model = ",
        SS.regression / SS.total)
  }

#' train.selection
#' Creation of a train variable, chan choose on within a range of a percentage.
#' @param percentage: percentage of the data that is going to be used for the training of the data
#' @param data: data itself where the training selection will be taken from.
#' @param seed: A fixed seed, that is taken to make the sample taken reproducable.
#' @param random.seed: determines if a random seed is even set.
#' 
#' @return train: vector of indexes that contain the training selection.
train.selection <- function(percentage, data, random.seed = TRUE, seed = as.numeric(1)){
  if(random.seed == TRUE){
    # Set a fixed seed to make the data reproducable.
    set.seed(seed)
  }
  
  # Generate a vector of indexes based on the sample and the seed.
  train <- sample(1:nrow(data), round(nrow(data)*percentage))
  return (train)
}

#' error.rate.prediction.trees
#' Predict the error rate of a tree.
#' @param tree.data: the tree that has been build.
#' @param dataset: dataset that has been used in the creation of the tree.
#' @param test.set: a list of indexes indicating the postions of the test dataset.
#' @param type.prediction: what kind of prediction has to be made in the predict().
#'
#' @return error rate of the prediction compared to the test set.
error.rate.prediction.trees <- function(tree.data, dataset, test.set, type.prediction = "class"){
  # Make a prediction using the tree data on the test data set.
  prediction.tree <- predict(tree.data, newdata = dataset[test.set,], type = type.prediction)
  
  # Get the y data out of the dataset, create a list. Otherwise you get errors in the return statement.
  test.data <- expr4T.filtered[test.set,]$tissue
  
  # Error test of the prediction against the test set.
  return( 1 - mean(prediction.tree == test.data) )
}

##Is this a function??
yhat.tree <-
  predict(tree.expr.clus, newdata = expr4T.filtered[-train.expr4T.data,], type = "class")
newdata.test <- expr4T.filtered[-train.expr4T.data,]$tissue
error.tree.clus <- 1 - mean(yhat.tree == newdata.test)
error.tree.clus


#' Cols
#' Creates a range of colors that can be used for plotting.
#' @param vec: is a vector that will be given colors ID's.
#'
#' @return Color ID's for the use in plotting
Cols <- function(vec) {
  cols <- rainbow(length(unique(vec)))
  return(cols[as.numeric(as.factor(vec))])
}

#' min.set.selection
#' Selection of the best size for trees and cross validation.
#' @param cvset: cross validated set.
#'
#' @return Minima set selection
min.set.selection <- function(cvset){ #rename
  min.dev <- which.min(cvset$dev)
  return (cvset$size[min.dev])
}

#' tissue.selection
#' Takes out two tissues
#' @param tissue1: first tissue name as a string.
#' @param tissue2: second tissue name as a string.
#' @param data: data as a data.frame.
#'
#' @return mydata: two tissue dataset as a data.frame
tissue.selection <- function(tissue1, tissue2, data = data.frame(expr4T.filtered)){
  mydata <- data[which(data$tissue == tissue1|data$tissue == tissue2),]
  mydata <- droplevels(mydata)
  print(dim(mydata))
  print(table(mydata$tissue))
  return (mydata)
}

#############
# Main code #
#############

#### WEEK 1 ####
#### Pre-processing ####  

# Load in the dataset.
# expr4T <- read.table("M:/Pattern Recognition/Project/expr4T.dat", sep="")
expr4T <- read.table("F:/Dropbox/Pattern Recognition/Week 1/Friday/expr4T.dat", sep="")
###
dim(expr4T) #Checking dimentions

####REDUCTION OF THE DATA SET####

nn <- ncol(expr4T) - 1 #Last column tissue label
m <- sapply(expr4T[,1:nn], mean)
s <- sapply(expr4T[,1:nn], sd)
sm <- s/m

# Filter genes based on minimum mean and stdev/mean
minsm <- mean(sm)
minm <- mean(m) 

m <- c(m, minm + 1)
sm <- c(sm, minsm + 1)

expr4T.filtered <- expr4T[,which(sm>minsm & m>minm)]
dim(expr4T.filtered)

#Remove variables that we are not going to use more
rm(m, minm, minsm, nn, s, sm)

# Check if these genes are in the dataset.
"ENSG00000271043.1_MTRNR2L2" %in% colnames(expr4T.filtered)
"ENSG00000229344.1_RP5.857K21.7" %in% colnames(expr4T.filtered)

####ANALYSIS OF THE GENE EXPRESSION LEVELS RELATIONSHIP####

#Creation of two tissue data set

tissue1 <- "brain_putamen"
tissue2 <- "brain_cortex"

mydat <- tissue.selection(tissue1, tissue2, data = expr4T.filtered)

#LINEAR REGReSSION MODELS

set.seed(1)
####Training set selection####
train.expr4T.data <-
  sample(1:nrow(expr4T.filtered), round(nrow(expr4T.filtered) * 0.35))

#Train set for two tissue data separated
train.mydat <- sample(1:nrow(mydat), round(nrow(mydat) * 0.35))

#### Performance of a linear model using all the data expected ENSG00000271043.1_MTRNR2L2 gene ####

lm.fit.1 <-
  lm(ENSG00000271043.1_MTRNR2L2 ~ .,
     data = expr4T.filtered,
     subset = train.expr4T.data)

lm.fit.1.two.tissues <-
  lm(ENSG00000271043.1_MTRNR2L2 ~ .,
     data = mydat[,1:70],
     subset = train.mydat)

# Check if there is evidence of non-linearity for ENSG00000271043.1_MTRNR2L2
plot(lm.fit.1)

# Check if there is evidence of non-linearity for ENSG00000271043.1_MTRNR2L2 fpr two tissue data
plot(lm.fit.1.two.tissues)

#### Approach for avoiding patterns in linearity plot(NOT SUCCESS) Remove? It is useless so maybe we can just omit it. ##### 
reduced <- c(names(lm.fit.1$fitted.values[which(lm.fit.1$fitted.values>500)]))
a<- expr4T.filtered[reduced,]

lm.fit.reduced <-
  lm(ENSG00000225972.1_MTND1P23 ~ .,
     data = a,
     subset = train.expr4T.data)
plot(lm.fit.reduced)

rm(reduced, a, lm.fit.reduced)

#### Checking predictions ####

PredictivePerformanceLm(
  y = "ENSG00000271043.1_MTRNR2L2",
  data.set = expr4T.filtered,
  training.data = train.expr4T.data,
  lm.training = lm.fit.1
)

# Check if there is collinearity, if true collinearity is present.
sqrt(vif(lm.fit.1)) > 2

# Data shows collinearity. Explains why the model is highly significant, while few individual predictors are significant.
# Therefore the results of the PredictivePerformanceLm are not trustworty

PredictivePerformanceLm(
  y = "ENSG00000271043.1_MTRNR2L2",
  data.set = mydat,
  training.data = train.mydat,
  lm.training = lm.fit.1.two.tissues
)
# Negative R-square shows that there is dependance of y on x is non-linear. Therefore  negative R-square.

#### Second gene of interest ENSG00000229344.1_RP5.857K21.7 ####

# Perform a linear model using all the terms in te model
lm.fit.2 <-
  lm(ENSG00000229344.1_RP5.857K21.7 ~ .,
     data = expr4T.filtered,
     subset = train.expr4T.data)
lm.fit.2.two.tissues <-
  lm(ENSG00000229344.1_RP5.857K21.7 ~ .,
     data = mydat[,1:70],
     subset = train.mydat)

# Check if there is evidence of non-linearity for ENSG00000271043.1_MTRNR2L2
plot(lm.fit.2)
# Check if there is evidence of non-linearity for ENSG00000271043.1_MTRNR2L2 in two tissue data set
plot(lm.fit.2.two.tissues)

PredictivePerformanceLm(
  y = "ENSG00000229344.1_RP5.857K21.7",
  data.set = expr4T.filtered,
  training.data = train.expr4T.data,
  lm.training = lm.fit.2
)
# Seems the trainings set explains quite a lot of the test set. So it seems like a good model.

PredictivePerformanceLm(
  y = "ENSG00000229344.1_RP5.857K21.7",
  data.set = mydat,
  training.data = train.mydat,
  lm.training = lm.fit.2.two.tissues
)

sqrt(vif(lm.fit.2)) > 2
#Again strong collinearity, results are not trustworthy.

####  OTHER TWO GENES ####
## ENSG00000198695.2_MT.ND6 was selected for being the one with the highest mean
## ENSG00000125144.9_MT1G lowest mean

# Perform a linear model using all the terms in te model
lm.fit.3 <-
  lm(ENSG00000125144.9_MT1G ~ .,
     data = expr4T.filtered,
     subset = train.expr4T.data)

lm.fit.3.two.tissues <-
  lm(ENSG00000125144.9_MT1G ~ .,
     data = mydat, #heRE WE USE THE WHOLE DATA SET, OTHERWISE THE GENE IS NOT FOUND
     subset = train.mydat)

# Check if there is evidence of non-linearity for ENSG00000271043.1_MTRNR2L2
plot(lm.fit.3)
# Check if there is evidence of non-linearity for ENSG00000271043.1_MTRNR2L2 in two tissue data set
plot(lm.fit.3.two.tissues)

PredictivePerformanceLm(
  y = "ENSG00000125144.9_MT1G",
  data.set = expr4T.filtered,
  training.data = train.expr4T.data,
  lm.training = lm.fit.3
)
# Seems the trainings set explains quite a lot of the test set. So it seems like a good model.

PredictivePerformanceLm(
  y = "ENSG00000125144.9_MT1G",
  data.set = mydat,
  training.data = train.mydat,
  lm.training = lm.fit.3.two.tissues
)

sqrt(vif(lm.fit.3)) > 2

#### SUMMARIES OF ALL LINEAR MODELS ####
summary(lm.fit.1) #14 SIGNIFICANT COEFF
#ENSG00000225972.1_MTND1P23                4.203e-02  1.210e-02   3.474 0.000595
#ENSG00000064787.8_BCAS1                   1.277e+00  8.954e-01   1.426 0.154862

summary(lm.fit.1.two.tissues) # 4
#ENSG00000225972.1_MTND1P23                  0.05533    0.03846   1.439   0.2237
#ENSG00000104419.10_NDRG1                   11.57953    3.43752   3.369   0.0281 *

summary(lm.fit.2) # 5
#ENSG00000125144.9_MT1G                   -5.575e+00  1.840e+00  -3.029 0.002681 **
#ENSG00000173267.9_SNCG                    6.329e+00  1.899e+00   3.332 0.000977 ***

summary(lm.fit.2.two.tissues) # 6
#ENSG00000234745.5_HLA.B                     58.3789    20.2977   2.876   0.0452 *
#ENSG00000237973.1_hsa.mir.6723               0.5232     0.1702   3.075   0.0371 *

summary(lm.fit.3) # 11
#tissuebrain_anteriorcortex                6.144e+01  1.379e+01   4.456 1.21e-05 ***
#ENSG00000125148.6_MT2A                    1.746e-01  4.950e-02   3.527 0.000491 ***

summary(lm.fit.3.two.tissues) # 0



# b) How many variables are significant? How do you decide on significance?
# Perform summary on all of them, but if the results show that themodel is highly significant, while few individual predictors are significant. Then collinearity

# Perform a log on everything, besides the tissues. 

###CONVERTION OF THE DATA ####
log.mydat <- log(mydat[-ncol(mydat)])
log.expr4T <- log(expr4T.filtered[-ncol(expr4T.filtered)])

is.finite.data.frame <- function(obj){
  sapply(obj,FUN = function(x) all(is.finite(x)))
}

# Check if there are any infs in the log-transformed datasets.
is.finite.data.frame(log.mydat)
is.finite.data.frame(log.expr4T)

# Works only on matrixes
log.mydat <- as.matrix(log.mydat)

# Replace any non finites with 0, this is needed to get the logs to work.
log.mydat[!is.finite(log.mydat)] <- 0
log.expr4T <- as.matrix(log.expr4T)
log.expr4T[!is.finite(log.expr4T)] <- 0

# Converted back to a data.frame for the lm's.
log.mydat <- as.data.frame(log.mydat)
log.expr4T <- as.data.frame(log.expr4T)


#### Perform a linear model on the log-transformed data. ####
lm.fit.1.log <-
  lm(ENSG00000271043.1_MTRNR2L2 ~ .,
     data = log.expr4T,
     subset = train.expr4T.data)

plot(lm.fit.1.log)

#Predicted values
test.pred <- predict(lm.fit.1.log, newdata = log.expr4T[-train.expr4T.data, ])
#Plot predicted against observed
plot(log.expr4T$ENSG00000271043.1_MTRNR2L2[-train.expr4T.data], test.pred,
     xlab = "Predicted values for ENSG00000271043.1_MTRNR2L2", 
     ylab = "Observed values for ENSG00000271043.1_MTRNR2L2")
abline(0,1, col = "red")

log.mydat <- log.mydat[,1:100] #again we select a small data set, otherwise it doesn't work
#Two tissue model
lm.fit.1.two.tissues.log <-
  lm(ENSG00000271043.1_MTRNR2L2 ~ .,
     data = log.mydat, #again we select a small data set, otherwise it doesn't work
     subset = train.mydat)

plot(lm.fit.1.two.tissues.log) 

#Predicted values
test.pred <- predict(lm.fit.1.two.tissues.log, newdata = log.mydat[-train.mydat, ])
#Plot predicted against observed
plot(log.mydat$ENSG00000271043.1_MTRNR2L2[-train.mydat], test.pred,
     xlab = "Predicted values for ENSG00000271043.1_MTRNR2L2", 
     ylab = "Observed values for ENSG00000271043.1_MTRNR2L2")
abline(0,1, col = "red")

#### Checking predictive performances ####

PredictivePerformanceLm(
  y = "ENSG00000271043.1_MTRNR2L2",
  data.set = log.expr4T,
  training.data = train.expr4T.data,
  lm.training = lm.fit.1.log
)

PredictivePerformanceLm(
  y = "ENSG00000271043.1_MTRNR2L2",
  data.set = log.mydat[,1:50],
  training.data = train.mydat,
  lm.training = lm.fit.1.two.tissues.log
)
#We can see that the R squared for two tissue data is worst.


#### Gene ENSG00000229344.1_RP5.857K21.7 model ####
lm.fit.2.log <-
  lm(ENSG00000229344.1_RP5.857K21.7 ~ .,
     data = log.expr4T,
     subset = train.expr4T.data)
lm.fit.2.two.tissues.log <-
  lm(ENSG00000229344.1_RP5.857K21.7 ~ .,
     data = log.mydat,
     subset = train.mydat)

plot(lm.fit.2.log)
plot(lm.fit.2.two.tissues.log)


test.pred <- predict(lm.fit.2.log, newdata = log.expr4T[-train.expr4T.data, ])


plot(log.expr4T$ENSG00000229344.1_RP5.857K21.7[-train.expr4T.data], test.pred,
     xlab = "Predicted values for ENSG00000229344.1_RP5.857K21.7", 
     ylab = "Observed values for ENSG00000229344.1_RP5.857K21.7")
abline(0,1, col = "red")

test.pred <- predict(lm.fit.2.two.tissues.log, newdata = log.mydat[-train.mydat, ])

plot(log.mydat$ENSG00000229344.1_RP5.857K21.7[-train.mydat], test.pred,
     xlab = "Predicted values for ENSG00000229344.1_RP5.857K21.7", 
     ylab = "Observed values for ENSG00000229344.1_RP5.857K21.7")
abline(0,1, col = "red")

#### Checking predictive performances ####
PredictivePerformanceLm(
  y = "ENSG00000229344.1_RP5.857K21.7",
  data.set = log.expr4T,
  training.data = train.expr4T.data,
  lm.training = lm.fit.2.log
)

PredictivePerformanceLm(
  y = "ENSG00000229344.1_RP5.857K21.7",
  data.set = log.mydat,
  training.data = train.mydat,
  lm.training = lm.fit.2.two.tissues.log
)











####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################



####  WEEK 2 ####




####  WEEK 3  ####

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

# Calculate the within groups sum of squares
wss <-
  (nrow(clustering.data) - 1) * sum(apply(clustering.data, 2, var))

# Loop trough all centers that have been used to create the clustering.
for (i in 1:ncenters) {
  wss[i] <- sum(kmeans(clustering.data,
                       centers = i)$withinss)
}

# Plot the SSE plot to see if clusters are present within the data.
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

## @ TODO: REMOVE THIS AFTER WEEK 2 HAS BEEN ADDED!!
p <- ncol(expr4T.filtered) - 1
p.2 <- p / 2
p.3 <- sqrt(p)


# Perform a randomforst on the pre-processed data of the clustering, predictors were taken from the interesting.genes vector.
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

# Perform a randomforst on the 'normal'data, predictors were taken from the all the genes.
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

# Performance of the randomforest 

# Create a tree that is based on the predictors of the pre-processed data.
tree.expr.clus <-
  tree(
    tissue ~ ENSG00000225972.1_MTND1P23 + ENSG00000225630.1_MTND2P28 + ENSG00000237973.1_hsa.mir.6723 + ENSG00000229344.1_RP5.857K21.7 + ENSG00000154146.8_NRGN + ENSG00000131095.7_GFAP + ENSG00000266844.1_RP11.862L9.3 + ENSG00000197971.10_MBP + ENSG00000123560.9_PLP1 + ENSG00000226958.1_CTD.2328D6.1 + ENSG00000198695.2_MT.ND6 + ENSG00000101405.3_OXT + ENSG00000101200.5_AVP,
    data = expr4T.filtered,
    subset = train.expr4T.data
  )

# @ TODO: Replace this repetative code as a function.
# Perform cross-validation and pruning on the tree
cv.tissue <- cv.tree(tree.expr.clus, FUN = prune.misclass)
best.set <- min.set.selection(cv.tissue)

prune.expr.clus <- prune.misclass(tree.expr.clus, best = best.set)



tree.all.genes <- tree(tissue ~ .,
                       data = expr4T.filtered,
                       subset = train.expr4T.data)



cv.tissue <- cv.tree(tree.all.genes, FUN = prune.misclass)
best.set <- min.set.selection(cv.tissue)

prune.expr.all.genes <- prune.misclass(tree.all.genes, best = best.set)

# Error rate unpruned tree of the clustering data.
error.rate.prediction.trees(tree.data = tree.expr.clus, dataset = expr4T.filtered, test.set = -train.expr4T.data)

# Error rate pruned tree of the clustering data.
error.rate.prediction.trees(tree.data = prune.expr.clus, dataset = expr4T.filtered, test.set = -train.expr4T.data)

# Error rate of unpruned all gene data.
error.rate.prediction.trees(tree.data = tree.all.genes, dataset = expr4T.filtered, test.set = -train.expr4T.data)

# Error rate of pruned all gene data.
error.rate.prediction.trees(tree.data = prune.expr.all.genes, dataset = expr4T.filtered, test.set = -train.expr4T.data)
