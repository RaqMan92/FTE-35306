
expr4T <- read.csv("M:/Pattern/Project/expr4T.dat", sep="")
expr4T <- read.csv("~/Desktop/expr4T.dat", sep="")


##Packages required##
install.packages("tree")
install.packages("randomForest")
install.packages("glmnet")
install.packages("gbm")
install.packages("e1071")

##Libraries##
library(tree)
library(randomForest)
library(glmnet)
library(e1071)
library(gbm)



#############
# FUNCTIONS #
#############

#SELECTION OF TISSUES FOR DATA SET
tissue.selection <- function(tissue1, tissue2){ #add data variable
  dat = data.frame(expr4T)
  mydata = dat[which(expr4T$tissue==tissue1|expr4T$tissue==tissue2),]
  mydata = droplevels(mydata)
  print(dim(mydata))
  print(table(mydata$tissue))
  return (mydata)
}

#CREATION OF TRAIN VARIABLE, CAN CHOOSE IN A RANGE OF PERCENTAGE
train.selection <- function(percentage, data, seed = as.numeric(1)){
  set.seed(seed)
  train = sample(1:nrow(data), round(nrow(data)*percentage))
  return (train)
}

#SELECTION OF BEST SIZE FOR TREES, CROSS VALIDATION
min.set.selection <- function(cvset){ #rename
  min.dev <- which.min(cvset$dev)
  return (cvset$size[min.dev])
  
}

###This function is not working!!Problems with - symbol in test.
yhat.function <- function(tree, dataset ,y ,test.set){
 
  yhat <- predict(tree, newdata = test.set , type = "class")
  test <- test.set$y
  error.test <- 1-mean(yhat==test)
  #print(table(yhat, dataset[-train,]$y))
  return(error.test)
}

# Not necessary

# #Random forest function
# rf.function <-function(train.x, train.y, test.x,test.y, p.value){
#   rf.model <- randomForest(train.X,
#                            train.Y,
#                            xtest = test.X,
#                            ytest = test.Y,
#                            mtry = p.value,
#                            ntree = 500)
# }
#   
# ##Not working, it's giving a 0% error. I don't know why.
# rf.p.selected <- function(y.value, dataset, train.value, best.p){
#   rf.model <- randomForest(
#     y.value ~ .,
#     data = dataset, 
#     subset = train.value,
#     mtry = best.p,
#     ntree = 500,
#     importance = TRUE
#   )
# }


#######################
## PREPROCESSING DATA #
#######################

#SELECTION OF SMALL DATASET WITH ONLY SELECTED TISSUES
tissue1 = "brain_amygdala"
tissue2 = "brain_putamen"
tissue.data <- tissue.selectionn(tissue1, tissue2)

##Creation of a newdata frame so we can add the column of tissues and have just a few variables.
#A small data set is selected due to we have too much variables.
tissue <- tissue.data$tissue
newdata = data.frame(tissue.data[,1:10], tissue)

#SELECTION OF TRAIN VALUES
#The function will take the determined percentage to get the training indexes. (It is rounded)

train <- train.selection(0.35, newdata) #Default seed = 1



#############
# MAIN CODE #
#############

##DECISION TREES
tree.expr <- tree(tissue~., 
                  data = newdata,
                  subset = train)

####If we increase the number of variables in the model the graph change. 
#So what we can do is to select the most significant genes and rebuild the tree.
summary(tree.expr)
plot(tree.expr)
text(tree.expr)

##Unpruned tree to compare with the pruned one.
yhat.tree <- predict(tree.expr, newdata = newdata[-train,], type = "class")
newdata.test <- newdata[-train,]$tissue
plot(yhat.tree, newdata.test)
error.tree <-1-mean(yhat.tree==newdata.test)
##confusing matrix
(table(yhat.tree, newdata.test))
##Nice graph! :)


#CROSS VALIDATION, prune tree
cv.tissue <- cv.tree(tree.expr, FUN=prune.misclass)
best.set <- min.set.selection(cv.tissue)

prune.expr <- prune.misclass(tree.expr, best = best.set)
summary(prune.expr)
tree.pred <- predict(prune.expr, newdata[-train,], type = "class")
table(tree.pred, newdata[-train,]$tissue)
error.prune.tree <- 1-mean(tree.pred==newdata[-train,]$tissue) #Test error rate

#Comparison bt two errors, prune is this case, is not better.
error.tree 
error.prune.tree


###More variables, the test error rates decreases.BUT too much variables also 
#leads to increase the error?
###Do graph to compare.

##PLOTTING THE TREES
plot(tree.expr)
text(tree.expr, pretty = 0)

plot(prune.expr)
text(prune.expr)


##REGRESSION TREES

tree.regression <- tree(ENSG00000271043.1_MTRNR2L2 ~ . - tissue,
                        data = tissue.data,
                        subset = train,
                        method = "anova")
summary (tree.regression)

##PREDICTION, CALCULATION OF MSE
yhat <- predict(tree.regression, newdata = tissue.data[-train,])
tree.test <- tissue.data[-train,"ENSG00000271043.1_MTRNR2L2"]
mean((yhat - tree.test)^2) ##MSE

##Predictions are not right!!! You cannot obtain the exactly same number as a prediction, 
#it's very difficult. Even if you round them.
1-mean(yhat== tissue.data[-train,]$ENSG00000271043.1_MTRNR2L2) #test error rate


#CROSS VALIDATION
cv.regression <- cv.tree(tree.regression, FUN = prune.tree)
best.set <- min.set.selection(cv.regression)
prune.expr.reg <- prune.tree(tree.regression, best = best.set)


# Make predictions based on this pruned tree.
yhat.prune <- predict(prune.expr.reg, newdata = tissue.data[-train,])
mean((yhat.prune - tree.test)^2)#MSE
1-mean(yhat.prune==tissue.data[-train,]$ENSG00000271043.1_MTRNR2L2)#test error rate

##I think the error is that the tree is taking only ONE variable!! Check it in summary
#Not pruned, same result... something must be wrong! Numbers do not fit.

##PLOTTING THE TREE
plot(tree.regression)
text(tree.regression, pretty = 0)

plot(prune.expr.reg)
text(prune.expr.reg)

#-----------------------------------------------------------------------
#################################
# SVM with two types of kernel  #

# Before we perform a SVM, we need to decide the parameters so use tune
set.seed(42)
tune.out.linear <- tune(svm,
  tissue ~ .,
  data = newdata ,
  kernel = "linear",
  ranges = list(cost = seq(0.1, 1.0, .1))
)
tune.out.linear


# Perform a SVM using a linear kernel on the classification problem.
    ##WHY 0.1???
svm.linear <- svm(tissue ~ .,
      data = newdata,
      kernel = "linear",
      cost = 0.1)

summary(svm.linear)

# Make a prediciton of the SVM result of the linear kernel.
ypred.linear <- predict(svm.linear, newdata)
table(predict = ypred.linear , truth = newdata$tissue)
1 - mean(ypred.linear ==  newdata$tissue) #test error rate

set.seed(42)
# Performing from c(1, 10, 20, 30, 40, 50, 60), lead to 40.
# c(40, 41, 42, 43, 44, 46, 48) leads to 43
# seq(42.7, 44, 0.1) leads to 43.2
tune.out.poly <- tune(
  svm,
  tissue ~ .,
  data = newdata ,
  kernel = "polynomial",
  ranges = list(cost = seq(1, 50, 1))
)
tune.out.poly

# Perform a SVM using a polynomial kernel on the classification problem.
svm.poly <-   svm(tissue ~ .,
                  data = newdata,
                  kernel = "polynomial",
                  cost = 43.2)
summary(svm.poly)
# Make a prediciton of the SVM result of the polynomial kernel.
ypred.poly <- predict(svm.linear, newdata)
table(predict = ypred.poly , truth = newdata$tissue)
1 - mean(ypred.poly ==  newdata$tissue) #test error rate



# Perform a SVM using a linear kernel on the REGRESSION problem.
svm.linear.r <-
  svm(
    ENSG00000225972.1_MTND1P23 ~ . - tissue,
    data = newdata,
    kernel = "linear",
    cost = 1
  )

# Perform a SVM using a polynomial kernel on the regression problem.
ypred.linear.r <- predict(svm.linear.r, newdata)
table(
  predict = round(ypred.linear.r, 1) , truth = round(newdata$ENSG00000225972.1_MTND1P23, 1)
)

# Error rate
1 - mean(round(ypred.linear.r, 1) == round(newdata$ENSG00000225972.1_MTND1P23, 1))
# MSE
mean((
  round(ypred.linear.r, 1) - round(newdata$ENSG00000225972.1_MTND1P23, 1)
) ^ 2)

# Perform a SVM using a polynomial kernel on the regression problem.
svm.poly.r <-   svm(
  ENSG00000225972.1_MTND1P23 ~ . - tissue,
  data = newdata,
  kernel = "polynomial",
  cost = 43.2,
)

# Perform a SVM using a polynomial kernel on the regression problem.
ypred.poly.r <- predict(svm.poly.r, newdata)
table(
  predict = round(ypred.poly.r, 1) , truth = round(newdata$ENSG00000225972.1_MTND1P23, 1)
)
1 - mean(round(ypred.poly.r, 1) == round(newdata$ENSG00000225972.1_MTND1P23, 1)) #test error rate
mean((round(ypred.poly.r, 1) - round(newdata$ENSG00000225972.1_MTND1P23, 1)) ^ 2) # MSE

##Compare the performances!!

#-----------------------------------------------------------------------
# Use different P's to see the effect on randomForest.
p <- ncol(newdata) - 1
p.2 <- p / 2
p.3 <- sqrt(p)

# Create different train and test sets for the data for the classification problem

train.X <- newdata[train,- 11]
test.X <- newdata[-train,- 11]
train.Y <- newdata[train, 11]
test.Y <- newdata[-train, 11]


# Run the randomforest, each of them has a different P. But all of them have the same ntree and datasets.


rf.tissue.p <- rf.function(train.X, train.Y, test.X, test.Y, p)

rf.tissue.p.2 <- rf.function(train.X, train.Y, test.X, test.Y, p.2)

rf.tissue.p.3 <- rf.function(train.X, train.Y, test.X, test.Y, p.3)

# See the mean error rate of each method of P.
mean(rf.tissue.p$test$err.rate)
mean(rf.tissue.p.2$test$err.rate) # Lowest error rate, so we use this.
mean(rf.tissue.p.3$test$err.rate)

# Perform a randomForest on the tissues to see which predictor is important for the tissues.

a <-rf.p.selected(tissue, newdata, train, p.2)#Not working, 0% errors. I leave the actual function as originally below

rf.tissues <- randomForest(
  tissue ~ .,
  data = newdata, #Only 10 variables.
  subset = train,
  mtry = p.2,
  ntree = 500,
  importance = TRUE
)

column.index <- grep("ENSG00000225972.1_MTND1P23", colnames(newdata))

# Remove the tissue data for the regression problem. 
mydat.regression.data <- newdata[,-11]

# Use different P's to see the effect on randomForest.
p.r <- ncol(mydat.regression.data) - 1
p.2.r <- p / 2
p.3.r <- sqrt(p)


# Create different train and test sets for the data for the classification problem
train.X <- mydat.regression.data[train,-1] #column.index is giving me problem, so I am just going to use the number
test.X <- mydat.regression.data[-train,-1]
train.Y <- mydat.regression.data[train, 1]
test.Y <- mydat.regression.data[-train, 1]

# Run the randomforest, each of them has a different P. But all of them have the same ntree and datasets.



rf.tissue.p.re <- rf.function(train.X, train.Y, test.X, test.Y, p.r)

rf.tissue.p.2.re <- rf.function(train.X, train.Y, test.X, test.Y, p.2.r)

rf.tissue.p.3.re <- rf.function(train.X, train.Y, test.X, test.Y, p.3.r)


# See the mean MSE of each method of P.
mean(rf.tissue.p.r$test$mse)
mean(rf.tissue.p.2.r$test$mse) 
mean(rf.tissue.p.3.r$test$mse)

rf.regression <- randomForest(
  ENSG00000225972.1_MTND1P23 ~ . - tissue,
  data = newdata,
  subset = train,
  mtry = p.2.r,
  ntree = 500,
  importance = TRUE
)

mean(rf.regression$mse)

##CONCLUSION?

# Method 3

boost.tissues <- gbm(
  tissue ~ .,
  data = newdata[train,],
  distribution = "gaussian",
  n.trees = 5000 ,
  interaction.depth = 4
)
# See the best five entries
summary(boost.tissues)[0:5,]

# Plot to see how the ENSG00000225972.1_MTND1P23 influences the data.
plot(boost.tissues ,i = "ENSG00000225972.1_MTND1P23")

# Produces numbers instead of tissues strings.
yhat.boost <- predict(boost.tissues,
                      newdata = newdata[-train,],
                      n.trees = 5000)
# # Calculate the test error rate (this is not the MSE).
1 - mean(yhat.boost == newdata[-train,]$ENSG00000225972.1_MTND1P23)
#Giving 100% error something must be wrong

#-----------------------------------------------------------------------
###ex2

#The glmnet() function has an alpha argument that determines what type
#of model is fit. If alpha=0 then a ridge regression model is fit, and if alpha=1
#then a lasso model is fit.

x <- model.matrix(tissue~.,newdata)
y <- newdata$tissue
grid =10^ seq (10,-2, length =100)

#Ridge regression, use of multinomial family due to the several classes of tissue we have
#Also, using lambda = grid we have a wide range of lambda to compare with (possibility of a graph!)
ridge.fit <- glmnet(x,y, alpha=0, lambda = grid, family = "binomial")
ridge.pred <-predict(ridge.fit, type = "coefficients")

#Plot of predictions of two tissues
plot(ridge.pred$brain_amygdala, ridge.pred$brain_hippocampus)
quantile(ridge.pred)
ridge.pred

#Here we do again a ridge regression but with a specific gene
x <- model.matrix(ENSG00000225972.1_MTND1P23~.-tissue, newdata)
y<- newdata$ENSG00000225972.1_MTND1P23
ridge.fit <- glmnet(x[train,],y[train], alpha=0, lambda = grid)
ridge.pred <- predict(ridge.fit, newx = x[-train,], s = 10) ##Changes s, which is lambda the MSE is different
mean((ridge.pred - y[-train])^2) ##test MSE

#The MSE increments with its value, bigger lambda, bigger MSE. If default values is used the MSE is very high (a lot)
##INCLUDE THIS IN THE REPORT!!! WITH A NICE GRAPH!!!

cv.out <- cv.glmnet(x[train,], y[train], alpha = 0, type.measure = "mse")
plot(cv.out)

ridge.pred = predict(ridge.fit, s = cv.out$lambda, newx = x[-train,])
mean((ridge.pred - y[-train])^2)

##We try lasso now
x <- model.matrix(tissue~., expr4T.filtered)
y <- expr4T.filtered$tissue
grid =10^ seq (10,-2, length =100)


###NOW A GLM TO FEATURE SELECTION

glm.fit<- glm(tissue~.,
              data = newdata,
              family = binomial,
              subset = train)

glm.fit$coefficients
summary(glm.fit)


##Jeroen has significant genes, we are going to do two random forest. Here with NOT SIGNIFICANT genes so we can compare.

rf.newdata <- randomForest(tissue~.,
                           data = newdata,
                           subset = train,
                           mtry=10,
                           importance = TRUE) 

yhat.rf <- predict(rf.newdata, newdata = newdata[-train,])
newdata.test <- newdata[-train, "tissue"]
1-mean(yhat.rf == newdata.test)

