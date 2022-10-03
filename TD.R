library(abind)
library(MASS)
library(glmnet)
library(rTensor)
#library(caret)
library(ranger)
#library(tidyverse)
#library(ggpubr)
library(abind)
library(liquidSVM)

args = commandArgs(trailingOnly = TRUE)
infile <- as.numeric(args[1])
data1=infile

infile <- as.numeric(args[2])
sweep_train_sample=infile

infile <- as.numeric(args[3])
neut_train_sample=infile

infile <- as.numeric(args[4])
sweep_test_sample=infile

infile <- as.numeric(args[5])
neut_test_sample=infile




#Reading Train data
added <- NULL
for (i in 1:as.numeric(sweep_train_sample))
{
  mean_sweep <- read.csv(paste0("./Data/CSV_files/sweep_train_align_",i,".csv"), header = TRUE,  row.names = 1)
  dim(mean_sweep)
  #mean_sweep<-as.matrix(mean_sweep)
  added<-abind(added, mean_sweep, along = 3)
  print(i)
}
dim(added)



for (i in 1:neut_train_sample)
{
  mean_sweep <- read.csv(paste0("./Data/CSV_files/neut_train_align_",i,".csv"), header = TRUE, row.names = 1)
  dim(mean_sweep)
  #mean_sweep<-as.matrix(mean_sweep)
  added<-abind(added, mean_sweep, along = 3)
  print(i)
}
dim(added)


add<-aperm(added, c(3,1,2))
dim(add)



#Reading Test data

addd <- NULL
for (i in 1:sweep_test_sample)
{
  i
  mean_sweep <- read.csv(paste0("./Data/CSV_files/sweep_test_align_",i,".csv"), header = TRUE, row.names = 1)
  dim(mean_sweep)
  #mean_sweep<-as.matrix(mean_sweep)
  addd<-abind(addd, mean_sweep, along = 3)
  print(i)
}
dim(addd)

for (i in 1:neut_test_sample)
{
  i
  mean_sweep <- read.csv(paste0("./Data/CSV_files/neut_test_align_",i,".csv"), header = TRUE, row.names = 1)
  dim(mean_sweep)
  #mean_sweep<-as.matrix(mean_sweep)
  addd<-abind(addd, mean_sweep, along = 3)
  print(i)
}
dim(addd)


ad<-aperm(addd, c(3,1,2))




# Center and scale tensor (see Rasmus Brow Tutorial on PARAFAC)
meansTrain <- matrix(NA, nrow = 64, ncol = 64) # One for each feature
sdsTrain <- c() # One for each observation
sdsTest <- c() # One for each observation

for(i in 1:(sweep_train_sample+neut_train_sample)) {
  sdsTrain[i] = sqrt(sum(add[i,,]^2))
  add[i,,] = add[i,,] / sdsTrain[i]
  
}

for(i in 1:sweep_test_sample+neut_test_sample) {
  
  sdsTest[i] = sqrt(sum(ad[i,,]^2))
  ad[i,,] = ad[i,,] / sdsTest[i]
}

for(j in 1:64) {
  for(k in 1:64) {
    meansTrain[j,k] = mean(add[,j,k])
   
    add[,j,k] = add[,j,k] - meansTrain[j,k]
    ad[,j,k] = ad[,j,k] - meansTrain[j,k]
  }
}
Ytrain <-data.frame(class = as.factor( c( rep("Sweep", sweep_train_sample), rep("Neutral", neut_train_sample) ) ) )
Y<-data.frame(class = as.factor( c( rep("Sweep", sweep_test_sample), rep("Neutral", neut_test_sample) ) ) )
nrow(Ytrain)
ntrees = 5000
nrandsplits = 10
alphas <- seq(0, 1, by = 0.1)
length(alphas)

# CP Tensor decomposition. 
max_rank_cp <-data1
deviances_cp <- c()
lambdas_cp <- c()
oob_errors_cp <- c()
cv_errors_cp_svm <- c()

  cp_decomp <- cp(as.tensor(add), num_components = data1)
  A <- cp_decomp$U[[1]]
write.csv(A,file="A.csv")
  B <- cp_decomp$U[[2]]
write.csv(B,file="B.csv")
  C <- cp_decomp$U[[3]]
write.csv(C,file="C.csv")

  Xtrain <- sweep( sweep( A, 2, colMeans(A) ), 2, FUN = "/", apply(A, 2, sd) )
  colnames(Xtrain) = paste("Component", seq(1, ncol(Xtrain)), sep = "")
 


 
  # ELASTIC-NET
  for(aIndex in 1:length(alphas)) {
    cvfit <- cv.glmnet(Xtrain, as.matrix(Ytrain), family = "binomial", alpha = alphas[aIndex], nfolds = 10, type.measure = "deviance")

deviances_cp[aIndex] = min(cvfit$cvm)
    
lambdas_cp[aIndex] = cvfit$lambda.min

  }


  # RANDOM FOREST
  trainData <- as.data.frame( cbind(Xtrain, Ytrain) )
 
  rf_trained <- ranger(class ~., data = trainData, num.trees = ntrees, probability = TRUE, oob.error = TRUE)
 
  oob_errors_cp[1] = rf_trained$prediction.error



dim(deviances_cp)
est_rank_cp <- data1
est_rank_cp_rf<-data1
est_rank_cp_svm<-data1

est_alpha_cp <- alphas[ which(deviances_cp == min(deviances_cp), arr.ind = TRUE)]
write.csv(est_alpha_cp,file="./Data/Results/alpha.csv")


est_lambda_cp <- lambdas_cp[which(deviances_cp == min(deviances_cp), arr.ind = TRUE)]
write.csv(est_lambda_cp,file="./Data/Results/lambda.csv")

est_alpha_cp
est_lambda_cp


# TEST FOR ELASTIC NET
Atrain_cp <- cp_decomp$U[[1]]
B <- cp_decomp$U[[2]]
C <- cp_decomp$U[[3]]

lambda_matrix_inv <- diag(1/cp_decomp$lambdas)
mode1_unfolded <- k_unfold(as.tensor(ad), m = 1)@data


Atest_cp <- mode1_unfolded %*% khatri_rao(C, B) %*% ginv( (t(C)%*%C) * (t(B)%*%B ) ) %*% lambda_matrix_inv
write.csv(Atest_cp,file="Atest_cp.csv")

Xtrain <- sweep( sweep( Atrain_cp, 2, colMeans(Atrain_cp) ), 2, FUN = "/", apply(Atrain_cp, 2, sd) )
fit_cp <- glmnet(Xtrain, as.matrix(Ytrain), family = "binomial", alpha = est_alpha_cp, lambda = est_lambda_cp)

X <- sweep( sweep( Atest_cp, 2, colMeans(Atrain_cp) ), 2, FUN = "/", apply(Atrain_cp, 2, sd) )
probs_est_cp <- predict(fit_cp, X, s = est_lambda_cp, type = "response")


probs_est_cp <- data.frame( Sweep = c(probs_est_cp),Neutral = 1-c(probs_est_cp)  )
write.csv(probs_est_cp,file="./Data/Results/Probs_EN.csv")

Yest_cp <- data.frame("Class" = ifelse(probs_est_cp$Sweep >= 0.5, "Sweep", "Neutral"))
write.csv(Yest_cp,file="./Data/Results/Class_EN.csv")



#TEST FOR RANDOM FOREST.



X <- sweep( sweep( Atest_cp, 2, colMeans(Atrain_cp) ), 2, FUN = "/", apply(Atrain_cp, 2, sd) )
colnames(X) = paste("Component", seq(1, ncol(X)), sep = "")
preds <- predict(rf_trained, X)
write.csv(preds,file="./Data/Results/Probs_RF.csv")
probs_est_cp_rf <- as.data.frame(preds$predictions)


Yest_cp_rf <- data.frame("Class" = ifelse(probs_est_cp_rf$Neutral < 0.5, "Sweep", "Neutral"))
write.csv(Yest_cp_rf,file="./Data/Results/Class_RF.csv")





#SVM
svm_trained <- mcSVM(x = Xtrain, y = Ytrain$class, mc_type = "OvA_hinge", folds = 10, scale = FALSE, do.select = TRUE, max_gamma=100000)

cv_errors_cp_svm[1] = mean(svm_trained$select_errors$val_error)


# TEST FOR SUPPORT VECTOR MACHINE


Xtrain <- sweep( sweep( Atrain_cp, 2, colMeans(Atrain_cp) ), 2, FUN = "/", apply(Atrain_cp, 2, sd) )
X <- sweep( sweep( Atest_cp, 2, colMeans(Atrain_cp) ), 2, FUN = "/", apply(Atrain_cp, 2, sd) )
preds <- c(predict(svm_trained, X))
pred_probs <- (preds + 1) / 2
write.csv(pred_probs,file="./Data/Results/Probs_SVM.csv")
probs_est_cp_svm <- data.frame("Sweep" = pred_probs, "Neutral" = 1-pred_probs)
Yest_cp_svm <- data.frame("Class" = ifelse(preds <= 0.5, "Neutral", "Sweep"))
write.csv(pred_probs,file="./Data/Results/Class_SVM.csv")





	

