#read in required packages
suppressMessages(library(tidyverse))
suppressMessages(library(janitor))
suppressMessages(library(here))
suppressMessages(library(arm))
suppressMessages(library(gbm))
suppressMessages(library(caret))
suppressMessages(library(e1071))
suppressMessages(library(RANN))
suppressMessages(library(caTools))

#read in data file containing the White Nose Syndrome status of each species;
#data file has been previously cleaned and missing data imputed
wns_clean <- read_csv(here("WNS_Projects","Susceptibility_Modeling","Data","Imputed_Dataset.csv")) %>%
  clean_names()

#remove the 'method' column
wns_clean$method <- NULL

#package gbm requires that the response variable, in this case 'wns_status'
#to be in an integer but make sure not to standardize it

#convert columns 'genus' and 'for_strat_value' to factors
wns_clean$genus <- as.factor(wns_clean$genus)
wns_clean$for_strat_value <- as.factor(wns_clean$for_strat_value)

#class separation in classification is when the coding of a variable or something as simply
#as the labels themselves perfectly allow for prediction of different classes and has nothing
#correlative or causal to do with the response variable; need to be vigilant for this in
#the variables

#need to scale the numeric variables so that are comparable in terms of mean and 
#deviation to the binary variables; to do this use the 'arm' package; first
#subset those columns that are numeric and remove the response variable from result
num_cols <- colnames(wns_clean)[sapply(wns_clean, is.numeric)]
num_cols <- num_cols[-grep(".*wns_st",num_cols)]

wns_clean[,num_cols] <- apply(wns_clean[,num_cols],2,rescale)

#going to create 'one-hot' encoding for the factor varialbles
GenusMM <- model.matrix(data = wns_clean, wns_status ~ genus - 1)
ForStratMM <- model.matrix(data = wns_clean, wns_status ~ for_strat_value - 1)

wns_clean <- cbind(wns_clean[, !names(wns_clean) %in% c("genus","for_strat_value")], GenusMM, ForStratMM)
View(wns_clean)

#need to create a stratified random sample for the k-fold cross validation
set.seed(10)
TrainSplits <- unlist(createDataPartition(wns_clean$wns_status, p = 0.80))

train_set <- wns_clean[TrainSplits,]
test_set <- wns_clean[-TrainSplits,]

#subset train_set to not include species identity as a covariate
sub_train <- train_set[,-1]

sub_train$wns_status <- as.factor(sub_train$wns_status)
sub_train <-as.data.frame(sub_train)


#build a parameter tuning grid; the 3 parameters that should be tuned for a 
#boosted regression tree are number of trees (ntree), interaction depth,
#and shrinkage(or learning rate)
gbmGrid <- expand.grid(interaction.depth = c(1,2,3),
                       n.trees = (seq(from = 1000, to =8000,100)),
                       shrinkage = c(0.1,0.01,0.001,0.0001,0.00001),
                       n.minobsinnode = c(2,3,5))
                       

#specify how the model training should take place when determining
#the optimal parameters; in this case use 10-fold cross validation
fitControl <- trainControl(
  method = "cv",
  number = 10)


#specify the training settings for the parameter tuning; would be interesting to see 
#how this changes with imputation of missing values
set.seed(825)
gbmTune <- train(wns_status ~ ., data = sub_train,  
                 method = "gbm", 
                 trControl = fitControl, 
                 #preProcess = c("bagImpute"),
                 verbose = FALSE, 
                 tuneGrid = gbmGrid,
                 #na.action = na.pass,
                 metric = "kappa")

#look at tuning results
tune_results <- (na.omit(gbmTune$results))
View(tune_results)

#look at the best tuning parameter settings as chosen using
#accuracy
gbmTune$bestTune

#the parameters chosen as the best parameters include:
#n.trees = 1000; interaction.depth = 2, learning rate = 0.001 and n.minobsinnode = 2

#look at the model summary for 'gbmTune'
summary(gbmTune)

 #write the model summary to dataframe
gbmSummary <- as.data.frame(summary(gbmTune))

#obtain predictions for the test data and then calculate AUC
gbmPredictions.test <- predict(gbmTune, newdata = test_set, type = "prob")
gbmPredictions.train <- predict(gbmTune, newdata = train_set, type = "prob")

full_train <- na.omit(train_set)
full_test <- na.omit(test_set)

#need to calculate AUC still, was getting auc of 1 which can't be right
colAUC(gbmPredictions.test,full_test$wns_status, plotROC = TRUE)
View(gbmPredictions.test)
