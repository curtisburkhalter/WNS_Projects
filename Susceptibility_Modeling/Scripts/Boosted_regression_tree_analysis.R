#read in required packages
suppressMessages(library(tidyverse))
suppressMessages(library(janitor))
suppressMessages(library(here))
suppressMessages(library(arm))
suppressMessages(library(gbm))
suppressMessages(library(caret))
suppressMessages(library(e1071))
suppressMessages(library(RANN))

#read in data file containing the White Nose Syndrome status of each North American species;
#data file has been previously cleaned
wns_clean <- read_csv(here("WNS_Projects","Susceptibility_Modeling","Data","Cleaned_WNS_file.csv")) %>%
  clean_names()

#the only values that wns_clean$quarries provides is either a 0 or an NA, as such
#it provides no information so remove it
wns_clean$quarries <- NULL

#subset the columns where the data is only of type integer, in this case it is the binary columns
int_vars <- colnames(wns_clean)[sapply(wns_clean,is.integer)]

#package gbm requires taht the response variable, in this case 'disease present'
#to be in an integer so remove it from int_vars
int_vars <- int_vars[-grep(".*disease",int_vars)]

#for the columns that were subsetted into the vector 'int_vars' convert to a factor using
#lapply function
wns_clean[,int_vars] <- lapply(wns_clean[,int_vars],as.factor)

#need to remove change in lambda from the dataset b/c it leads to 
#class separation
wns_clean <- wns_clean[,-grep(".*change",colnames(wns_clean))]

#need to scale the numeric variables so that are comparable in terms of mean and 
#deviation to the binary variables; to do this use the 'arm' package; first
#subset those columns that are numeric and remove the response variable from result
num_cols <- colnames(wns_clean)[sapply(wns_clean, is.numeric)]
num_cols <- num_cols[-grep(".*dise",num_cols)]

wns_clean[,num_cols] <- apply(wns_clean[,num_cols],2,rescale)

#need to create a stratified random sample for the k-fold cross validation
set.seed(10)
TrainSplits <- unlist(createDataPartition(wns_clean$disease_present, p = 0.80))

train_set <- wns_clean[TrainSplits,]
test_set <- wns_clean[-TrainSplits,]

#subset train_set to not include species identity as a covariate
sub_train <- train_set[,-1]

sub_train$disease_present <- as.factor(sub_train$disease_present)
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
gbmTune <- train(disease_present ~ ., data = sub_train,  
                 method = "gbm", 
                 trControl = fitControl, 
                 preProcess = c("bagImpute"),
                 verbose = FALSE, 
                 tuneGrid = gbmGrid,
                 na.action = na.pass,
                 metric = "kappa")

#look at tuning results
View(na.omit(gbmTune$results))

#look at the best tuning parameter settings as chosen using
#accuracy
gbmTune$bestTune
