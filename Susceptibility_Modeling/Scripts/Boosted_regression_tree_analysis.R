#read in required packages
suppressMessages(library(tidyverse))
suppressMessages(library(janitor))
suppressMessages(library(here))
suppressMessages(library(arm))
suppressMessages(library(gbm))
suppressMessages(library(caret))

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
