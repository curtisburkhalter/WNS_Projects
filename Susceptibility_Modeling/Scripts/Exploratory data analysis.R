#read in required packages
suppressMessages(library(tidyverse))
suppressMessages(library(janitor))
suppressMessages(library(here))

#read in data file containing the White Nose Syndrome status of each North American species;
#data file has been previously cleaned
wns_clean <- read_csv(here("WNS_Projects","Susceptibility_Modeling","Data","Cleaned_WNS_file.csv")) %>%
  clean_names()

#subset the columns where the data is only of type integer, in this case it is the binary columns
int_vars <- colnames(wns_clean)[sapply(wns_clean,is.integer)]

#for the columns that were subsetted into the vector 'int_vars' convert to a factor using
#lapply function
wns_clean[,int_vars] <- lapply(wns_clean[,int_vars],as.factor)
