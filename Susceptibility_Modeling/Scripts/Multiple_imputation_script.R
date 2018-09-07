#read in required packages
suppressMessages(library(tidyverse))
suppressMessages(library(janitor))
suppressMessages(library(here))
suppressMessages(library(mice))
suppressMessages(library(VIM))

#read in updated cleaned up dataset
clean <- read_csv(here("WNS_Projects","Susceptibility_Modeling","Data","updated_cleaned_WNS_file09042018.csv")) %>%
  clean_names()

#look at patterns of missing data using VIM::aggr function
aggr_plot <- aggr(clean[,4:19], col=c('navyblue','red'), numbers=TRUE, sortVars=TRUE, labels=abbreviate(names(clean)[4:19], minlength = 10), cex.axis=.7, gap=3, ylab=c("Histogram of missing data","Pattern"))

#perform the imputation step using a few different methods
ImpData <- mice(clean, meth = "rf", ntree = 4)

#get completed data
completed <- complete(ImpData)

