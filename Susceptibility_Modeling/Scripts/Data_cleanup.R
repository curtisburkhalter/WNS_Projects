#read in required packages
suppressMessages(library(tidyverse))
suppressMessages(library(janitor))
suppressMessages(library(here))

#read in data file containing the White Nose Syndrome status of each North American species;
#data file also includes associated habitat and behavioral data

wns_raw <- read_csv(here("WNS_Projects","Susceptibility_Modeling","Data","NABats_WNS.csv")) %>%
  clean_names()

#there are a bunch of extra rows and columns so remove those
wns_raw <- wns_raw[1:48,1:19]

#need to change some column names
colnames(wns_raw)[1] <- "species"
colnames(wns_raw)[grep(".*min_temp",colnames(wns_raw))] <- "min_temp"
colnames(wns_raw)[grep(".*develop",colnames(wns_raw))] <- "disease_present"
colnames(wns_raw)[grep(".*percent",colnames(wns_raw))] <- "change_lambda"

#need to change any row entries that are either 'unk' or 'NULL' to an NA
wns_raw[wns_raw == "unk"|wns_raw == "Null"] <- NA

#change 'disease_present' to factor
wns_raw$disease_present <- as.factor(wns_raw$disease_present)

#write function that checks if a column is binary or not; because the binary columns were
#read in as character you have to use the numbers in quotes as opposed to 'x in 0:1' as if 
#they were numeric
bin_check <- function(x) {
  all(x == "0"|x=="1" | is.na(x)) 
}

#apply the 'bin_check' function to the columns of wns_raw and subset the column names
#where TRUE
bin_vars <- colnames(wns_raw)[apply(wns_raw,2,bin_check)]

#using the subsetted column names in 'bin_var' convert all the selected columns to factors
wns_raw[,bin_vars] <- lapply(wns_raw[,bin_vars], as.factor)

#I need to convert min_temp, max_temp and change_lambda to numeric;
#first build function that checks to see if a column is of type character
#and if any of the column entries contain digits
 num_check <- function(x) {
   is.character(x) & any(grepl("\\d",x))
 }
 
#using the 'num_check' function subset the column names that satisfy the conditions
#laid out in the funciton; when applying the num_check function we have to use
#sapply versus apply. This is because when using is.character(x) the dataframe is
#coerced to a matrix and the factor variables are then treated as character which 
#in turn returns any column that is a factor, which in this case is most of the covariates
num_vars <- colnames(wns_raw)[sapply(wns_raw, num_check)]
 
wns_raw[,num_vars] <- lapply(wns_raw[,num_vars], as.numeric)

wns_cleaned <- wns_raw

write_csv(wns_cleaned,here("WNS_Projects","Susceptibility_Modeling","Data","Cleaned_WNS_file.csv"),col_names = TRUE)
