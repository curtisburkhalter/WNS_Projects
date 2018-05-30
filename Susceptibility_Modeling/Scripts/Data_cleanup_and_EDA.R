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
