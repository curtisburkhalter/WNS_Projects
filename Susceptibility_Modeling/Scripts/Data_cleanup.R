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


################################################################
################################################################
#READ IN UPDATED WORLDWIDE DATA SHEET
################################################################
################################################################

#read in the worldwide species list with their corresponding
#WNS status
wns_ww <- read_csv(here("WNS_Projects","Susceptibility_Modeling","Data","worldwide_WNS_status.csv")) %>%
  clean_names()

#remove last column; contains no information
wns_ww <- wns_ww[,1:3]

#provide some column names
colnames(wns_ww) <- c("index", "spp", "wns_status")

#set 'wns_status' as factor
wns_ww$wns_status <- factor(wns_ww$wns_status)

#read in updated covariate data and join to wns_ww
wns_all <- read_csv(here("WNS_Projects","Susceptibility_Modeling","Data","worldwide_WNS_covariates.csv")) %>%
  clean_names() %>%
  arrange(pan) %>%
  select(-label) %>%
  left_join(., wns_ww[,2:3], by = c("pan" = "spp"))

#obtain a quick summary of wns_status column
#mostly unknowns (~97%)
summary(wns_all$wns_status)

ggplot(wns_all, aes(x=wns_status, fill = wns_status)) +
  geom_bar() +
  xlab("White nose syndrome status") +
  ylab("Count") +
  labs(fill = "WNS Status") +
  scale_fill_manual(values = hue_pal()(2), breaks = c(0,1,NA),labels = c("Uninfected", "Infected", "Unknown"),na.value = "black") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

#subset worldwide data where WNS status is not an NA
sub_ww <- wns_all[!is.na(wns_all$wns_status),]

#count NAs by column for the data where status is not an NA
NAs_bycol <- as.data.frame(t(sub_ww %>%
            summarise_all(funs(sum(is.na(.))))))

#attach a new column name 'count' which is the count of NAs for the 
#variable across all 34 records of sub_ww
colnames(NAs_bycol)[1] <- "count"

#attach a new column name 'covariate' which is the name of the
#covariates using the row names already attached
NAs_bycol$covariate <- row.names(NAs_bycol)

#delete the row names of 'NAs_bycol'
row.names(NAs_bycol) <- NULL

#determine the percentage of covariate data that is an NA for
#each covariate
NAs_bycol$pctNA <- round((NAs_bycol$count/nrow(sub_ww))*100,digits = 1)

#determine how many rows of NAs_bycol have 0% covariate data missing
nrow(NAs_bycol[NAs_bycol$pctNA < 0.0001,])
summary(NAs_bycol)

#begin the process of cleaning up the covariate data

#after looking the diet and activity variables, all of the species
#are 100% nectar feeding and 100%nocturnal so we can't use any of the
#diet or activity variables
sub_ww <- sub_ww[,-grep(".*diet_i|.*diet_n|.*diet_v|.*diet_s|.*diet_f|.*diet_p|.*activity_crep|.*activity_noc|.*activity_diur",colnames(sub_ww))]

#all of the species that we have records for in terms of their WNS status
#are in the family Vespertilionidae so we can remove all of the family
#names because those are useless now due to lack of variation
sub_ww <- sub_ww[,-grep(".*idae",colnames(sub_ww))]

nrow(NAs_bycol[NAs_bycol$pctNA >= 50,])
