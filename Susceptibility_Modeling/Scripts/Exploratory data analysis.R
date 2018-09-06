#read in required packages
suppressMessages(library(tidyverse))
suppressMessages(library(janitor))
suppressMessages(library(here))
suppressMessages(library(scales))

#read in data file containing the White Nose Syndrome status of each North American species;
#data file has been previously cleaned
wns_clean <- read_csv(here("WNS_Projects","Susceptibility_Modeling","Data","updated_cleaned_WNS_file09042018.csv")) %>%
  clean_names()

######################################################################################################
#THIS SECTION HERE IS REALLY GOOD IF THERE ARE LOTS OF CATEGORICAL VARIABLES BUT THE UPDATED COVARIATE
#SET DOES NOT INCLUDE MANY SO I"M GOING TO COMMENT THIS SECTION OUT

# #subset the columns where the data is only of type integer, in this case it is the binary columns
# int_vars <- colnames(wns_clean)[sapply(wns_clean,is.integer)]
# 
# #for the columns that were subsetted into the vector 'int_vars' convert to a factor using
# #lapply function
# wns_clean[,int_vars] <- lapply(wns_clean[,int_vars],as.factor)
# 
# #subset the columns in wns_clean that are factors
# wns_factors <- colnames(wns_clean)[sapply(wns_clean,is.factor)]
# 
# #use the function base::expand.grid to create all combinations of each factor name, write as dataframe
# #and ensure that the variables are no longer factors by using the argument 'stringsAsFactors = F'
# fac_combs <- data.frame(lapply(expand.grid(wns_factors,wns_factors),as.character),stringsAsFactors = F)
# 
# #remove the rows where column entries are equal
# fac_combs <- fac_combs[!(fac_combs$Var1 == fac_combs$Var2),]
# 
# #sort the rows of fac_combs to be in alphabetical order and transpose the result so that 
# #the result is a long columnar dataframe
# fac_combs.sort <- as.data.frame(t(apply(fac_combs,1,sort)))
# 
# #remove duplicated rows from the sorted dataframe
# unique_combs <- fac_combs.sort[!duplicated(fac_combs.sort[1:2]),]
# 
# #write function that builds 2x2 contigency table for each of the pairwise comparisons of factors
# #stored in each row of unique_combs; each contingency table is stored in list item
# wns_Ftable <- apply(unique_combs,1, function (x) as.data.frame(table(wns_clean[,x])))
# 
# #write function that adds a column called 'Combo' to each dataframe stored within
# #wns_Ftable
# add_df_col <- function (x) {
#  cbind(x, "Combo" = c(paste("Neither", paste(colnames(x)[1],colnames(x)[2],sep=" or "),sep = ":"),
#                   paste(colnames(x)[1], "Only",sep=" "),
#                   paste(colnames(x)[2], "Only",sep=" "), 
#                   paste("Both", paste(colnames(x)[1],colnames(x)[2],sep=" or "),sep = ":"))
#   )}
# 
# wns_Ftable <- lapply(wns_Ftable, add_df_col)
# new_names <- c("Var1","Var2","Freq","Combo")
# 
# #RENAME THE COLUMNS FOR EACH DATAFRAME IN THE LIST
# wns_Ftable <- lapply(wns_Ftable, setNames, new_names)
# 
# #rbind all of the separate dataframes
# wns_ConTbl_df <- data.frame(do.call(rbind,wns_Ftable))
# 
# #add comparison index to each grouping in wns_ConTbldf
# wns_ConTbl_df$Index <- rep(1:nrow(unique_combs),each = nrow(wns_ConTbl_df)/nrow(unique_combs))
# 
# #add a 'ComboCode' column to wns_ConTbl_df
# wns_ConTbl_df$ComboCode <- paste(wns_ConTbl_df$Var1,wns_ConTbl_df$Var2,sep="")
# 
# #add the actual variable names to the dataframe 'wns_ConTbl_df'
# wns_ConTbl_df$Var1Name <- rep(as.character(unique_combs[,1]),each=4)
# wns_ConTbl_df$Var2Name <- rep(as.character(unique_combs[,2]),each=4)
# head(wns_ConTbl_df)
# 
# #create labels that will be used in the plotting
# labels = paste(abbreviate(gsub("_"," ",wns_ConTbl_df$Var1Name), min = 3),abbreviate(gsub("_"," ",wns_ConTbl_df$Var2Name), min = 3),sep = "/")
# panel_labels <- labels[!duplicated(labels)]
# 
# #create facetted plot showing frequency of binary variables
# #from 2x2 contingency tables
# ggplot(wns_ConTbl_df, aes(x = ComboCode, y = Freq, fill = ComboCode)) +
#   geom_bar(stat = "Identity") +
#   facet_wrap(~factor(Index,labels = panel_labels),scales = "free_y") +
#   guides(fill = FALSE) +
#   xlab("Combination Code") +
#   ylab("Frequency")
# 
# ggsave("Contingency Table Plots for Binary Variables.png",device = "png", path = here("WNS_Projects","Susceptibility_Modeling","EDA_Figures"), dpi = 400)
#######################################################################################################################################

#convert WNS status to factor
wns_clean$wns_status <- as.factor(wns_clean$wns_status)

#count the number of 0's and 1's indicating WNS status for each genus and plot
wns_clean %>%
  group_by(msw05_genus, wns_status) %>%
  summarize(count = n()) %>%
ggplot(., aes(fill=wns_status, y=count, x=msw05_genus)) + 
  geom_bar(position="dodge", stat="identity") +
  xlab("Genus") +
  ylab("Count") +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_discrete(name = "WNS Status", breaks = c("0","1"), labels = c("WN Negative", "WN Positive"))

ggsave("Count of WNS status for each genus.png",device = "png", path = here("WNS_Projects","Susceptibility_Modeling","EDA_Figures"), dpi = 400)

#count the number of 0's and 1's indicating WNS status for each foraging strategy type and plot
wns_clean %>%
  group_by(for_strat_value, wns_status) %>%
  summarize(count = n()) %>%
  ggplot(., aes(fill=wns_status, y=count, x=for_strat_value)) + 
  geom_bar(position="dodge", stat="identity") +
  xlab("Foraging Strategy Type") +
  ylab("Count") +
  scale_fill_discrete(name = "WNS Status", breaks = c("0","1"), labels = c("WN Negative", "WN Positive"))

ggsave("Count of WNS status for each foraging strategy type.png",device = "png", path = here("WNS_Projects","Susceptibility_Modeling","EDA_Figures"), dpi = 400)

#get logical vector indicating which columns are numeric in type
nums <- unlist(lapply(wns_clean, is.numeric))  

#subset column names to only include the numerics
num_cols <- colnames(wns_clean[nums])

#write function that creates a boxplot for each continuous variable against the response variable 'wns_status'
#and writes it to disk
plot_data_column = function (data, column) {
  ggplot(data = wns_clean, aes_string(x = "wns_status", y = column, fill = "wns_status")) +
  geom_boxplot(color = "black") +
  xlab("wns_status")
  
  ggsave(paste("Relationship between ", column, " and WNS status.png", sep = ""),device = "png", path = here("WNS_Projects","Susceptibility_Modeling","EDA_Figures"), dpi = 400)
  
}
#apply the boxplot function to only the numeric columns
fig_list <- lapply(num_cols, plot_data_column, data = wns_clean)
fig_list

#rewrite the 'plot_data_column' function but remove the 'ggsave' portion so that it prints the plots
#to the plot viewer window
plot_data_column_show = function (data, column) {
  ggplot(data = wns_clean, aes_string(x = "wns_status", y = column, fill = "wns_status")) +
    geom_boxplot(color = "black") +
    xlab("wns_status")
  
}
#apply the boxplot function to only the numeric columns
fig_show <- lapply(num_cols, plot_data_column_show, data = wns_clean)
fig_show


