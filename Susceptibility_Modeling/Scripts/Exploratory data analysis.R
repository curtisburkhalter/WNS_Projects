#read in required packages
suppressMessages(library(tidyverse))
suppressMessages(library(janitor))
suppressMessages(library(here))

#read in data file containing the White Nose Syndrome status of each North American species;
#data file has been previously cleaned
wns_clean <- read_csv(here("WNS_Projects","Susceptibility_Modeling","Data","Cleaned_WNS_file.csv")) %>%
  clean_names()

#the only values that wns_clean$quarries provides is either a 0 or an NA, as such
#it provides no information so remove it
wns_clean$quarries <- NULL

#subset the columns where the data is only of type integer, in this case it is the binary columns
int_vars <- colnames(wns_clean)[sapply(wns_clean,is.integer)]

#for the columns that were subsetted into the vector 'int_vars' convert to a factor using
#lapply function
wns_clean[,int_vars] <- lapply(wns_clean[,int_vars],as.factor)

#subset the columns in wns_clean that are factors
wns_factors <- colnames(wns_clean)[sapply(wns_clean,is.factor)]

#use the function base::expand.grid to create all combinations of each factor name, write as dataframe
#and ensure that the variables are no longer factors by using the argument 'stringsAsFactors = F'
fac_combs <- data.frame(lapply(expand.grid(wns_factors,wns_factors),as.character),stringsAsFactors = F)

#remove the rows where column entries are equal
fac_combs <- fac_combs[!(fac_combs$Var1 == fac_combs$Var2),]

#sort the rows of fac_combs to be in alphabetical order and transpose the result so that 
#the result is a long columnar dataframe
fac_combs.sort <- as.data.frame(t(apply(fac_combs,1,sort)))

#remove duplicated rows from the sorted dataframe
unique_combs <- fac_combs.sort[!duplicated(fac_combs.sort[1:2]),]

#write function that builds 2x2 contigency table for each of the pairwise comparisons of factors
#stored in each row of unique_combs; each contingency table is stored in list item
wns_Ftable <- apply(unique_combs,1, function (x) as.data.frame(table(wns_clean[,x])))

#write function that adds a column called 'Combo' to each dataframe stored within
#wns_Ftable
add_df_col <- function (x) {
 cbind(x, "Combo" = c(paste("Neither", paste(colnames(x)[1],colnames(x)[2],sep=" or "),sep = ":"),
                  paste(colnames(x)[1], "Only",sep=" "),
                  paste(colnames(x)[2], "Only",sep=" "), 
                  paste("Both", paste(colnames(x)[1],colnames(x)[2],sep=" or "),sep = ":"))
  )}

wns_Ftable <- lapply(wns_Ftable, add_df_col)
new_names <- c("Var1","Var2","Freq","Combo")

#RENAME THE COLUMNS FOR EACH DATAFRAME IN THE LIST
wns_Ftable <- lapply(wns_Ftable, setNames, new_names)

#rbind all of the separate dataframes
wns_ConTbl_df <- data.frame(do.call(rbind,wns_Ftable))

#add comparison index to each grouping in wns_ConTbldf
wns_ConTbl_df$Index <- rep(1:nrow(unique_combs),each = nrow(wns_ConTbl_df)/nrow(unique_combs))

#add a 'ComboCode' column to wns_ConTbl_df
wns_ConTbl_df$ComboCode <- paste(wns_ConTbl_df$Var1,wns_ConTbl_df$Var2,sep="")

#add the actual variable names to the dataframe 'wns_ConTbl_df'
wns_ConTbl_df$Var1Name <- rep(as.character(unique_combs[,1]),each=4)
wns_ConTbl_df$Var2Name <- rep(as.character(unique_combs[,2]),each=4)
head(wns_ConTbl_df)

#create labels that will be used in the plotting
labels = paste(abbreviate(gsub("_"," ",wns_ConTbl_df$Var1Name), min = 3),abbreviate(gsub("_"," ",wns_ConTbl_df$Var2Name), min = 3),sep = "/")
panel_labels <- labels[!duplicated(labels)]

#create facetted plot showing frequency of binary variables
#from 2x2 contingency tables
ggplot(wns_ConTbl_df, aes(x = ComboCode, y = Freq, fill = ComboCode)) +
  geom_bar(stat = "Identity") +
  facet_wrap(~factor(Index,labels = panel_labels),scales = "free_y") +
  guides(fill = FALSE) +
  xlab("Combination Code") +
  ylab("Frequency")

ggsave("Contingency Table Plots for Binary Variables.png",device = "png", path = here("WNS_Projects","Susceptibility_Modeling","EDA_Figures"), dpi = 400)

#create plots of the various continous variables to look
#for associations
ggplot(wns_clean,aes(x = min_temp, y = max_temp)) +
  geom_point(color = "blue") +
  geom_smooth(stat = "smooth",color = "red") +
  xlab("Minimum Temperature (C)") +
  ylab("Maximum Temperature (C)")

ggsave("Relationship between min temp and max temp.png",device = "png", path = here("WNS_Projects","Susceptibility_Modeling","EDA_Figures"), dpi = 400)

#ignore the warnings below, there are issues with estimating
#the prediction intervals where there are so many NAs
#but this intended only to give us a very rough idea of
#association between these variables
ggplot(wns_clean,aes(x = min_temp, y = change_lambda)) +
  geom_point(color = "blue") +
  geom_smooth(stat = "smooth",color = "red") +
  xlab("Minimum Temperature (C)") +
  ylab("Change in lambda")

ggsave("Relationship between min temp and change in lambda.png",device = "png", path = here("WNS_Projects","Susceptibility_Modeling","EDA_Figures"), dpi = 400)

ggplot(wns_clean,aes(x = max_temp, y = change_lambda)) +
  geom_point(color = "blue") +
  geom_smooth(stat = "smooth",color = "red") +
  xlab("Maximum Temperature (C)") +
  ylab("Change in lambda")

ggsave("Relationship between max temp and change in lambda.png",device = "png", path = here("WNS_Projects","Susceptibility_Modeling","EDA_Figures"), dpi = 400)

cor(wns_clean[,grep(".*temp|.*change",colnames(wns_clean))], use = "pairwise.complete.obs")


con_tbl_plotdf <- data.frame(Variable = unique(wns_ConTbl_df$Var1Name), Abbrev = unique(abbreviate(gsub("_"," ",wns_ConTbl_df$Var1Name),min = 3)))

write_csv(con_tbl_plotdf, path = here("WNS_Projects","Susceptibility_Modeling","Contingency_Table_Abbreviations.csv"), col_names = TRUE)
