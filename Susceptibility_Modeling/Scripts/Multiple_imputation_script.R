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

#perform the imputation step using a few different methods; first using random forest
ImpData_RF <- mice(clean, meth = "rf", ntree = 4)

#get completed data
completed_RF <- complete(ImpData_RF)

#impute using Bayesian linear regression
ImpData_Bayes <- mice(clean, meth = "norm")

#get completed data
completed_Bayes <- complete(ImpData_Bayes)

#impute using arithmetic mean
ImpData_mean <- mice(clean, meth = "mean")

#get completed data
completed_mean <- complete(ImpData_mean)

#combine the different datasets completed by different imputation techniques 
#into one dataset and plot for differences

completed_mean$method = rep("mean", times = nrow(completed_mean))
completed_Bayes$method = rep("Bayes", times = nrow(completed_Bayes))
completed_RF$method = rep("RF", times = nrow(completed_RF))

combined_Imp <- rbind(completed_mean,completed_Bayes,completed_RF)


#get logical vector indicating which columns are numeric in type
nums <- unlist(lapply(combined_Imp, is.numeric))  

#subset column names to only include the numerics
num_cols <- colnames(combined_Imp[nums])

#write function that plots the distribution for each variable comparing across methods
plot_data_column_show = function (data, column) {
  ggplot(data = combined_Imp, aes_string(x = "method", y = column, fill = "method")) +
    geom_boxplot(color = "black") +
    xlab("method")
  
}
#apply the boxplot function to only the numeric columns and print to 'Plots' window
fig_show <- lapply(num_cols, plot_data_column_show, data = combined_Imp)
fig_show

#I'm going to use the random forest method for imputation because it tends to produce fewer
#outlier variables in comparison to the other methods

#write the 'completed_RF' dataframe to file
write_csv(completed_RF, path = here("WNS_Projects","Susceptibility_Modeling","Data","Imputed_Dataset.csv"), col_names = TRUE)
