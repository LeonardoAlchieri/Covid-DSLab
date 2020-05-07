suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(forecast))
suppressPackageStartupMessages(library(lmtest))
suppressPackageStartupMessages(library(strucchange))
suppressPackageStartupMessages(library(reshape))
suppressPackageStartupMessages(library(tseries))
suppressPackageStartupMessages(library(lubridate))
suppressPackageStartupMessages(library(FitAR))
suppressPackageStartupMessages(library(propagate))
# The last package is used for uncertainty propagation.
# It allows for quick 2nd degree Taylor approximation 
# or Monte-Carlo evaluation

# Read from a pre-cleaned file. 
Lombardia_tot <- read.csv("Lombardia.csv", colClasses = c('Date', 'integer', 
                                                      'character', 'integer')) %>% 
            dplyr::group_by(DATA) %>% 
            dplyr::summarise(DECESSI = sum(DECESSI))
# 
# This cleaning is meant to speed up things, and allow
# to have day-deaths 

# Define two variable for plotting purposes
Lombardia_tot["Month.Day"] = as.Date(format(Lombardia_tot$DATA, format="%m-%d"), 
                                 format="%m-%d")
Lombardia_tot["Year"] =format(Lombardia_tot$DATA, format="%Y")

# Plot and save to pdf the deaths for all data 
pdf("images/lombardy_deaths_2015_2020.pdf")
ggplot(data=Lombardia_tot, aes(x=Month.Day, y=DECESSI, color=Year)) + 
geom_line() + theme(plot.margin = margin(3,0,3,0, "cm")) + 
labs(title="Deaths in Lombardy, first 5 months, 2015-2019", 
     y='Deaths', x='Day')
dev.off()
ggplot(data=Lombardia_tot, aes(x=Month.Day, y=DECESSI, color=Year)) + 
geom_line() + theme(plot.margin = margin(3,0,3,0, "cm")) + 
labs(title="Deaths in Lombardy, first 5 months, 2015-2019", 
     y='Deaths', x='Day')

# Read from a pre-cleaned file. 
Lombardia <- read.csv("Lombardia.csv", colClasses = c('Date', 'integer', 
                                                      'character', 'integer')) %>% 
            dplyr::filter(DATA < "2020-01-01") %>%
            dplyr::group_by(DATA) %>% 
            dplyr::summarise(DECESSI = sum(DECESSI))
# 
# This cleaning is meant to speed up things, and allow
# to have day-deaths 

# Define two variable for plotting purposes
Lombardia["Month.Day"] = as.Date(format(Lombardia$DATA, format="%m-%d"), 
                                 format="%m-%d")
Lombardia["Year"] =format(Lombardia$DATA, format="%Y")

# Plot and save to pdf the deaths for 2015-2019
pdf("images/lombardy_deaths_2015_2019.pdf")
ggplot(data=Lombardia, aes(x=Month.Day, y=DECESSI, color=Year)) + 
geom_line() + theme(plot.margin = margin(3,0,3,0, "cm")) + 
labs(title="Deaths in Lombardy, first 5 months, 2015-2019", 
     y='Deaths', x='Day')
dev.off()

ggplot(data=Lombardia, aes(x=Month.Day, y=DECESSI, color=Year)) + 
geom_line() + theme(plot.margin = margin(3,0,3,0, "cm")) + 
labs(title="Deaths in Lombardy, first 5 months, 2015-2019", y='Deaths', x='Day')

# We create a time series with only the deaths.
# We do not take into account the date, and out
# each year next to the other
tsLomb = ts(Lombardia$DECESSI)

# This model is the best one we found
# through trial and error
modLomb = Arima(tsLomb, order=c(3,0,2))
summary(modLomb)

coeftest(modLomb)
# All of the coeffients are significant

pdf("images/ACF.pdf")
ggAcf(modLomb$residuals, lag.max = 200) + labs(title="ACF for Lombardy Arima model")
dev.off()
# We can notice a small problem around lag 120, 
# but we can infer it's in the 5% error.
# We strenthen that our model is not seasonal
ggAcf(modLomb$residuals, lag.max = 200) + labs(title="ACF for Lombardy Arima model")

pdf("images/check_residuals.pdf")
checkresiduals(modLomb, test=FALSE)
dev.off()
checkresiduals(modLomb)
# Normality test is passed, thus our model 
# is decent

LjungBoxTest(residuals(modLomb), k = 1, lag.max = 20)

#Normality test
jarque.bera.test(residuals(modLomb))
# According to this test it is not normal
# But, looking at the graph above, we can infer that 
# it is most likely due to a few outliers

mean(tsLomb)

pdf("images/fit.pdf")
plot(tsLomb, col='blue', lty=1, main="Train vs Fitted", ylab="Deaths", xlab="Day since Jan 01, 2015")
lines(modLomb$fitted, lwd=5)
dev.off()
plot(tsLomb, col='blue', lty=1, main="Train vs Fitted", ylab="Deaths", xlab="Day since Jan 01, 2015")
lines(modLomb$fitted, lwd=5)
# The fitting of the model is very accurate

N = 67
futForecast = forecast(modLomb,h=N, level=c(99))
pdf("images/forecast.pdf")
plot(futForecast, type='l', main=paste("Forecast for first", N, "days of 2020"),
    xlab='Day since Jan 01, 2015', ylab='Deaths') 
dev.off()
plot(futForecast, type='l', main=paste("Forecast for first", N, "days of 2020"),
    xlab='Day since Jan 01, 2015', ylab='Deaths') 

# Load the data again but consider only for 2020
Lombardia2020 <- read.csv("Lombardia.csv", colClasses = c('Date', 'integer', 
                                                          'character', 'integer')) %>% 
            dplyr::filter(DATA > "2020-01-01" & DATA < "2020-03-27") %>%
            dplyr::group_by(DATA) %>% 
            dplyr::summarise(DECESSI = sum(DECESSI))
# I know, loading twice the same dataset is bad pratice,
# but it's fast so there is no big downside

# Subtract forecast value to real numbers
# We select a 90% confidence interval
futurVal <- forecast(modLomb,h=length(Lombardia2020$DECESSI), level=c(99))
futDf = data.frame(futurVal)

Lombardia2020 = cbind(Lombardia2020,futDf)

Lombardia2020["COVID_DEC"] = Lombardia2020$DECESSI - Lombardia2020$Point.Forecast
Lombardia2020["COVID_DEC_HI"] = Lombardia2020$DECESSI - Lombardia2020$Lo.99
Lombardia2020["COVID_DEC_LO"] = Lombardia2020$DECESSI - Lombardia2020$Hi.99

pdf("images/lombardy_clean.pdf")
ggplot(Lombardia2020, aes(DATA, COVID_DEC, ymin = COVID_DEC_LO, ymax = COVID_DEC_HI)) +
geom_line() + geom_ribbon(aes(ymin=COVID_DEC_LO, ymax=COVID_DEC_HI), alpha=0.2) + 
labs(title="Deaths in Lombardy, scaled to forecast values. (2020)", y='Deaths', x='Day (2020)')
dev.off()

ggplot(Lombardia2020, aes(DATA, COVID_DEC, ymin = COVID_DEC_LO, ymax = COVID_DEC_HI)) +
geom_line() + geom_ribbon(aes(ymin=COVID_DEC_LO, ymax=COVID_DEC_HI), alpha=0.2) + 
labs(title="Deaths in Lombardy, scaled to forecast values. (2020)", y='Deaths', x='Day (2020)')

# Read recorder covid-19 deaths
lomb <- read.csv("regione_giorno.csv") %>%
  dplyr::filter(denominazione_regione == "Lombardia") %>%
  dplyr::select(c(data, deceduti )) %>%
  arrange(data)
lomb$data = as.Date(lomb$data,format="%Y-%m-%d")

plot(x=lomb$data, y=lomb$deceduti, type='l')

# The data is cumulative, so I need to do 
# running differences
lomb["diff"] = c(6, diff(lomb$deceduti))
# Change the name to fit the schema of the analysis
names(lomb)[names(lomb) == "data"] <- "DATA"
# Merge the two dataset
final_df = merge(Lombardia2020, lomb, by="DATA")

pdf("images/confront.pdf")
ggplot(final_df, mapping=aes(x=DATA, ymin = COVID_DEC_LO, ymax = COVID_DEC_HI)) +
geom_line(aes(y=COVID_DEC, colour='1')) + 
geom_ribbon(aes(ymin=COVID_DEC_LO, ymax=COVID_DEC_HI), alpha=0.2) + 
geom_line(mapping=aes(y=diff, colour='2')) + 
labs(title="Comparison between estimated and reported covid-19 deaths in Lombardy", y='Deaths', x='Day (2020)') +
scale_color_discrete(name = "", labels = c("Estimated", "Reported"))
dev.off()

ggplot(final_df, mapping=aes(x=DATA, ymin = COVID_DEC_LO, ymax = COVID_DEC_HI)) +
geom_line(aes(y=COVID_DEC, colour='1')) + 
geom_ribbon(aes(ymin=COVID_DEC_LO, ymax=COVID_DEC_HI), alpha=0.2) + 
geom_line(mapping=aes(y=diff, colour='2')) + 
labs(title="Comparison between estimated and reported covid-19 deaths in Lombardy", y='Deaths', x='Day (2020)') +
scale_color_discrete(name = "", labels = c("Estimated", "Reported"))

final_df[c('Point.Forecast', 'DECESSI', 'Lo.99', 'Hi.99', 'deceduti')] = NULL 

names(final_df)[names(final_df) == "COVID_DEC"] = "Estimated.Deceased"

names(final_df)[names(final_df) == "COVID_DEC_LO"] = "Estimated.Deceased.Lower"
names(final_df)[names(final_df) == "COVID_DEC_HI"] = "Estimated.Deceased.Higher"

names(final_df)[names(final_df) == "diff"] = "Reported.Deceased"

# Evaluate the deviances for each point
radii = final_df$Estimated.Deceased.Higher - final_df$Estimated.Deceased
# Create a joiend dataset with the estimated values
# and their deviances
df_for_propagation = cbind(final_df["Estimated.Deceased"], radii)
# For the library we use, I have to take the transpose
df_for_propagation = t(df_for_propagation)
# We have to assign names to the columns
names_arr = array()
for (i in seq(1,ncol(df_for_propagation))){
    names_arr = c(names_arr, paste0("x",i))
}
names_arr = names_arr[-1]
colnames(df_for_propagation) = names_arr
# The library works with expression, thus 
# I have to create one
my_sum = ""
for (name in colnames(df_for_propagation)){
    if(name == "x32"){
        my_sum = paste0(my_sum, name, sep="")
    }
    else{
        my_sum = paste0(my_sum, name, sep=" + ")
    }

}
EXPR1 = parse(text=my_sum)

# Evaluate propagation of uncertainty with 
# 2nd-degree Taylor approximation and 
# Monte-Carlo simulation
RES1 <- propagate(expr = EXPR1, data = df_for_propagation, type = "stat", 
                  do.sim = TRUE, verbose = TRUE, cov=TRUE,
                  nsim = 1e7)

summary(RES1)

# I report here the final results.
# They have been written by hand for laziness
paste("For the period between", final_df[1,]$DATA, " and", final_df[nrow(final_df),]$DATA,
     "we estimated the number of deceased due to covid-19.")  %>% cat()
paste("\nWe found that the total number of deaths is between", 9305 - 453, "and ", 9305 + 453, ".
        This results have been estimated with a 99% confidence interval.") %>% cat()
paste("\nBy comparison, the total number of reported cases in the same time period is:", sum(final_df$Reported.Deceased)) %>% 
cat()
paste(".\nThere is a ",as.integer((9305 - 453)/sum(final_df$Reported.Deceased)*100),"% to ", as.integer((9305 + 453)/sum(final_df$Reported.Deceased)*100), 
"% increase of covid-19 related deaths, according to our calculations.") %>% cat()
