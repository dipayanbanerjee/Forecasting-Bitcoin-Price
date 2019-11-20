rm(list = ls())

setwd("~/Documents/Sem-4/Time Series Analysis")
getwd()

#-------------------------------------------------------#

library(TSA)
library(fUnitRoots)
library(forecast)
library(CombMSC)
library(lmtest)
library(fGarch)
library(rugarch)

Bitcoin<- read.csv("Bitcoin_Historical_Price.csv",sep = ",")
View(Bitcoin)

Bitcoin$Close <- as.numeric(as.character(gsub(",","",Bitcoin$Close)))
BitcoinTs<- ts(as.vector(Bitcoin$Close),start = c(2013,117),frequency = 365.25)

par(mfrow=c(1,1))
plot(BitcoinTs, ylab='Bitcoin Closing Rate',xlab='Time',main='Bitcoin closing rate 2013-2018',type='l')
grid(nx=100,ny=100,col='LightGray',lty=2,lwd = par('lwd'))

# Trend -Yes
# Changing variance-Yes
# Intervention - Yes
# No Seasonality
# Moving Average Behaviour

residual.analysis <- function(model, std = TRUE,start = 2, class = c("ARIMA","GARCH","ARMA-GARCH")[1]){
  # If you have an output from arima() function use class = "ARIMA"
  # If you have an output from garch() function use class = "GARCH"
  # If you have an output from ugarchfit() function use class = "ARMA-GARCH"
  library(TSA)
  library(FitAR)
  if (class == "ARIMA"){
    if (std == TRUE){
      res.model = rstandard(model)
    }else{
      res.model = residuals(model)
    }
  }else if (class == "GARCH"){
    res.model = model$residuals[start:model$n.used]
  }else if (class == "ARMA-GARCH"){
    res.model = model@fit$residuals
  }else {
    stop("The argument 'class' must be either 'ARIMA' or 'GARCH' ")
  }
  par(mfrow=c(3,2))
  plot(res.model,type='o',ylab='Standardised residuals', main="Time series plot of standardised residuals")
  abline(h=0)
  hist(res.model,main="Histogram of standardised residuals")
  acf(res.model,main="ACF of standardised residuals")
  pacf(res.model,main="PACF of standardised residuals")
  qqnorm(res.model,main="QQ plot of standardised residuals")
  qqline(res.model, col = 2)
  print(shapiro.test(res.model))
  k=0
  LBQPlot(res.model, lag.max = 30, StartLag = k + 1, k = 0, SquaredQ = FALSE)
}


ar(diff(BitcoinTs))
adfTest(BitcoinTs,lags=33,title=NULL,description=NULL)

diff1=diff(log(BitcoinTs))
ar(diff(diff1))
adfTest(diff1,lags=32,title = NULL,description = NULL)

acf(diff1,ci.type='ma')
pacf(diff1)
eacf(diff1)


diff2=diff(diff(log(BitcoinTs)))
ar(diff(diff2))
adfTest(diff2,lags=29,title = NULL,description = NULL)

par(mfrow=c(1,2))
acf(diff2,ci.type='ma')
pacf(diff2)
eacf(diff2)

# ARIMA(4,2,5),ARIMA(5,2,5),ARIMA(5,2,6)

#transformation=BoxCox.ar(diff2+abs(min(diff2))+0.01)
#transformation$ci
#lambda=1.35
#BC.Bitcoin=((diff1+abs(min(diff1))+0.01)^(lambda)-1)/lambda
par(mfrow=c(1,1))
qqnorm(diff2)
qqline(diff2,col=2)
plot(diff2,type='o')

source('TSHandy.r')

modelList <- list(c(4,2,5), c(5,2,5), c(5,2,6),c(4,2,6))
modelEstimation <- myCandidate(log(BitcoinTs), orderList = modelList, methodType = "ML")
# Nearly all models have all coefficients significant.
modelEstimation$IC
# The smalles AIC, AICc and BIC come from ARIMA(5,2,6) model which is 3rd model in the modelEstimation object
# But after looking at coefficients (4,2,5) model looks more appropriate.

coeftest(modelEstimation$model[[1]])
coeftest(modelEstimation$model[[2]])

residual.analysis(modelEstimation$model[[1]], std = TRUE,start = 1)

  m425_residuals = modelEstimation$model[[1]]$residuals

abs.res = abs(m425_residuals)
sq.res = m425_residuals^2

par(mfrow=c(1,2))
acf(abs.res, ci.type="ma",main="The sample ACF plot for absolute residual series")
pacf(abs.res, main="The sample PACF plot for absolute residual series")
eacf(abs.res)

# From the EACF, we can identify ARMA(2,2), ARMA(2,3) models for absolute residual series. So the corresponding 
# tentative GARCH models are GARCH(2,2), GARCH(3,2)


par(mfrow=c(1,2))
acf(sq.res, ci.type="ma",main="The sample ACF plot for absolute residual series")
pacf(sq.res, main="The sample PACF plot for absolute residual series")
eacf(sq.res)

# From the EACF, we can identify ARMA(4,4) models for absolute residual series. So the corresponding 
# tentative GARCH models are GARCH(4,4)

model0 <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 2)), 
                     mean.model = list(armaOrder = c(4, 5), include.mean = FALSE), 
                     distribution.model = "norm")
m.11 <- ugarchfit(spec = model0, data = diff2, out.sample = 100)
m.11
plot(m.11)

model1<-ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(2, 2)), 
                   mean.model = list(armaOrder = c(4, 5), include.mean = FALSE), 
                   distribution.model = "norm")
m.22<-ugarchfit(spec = model1, data = diff2, out.sample = 100)
# I'm fitting ARMA model here not ARIMA!
# Therefore I need to send the differenced and transformed series to ugarchfit().
m.22  # AIC = -3.7622
plot(m.22)
9# We display residual plots with selection of 8, 9, 10 and 11. 


model2<-ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(3, 2)), 
                    mean.model = list(armaOrder = c(4, 5), include.mean = FALSE), 
                    distribution.model = "norm")
m.32<-ugarchfit(spec = model2, data = diff2, out.sample = 100)
# I'm fitting ARMA model here not ARIMA!
# Therefore I need to send the differenced and transformed series to ugarchfit().
m.32  # AIC = -3.7629
plot(m.32)
# We display residual plots with selection of 8, 9, 10 and 11. 

model3<-ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(4, 4)), 
                    mean.model = list(armaOrder = c(4, 5), include.mean = FALSE), 
                    distribution.model = "norm")
m.44<-ugarchfit(spec = model3, data = diff2, out.sample = 100)
# I'm fitting ARMA model here not ARIMA!
# Therefore I need to send the differenced and transformed series to ugarchfit().
m.44  # AIC = -3.7593
plot(m.44)
# We display residual plots with selection of 8, 9, 10 and 11. 

forc.32 = ugarchforecast(m.32, data = diff2, n.ahead = 10, n.roll = 10)
plot(forc.32, which = "all")
forc.32

forc.44 = ugarchforecast(m.44, data = diff2, n.ahead = 10, n.roll = 10)
plot(forc.44,which="all")
for.44









