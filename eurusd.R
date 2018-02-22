# Clean environment
rm(list=ls())
gc()
options(java.parameters = "-Xmx8000m")
# Set the directory
setwd('/mnt/c/Users/cyril/Desktop/Master application')

## GLOBAL
# Switch to old/new dataset or brand new with ''
NEW = FALSE
TARGETS = 'EUR.USD'
ALL.TARGET = TRUE
## END GLOBAL

if(is.logical(NEW)) {
  if(NEW==FALSE) {
    NEW = '#o'
  } else {
    NEW = '#new'
  }
}

# Load directories
settings <- NULL
settings$DATA_PATH = 'data/'
settings$PACKAGES_PATH = 'packages/'

# Install packages
packages.name <- list.files(settings$PACKAGES_PATH )
packages.path <- paste0(settings$PACKAGES_PATH,packages.name)
#install.packages(packages.path, dependencies=T)

# add caret
library(mlr)
library(extraTrees)
# Load packages
packages.name <- strsplit(packages.name,"_")
packages.name <- c(do.call(rbind,packages.name)[,1])
lapply(packages.name, require, character.only = TRUE)
library(parallel)

# Register parallel background
if(Sys.info()['sysname'] == "Windows") {
  library(doParallel)
  registerDoParallel(cores = detectCores())
} else {
  library(doMC)
  registerDoMC(cores = detectCores())
}

# Load functions
source('functions.R')

all.curr <- downloadData('20:00','10 Y')
returns <- NULL
ratios <- NULL
if(ALL.TARGET==TRUE) TARGETS = names(all.curr)

for(TARGET in TARGETS) {
    # Get the predictor value 
    variables <- lag(diff(log(all.curr[,TARGET])),-1)
    
    # Get the lags
    greps <- paste(strsplit(TARGET,'[.]')[[1]],collapse = '|')
    variables <- cbind(variables,getLags(all.curr[,grep(greps,colnames(all.curr))],4))
    
    # Get arima garch and arima Monte Carlo
    garch <- arima.garch(all.curr)
    mc <- arima.mc(all.curr)
    
    if(NEW != '#new') {
      garch <- format.garch(garch)
      mc <- format.mc(mc)
    }
    
    # Merge variables
    lookback = 500
    variables <- data.frame(variables[(lookback+1):nrow(variables),],garch,mc)
    variables <- selectVariables(variables,greps)
    
    # Get past Returns
    pastReturns <- getPastReturns(variables[,-1])
    
    # Merge variables
    variables <-cbind(variables,pastReturns)
    garch5 <-garch[,grep('argarch5',colnames(garch))]
    garch5 <-garch5[,grep(greps,colnames(garch5))]
    variables <- cbind(variables,garch5)
    
    # Preprocessing
    variables[is.na(variables[,1]),1] = 0
    variables <- as.data.frame(variables)
    
    # Divide into training, validation and testing set
    dividing(variables,cutoff=0.81,validation=F)
    
    # Train Random Forest
    train <- train.rf(train.set,1)
    extra <- train.extra(train.set)
    
    # Rank Important Features 
    important.Variables <- importantVariables(train,ncol(train.set)-1)
    important.Variables.extra <- extra
    
    # Retrain Random Forest with optimal parameters
    train <- train.rf(train.set[,c(TARGET,important.Variables[1:30])],10)
    extra <- extraTrees(train.set[,important.Variables.extra],train.set[,TARGET])
    
    # Compute in sample predictions for Cutting off risk
    in.sample.preds <- computeInSamplePreds(train.set[,important.Variables[1:30]],train,10)
    in.sample.preds.extra <- predict(extra,train.set[,important.Variables.extra])
    
    # Test Models
    preds.rf <- test.rf(test.set,train,10)
    preds.extra <- predict(extra,test.set[,important.Variables.extra])
    
    # Cut off risk points beyond 45% or 55% accuracy compared to in sample predictions, make a trade
    preds.rf <- cutOffRisk(in.sample.preds,preds.rf,10,12)
    preds.extra <- cutOffRisk(in.sample.preds.extra,preds.extra,10,12)
    preds = ifelse(preds.rf > 0 & preds.extra > 0,1,
                   ifelse(preds.rf < 0 & preds.extra < 0,-1,0) )
    
    # Compute Performance
    returns<- cbind(returns,xts(test.set[,1]*preds,order.by=as.POSIXct(row.names(variables)[(nrow(variables)-length(preds)+1):(nrow(variables))])))
    ratios <- rbind(ratios,performance.xts(test.set,preds))
    print(ratios)
  }
ratiostmp <- cbind(TARGETS,ratios)
ratiostmp
total.returns = apply(returns,1,sum)
# Plot results
charts.PerformanceSummary(total.returns)
SharpeRatio.annualized(total.returns)
CalmarRatio(total.returns)
Return.annualized(total.returns)
