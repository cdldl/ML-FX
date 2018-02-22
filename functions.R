arima.garch <- function(data) {
  
  if(file.exists(paste0(settings$DATA_PATH,NEW,'results arimagarch units.rds'))) {
    return(readRDS(paste0(settings$DATA_PATH,NEW,'results arimagarch units.rds')))
  }
  
  ret <- logret(data)
  grid.param <- expand.grid(names(data),500:nrow(ret)) 
  
  
  fe <- foreach(param = iter(grid.param, by = "row"), 
                .verbose = TRUE, .errorhandling = "pass",  
                .multicombine = TRUE, .maxcombine = max(2, nrow(grid.param)),
                .export=c("selectArima","xts","arima","predict","arpredict",'ugarchspec'))
  fe$args <- fe$args[1]
  fe$argnames <- fe$argnames[1]
  
  results <- fe %dopar% {
    #param = grid.param[2094,]
    data.tmp <- ret[,grep(param[,1],names(ret))]  
    
    preds <- c()
    T=param[,2]
    tmp.ret = data.tmp[(T-499):T]
    ar = selectArima(coredata(tmp.ret))
    
    spec = ugarchspec(
      variance.model=list(garchOrder=c(1,1)),
      mean.model=list(armaOrder=c(ar[1], ar[2]), include.mean=T),
      distribution.model="std"
    )
    fit = tryCatch(
      ugarchfit(
        spec, tmp.ret, solver = 'hybrid'
      ), error=function(e) e, warning=function(w) w
    )
    
    if(is(fit, "warning")) {
      preds <- c(0,0)
    } else {
      fore = ugarchforecast(fit, n.ahead=5)
      ind = fore@forecast$seriesFor
      preds = ind[1]
      ind = mean(ind,na.rm=T)
      preds = c(preds,ind)
    } 
    
    
    return(preds)
  }
  names(results) <-  apply(grid.param,1,paste,collapse=",")
  saveRDS(results,paste0(settings$DATA_PATH,NEW,'results arimagarch units.rds'))
  results
}

arima.mc <- function(data) {
  if(file.exists(paste0(settings$DATA_PATH,NEW,'results arimaMcFx.rds'))) {
    return(readRDS(paste0(settings$DATA_PATH,NEW,'results arimaMcFx.rds')))
  }
  ret <- logret(data)
  grid.param <- expand.grid(names(data))
  fe <- foreach(param = iter(grid.param, by = "row"), 
                .verbose = TRUE, .errorhandling = "pass",  
                .multicombine = TRUE, .maxcombine = max(2, nrow(grid.param)),
                .export=c("selectArima","xts","arima","simulate","monteCarlo","weekdays"))
  fe$args <- fe$args[1]
  fe$argnames <- fe$argnames[1]
  
  results <- fe %dopar% {
    #param = grid.param[1,]
    
    data.tmp <- ret[,grep(param,names(ret))]  
    
    #T=param[,2]
    beg=500
    preds <- numeric(length(beg:nrow(ret)))
    for(T in beg:nrow(ret)) {
      #T=1  
      tmp.ret = data.tmp[(T-499):T]
      p<-monteCarlo(coredata(tmp.ret))
      preds[T-499] = p
    }
    return(preds)
  }
  names(results) <-  apply(grid.param,1,paste,collapse=",")
  saveRDS(results,paste0(settings$DATA_PATH,NEW,'results arimaMcFx.rds'))
  results
}

# Compute insample preds
computeInSamplePreds <- function(data,rfs,number=1) {
  in.sample.preds <- NULL
  
  for(i in 1:number) in.sample.preds =cbind(in.sample.preds,predict(rfs[[i]],as.data.frame(data))$predictions)
  in.sample.preds <- apply(in.sample.preds,1,mean)
  in.sample.preds
}
# Cut off probabilities based on past proba - Triggering buy and sell orders
cutOffRisk <- function(in.sample.preds,preds,quant1,quant2) {
  preds <- ifelse(preds>quantile(in.sample.preds,seq(0,1,0.05))[quant2]
                  ,1,ifelse(preds<quantile(in.sample.preds,seq(0,1,0.05))[quant1]
                            ,-1,0))
  preds
}
# Create training, validation and testing dataset
dividing <- function(variables,cutoff=0.8,validation=T) {
  train.set.index <- 1:round(nrow(variables)*cutoff)
  
  if(validation == TRUE) {
    train.set.index <- 1:round(nrow(variables)*(cutoff-0.1))
    val.set.index <- (round(nrow(variables)*(cutoff-0.1)+1)):round(nrow(variables)*cutoff)
    val.set <<- variables[val.set.index,]
  }
  
  test.set.index <- (round(nrow(variables)*cutoff+1)):nrow(variables)
  train.set <<- variables[train.set.index,]
  test.set <<- variables[test.set.index,]
}
downloadData <- function(time,dur) {
  if(file.exists(paste0(settings$DATA_PATH,NEW,'all curr daily tws 8pm.rds'))) {
    return(readRDS(paste0(settings$DATA_PATH,NEW,'all curr daily tws 8pm.rds')))
  } 
  
  tws <- twsConnect()
  ccy1 <- reqContractDetails(tws, twsCurrency("USD", "JPY"))[[1]]$contract
  ccy2 <- reqContractDetails(tws, twsCurrency("USD", "CAD"))[[1]]$contract
  ccy3 <- reqContractDetails(tws, twsCurrency("USD", "CHF"))[[1]]$contract
  ccy4 <- reqContractDetails(tws, twsCurrency("EUR", "USD"))[[1]]$contract
  ccy5 <- reqContractDetails(tws, twsCurrency("GBP", "USD"))[[1]]$contract
  ccy6 <- reqContractDetails(tws, twsCurrency("AUD", "USD"))[[1]]$contract
  
  x1 <- reqHistoricalData(tws, ccy1,barSize = "1 hour",whatToShow='MIDPOINT',duration=dur)
  x2 <- reqHistoricalData(tws, ccy2,barSize = "1 hour",whatToShow='MIDPOINT',duration=dur)
  x3 <- reqHistoricalData(tws, ccy3,barSize = "1 hour",whatToShow='MIDPOINT',duration=dur)
  x4 <- reqHistoricalData(tws, ccy4,barSize = "1 hour",whatToShow='MIDPOINT',duration=dur)
  x5 <- reqHistoricalData(tws, ccy5,barSize = "1 hour",whatToShow='MIDPOINT',duration=dur)
  x6 <- reqHistoricalData(tws, ccy6,barSize = "1 hour",whatToShow='MIDPOINT',duration=dur)
  
  
  x7 <- x1[which(format(index(x1),"%H:%M")==time),]
  x8 <- x2[which(format(index(x2),"%H:%M")==time),]
  x9 <- x3[which(format(index(x3),"%H:%M")==time),]
  x10 <- x4[which(format(index(x4),"%H:%M")==time),]
  x11 <- x5[which(format(index(x5),"%H:%M")==time),]
  x12 <- x6[which(format(index(x6),"%H:%M")==time),]
  all <- cbind(x7[,4],x8[,4],x9[,4],x10[,4],x11[,4],x12[,4])
  
  names(all) <- c("USD.JPY","USD.CAD","USD.CHF","EUR.USD","GBP.USD","AUD.USD")
  tmp <- all
  
  tmp$EUR.GBP <- tmp$EUR.USD * (1/tmp$GBP.USD)
  
  tmp$EUR.AUD <- tmp$EUR.USD * (1/tmp$AUD.USD)
  tmp$GBP.AUD <- tmp$GBP.USD * (1/tmp$AUD.USD)
  
  tmp$EUR.JPY <- tmp$EUR.USD * tmp$USD.JPY
  tmp$GBP.JPY <- tmp$GBP.USD * tmp$USD.JPY
  tmp$AUD.JPY <- tmp$AUD.USD * tmp$USD.JPY
  
  tmp$EUR.CHF <- tmp$EUR.USD * tmp$USD.CHF
  tmp$GBP.CHF <- tmp$GBP.USD * tmp$USD.CHF
  tmp$AUD.CHF <- tmp$AUD.USD * tmp$USD.CHF
  tmp$CHF.JPY <- tmp$USD.JPY * (1/tmp$USD.CHF)
  
  tmp$EUR.CAD <- tmp$EUR.USD * tmp$USD.CAD
  tmp$GBP.CAD <- tmp$GBP.USD * tmp$USD.CAD
  tmp$AUD.CAD <- tmp$AUD.USD * tmp$USD.CAD
  tmp$CAD.CHF <- (1/tmp$USD.CAD) * tmp$USD.CHF
  tmp$CAD.JPY <- (1/tmp$USD.CAD) * tmp$USD.JPY
  dir = readRDS(paste0(settings$DATA_PATH,NEW,'all curr daily tws 8pm.rds'))
  tmp <- tmp[,colnames(dir)]
  saveRDS(tmp,paste0(settings$DATA_PATH,NEW,'all curr daily tws 8pm.rds'))
  tmp
}
format.garch <- function(garch) {
  data <- garch
  curr <- unique(substr(names(data),1,7))
  for(i in 1:length(data)) {
    if(length(grep('error',class(data[[i]])))!=0)   data[[i]] = c(0,0)
  }
  data <- do.call(rbind,data)
  
  all <- list()
  for(i in curr) all[[i]] = data[grep(i,row.names(data)),]
  all <- do.call(cbind,all)
  name <- c()
  for(i in 1:length(curr)) name = c(name,paste0('argarch1',curr[i]),paste0('argarch5',curr[i]))
  colnames(all) = name
  all <- as.data.frame(all)
  all
}
format.mc <- function(data) {
  do.call(cbind,data)
}
getLags <- function(var,number.of.lags) {
  test.ret <- list()
  ind =1 
  for(j in 1:ncol(var)) for(k in 0:number.of.lags)  {
    test.ret[[ind]] <- lag(diff(log(var[,j])),k) #diff(log(
    names(test.ret[[ind]]) = paste0(colnames(var[,j]),k) 
    ind = ind+1
  }
  for(j in 1:length(test.ret)) test.ret[[j]] <- coredata(test.ret[[j]])
  test.ret <- do.call(cbind,test.ret)
  test.ret[is.na(test.ret)] <- 0
  test.ret
}
# Get past returns variable
getPastReturns <- function(all) {
  sum3 <- rollapply(all[,which(2:57 %% 5   == 0)-3],3,sum,align='right',fill=0)
  colnames(sum3) = paste0(colnames(sum3),'sum3')
  sum5<- rollapply(all[,which(2:57 %% 5   == 0)-3],5,sum,align='right',fill=0)
  colnames(sum5) = paste0(colnames(sum5),'sum5')
  sum15<- rollapply(all[,which(2:57 %% 5   == 0)-3],15,sum,align='right',fill=0)
  colnames(sum15) = paste0(colnames(sum15),'sum15')
  sum30<- rollapply(all[,which(2:57 %% 5   == 0)-3],30,sum,align='right',fill=0)
  colnames(sum30) = paste0(colnames(sum30),'sum30')
  pastReturns <- cbind(sum3,sum5,sum15,sum30)
  pastReturns
}
# Get importance of variables
importantVariables <- function(train,number) {
  imps <- sapply(train,function(rf) importance(rf))
  imp <- apply(imps,1,sum)
  names.imp <-names(sort(imp)[1:number])  
  names.imp
}
monteCarlo <- function(ts.ts, N=1000) {
  pq<-selectArima(ts.ts);
  p=pq[1];
  q=pq[2];
  fit<-arima(ts.ts, c(p, 0, q));
  ts.sim<-ts.ts;
  n=length(ts.sim)
  
  # We simulate different trajectories
  
  S <- vector(1000,mode='numeric');
  h <- 1;
  #doPar: parallelized it 
  for(i in 1:N){
    S[i] <- as.vector(simulate(fit,h,future=TRUE,bootstrap = TRUE, innov=NULL));
  }
  # Double check this statement
  S[is.na(S)] = 0;
  
  p <- sum(S>0)/N;
  return(p);
}
# Performance Analytics
performance <- function(test,preds) {
  # Get Stats
  returns<- test[,1]*ifelse(preds==0,0,ifelse(preds>0,1,-1))
  # Area under ROC
  auc <- auc(ifelse(test[which(preds!=0),1]>0,1,-1),ifelse(preds[which(preds!=0)]>0,1,-1))
  #rmse <- rmse(test.all[which(preds!=0),1],preds[which(preds!=0)])
  
  #sharpe <- SharpeRatio.annualized(returns) 
  sharpe <- (mean(returns)*sqrt(252)) / sd(returns)
  names(sharpe) = 'sharpe'
  
  # Maximum drawdown ratio
  
  calmar <- maxDrawdown(returns)
  names(calmar) = 'drawdown'
  
  yearly.return <- sum(returns)
  names(yearly.return) = 'yearly.return'
  
  win.rate <- length(which(returns[which(preds!=0)]>0))/length(returns[which(preds!=0)])
  
  # Number of trades over the period
  n.trades <- length(returns[which(preds!=0)])
  
  # Pips per trade
  bps.trade = mean(returns[which(preds!=0)])*10000
  
  ratios = data.frame(auc,sharpe=sharpe,maxDrawdown=calmar,return=yearly.return,ntrades=n.trades,win.rate,bps.trade) #rmse
  row.names(ratios) = NULL
  ratios
} 
performance.xts <- function(test,preds) {
  # Get Stats
  returns<- xts(test[,1]*ifelse(preds==0,0,ifelse(preds>0,1,-1)),order.by=as.POSIXct(row.names(test)))
  # Area under ROC
  auc <- auc(ifelse(test[which(preds!=0),1]>0,1,-1),ifelse(preds[which(preds!=0)]>0,1,-1))
  #rmse <- rmse(test.all[which(preds!=0),1],preds[which(preds!=0)])
  
  sharpe <- SharpeRatio.annualized(returns) 
  colnames(sharpe) = 'sharpe'
  
  # Maximum drawdown ratio
  calmar <- CalmarRatio(returns)
  colnames(calmar) = 'calmar'
  
  yearly.return <- Return.annualized(returns)
  colnames(yearly.return) = 'yearly.return'
  
  win.rate <- length(which(returns[which(preds!=0)]>0))/length(returns[which(preds!=0)])
  
  # Number of trades over the period
  n.trades <- length(returns[which(preds!=0)])
  
  # Pips per trade
  bps.trade = mean(returns[which(preds!=0)])*10000
  
  ratios = data.frame(auc,sharpe=sharpe,calmar=calmar,yearly.return=yearly.return,ntrades=n.trades,win.rate,bps.trade) #rmse
  row.names(ratios) = NULL
  ratios
} 
selectArima <- function(timeseries.ts) {
  final.aic <- Inf
  final.order <- c(0,0,0)
  for (p in 0:5) for (q in 0:5) {
    if ( p == 0 && q == 0) {
      next
    }
    
    arimaFit = tryCatch( arima(timeseries.ts, order=c(p, 0, q)),
                         error=function( err ) FALSE,
                         warning=function( err ) FALSE )
    
    if( !is.logical( arimaFit ) ) {
      current.aic <- AIC(arimaFit)
      if (current.aic < final.aic) {
        final.aic <- current.aic
        final.order <- c(p, 0, q)
      }
    } else {
      next
    }
  }
  return(c(final.order[1],final.order[3]))
}
selectVariables <- function(all,greps) {
  all <- all[,-grep('argarch5',colnames(all))]
  all <- all[,grep(greps,colnames(all))]
  all
}
train.extra <- function(train.set) {
  
  regr.task = makeRegrTask(id = "x", data = train.set, target = TARGET)
  lrn = makeFilterWrapper(learner = "regr.extraTrees", fw.method = "rank.correlation",fw.perc =0.3)
  rdesc = makeResampleDesc("CV", iters = 5)
  r = resample(learner = lrn, task = regr.task, resampling = rdesc, show.info = TRUE, models = TRUE)
  sfeats = sapply(r$models, getFilteredFeatures)
  len = length(sfeats)
  min = min(45,len)
  names <- names(rev(sort(table(sfeats))))[1:min]
  names
}
train.rf <- function(data,number=1) {
  data <- as.data.frame(data)
  formula <- as.formula(as.formula(paste0(TARGET,'~',paste0(names(data[,-1]),collapse='+'))))
  system.time({
    rfs <- list()
  for(i in 1:number) {
    #data<-data[sample(nrow(data)),]
    rfs[[i]] <- ranger(formula,data=data,write.forest = TRUE,importance='impurity',min.node.size=4)
  }})
  rfs
}
test.rf <- function(data,rfs,number=1) {
  yo.reg <- NULL
  for(i in 1:number) yo.reg <- cbind(yo.reg,predict(rfs[[i]],as.data.frame(data))$predictions)
  preds <- apply(yo.reg,1,mean)
  preds
}


