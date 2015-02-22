fun <- function(symbol = "CAT", daysLeft=90, rfr=0.02, callStrike=100, putStrike=100, strategy="straddle", trials=100)
{
  #2 packages that are required to run this
  require(quantmod)
  require(ggplot2)
  
  #pull in stock prices and calculate returns, mean, sd, annual sd
  stock <- getSymbols(symbol,auto.assign=FALSE,warnings=FALSE,verbose=FALSE)
  sPrice <- st[,6]
  obs <- nrow(sPrice)
  adjPrice <- sPrice[(obs-252):obs,]
  returns <- na.omit(dailyReturn(adjPrice,leading=FALSE))
  meanRet <- mean(returns,na.omit=TRUE)
  volRet <- sd(returns)
  annualVol <- volRet*sqrt(252)
  lastPrice <- as.numeric(adjPrice[nrow(adjPrice)])
  
  #tSteps is days to maturity
  TTM <- daysLeft/252
  
  #random walk future stock price estimation
  priceEstimates <- as.numeric(lastPrice)*(exp((meanRet-0.5*annualVol^2)*TTM+annualVol*sqrt(TTM)*rnorm(trials,0,1))) #estimated prices at option expiration
  
  #pull in dividends and calculate dividend yield (DY)
  divs <- getDividends(symbol,auto.assign=FALSE,warnings=FALSE,verbose=FALSE)
  lastDiv <- as.numeric(divs[nrow(divs),1])
  if(is.null(lastDiv) || is.na(lastDiv))
  {
    DY <- 0
  }
  else
  {
    annualDiv <- 4*lastDiv
    DY <- annualDiv/lastPrice
  }
  
  #DY <- 0
  
  #use Black Scholes Merton to price calls and puts
  #ex: sapply(c(100,105,110),FUN=BSM,s0=90,vol=0.25,rfr=0.02,tiy=0.5,DY=0)
  
  if(strategy=="straddle")
  {
    callPrice <- BSM(callStrike,s0=lastPrice,vol=annualVol,rfr=rfr,tiy=TTM,DY=DY)$callPrice
    putPrice <- BSM(putStrike,s0=lastPrice,vol=annualVol,rfr=rfr,tiy=TTM,DY=DY)$putPrice
    longCall <- pmax(rep(0,trials),priceEstimates-callStrike) -callPrice
    longPut <- pmax(rep(0,trials),putStrike-priceEstimates) -putPrice
    longPayoff <- longCall + longPut
    shortPayoff <- -longPayoff
  }
  summary(cbind(longPayoff,shortPayoff))
}

BSM <- function(strike,s0,vol,rfr,tiy,DY)
{
  #Merton 
  d1M <- (log(s0/strike)+(rfr-DY+0.5*(vol*vol))*tiy)/(vol*sqrt(tiy))
  d2M <- d1M-(vol*sqrt(tiy))
  CallM <- exp(-DY*tiy)*s0*pnorm(d1M)-strike*exp(-rfr*tiy)*pnorm(d2M)
  PutM <- -exp(-DY*tiy)*s0*pnorm(-d1M)+strike*exp(-rfr*tiy)*pnorm(-d2M)
  list(callPrice = CallM, putPrice = PutM)
}