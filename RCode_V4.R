#PS7 solution
BSprice<-function(pc, S, k, vol, d, r, t)
{
  #pc  put/call indicator call=1, put=-1
  #S   Stock price at 0
  #K   strike
  #vol volatility
  #d   dividend yield
  #r   riskless rate
  #t   time to maturity
  
  
  d1 = (log(S / k) + t * (r - d + (vol ^ 2) / 2)) / (vol * sqrt(t))
  d2 = d1 - vol * sqrt(t)
  
  BSprice = pc * exp(-d * t) * S * 
    pnorm(pc * d1) - pc * k * exp(-r * t) * pnorm(pc * d2)
  return(BSprice)
}

load(file = "fm408_exam_data.RData")

volvec<-vector()
ttmvec<-vector()
moneynessvec<-vector()
length(volvec) = nrow(raw)
length(ttmvec) = nrow(raw)
length(moneynessvec) = nrow(raw)

for (i in 1:nrow(raw)) {
  f<-function (x) BSprice(1, raw$Forward[i]*exp(-(raw$InterestRate[i]-raw$DividendYield[i])*(raw$Term[i]/12)), raw$Strike[i], x, raw$DividendYield[i], raw$InterestRate[i], raw$Term[i]/12) - raw$CallPrice[i]
  
  voli<-uniroot(f, lower=0, upper=3)[[1]]
  ttm<-raw$Term[i]/12
  moneyness<-log(raw$Strike[i]/raw$Forward[i])
  
  #assign results to vectors
  volvec[i]<-voli
  ttmvec[i]<-ttm
  moneynessvec[i]<-moneyness
}

raw["ImpliedVol"]<-volvec
raw["TTM"]<-ttmvec
raw["Moneyness"]<-moneynessvec

#install.packages("plotly")
library("plotly")

plotdata <- raw[raw[1] == "2006-01-31", ]

p<-plot_ly(plotdata, x=plotdata$TTM, y=plotdata$Moneyness, z=plotdata$ImpliedVol) %>%
           add_markers() %>%
           layout(scene = list(xaxis = list(title = 'Time to Maturity'),
                               yaxis = list(title = 'Moneyness [log(K/F)]'),
                               zaxis = list(title = 'Implied Volatility')))

p

Sys.setenv("plotly_username"="tneuber")
Sys.setenv("plotly_api_key"="7CknaAatVziORktIj116")
#chart_link = plotly_POST(p, filename="ImpliedVol")
#chart_link

## p: vector of svi parameters in the following order
##    (a, b, rho, m, sigma)
## k: log moneyness (log(K/F))
## maturity: options' time to maturity in years
## v: observed implied volatility
vfunc=function(p,k, maturity, v){
  a=p[1]
  b=p[2]
  rho=p[3]
  m=p[4]
  sigma2=p[5]^2
  
  kk=(k-m)^2+sigma2
  totalIV=a+b*(rho * (k-m)+sqrt(kk))
  
  if(!is.numeric(totalIV)) return(1000)
  if(min(totalIV)<0) return(1000)
  
  ## students should add restrictions here
  # e.g. if(some conditions give TRUE/FALSE as output) return(1000)
  # restrictions from task description
  if(!is.numeric(a)) return(1000)
  if(b<0) return(1000)
  if(abs(rho)>=1) return(1000)
  if(!is.numeric(m)) return(1000) 
  if(p[5]<=0) return(1000) #p[5] = sigma
  if((a + b*p[5]*sqrt(1-rho^2)) < 0) return(1000)
  
  
  # sum of squares
  res=sum((totalIV-maturity*((v)^2))^2)
  
  return(res)
}

## parameters for optimisation function DEoptim
## l: vector of lower bound of parameters
## u: vector of upper bound of parameters
## itermax: number of iterations
## VTR: minimum sum of squares to be reached

library("DEoptim")

l <- c(-20,0,-0.99,-20 ,1e-5)
u <- c(20,50,0.99,5,20)


meanvec<-vector()
#annualised
varvec<-vector()
termvec<-vector()
#c(1, 3, 6, 12, 24, 36, 48, 60, 84, 120)
for(term in c(1, 3, 6, 12, 24, 36, 48, 60, 84, 120)) {
  
  #TODO: make selection variable baseds
  optionsToFit<-raw[raw$ValuationDate == "2006-01-31", ]
  optionsToFit<-optionsToFit[optionsToFit$Term == term, ]
  
  #TODO: no extrapolation?!
  lIntegrBound <- optionsToFit$Moneyness[1]
  uIntegrBound <- optionsToFit$Moneyness[nrow(optionsToFit)]
  
  #initialize negative
  var<--1
  
  while(var < 0) {
    func = function(p) {
      sum = 0
      for (i in 1:nrow(optionsToFit)) {
        sum =+ vfunc(p, optionsToFit$Moneyness[i], optionsToFit$TTM[i], optionsToFit$ImpliedVol[i])
      }
      return(sum)
    }
    
    #Then write the optimization function applied to vfunc
    outDEoptim <- DEoptim(func, l, u, DEoptim.control(VTR = 1e-4, itermax =  1e4))
    summary(outDEoptim)
    
    ## p: vector of svi parameters in the following order
    ##    (a, b, rho, m, sigma)
    ## k: log moneyness (log(K/F))
    ## maturity: options' time to maturity in years
    implVol = function(p, k, maturity) {
      a=p[1]
      b=p[2]
      rho=p[3]
      m=p[4]
      sigma2=p[5]^2
      
      kk=(k-m)^2+sigma2
      totalIV=a+b*(rho * (k-m)+sqrt(kk))
      return(sqrt(totalIV/maturity))
    }
    
    bestFit = function(k) {implVol(outDEoptim$optim$bestmem, k, optionsToFit$TTM[1])}
    
    curve(bestFit, from=-1, to=3, xlab="Moneyness [log(K/F)]", ylab="Implied Vol")
    
    secondDerivative = function(f, x, delta) {
      res = (f(x - delta) - 2*f(x) + f(x + delta)) / (delta^2)
      return (res)
    }
    
    #TODO: Check second argument (discounting?!)
    C = function(K) {BSprice(1, optionsToFit$Forward[1]*exp(-(optionsToFit$InterestRate[1] - optionsToFit$DividendYield[1])*optionsToFit$TTM[1]), K, bestFit(log(K/optionsToFit$Forward[1])), optionsToFit$DividendYield[1], optionsToFit$InterestRate[1], optionsToFit$TTM[1])}
    
    implDensity = function(K) {
      Dt = exp(-optionsToFit$TTM[1]*(optionsToFit$InterestRate[1] - optionsToFit$DividendYield[1]))
      secDeriv = secondDerivative(C, K, 0.5)
      return(secDeriv / Dt)
    }
    
    curve(implDensity, from=500, to=3000, xlab="K", ylab="Implied Density")
    
    returnImplDensityUnscaled = function(r) {return(implDensity(exp(r)*optionsToFit$Forward[1]))} 
    returnImplDensity = function(r) {return(returnImplDensityUnscaled(r)/integrate(returnImplDensityUnscaled, lower=lIntegrBound,upper=uIntegrBound)$value)}
    
    curve(returnImplDensity, from=-1, to=1, xlab="R", ylab="Implied Density")
    
    mean<-integrate(function(x) {return(x*returnImplDensity(x))}, lIntegrBound, uIntegrBound)$value
    var<-(integrate(function(x) {return((x-mean)^2*returnImplDensity(x))}, lIntegrBound, uIntegrBound)$value)/(term/12)
  }
  
  meanvec<-c(meanvec, mean)
  varvec<-c(varvec, var)
  termvec<-c(termvec, term)
  
}

varvec  

integrate(returnImplDensityUnscaled, -3, 3)$value/integrate(returnImplDensityUnscaled, lower=-5,upper=Inf)$value
integrate(returnImplDensity, -5, 5)$value

implDensity(150)

BSprice(1, raw[[6]][option]*exp(-(raw[[4]][option]-raw[[5]][option])*(raw[[2]][option]/12)), raw$Strike[option], raw$ImpliedVol[option], raw$DividendYield[option], raw$InterestRate[option], raw$TTM[option])

secondDerivative(C, 2000, 0.10)
curve(C, from=1900, to=2100, xlab="K", ylab="C")
C(1670.3135)
bestFit(1500)
bestFit(2000)
bestFit(2500)
raw$Strike[option]
C(raw$Strike[option])
C(raw[[6]][option]*exp(-(raw[[4]][option]-raw[[5]][option])*(raw[[2]][option]/12)))
bestFit(raw$Strike[option])
bestFit(log(raw$Strike[option]/(raw[[6]][option]*exp(-(raw[[4]][option]-raw[[5]][option])*(raw[[2]][option]/12)))))