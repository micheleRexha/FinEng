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
  if(b<0) return(1000)
  if(abs(p)>=1) return(1000)
  if(!is.numeric(m)) return(1000) #TODO: check earlier?
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

#Then write the optimization function applied to vfunc

##typical criteria would be 
itermax <- 1e4
VTR <- 1e-4
l <- c(-20,0,-0.99,-20 ,1e-5)
u <- c(20,50,0.99,5,20)

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


#PS7 solution
BSvega<-function(pc, S, k, vol, d, r, t)
{
  #pc  put/call indicator call=1, put=-1
  #S   Stock price at 0
  #K   strike
  #vol volatility
  #d   dividend yield
  #r   riskless rate
  #t   time to maturity
  
  d1 = (log(S / k) + t * (r - d + (vol ^ 2) / 2)) / (vol * sqrt(t))
  
  BSvega = exp(-d * t) * S * sqrt(t) * exp((-d1 ^ 2) / 2) / (sqrt(2 * pi))
  return(BSvega)
}


#PS7 solution
BSvol<-function(pc, S, k, price, d, r, t, start = 0.2)
{
  #pc    put/call indicator call=1, put=-1
  #S     Stock price at 0
  #K     strike
  #price option premium
  #d     dividend yield
  #r     riskless rate
  #t     time to maturity
  #start starting value for vol, optional, by default=0.2
  
  voli = start
  pricei = BSprice(pc, S, k, voli, d, r, t)
  vegai = BSvega(pc, S, k, voli, d, r, t)
  i = 0

  while(abs(price - pricei) > 0.000001) 
  {
    i =+ 1
    voli<-voli + (price - pricei) / vegai
    
    #TODO: catch vol divergence
    #if(voli > 10^3) return(NA)
    #if(i > 10^4) return(NA)
    print(voli)
    
    pricei<-BSprice(pc, S, k, voli, d, r, t)
    vegai<-BSvega(pc, S, k, voli, d, r, t)
  }
  
  BSvol = voli
  return(BSvol)
}

load(file = "fm408_exam_data.RData")

#allocate vectors
volvec<-vector()
ttmvec<-vector()
moneynessvec<-vector()
length(volvec) = nrow(raw)
length(ttmvec) = nrow(raw)
length(moneynessvec) = nrow(raw)


#iterate over calls, calculate implied vols
for (i in 1:3360) {
  #TODO: quick fix: discard calls that are far ITM
  #if(raw[[7]][i] <= 20) {
  #  voli<-BSvol(1, raw[[6]][i], raw[[3]][i], raw[[7]][i], raw[[5]][i], raw[[4]][i], raw[[2]][i]/12)
  #} else voli<-NA
  voli<-BSvol(1, raw[[6]][i]*exp(-(raw[[4]][i]-raw[[5]][i])*(raw[[2]][i]/12)), raw[[3]][i], raw[[7]][i], raw[[5]][i], raw[[4]][i], raw[[2]][i]/12)
  
  ttm<-raw[[2]][i]/12
  moneyness<-raw[[6]][i]/raw[[3]][i]
  
  #assign results to vectors
  volvec[i]<-voli
  ttmvec[i]<-ttm
  moneynessvec[i]<-moneyness
  
  print(voli)
  flush.console()
}

plotvalue = data.frame(ttmvec, moneynessvec, volvec)  

#install.packages("plotly")
library("plotly")


volvec<-vector()
ttmvec<-vector()
moneynessvec<-vector()
length(volvec) = nrow(raw)
length(ttmvec) = nrow(raw)
length(moneynessvec) = nrow(raw)

for (i in 1:nrow(raw)) {
  f<-function (x) BSprice(1, raw[[6]][i]*exp(-(raw[[4]][i]-raw[[5]][i])*(raw[[2]][i]/12)), raw[[3]][i], x, raw[[5]][i], raw[[4]][i], raw[[2]][i]/12) - raw[[7]][i]
  
  voli<-uniroot(f, lower=0, upper=10)[[1]]
  ttm<-raw[[2]][i]/12
  moneyness<-raw[[6]][i]/raw[[3]][i]
  
  #assign results to vectors
  volvec[i]<-voli
  ttmvec[i]<-ttm
  moneynessvec[i]<-moneyness
}

for (i in 1:1) {
  f<-function (x) BSprice(1, raw[[6]][i]*exp(-(raw[[4]][i]-raw[[5]][i])*(raw[[2]][i]/12)), raw[[3]][i], x, raw[[5]][i], raw[[4]][i], raw[[2]][i]/12) - raw[[7]][i]
  #f<-function (x) BSprice(1, raw[[6]][i], raw[[3]][i], x, raw[[5]][i], raw[[4]][i], raw[[2]][i]/12) - raw[[7]][i]
  
  voli<-uniroot(f, lower=0, upper=10)[[1]][1]
  ttm<-raw[[2]][i]/12
  moneyness<-raw[[6]][i]/raw[[3]][i]
  print(voli)
  
  #assign results to vectors
  volvec[i]<-voli
  ttmvec[i]<-ttm
  moneynessvec[i]<-moneyness
}

plotvalue = data.frame(ttmvec, moneynessvec, volvec)  





BSprice(1, raw[[6]][1]*exp(-(raw[[4]][1]-raw[[5]][1])*(raw[[2]][1]/12)), raw[[3]][1], 1.104639, raw[[5]][1], raw[[4]][1], raw[[2]][1]/12)

#BSvol(1, raw[[6]][2], raw[[3]][2], raw[[7]][2], raw[[5]][2]/12, raw[[4]][2]/12, raw[[2]][2]/12)
#raw[[7]][1]
#BSprice(1, raw[[6]][1]*exp(-(raw[[4]][1]-raw[[5]][1])*(raw[[2]][1]/12)), raw[[3]][1], 1.1, raw[[5]][1], raw[[4]][1], raw[[2]][1]/12)
#BSprice(1, raw[[6]][3], raw[[3]][3], 0.1000, raw[[5]][3], raw[[4]][3], raw[[2]][3]/12)
#BSprice(1, raw[[6]][3], raw[[3]][3], 0.05000, raw[[5]][3], raw[[4]][3], raw[[2]][3]/12)
#BSprice(1, raw[[6]][3], raw[[3]][3], 0.0000, raw[[5]][3], raw[[4]][3], raw[[2]][3]/12)
#raw[[6]][1]*exp(-(raw[[4]][1]-raw[[5]][1])*(raw[[2]][1]/12))
#BSvol(1, raw[[6]][1]*exp(-(raw[[4]][1]-raw[[5]][1])*(raw[[2]][1]/12)), raw[[3]][1], raw[[7]][1], raw[[5]][1], raw[[4]][1], raw[[2]][1]/12)