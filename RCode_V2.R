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
  f<-function (x) BSprice(1, raw[[6]][i]*exp(-(raw[[4]][i]-raw[[5]][i])*(raw[[2]][i]/12)), raw[[3]][i], x, raw[[5]][i], raw[[4]][i], raw[[2]][i]/12) - raw[[7]][i]

  voli<-uniroot(f, lower=0, upper=10)[[1]]
  ttm<-raw[[2]][i]/12
  moneyness<-raw[[6]][i]/raw[[3]][i]
  
  #assign results to vectors
  volvec[i]<-voli
  ttmvec[i]<-ttm
  moneynessvec[i]<-moneyness
}

raw["ImpliedVol"]<-volvec
raw["TTM"]<-ttmvec
raw["Moneyness"]<-moneynessvec
plotvalue=data.frame(ttmvec, moneynessvec, volvec)

#install.packages("plotly")
library("plotly")

p<-plot_ly(raw, x=ttmvec, y=moneynessvec, z=volvec) %>%
           add_markers() %>%
           layout(scene = list(xaxis = list(title = 'Time to Maturity'),
                               yaxis = list(title = 'Moneyness [S/K]'),
                               zaxis = list(title = 'Implied Volatility')))

p
chart_link = plotly_POST(p, filename="ImpliedVol")
chart_link

p2<-plot_ly(plotvalue, x=plotvalue[[1]], y=plotvalue[[2]], z=plotvalue[[3]]) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'Time to Maturity'),
                      yaxis = list(title = 'Moneyness [S/K]'),
                      zaxis = list(title = 'Implied Volatility')))

p2

Sys.setenv("plotly_username"="tneuber")
Sys.setenv("plotly_api_key"="7CknaAatVziORktIj116")
plotly_POST(p2, filename = "ImpliedVol")

raw[raw[1] == "2006-01-31", ]
