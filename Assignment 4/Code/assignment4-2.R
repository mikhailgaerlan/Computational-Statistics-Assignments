par(family = 'serif')
setwd("/Users/mikhailgaerlan/Box Sync/Education/UC Davis/2016-2017 Spring/STA 243 Computational Statistics/Assignments/Assignment 4")
rm(list=ls())

nu = c(0.1,1,10)

for(i in nu){
  m = 10000
  x = rnorm(m,1.5,i)
  f = 0*x
  for(j in 1:m){
    if(x[j]>=1&&x[j]<=2){
      f[j] = 1
    }
  }
  g = dnorm(x,1.5,i)
  h = dnorm(x,0,1)
  w = f/g
  mu = mean(h*w)
  var = sum((h*w-mu)^2)/(m-1)
  hist(h*w)
  print(c(mu,var))
}