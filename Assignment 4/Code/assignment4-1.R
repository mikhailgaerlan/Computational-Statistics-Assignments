par(family = 'serif')
setwd("/Users/mikhailgaerlan/Box Sync/Education/UC Davis/2016-2017 Spring/STA 243 Computational Statistics/Assignments/Assignment 4")
rm(list=ls())

n = 1000
x = runif(n)
h = x^2
mu = mean(h)
var = sum((h-mu)^2)/(n-1)
print(c(mu,var))

n = 1000
x = cbind(runif(n,-2,2),runif(n,0,1))
h = 4*x[,1]^2*cos(x[,1]*x[,2])
mu = mean(h)
var = sum((h-mu)^2)/(n-1)
print(c(mu,var))

n = 1000
x = rexp(n,1/4)
h = 3*x^4*exp(-x^3/4+x/4)
mu = mean(h)
var = sum((h-mu)^2)/(n-1)
print(c(mu,var))