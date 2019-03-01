par(family = 'serif')
setwd("/Users/mikhailgaerlan/Box Sync/Education/UC Davis/2016-2017 Spring/STA 243 Computational Statistics/Assignments/Assignment 3")
rm(list=ls())

n = 5000
u = runif(n)
x = -log(1-(1-exp(-2))*u)

h = hist(x,plot=F)
hcum = h
hcum$counts=cumsum(hcum$counts)
plot(hcum)
plot(h,add=T,col='grey')

d = density(x)
lines(d$x,d$y*length(x)*diff(h$breaks)[1],lwd=2)
lines(d$x,cumsum(d$y)/max(cumsum(d$y))*length(x),lwd=2)