par(family = 'serif')
setwd("/Users/mikhailgaerlan/Box Sync/Education/UC Davis/2016-2017 Spring/STA 243 Computational Statistics/Assignments/Assignment 3")
rm(list=ls())

alpha1 = 1
alpha2 = pi/2

q = function(x){
  return(exp(-x)/(1+x^2))
}
g1 = function(x){
  return(exp(-x))
}
g2 = function(x){
  return(2/(pi*(1+x^2)))
}
sampleg1 = function(n){
  x = 0*(1:n)
  for(i in 1:n){
    u = runif(1)
    x[i] = -log(1-u)
  }
  return(x)
}
sampleg2 = function(n){
  x = 0*(1:n)
  for(i in 1:n){
    u = runif(1)
    x[i] = tan(pi*u/2)
  }
  return(x)
}
samplef = function(n,o){
  x2 = 0*(1:n)
  j = 0
  for (i in 1:n){
    test = T
    while (test){
      if (o == 1){
        test2 = T
        while (test2){
          x = sampleg1(1)
          if (x < 5){
            test2 = F
          }
        }
        j = j + 1
        if (runif(1) < q(x)/(alpha1*g1(x))){
          x2[i] = x
          test = F
        }
      } else {
        test2 = T
        while (test2){
          x = sampleg2(1)
          if (x < 5){
            test2 = F
          }
        }
        j = j + 1
        if (runif(1) < q(x)/(alpha2*g2(x))){
          x2[i] = x
          test = F
        }
      }
    }
  }
  print(j)
  return(x2)
}

x = samplef(5000,1)
h = hist(x,plot=F)
hcum = h
hcum$counts=cumsum(hcum$counts)
plot(hcum)
plot(h,add=T,col='grey')
d = density(x)
lines(d$x,d$y*length(x)*diff(h$breaks)[1],lwd=2)
lines(d$x,cumsum(d$y)/max(cumsum(d$y))*length(x),lwd=2)

x = samplef(5000,2)
h = hist(x,plot=F)
hcum = h
hcum$counts=cumsum(hcum$counts)
plot(hcum)
plot(h,add=T,col='grey')
d = density(x)
lines(d$x,d$y*length(x)*diff(h$breaks)[1],lwd=2)
lines(d$x,cumsum(d$y)/max(cumsum(d$y))*length(x),lwd=2)