par(family = 'serif')
setwd("/Users/mikhailgaerlan/Box Sync/Education/UC Davis/2016-2017 Spring/STA 243 Computational Statistics/Assignments/Assignment 4")
rm(list=ls())

options(digits=6)

betas = c(0.5,1,2)
m = 100  #Number of experiments
n = 1000 #Number of samples
x = array(0,c(length(betas),m,n))
theta1 = 1.5
theta2 = 2
for(i in 1:length(betas)){
  means = 1:m
  means1 = 1:m
  for (j in 1:m){
    f = function(z){
      return((z^(-3/2))*exp(-theta1*z-theta2/z+2*sqrt(theta1*theta2)+log(sqrt(2*theta2))))
    }
    
    g = function(x,y){
      alpha = y
      return(dgamma(x,alpha,betas[i]))
    }
    
    r = function(x,y){
      return(min(c(f(y)*g(x,y)/(f(x)*g(y,x)),1)))
    }
    
    x[i,j,1] = 1
    for(k in 2:n){
      y = rgamma(1,x[i,j,k-1],betas[i])
      if(runif(1) <= r(x[i,j,k-1],y)){
        x[i,j,k] = y
      }else{
        x[i,j,k]=x[i,j,k-1]
      }
    }
    means[j] = mean(x[i,j,])
    means1[j] = mean(1/x[i,j,])
    #print(c(betas[i],mu,vr,sqrt(2/1.5),mu1,vr1,sqrt(1.5/2)+1/(2*2)))
  }
  hist(x[i,1,])
  print(c(betas[i],mean(means),var(means),sqrt(theta2/theta1),mean(means1),var(means1),sqrt(theta1/theta2)+1/(2*theta2)))
}