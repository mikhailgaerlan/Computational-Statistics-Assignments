par(family = 'serif')
setwd("/Users/mikhailgaerlan/Box Sync/Education/UC Davis/2016-2017 Spring/STA 243 Computational Statistics/Assignments/Assignment 4")
rm(list=ls())
set.seed(0517)

phi = function(u){
  return(exp(-u^2/2)/sqrt(2*pi))
}

f = function(x){
  return(1.5*phi((x-0.35)/0.15)-phi((x-0.8)/0.04))
}

smoothingmatrix = function(x,y){
  
}

aicc = function(lambda,x,y){
  X = array(0,c(n,p+1+k))
  for(i in 1:(p+1)){
    for(j in 1:n){
      X[j,i] = x[j]^(i-1)
    }
  }
  for(i in 1:k){
    for(j in 1:n){
      X[j,p+1+i] = (max(c((x[j]-knots[i]),0)))^p
    }
  }
  
  D = diag(c(0*(1:(p+1)),0*(1:k)+1))
  hlambda = X %*% (solve((t(X) %*% X)+lambda*D) %*% t(X))
  fhatlambda = hlambda %*% y
  
  n = length(x)
  tr = sum(diag(hlambda))
  norm = sum((y-fhatlambda)^2)
  return(log(norm)+2*(tr+1)/(n-tr-2))
}

# Parameter Values
p = 3
n = 200
j = 1
k = 30

# Data generation
x = (1:200-0.5)/n
knots = min(x) + 1:k*(max(x)-min(x))/(k+1)
sigma = 0.02 + 0.04*(j-1)^2
epsilon = rnorm(n,0,1)
y = f(x) + sigma*epsilon
plot(x,y,pch=".")

# Spline Regression
X = array(0,c(n,p+1+k))
for(i in 1:(p+1)){
  for(j in 1:n){
    X[j,i] = x[j]^(i-1)
  }
}
for(i in 1:k){
  for(j in 1:n){
    X[j,p+1+i] = (max(c((x[j]-knots[i]),0)))^p
  }
}
fhat = (X %*% (solve(t(X) %*% X) %*% t(X))) %*% y
print(aicc(0,x,y))
lines(x,fhat,col="blue")

# Penalized Spline Regression
lambda = 0
minlambda = nlm(aicc(lambda,x,y),c(lambda))
lambda = minlambda$estimate
D = diag(c(0*(1:(p+1)),0*(1:k)+1))
hlambda = X %*% (solve((t(X) %*% X)+lambda*D) %*% t(X))
fhatlambda = hlambda %*% y
print(aicc(lambda,x,y))
lines(x,fhatlambda,col="red")

D = diag(c(0*(1:(p+1)),0*(1:k)+1))
hlambda = X %*% (solve((t(X) %*% X)+0.000000000001*D) %*% t(X))
fhatlambda = hlambda %*% y
print(aicc(0.000000000001,x,y))
lines(x,fhatlambda,col="green")
