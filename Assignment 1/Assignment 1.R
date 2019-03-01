par(family = 'serif')
setwd("/Users/mikhailgaerlan/Box Sync/Education/UC Davis/2016-2017 Spring/STA 243 Computational Statistics/Assignments/Assignment 1")

#=========
#    1
#=========
#----------
#   1(c)
#----------
sample = c(-13.87,-2.53,-2.44,-2.40,-1.75,-1.34,-1.05,-0.23,-0.07,0.27,1.77,2.76,3.29,3.47,3.71,3.80,4.24,4.53,43.21,56.75)
n = length(sample)
thetas = seq(min(sample),max(sample),0.01)
loglikelihood = thetas*0
for (x in sample)
  loglikelihood = loglikelihood - log(1+(thetas-x)^2)
loglikelihood = loglikelihood - n*log(pi)

plot(loglikelihood ~ thetas,type='l',xlab=expression(theta),ylab=expression(paste(italic("l"),"(",theta,")")),main="Log Likelihood Function")

#----------
#   1(d)
#----------
starting = c(-11,-1,0,1.4,4.1,4.8,7,8,38)

print("1d")
maxiter = 200
err = 10^(-5)
for (theta0 in starting){
  theta = theta0
  for (t in 1:maxiter){
    oldtheta = theta
    lptheta = 0
    for (x in sample){
      lptheta = lptheta - 2*(theta-x)/(1+(theta-x)^2)
    }
    lpptheta = 0
    for (x in sample){
      lpptheta = lpptheta - 2*(1-(theta-x)^2)/(1+(theta-x)^2)^2
    }
    h = -lptheta/lpptheta
    theta = theta + h
    experr = abs(theta - oldtheta)
    if (experr < err) {
      break
    }
  }
  bigltheta = 1
  for (x in sample){
    bigltheta = bigltheta * 1/(pi*(1+(x-theta)^2))
  }
  print(sprintf("theta = %12.5e, L(theta) = %12.5e",theta,bigltheta))
}

#----------
#   1(e)
#----------
print("1e")
for (theta0 in starting){
  theta = theta0
  
  #Fisher-Scoring
  maxiter = 100
  err = 10^(-2)
  for (t in 1:maxiter){
    oldtheta = theta
    lptheta = 0
    for (x in sample){
      lptheta = lptheta - 2*(theta-x)/(1+(theta-x)^2)
    }
    h = 2*lptheta/n
    theta = theta + h
    experr = abs(theta - oldtheta)
    if (experr < err) {
      break
    }
  }
  
  #Newton-Raphson
  maxiter = 200
  err = 10^(-5)
  for (t in 1:maxiter){
    oldtheta = theta
    lptheta = 0
    for (x in sample){
      lptheta = lptheta - 2*(theta-x)/(1+(theta-x)^2)
    }
    lpptheta = 0
    for (x in sample){
      lpptheta = lpptheta - 2*(1-(theta-x)^2)/(1+(theta-x)^2)^2
    }
    h = -lptheta/lpptheta
    theta = theta + h
    experr = abs(theta - oldtheta)
    if (experr < err) {
      break
    }
  }
  bigltheta = 1
  for (x in sample){
    bigltheta = bigltheta * 1/(pi*(1+(x-theta)^2))
  }
  print(sprintf("theta = %12.5e, L(theta) = %12.5e",theta,bigltheta))
}

rm(list=ls())
#==========
#    2
#==========
#----------
#   2(a)
#----------
sample = c(0.52,1.96,2.22,2.28,2.28,2.46,2.50,2.53,2.54,2.99,3.47,3.53,3.70,3.88,3.91,4.04,4.06,4.82,4.85,5.46)
thetas = seq(-pi,pi,0.01)
loglikelihood = thetas*0
for (x in sample){
  loglikelihood = loglikelihood + log((1-cos(x-thetas))/(2*pi))
}
plot(loglikelihood ~ thetas,type='l',xlab=expression(theta),ylab=expression(paste(italic("l"),"(",theta,")")),main="Log Likelihood Function")

#----------
#   2(c)
#----------
print("2c")
maxiter = 200
err = 10^(-5)
theta = pi
for (t in 1:maxiter){
  oldtheta = theta
  lptheta = 0
  for (x in sample){
    lptheta = lptheta - sin(x-theta)/(1-cos(x-theta))
  }
  lpptheta = 0
  for (x in sample){
    lpptheta = lpptheta +  1/(cos(x-theta)-1)
  }
  h = -lptheta/lpptheta
  theta = theta + h
  experr = abs(theta - oldtheta)
  if (experr < err) {
    break
  }
}
bigltheta = 1
for (x in sample){
  bigltheta = bigltheta * ((1-cos(x-theta))/(2*pi))
}
print(sprintf("theta = %12.5e, L(theta) = %12.5e",theta,bigltheta))

#----------
#   2(d)
#----------
starting = c(-2.7,2.7)

print("2d")
maxiter = 200
err = 10^(-5)
for (theta0 in starting){
  theta = theta0
  for (t in 1:maxiter){
    oldtheta = theta
    lptheta = 0
    for (x in sample){
      lptheta = lptheta - sin(x-theta)/(1-cos(x-theta))
    }
    lpptheta = 0
    for (x in sample){
      lpptheta = lpptheta +  1/(cos(x-theta)-1)
    }
    h = -lptheta/lpptheta
    theta = theta + h
    experr = abs(theta - oldtheta)
    if (experr < err) {
      break
    }
  }
  bigltheta = 1
  for (x in sample){
    bigltheta = bigltheta * ((1-cos(x-theta))/(2*pi))
  }
  print(sprintf("theta = %12.5e, L(theta) = %12.5e",theta,bigltheta))
}
#----------
#   2(e)
#----------
starting = seq(-pi,pi,pi/100)
ending = starting*0
bigltheta = 0*ending+1

print("2e")
maxiter = 300
err = 10^(-5)
for (i in 1:length(starting)){
  ending[i] = starting[i]
  for (t in 1:maxiter){
    oldtheta = ending[i]
    lptheta = 0
    for (x in sample){
      lptheta = lptheta - sin(x-ending[i])/(1-cos(x-ending[i]))
    }
    lpptheta = 0
    for (x in sample){
      lpptheta = lpptheta +  1/(cos(x-ending[i])-1)
    }
    h = -lptheta/lpptheta
    ending[i] = ending[i] + h
    experr = abs(ending[i] - oldtheta)
    if (experr < err) {
      break
    }
  }
  for (x in sample){
    bigltheta[i] = bigltheta[i] * ((1-cos(x-ending[i]))/(2*pi))
  }
#  print(sprintf("i = %3d theta = %8.5f, L(theta) = %12.5e",i,ending[i],bigltheta[i]))
}

i = 1
while (i < length(starting)){
  for (j in i:length(starting)){
    if (sprintf("%8.5f",ending[j]) != sprintf("%8.5f",ending[j+1])){
      print(sprintf("a = %8.5f, b = %8.5f, theta = %8.5f, MLE = %12.5e",starting[i],starting[j],ending[i],bigltheta[i]))
      i = j
      break
    }
  }
  i = i+1
}

rm(list=ls())
#==========
#    3
#==========
#----------
#   3(a)
#----------
print("3a")
x = c(0.02,0.06,0.11,0.22,0.56,1.10)
y1 = c(47,97,123,152,191,200)

xm = 1/x
ym = 1/y1

model = lm(formula=ym~xm)
theta1 = c(1/coef(model)[1],coef(model)[2]/coef(model)[1])
print(sprintf("theta_1 = %8.5f, theta_2 = %8.5f",theta1[1],theta1[2]))
xtheo = seq(min(x),max(x),0.01)
ytheo = theta1[1]*xtheo/(xtheo+theta1[2])
plot(x,y1,xlab="substrate concentration",ylab="velocity")
lines(xtheo,ytheo)

y2 = c(76,107,139,159,201,207)

xm = 1/x
ym = 1/y2

model = lm(formula=ym~xm)
theta2 = c(1/coef(model)[1],coef(model)[2]/coef(model)[1])
print(sprintf("theta_1 = %8.5f, theta_2 = %8.5f",theta2[1],theta2[2]))
xtheo = seq(min(x),max(x),0.01)
ytheo = theta2[1]*xtheo/(xtheo+theta2[2])
plot(x,y2,xlab="substrate concentration",ylab="velocity")
lines(xtheo,ytheo)

#----------
#   3(b)
#----------
print("3b")
theta = theta1
maxiter = 400
err = 10^(-5)
for (t in 1:maxiter){
  dRSSdtheta1 = 0
  for (i in 1:length(x)){
    dRSSdtheta1 = dRSSdtheta1 - 2*(y1[i]-theta[1]*x[i]/(x[i]+theta[2]))*(x[i]/(x[i]+theta[2]))
  }
  dRSSdtheta2 = 0
  for (i in 1:length(x)){
    dRSSdtheta2 = dRSSdtheta2 + 2*(y1[i]-theta[1]*x[i]/(x[i]+theta[2]))*(theta[1]*x[i]/(x[i]+theta[2])^2)
  }
  d2RSSdtheta12 = 0
  for (i in 1:length(x)){
    d2RSSdtheta12 = d2RSSdtheta12 + 2*(x[i]/(x[i]+theta[2]))^2
  }
  d2RSSdtheta22 = 0
  for (i in 1:length(x)){
    d2RSSdtheta22 = d2RSSdtheta22 + 2*((theta[1]^2*x[i]^2)/(x[i]+theta[2])^4-(2*theta[1]*x[i]/(x[i]+theta[2])^3)*(y1[i]-theta[1]*x[i]/(x[i]+theta[2])))
  }
  d2RSSdtheta1dtheta2 = 0
  for (i in 1:length(x)){
    d2RSSdtheta1dtheta2 = d2RSSdtheta1dtheta2 - 2*(-x[i]*y1[i]/(x[i]+theta[2])^2+2*theta[1]*x[i]^2/(x[i]+theta[2])^3)
  }
  gradRSS = c(dRSSdtheta1,dRSSdtheta2)
  hessRSS = matrix(c(d2RSSdtheta12,d2RSSdtheta1dtheta2,d2RSSdtheta1dtheta2,d2RSSdtheta22),nrow=2,ncol=2)
  theta = theta - solve(hessRSS) %*% gradRSS
}
print(sprintf("theta_1 = %8.5f, theta_2 = %8.5f",theta[1],theta[2]))
plot(x,y1,xlab="substrate concentration",ylab="velocity")
xtheo = seq(min(x),max(x),0.01)
ytheo = theta[1]*xtheo/(xtheo+theta[2])
lines(xtheo,ytheo)

theta = theta2
maxiter = 400
err = 10^(-5)
for (t in 1:maxiter){
  dRSSdtheta1 = 0
  for (i in 1:length(x)){
    dRSSdtheta1 = dRSSdtheta1 - 2*(y2[i]-theta[1]*x[i]/(x[i]+theta[2]))*(x[i]/(x[i]+theta[2]))
  }
  dRSSdtheta2 = 0
  for (i in 1:length(x)){
    dRSSdtheta2 = dRSSdtheta2 + 2*(y2[i]-theta[1]*x[i]/(x[i]+theta[2]))*(theta[1]*x[i]/(x[i]+theta[2])^2)
  }
  d2RSSdtheta12 = 0
  for (i in 1:length(x)){
    d2RSSdtheta12 = d2RSSdtheta12 + 2*(x[i]/(x[i]+theta[2]))^2
  }
  d2RSSdtheta22 = 0
  for (i in 1:length(x)){
    d2RSSdtheta22 = d2RSSdtheta22 + 2*((theta[1]^2*x[i]^2)/(x[i]+theta[2])^4-(2*theta[1]*x[i]/(x[i]+theta[2])^3)*(y2[i]-theta[1]*x[i]/(x[i]+theta[2])))
  }
  d2RSSdtheta1dtheta2 = 0
  for (i in 1:length(x)){
    d2RSSdtheta1dtheta2 = d2RSSdtheta1dtheta2 - 2*(-x[i]*y2[i]/(x[i]+theta[2])^2+2*theta[1]*x[i]^2/(x[i]+theta[2])^3)
  }
  gradRSS = c(dRSSdtheta1,dRSSdtheta2)
  hessRSS = matrix(c(d2RSSdtheta12,d2RSSdtheta1dtheta2,d2RSSdtheta1dtheta2,d2RSSdtheta22),nrow=2,ncol=2)
  theta = theta - solve(hessRSS) %*% gradRSS
}
print(sprintf("theta_1 = %8.5f, theta_2 = %8.5f",theta[1],theta[2]))
plot(x,y2,xlab="substrate concentration",ylab="velocity")
xtheo = seq(min(x),max(x),0.01)
ytheo = theta[1]*xtheo/(xtheo+theta[2])
lines(xtheo,ytheo)

#----------
#   3(c)
#----------
print("3c")
theta = theta1
maxiter = 600
err = 10^(-5)
for (t in 1:maxiter){
  dRSSdtheta1 = 0
  for (i in 1:length(x)){
    dRSSdtheta1 = dRSSdtheta1 + 2*(y1[i]-theta[1]*x[i]/(x[i]+theta[2]))*(x[i]/(x[i]+theta[2]))
  }
  dRSSdtheta2 = 0
  for (i in 1:length(x)){
    dRSSdtheta2 = dRSSdtheta2 - 2*(y1[i]-theta[1]*x[i]/(x[i]+theta[2]))*(theta[1]*x[i]/(x[i]+theta[2])^2)
  }
  gradRSS = c(dRSSdtheta1,dRSSdtheta2)
  alphat = 0.01
  theta = theta + alphat*gradRSS
}
print(sprintf("theta_1 = %8.5f, theta_2 = %8.5f",theta[1],theta[2]))
plot(x,y1,xlab="substrate concentration",ylab="velocity")
xtheo = seq(min(x),max(x),0.01)
ytheo = theta[1]*xtheo/(xtheo+theta[2])
lines(xtheo,ytheo)

theta = theta2
maxiter = 400
err = 10^(-5)
for (t in 1:maxiter){
  dRSSdtheta1 = 0
  for (i in 1:length(x)){
    dRSSdtheta1 = dRSSdtheta1 - 2*(y2[i]-theta[1]*x[i]/(x[i]+theta[2]))*(x[i]/(x[i]+theta[2]))
  }
  dRSSdtheta2 = 0
  for (i in 1:length(x)){
    dRSSdtheta2 = dRSSdtheta2 + 2*(y2[i]-theta[1]*x[i]/(x[i]+theta[2]))*(theta[1]*x[i]/(x[i]+theta[2])^2)
  }
  gradRSS = c(dRSSdtheta1,dRSSdtheta2)
  alphat = alphat = 0.01
  theta = theta - alphat*gradRSS
}
print(sprintf("theta_1 = %8.5f, theta_2 = %8.5f",theta[1],theta[2]))
plot(x,y2,xlab="substrate concentration",ylab="velocity")
xtheo = seq(min(x),max(x),0.01)
ytheo = theta[1]*xtheo/(xtheo+theta[2])
lines(xtheo,ytheo)

#----------
#   3(d)
#----------
print("3d")
theta = theta1
maxiter = 400
err = 10^(-5)
for (t in 1:maxiter){
  amat = matrix(data=0,nrow=length(x),ncol=2)
  for (i in 1:length(x)){
    amat[i,] = c(x[i]/(x[i]+theta[2]),-theta[1]*x[i]/(x[i]+theta[2])^2)
  }
  zvec = x*0
  for (i in 1:length(x)){
    zvec[i] = y1[i]-theta[1]*x[i]/(x[i]+theta[2])
  }
  theta = theta + solve(t(amat) %*% amat) %*% (t(amat) %*% zvec)
}
print(sprintf("theta_1 = %8.5f, theta_2 = %8.5f",theta[1],theta[2]))
plot(x,y1,xlab="substrate concentration",ylab="velocity")
xtheo = seq(min(x),max(x),0.01)
ytheo = theta[1]*xtheo/(xtheo+theta[2])
lines(xtheo,ytheo)

theta = theta2
maxiter = 400
err = 10^(-5)
for (t in 1:maxiter){
  amat = matrix(data=0,nrow=length(x),ncol=2)
  for (i in 1:length(x)){
    amat[i,] = c(x[i]/(x[i]+theta[2]),-theta[1]*x[i]/(x[i]+theta[2])^2)
  }
  zvec = x*0
  for (i in 1:length(x)){
    zvec[i] = y2[i]-theta[1]*x[i]/(x[i]+theta[2])
  }
  theta = theta + solve(t(amat) %*% amat) %*% (t(amat) %*% zvec)
}
print(sprintf("theta_1 = %8.5f, theta_2 = %8.5f",theta[1],theta[2]))
plot(x,y2,xlab="substrate concentration",ylab="velocity")
xtheo = seq(min(x),max(x),0.01)
ytheo = theta[1]*xtheo/(xtheo+theta[2])
lines(xtheo,ytheo)