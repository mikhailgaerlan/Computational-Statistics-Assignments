require("LaplacesDemon")
require("extrafont")
par(family = 'Times New Roman')
setwd("/Users/mikhailgaerlan/Box Sync/Education/UC Davis/2016-2017 Spring/STA 243 Computational Statistics/Assignments/Assignment 6/Code")
rm(list=ls())
set.seed(0609)

n = 50
theta = 3
B = 5000

X = runif(n,0,theta)
thetarep = array(0,c(B))

#Parametric Bootstrap
for(i in 1:B){
    xsample = runif(n,0,theta)
    thetarep[i] = max(xsample)
}
p1 = hist(thetarep)
print(theta*n/(n+1))
print(n*theta^2/((n+2)*(n+1)^2))
thetamean = mean(thetarep)
print(thetamean)
print(sum((thetarep-thetamean)^2)/(B-1))

pdf(file = "../Graphs/parametricboot.pdf",width=4,height=3.5,family="Times New Roman")
plot(p1,col=rgb(1,0,0,1/4),xlim=c(2.4,3),main="",xlab="")
dev.off()


#Nonparametric Bootstrap
for(i in 1:B){
    xsample = sample(X,n,T)
    thetarep[i] = max(xsample)
}
p2 = hist(thetarep)
thetamean = mean(thetarep)
print(thetamean)
print(sum((thetarep-thetamean)^2)/(B-1))

pdf(file = "../Graphs/nonparametricboot.pdf",width=4,height=3.5,family="Times New Roman")
plot(p2,col=rgb(0,0,1,1/4),xlim=c(2.4,3),main="",xlab="")
dev.off()

#Distribution
x = seq(2.4,3,0.001)
f = n*x^(n-1)/theta^n
pdf(file = "../Graphs/distribution.pdf",width=4,height=3.5,family="Times New Roman")
plot(x,f,type='l',xlim=c(2.4,3),xlab=expression(italic(x)),ylab=expression(italic(f(x))))
dev.off()
