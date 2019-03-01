par(family = 'serif')
setwd("/Users/mikhailgaerlan/Box Sync/Education/UC Davis/2016-2017 Spring/STA 243 Computational Statistics/Assignments/Assignment 4")
rm(list=ls())

#=======================
#       3 (a)
#=======================
n = 1500
u = runif(n)

h = 1/(1+u)
Imc = mean(h)
Imcvar = sum((h-Imc)^2)/n
print(c(Imc,Imcvar,log(2)))

c = 1+u
thetamc = mean(c)
thetavar = sum((c-thetamc)^2)/(n-1)

covar = sum((h-Imc)*(c-thetamc))/(n-1)
b = covar/thetavar
print(b)
Icv = mean(h)-b*(mean(c)-3/2)
Icvvar = Imcvar+b^2*thetavar-2*b*covar
print(c(Icv,Icvvar,log(2)))

par(family = 'serif')
setwd("/Users/mikhailgaerlan/Box Sync/Education/UC Davis/2016-2017 Spring/STA 243 Computational Statistics/Assignments/Assignment 4")
rm(list=ls())

#=======================
#       3 (d)
#=======================
n = 1500

u = runif(n)
h = 1/(1+u)
Imc = mean(h)
var = sum((h-Imc)^2)/n

c1 = 1+u
theta1 = mean(c1)
var1 = sum((c1-theta1)^2)/(n-1)

c2 = 1-u/2
theta2 = mean(c2)
var2 = sum((c2-theta2)^2)/(n-1)

covar1 = sum((h-Imc)*(c1-theta1))/(n-1)
covar2 = sum((h-Imc)*(c2-theta2))/(n-1)
covar12 = sum((c1-theta1)*(c2-theta2))/(n-1)
b1 = (covar12*covar2+covar1*var2)/(var1*var2-covar12^2)
b2 = (covar12*covar1+covar2*var1)/(var1*var2-covar12^2)
Icv = mean(h)-b1*(mean(c1)-(3/2))-b2*(mean(c2)-(3/4))
Icvvar = var+b1^2*var1+b2^2*var2-2*b1*covar1-2*b2*covar2-2*b1*b2*covar12
print(c(b1,b2))
print(c(Icv,Icvvar,log(2)))