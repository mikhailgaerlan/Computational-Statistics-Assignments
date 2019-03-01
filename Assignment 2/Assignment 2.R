par(family = 'serif')
setwd("/Users/mikhailgaerlan/Box Sync/Education/UC Davis/2016-2017 Spring/STA 243 Computational Statistics/Assignments/Assignment 2")
rm(list=ls())

#=========
#    1
#=========
#Distance Between Cities
distances = matrix(
   #1,2,3,4,5,6,7,8,9,0,1,2,3,4,5
   #A,B,C,D,E,F,G,H,I,J,K,L,M,N,O
  c(0,1,2,4,9,8,3,2,1,5,7,1,2,9,3,#A,1
    1,0,5,3,7,2,5,1,3,4,6,6,6,1,9,#B,2
    2,5,0,6,1,4,7,7,1,6,5,9,1,3,4,#C,3
    4,3,6,0,5,2,1,6,5,4,2,1,2,1,3,#D,4
    9,7,1,5,0,9,1,1,2,1,3,6,8,2,5,#E,5
    8,2,4,2,9,0,3,5,4,7,8,3,1,2,5,#F,6
    3,5,7,1,1,3,0,2,6,1,7,9,5,1,4,#G,7
    2,1,7,6,1,5,2,0,9,4,2,1,1,7,8,#H,8
    1,3,1,5,2,4,6,9,0,3,3,5,1,6,4,#I,9
    5,4,6,4,1,7,1,4,3,0,9,1,8,5,2,#J,10
    7,6,5,2,3,8,7,2,3,9,0,2,1,8,1,#K,11
    1,6,9,1,6,3,9,1,5,1,2,0,5,4,3,#L,12
    2,6,1,2,8,1,5,1,1,8,1,5,0,9,6,#M,13
    9,1,3,1,2,2,1,7,6,5,8,4,9,0,7,#N,14
    3,9,4,3,5,5,4,8,4,2,1,3,6,7,0 #O,15
  ),
  nrow=15,
  ncol=15,
  byrow=TRUE
)

#Objective Function
obj = function(theta){
  result = distances[1,theta[1]]
  for (i in 1:(length(theta)-1)){
    result = result + distances[theta[i],theta[i+1]]
  }
  result = result + distances[1,theta[length(theta)]]
  return(result)
}
alpha = function(tau,p){
  result = p*tau
  return(result)
}
beta = function(m){
  result = 100
  return(result)
}
neighborhoodsample = function(theta){
  switch = sample(1:length(theta),2)
  newtheta = theta
  newtheta[switch[1]] = theta[switch[2]]
  newtheta[switch[2]] = theta[switch[1]]
  return(newtheta)
}

#Simulated Annealing Parameters
for (p in c(0.999,0.97,0.9)){
  for (tau in c(400,200,100)){
    printl = TRUE
    writel = TRUE
    maxiter = 10000
    prob = 10^(-50)
    resultsname = paste("resultsp",toString(p*1000),"t",toString(tau),".csv",sep = "")
    allname = paste("allp",toString(p*1000),"t",toString(tau),".csv",sep = "")
    
    theta = sample(2:15,14)
    #theta = c(9,3,5,10,15,11,13,6,14,7,4,12,2,8)
    #theta = 2:15
    if (writel){
      write(sprintf("%d %d %d %s",obj(theta),0,0,paste(theta,collapse = " ")),file = resultsname,append = FALSE)
      write(sprintf("%d %d %d %s",obj(theta),0,0,paste(theta,collapse = " ")),file = allname,append = FALSE)
    }
    #Begin Simulated Annealing
    break_outer = FALSE
    for (j in 1:maxiter){
      for (m in 1:beta(j)){
        newtheta = neighborhoodsample(theta)
        if (printl) print(sprintf("%2s %3d",paste(newtheta,collapse=" "),obj(newtheta)))
        if (writel) write(sprintf("%d %d %d %s",obj(newtheta),j-1,m-1,paste(newtheta,collapse=" ")),file = allname,append = TRUE)
        delta = obj(newtheta) - obj(theta)
        if (exp(-delta/tau) < prob){
          break_outer = TRUE
          break
        }
        if (delta <= 0){
          theta = newtheta
          if (writel) write(sprintf("%d %d %d %s",obj(theta),j-1,m-1,paste(theta,collapse=" ")),file = resultsname,append = TRUE)
        } else if (runif(1) < (exp(-delta/tau))){
          theta = newtheta
          if (writel) write(sprintf("%d %d %d %s",obj(theta),j-1,m-1,paste(theta,collapse=" ")),file = resultsname,append = TRUE)
        } else {
          if (writel) write(sprintf("%d %d %d %s",obj(theta),j-1,m-1,paste(theta,collapse=" ")),file = resultsname,append = TRUE)
        }
      }
      if (break_outer) break
      tau = alpha(tau,p)
    }
    print(sprintf("%2s %3d",paste(theta,collapse=" "),obj(theta)))
  }
}
#print(obj(c(9,3,5,10,15,11,13,6,14,7,4,12,8,2)))
