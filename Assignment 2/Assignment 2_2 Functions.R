#=======================
# Program functions
#=======================
getparams = function(chromo){
  pieces = sum(chromo)+1
  breaks = 0*1:(pieces-1)
  breakindex = 0*1:(pieces-1)
  heights = 0*1:pieces
  j = 1
  for (n in 1:length(chromo)){
    if (chromo[n] == 1){
      breaks[j] = x[n+1]
      breakindex[j] = n+1
      if (j == 1){
        heights[j] = mean(y[1:n])
      } else{
        heights[j] = mean(y[(breakindex[j-1]+1):breakindex[j]])
      }
      sum = y[n]
      j = j + 1
    }
  }
  heights[pieces] = mean(y[(breakindex[pieces-1]):length(x)])
  return(c(pieces,breaks,heights,breakindex))
}

modelfunction = function(x,t,h){
  temp = 0
  for (j in 1:length(h)){
    if (j == 1){
      if ((x >= 0)&&(x < t[j])){
        temp = h[j]
        break
      }
    } else if (j == length(h)){
      if ((x >= t[j-1])&&(x <= 1)){
        temp = h[j]
        break
      }
    } else {
      if ((x >= t[j-1])&&(x < t[j])){
        temp = h[j]
        break
      }
    }
  }
  return(temp)
}

mdl = function(chromo){
  params = getparams(chromo)
  limits = params[2:(params[1])]
  heights = params[(params[1]+1):(2*params[1])]
  indices = params[(2*params[1]+1):length(params)]
  nj = 1:params[1]
  nj[1] = indices[1]-1
  for (i in 2:(params[1]-1)){
    nj[i] = indices[i]-indices[i-1]
  }
  nj[params[1]] = n+1-indices[params[1]-1]
  ymodel = getmodel(chromo,limits,heights)
  return(params[1]*log(n)+(1/2)*sum(log(nj))+(n/2)*log((1/n)*sum((y-ymodel)^2)))
}

aic = function(chromo){
  params = getparams(chromo)
  limits = params[2:(params[1])]
  heights = params[(params[1]+1):(2*params[1])]
  ymodel = getmodel(chromo,limits,heights)
  return(n*log((1/n)*sum(y-ymodel)^2)+2*params[1]*log(n))
}

getmodel = function(chromo,limits,heights){
  ymodel = 0*x
  for (i in 1:length(x)){
    ymodel[i] = modelfunction(x[i],limits,heights)
  }
  return(ymodel)
}

rank = function(population,fitness){
  rankings = 0*fitness
  sortedfitness = sort(fitness,decreasing=TRUE)
  for (i in 1:dim(population)[1]){
    for (j in 1:dim(population)[1]){
      if (sortedfitness[i]==fitness[j]){
        rankings[j] = i
        break
      }
    }
  }
  return(rankings)
}