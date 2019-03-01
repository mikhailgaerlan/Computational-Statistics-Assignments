par(family = 'serif')
setwd("/Users/mikhailgaerlan/Box Sync/Education/UC Davis/2016-2017 Spring/STA 243 Computational Statistics/Assignments/Assignment 4")
rm(list=ls())

#=======================
#       4 (a)
#=======================
require("Rlab")

n = 100
lambda = 2
p = 0.3
y = rpois(n,lambda)
r = rbern(n,p)
x = y*r
hist(x)

