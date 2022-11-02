library(future) # allows parralel processing in greed()
library(Matrix) # sparse matrix
library(ggplot2) # ploting and data 
library(greed)
library(dplyr)
library(ggpubr)

sim_sbm = function(n_nodes, alpha, pi){
  Q = length(alpha)
  Z = matrix(sample(1:Q, n_nodes, replace=TRUE, prob=alpha))
  test1 = Vectorize(function(x) pi[,x])
  test2 = function(x) Vectorize(function(y) x[y])(Z)
  X = Matrix(0,n_nodes,n_nodes)
  X = X + ((runif(n_nodes^2)+diag(1,n_nodes) <= apply(test1(Z), 2, test2)) + 0)
  X
}

K=6
Ks=3
#mu=bdiag(lapply(1:(K/Ks)))
sim_sbm(3, rep(1/3,3), matrix(rep(1/3,9),ncol=3,nrow=3))
