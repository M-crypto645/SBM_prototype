library(future) # allows parralel processing in greed()
library(Matrix) # sparse matrix
library(ggplot2) # ploting and data 
library(greed)
library(dplyr)
library(ggpubr)
library(FixedPoint)

sim_sbm = function(n_nodes, alpha, pi){
  Q = length(alpha)
  Z = matrix(sample(1:Q, n_nodes, replace=TRUE, prob=alpha))
  select_col = Vectorize(function(x) pi[,x])
  select_row = function(x) Vectorize(function(y) x[y])(Z)
  X = Matrix(0,n_nodes,n_nodes)
  pi = apply(select_col(Z), 2, select_row)
  X = X + ((runif(n_nodes^2)+diag(1,n_nodes) <= pi) + 0)
  list(Z, X, pi)
}

#NO SELF LOOPS

E_step = function(X, tau){
  Q = dim(tau)[2]
  pi = matrix(numeric(Q^2),Q,Q)
  n_nodes = dim(tau)[1]
  Q = dim(tau)[2]
  alpha = sapply(1:Q, function(q) mean(tau[,q]))
  for (l in 1:Q){
    for (q in 1:Q){
      mult_matrix = apply(matrix(1:n_nodes), 1, function(i) tau[i,q]) %*% 
        t(apply(matrix(1:n_nodes), 1, function(i) tau[i,l]))
      diag(mult_matrix) = 0
      pi[q,l] = sum(mult_matrix*X)/sum(mult_matrix)
    }
  }
  list(alpha, pi)
}

M_step = function(X, alpha, pi, tau_0){
  n_nodes = dim(X)[1]
  Q = length(alpha)
  #tau_0=matrix(rep(1/3,n_nodes*Q),nrow=n_nodes,ncol=Q)
  fixed_point_equation = function(tau_vec){
    for (q in 1:Q){
      res1 = matrix(tau_vec, n_nodes, Q)
      res = matrix(tau_vec, n_nodes, Q)
      sol = matrix(numeric(n_nodes))
      for(l in 1:Q){
        test333 = (X * log(pi[q,l]) + (1-X) * log(1-pi[q,l]))
        diag(test333)=0
        sol = sol + test333 %*% res1[,l]
      }
      res[,q] = (alpha[q] * exp(sol))[,1]
    }
    c(t(apply(res, 1, function(x) x/sum(x))))
  }
  matrix(FixedPoint(fixed_point_equation, c(tau_0), Method="Aitken")$FixedPoint,
         ncol=Q,nrow=n_nodes)
}

# future::plan("multisession", workers=5)
# library(foreach)
# library(doParallel)
# n_nodes=n_nodes
# alpha=alpha
# Q=Q
# log_b = log_b
# pi=pi
# X=X
# n_workers= 9
# cl <- makePSOCKcluster(n_workers)
# clusterSetRNGStream(cl)
# registerDoParallel(cl, cores = n_workers)

N  <- 400           # Number of node
K  <- 6             # Number of cluster
pi <- rep(1/K,K)    # Clusters proportions 
lambda   <- 0.1     # Building the connectivity matrix template
lambda_o <- 0.01
Ks <- 3
mu <- bdiag(lapply(1:(K/Ks), function(k){
  matrix(lambda_o,Ks,Ks)+diag(rep(lambda,Ks))}))+0.001
Z_X = sim_sbm(N, pi, mu)
X = Z_X[[2]]

Z = Z_X[[1]]
Z_0 = matrix(kmeans(X, 6)$cl)
tau_0=t(apply(Z_0,1, function(j) c(rep(0,j-1),1,rep(0,6-j))))
#tau_0=matrix(runif(N*K),nrow=N,ncol=K)
#tau_0 = t(apply(tau_0, 1, function(x) x/sum(x)))

alpha_pi_0 = E_step(X, tau_0)
alpha_0 = alpha_pi_0[[1]]
pi_0 = alpha_pi_0[[2]]
pi_0
tau_0=M_step(X, alpha_0, pi_0, tau_0)
sum(abs(pi_0-mu))
table(apply(tau_0, 1, which.max), Z)
