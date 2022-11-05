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
      if (sum(mult_matrix) == 0){
        pi[q,l] = 0
      } else {
        pi[q,l] = sum(mult_matrix*X)/sum(mult_matrix)
      }
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
        if (pi[q,l] != 0){
          test333 = (X * log(pi[q,l]) + (1-X) * log(1-pi[q,l]))
          diag(test333)=0
          sol = sol + test333 %*% res1[,l]
        }
      }
      res[,q] = (alpha[q] * exp(sol))[,1]
    }
    c(t(apply(res, 1, function(x) x/sum(x))))
  }
  matrix(FixedPoint(fixed_point_equation, c(tau_0), Method="Aitken")$FixedPoint,
         ncol=Q,nrow=n_nodes)
}

sim_ICL = function(X, Z_tilde, Q){
  print(Q)
  n_nodes = dim(X)[1]
  alpha_comp = 0
  for(q in 1:Q){
    num_Z_tilde_equ_q = sum(Z_tilde == q)
    if(num_Z_tilde_equ_q > 0){
      alpha_q = num_Z_tilde_equ_q / n_nodes
      alpha_comp = alpha_comp + num_Z_tilde_equ_q * log(alpha_q)
    }
  }
  pi_comp = 0
  for(q in 1:Q){
    for(l in 1:Q){
      num_Z_tilde_equ_q_l = (Z==q) %*% t(Z==l)
      diag(num_Z_tilde_equ_q_l) = 0
      if (sum(num_Z_tilde_equ_q_l) > 0) {
        m_1 = sum(num_Z_tilde_equ_q_l*X)
        m_2 = sum(num_Z_tilde_equ_q_l*(1-X))
        if (m_1 == 0){
          pi_q_l = 0
          pi_comp = pi_comp + (1/2)*sum(num_Z_tilde_equ_q_l*((1-X)*log(1-pi_q_l)))
        } else if (m_2 == 0){
          pi_q_l = 1
          pi_comp = pi_comp + (1/2)*sum(num_Z_tilde_equ_q_l*(X*log(pi_q_l)))
        } else{
          pi_q_l = m_1/(m_1+m_2)
          pi_comp = pi_comp + (1/2)*sum(num_Z_tilde_equ_q_l*(X*log(pi_q_l) + (1-X)*log(1-pi_q_l)))
        }
      } else{
        print("FAIL")
      }
    }
  }
  print("picomp")
  print(pi_comp)
  print("alpha")
  print(alpha_comp)
  print("REST")
  print(alpha_comp + pi_comp - (1/2)*Q*(Q+1)*log(n_nodes*(n_nodes-1)/2)+(Q-1)*log(n_nodes))
  alpha_comp + pi_comp - (1/2)*Q*(Q+1)*log(n_nodes*(n_nodes-1)/2)+(Q-1)*log(n_nodes)
}

Daudin_Var_Bayes = function(X, Q, epsilon=1e-6){
  n_nodes = dim(X)[1]
  abs_error_0 = 1
  abs_error_1 = 0
  Z_0 = matrix(kmeans(X, Q)$cl)
  tau_0=t(apply(Z_0, 1, function(j) c(rep(0,j-1), 1, rep(0,Q-j))))
  loop=0
  while(sum(abs_error_0-abs_error_1) > epsilon && loop <= 20){
    alpha_pi_0 = E_step(X, tau_0)
    alpha_0 = alpha_pi_0[[1]]
    pi_0 = alpha_pi_0[[2]]
    tau_0 = M_step(X, alpha_0, pi_0, tau_0)
    abs_error_1 = abs_error_0
    abs_error_0 = c(alpha_0, pi_0)
    loop = loop+1
  }
  if (loop==0){
    print("loop error")
  }
  Z_tilde = apply(tau_0, 1, which.max)
  Z_tilde
}

Daudin_Var_Bayes_model_selection = function(X, Q=20, epsilon=1e-6){
  Q_dir = Q
  ICL_prot_1 = sim_ICL(X, Daudin_Var_Bayes(X, Q_dir, epsilon), Q_dir)
  dir = 1
  ICL_prot_2 = sim_ICL(X, Daudin_Var_Bayes(X, Q_dir+dir, epsilon), Q_dir + dir)
  if (ICL_prot_2 <= ICL_prot_1){
    dir = -1
    ICL_prot_2 = sim_ICL(X, Daudin_Var_Bayes(X, Q_dir+dir, epsilon), Q_dir + dir) 
  }
  backup=0
  while(backup<=5 && ICL_prot_2 < -6000){
    if(ICL_prot_2 <= ICL_prot_1){
      ICL_prot_2 = sim_ICL(X, Daudin_Var_Bayes(X, Q_dir + dir, epsilon), Q_dir + dir)
      backup=backup+1
    } else{
    Q_dir = Q_dir + dir
    print(Q_dir)
    ICL_prot_1 = ICL_prot_2
    ICL_prot_2 = sim_ICL(X, Daudin_Var_Bayes(X, Q_dir + dir, epsilon), Q_dir + dir)
    backup=0
    }
  }
  print(Daudin_Var_Bayes(X, Q_dir + dir, epsilon))
  print(Q_dir + dir)
  Z_tilde = Daudin_Var_Bayes(X, Q_dir, epsilon)
  Z_tilde
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
K  <-  6           # Number of cluster
pi <- rep(1/K,K)    # Clusters proportions 
lambda   <- 0.1     # Building the connectivity matrix template
lambda_o <- 0.01
Ks <- 3
mu <- bdiag(lapply(1:(K/Ks), function(k){
  matrix(lambda_o,Ks,Ks)+diag(rep(lambda,Ks))}))+0.001
Z_X = sim_sbm(N, pi, mu)
X = Z_X[[2]]

Z = Z_X[[1]]
Z_0 = matrix(kmeans(X, K)$cl)
tau_0=t(apply(Z_0,1, function(j) c(rep(0,j-1),1,rep(0,K-j))))
#tau_0=matrix(runif(N*K),nrow=N,ncol=K)
#tau_0 = t(apply(tau_0, 1, function(x) x/sum(x)))

alpha_pi_0 = E_step(X, tau_0)
alpha_0 = alpha_pi_0[[1]]
pi_0 = alpha_pi_0[[2]]
pi_0
tau_0=M_step(X, alpha_0, pi_0, tau_0)
sum(abs(pi_0-mu))
Z_tilde = apply(tau_0, 1, which.max)
table(Z_tilde, Z)
#Daudin_Var_Bayes_model_selection(X)
d_test = Daudin_Var_Bayes_model_selection(X, Q=20)
table(d_test, Z)
sim_ICL(X, Daudin_Var_Bayes(X, 2), 2)
ICL(sol)
