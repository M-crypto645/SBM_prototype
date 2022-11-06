library(Matrix) # sparse matrix
library(FixedPoint)
library(cluster)

sim_sbm = function(n_nodes, alpha, pi){
  Q = length(alpha)
  Z = matrix(sample(1:Q, n_nodes, replace=TRUE, prob=alpha))
  select_col = Vectorize(function(x) pi[,x])
  select_row = function(x) Vectorize(function(y) x[y])(Z)
  pi_formated = apply(select_col(Z), 2, select_row)
  X = Matrix(0,n_nodes,n_nodes)
  X = X + (runif(n_nodes^2)+diag(1,n_nodes) <= pi_formated) + 0
  sbm = list(cl = c(Z), x = X, mu = pi, K = Q, N = n_nodes)
  class(sbm) <- "Stochastic Block Model (Daudin 2008)"
  sbm
}

E_step = function(X, tau){
  Q = dim(tau)[2]
  pi = matrix(numeric(Q^2),Q,Q)
  n_nodes = dim(tau)[1]
  alpha = sapply(1:Q, function(q) mean(tau[,q]))
  for (l in 1:Q){
    for (q in 1:Q){
      mult_matrix = apply(matrix(1:n_nodes), 1, function(i) tau[i,q]) %*% 
        t(apply(matrix(1:n_nodes), 1, function(i) tau[i,l]))
      diag(mult_matrix) = 0
      if (sum(mult_matrix) == 0){
        pi[q,l] = 0 #Doesnt really matter
      } else {
        pi[q,l] = sum(mult_matrix*X)/sum(mult_matrix)
      }
    }
  }
  res = list(pi=alpha, mu=pi + .0001)
  class(res) = "E_step (Daudin 2008)"
  res
}


M_step = function(X, alpha, pi, tau_0){
  n_nodes = dim(X)[1]
  Q = length(alpha)
  fixed_point_equation = function(tau_vec){
    res = matrix(tau_vec, n_nodes, Q)
    res1 = matrix(tau_vec, n_nodes, Q)
    for (q in 1:Q){
      sol = matrix(numeric(n_nodes))
      for(l in 1:Q){
        mult_matrix = X
        mult_matrix[X==1] = log(pi[q,l]) 
        mult_matrix[X==0] = log(1-pi[q,l])
        diag(mult_matrix)=0
        add_matrix = mult_matrix %*% res1[,l]
        add_matrix[is.na(add_matrix)] = 0
        sol = sol + add_matrix
      }
      res[,q] = (alpha[q] * exp(sol))[,1]
    }
    #TODO: What if one row is 0???
    c(t(apply(res, 1, function(x) x/sum(x))))
  }
  res=matrix(FixedPoint(fixed_point_equation, c(tau_0), Method="RRE", 
                        ConvergenceMetricThreshold=1e-3, MaxIter=5)$FixedPoint,
             ncol=Q,nrow=n_nodes)
  res
}

sim_ICL = function(X, Z_tilde, Q){
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
      num_Z_tilde_equ_q_l = Matrix(Z_tilde==q) %*% t(Matrix(Z_tilde==l))
      diag(num_Z_tilde_equ_q_l) = 0
      m_1 = sum(num_Z_tilde_equ_q_l*X)
      m_2 = sum(num_Z_tilde_equ_q_l*(1-X))
      if (m_1 > 0 && m_2 > 0) {
        pi_comp = pi_comp + (1/2)*(m_1*log(m_1/(m_1+m_2)) + m_2*log(m_2/(m_1+m_2)))
      }
    }
  }
  ICL_res = alpha_comp + pi_comp - 
    (1/2)*(Q*(Q+1)*log(n_nodes*(n_nodes-1)/2)/2+(Q-1)*log(n_nodes))
  print(paste("ICL: ", as.character(ICL_res), "for ", as.character(Q), "clusters"))
  ICL_res
}

Daudin_Var_Bayes = function(X, Q, epsilon=1e-3, in_cluster=FALSE){
  n_nodes = dim(X)[1]
  abs_error_0 = 1
  abs_error_1 = 0
  if (isFALSE(in_cluster) == FALSE){
    Z_0 = matrix(cutree(in_cluster, k=Q))
  } else{
    Z_0 = matrix(cutree(agnes(X, method = "ward"), k=Q))
  }
  tau_0=t(apply(Z_0, 1, function(j) c(rep(0,j-1), 1, rep(0,Q-j))))
  loop_counter=0
  while(sum(abs(abs_error_0-abs_error_1)) > epsilon){
    e_step = E_step(X, tau_0)
    tau_0_backup = tau_0
    tau_0 = M_step(X, e_step$pi, e_step$mu, tau_0)
    abs_error_1 = abs_error_0
    abs_error_0 = c(e_step$pi, e_step$mu)
    loop_counter = loop_counter + 1
    if(loop_counter == 5 || sum(is.na(tau_0)) > 0){
      tau_0 = tau_0_backup
      abs_error_0=0
      abs_error_1=0
    }
  }
  apply(tau_0, 1, which.max)
  #TODO RETURN CLASS
}

Daudin_Var_Bayes_model_selection = function(X, Q=20, epsilon=1e-6){
  Q_dir = Q
  in_cluster = agnes(X, method = "ward")
  ICL_prot_1 = sim_ICL(X, Daudin_Var_Bayes(X, Q_dir, epsilon, in_cluster), Q_dir)
  dir = 1
  ICL_prot_2 = sim_ICL(X, Daudin_Var_Bayes(X, Q_dir+dir, epsilon, in_cluster), Q_dir + dir)
  if (ICL_prot_2 <= ICL_prot_1){
    dir = -1
    ICL_prot_2 = sim_ICL(X, Daudin_Var_Bayes(X, Q_dir+dir, epsilon, in_cluster), Q_dir + dir) 
  }
  while(ICL_prot_2 > ICL_prot_1){
    Q_dir = Q_dir + dir
    ICL_prot_1 = ICL_prot_2
    ICL_prot_2 = sim_ICL(X, Daudin_Var_Bayes(X, Q_dir + dir, epsilon, in_cluster), Q_dir + dir)
  }
  Daudin_Var_Bayes(X, Q_dir, epsilon)
  #TODO S4-class
}

# N0  <- 400           # Number of node
# K0  <-  3           # Number of cluster
# pi0 <- rep(1/K0,K0)    # Clusters proportions
# lambda   <- 0.1     # Building the connectivity matrix template
# lambda_o <- 0.01
# Ks <- 3
# mu0 <- bdiag(lapply(1:(K0/Ks), function(k){
#   matrix(lambda_o,Ks,Ks)+diag(rep(lambda,Ks))}))+0.001
# sbm = sim_sbm(N0, pi0, mu0)
# 
# d_test = Daudin_Var_Bayes_model_selection(sbm$x, Q=6)
# table(d_test, sbm$cl)
# d_test_2 = Daudin_Var_Bayes(sbm$x, sbm$K)
# sim_ICL(sbm$x, d_test_2, 9)
# table(d_test_2, sbm$cl)

#bug_pi = savepoint[[2]]
#bug_mu = savepoint[[3]]
#bug_tau = savepoint_tau0
#M_step(sbm$x, bug_pi, bug_mu, bug_tau)

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

