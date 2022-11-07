#' @title Simulating a Stochastic Block Model
#'
#' @description Simulates an SBM
#' @param n_nodes number of nodes
#' @param alpha vector of class-probabilities
#' @param pi matriX of edge-probabilities
#' @return a \code{dgCMatrix} object representing the adjacency matrix of X
#' @examples
#' sim_sbm(6, rep(1/6,6), bdiag(lapply(1:(K/Ks), function(k){
#' matrix(lambda_o,Ks,Ks)+diag(rep(lambda,Ks))}))+0.001)
sim_sbm = function(n_nodes, alpha, pi){
  Q = length(alpha)
  Z = matrix(sample(1:Q, n_nodes, replace=TRUE, prob=alpha))
  select_col = Vectorize(function(x) pi[,x])
  select_row = function(x) Vectorize(function(y) x[y])(Z)
  pi_formated = apply(select_col(Z), 2, select_row)
  X = Matrix(0,n_nodes,n_nodes)
  X = X + (runif(n_nodes^2)+diag(1,n_nodes) <= pi_formated) + 0
  sbm = list(cl = c(Z), x = X, mu = pi, K = Q, N = n_nodes)
  class(sbm) <- "SBM"
  sbm
}

#' @title E-step of the Variational EM-Algorithm
#'
#' @description Calculates the E-step of the Variaationaal EM-Algorithm
#' @param sbm An \code{adj_matrix_type}-object for X
#' @param tau A stochastic matrix to estimate the posterior class probabilities
#' @return a \code{Var_E_step} object containing estimates of pi and alpha
var_e_step = function(sbm, tau){
  if(class(sbm) == "SBM"){
    X = sbm$x
  } else{
    X = sbm
  }
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
  class(res) = "Var_E_step"
  res
}

#' @title M-step of the Variational EM-Algorithm
#'
#' @description Calculates the M-step of the Variaationaal EM-Algorithm
#' @param sbm An \code{adj_matrix_type}-object for X
#' @param alpha vector of class-probabilities
#' @param pi matriX of edge-probabilities
#' @param tau_0 A stochastic matrix to estimate the posterior class probabilities
#' (used as initial value for the fixed-point-solver)
#' @param epsilon=1e-3 Convergence threshold for the fixed-point-solver
#' @param max_iter=1e-3 Maximal number of iterations for the fixed-point-solver
#' @return a \code{Var_M_step} object containing an estimate of tau
var_m_step = function(sbm, alpha, pi, tau_0, epsilon=1e-3, max_iter=5){
  if(class(sbm) == "SBM"){
    X = sbm$x
  } else{
    X = sbm
  }
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
    #What if one row is 0???
    c(t(apply(res, 1, function(x) x/sum(x))))
  }
  fixed_point = FixedPoint(fixed_point_equation, c(tau_0), Method="RRE", 
                   ConvergenceMetricThreshold=epsilon, MaxIter = max_iter)
  res=list(fixed_point_epsilon=epsilon, method="RRE", fpevals = fixed_point$fpevals,
           tau_0=matrix(fixed_point$FixedPoint, ncol=Q,nrow=n_nodes))
  class(res) <- "Var_M_step"
  res
}

#' @title The ICL (Daudin, 2008)
#'
#' @description Simulates the ICL from Daudin's 2008 paper
#' @param sbm An \code{adj_matrix_type}-object for X
#' @param Z_tilde A clustering of X
#' @return a number (the ICL)
sim_ICL = function(sbm, Z_tilde, Q){
  if(class(sbm) == "SBM"){
    X = sbm$x
  } else{
    X = sbm
  }
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

#' @title The Variational EM-Algorithm
#'
#' @description Cycles between the E- and M-step of the Variational EM-Algorithm
#' until a convergences threshold is passed
#' @param sbm An \code{adj_matrix_type}-object for X
#' @param Q number of clusters
#' @param epsilon=1e-3 convergence threshold for the algorithm
#' @param in_cluster=FALSE A prior clustering of X. Initialized if FALSE
#' @param max_iter Maximal number of iterations for the algorithm
#' @param conv_threshold_m_step=1e-3 Convergence threshold for the fixed-point-solver of the m-step
#' @param max_iter_m_step=1e-3 Maximal number of iterations for the fixed-point-solver of the m-step
#' @return a \code{Var_Bayes_output} object all relevant informations + estimates
var_bayes = function(sbm, Q, epsilon=1e-3, in_cluster=FALSE, max_iter=5,
                     conv_threshold_m_step=1e-3, max_iter_m_step=5){
  if(class(sbm) == "SBM"){
    X = sbm$x
  } else{
    X = sbm
  }
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
    e_step = var_e_step(X, tau_0)
    m_step = var_m_step(X, e_step$pi, e_step$mu, tau_0, 
                        max_iter=max_iter_m_step, 
                        epsilon=conv_threshold_m_step)
    tau_0 = m_step$tau_0
    abs_error_1 = abs_error_0
    abs_error_0 = c(e_step$pi, e_step$mu)
    loop_counter = loop_counter + 1
    if(loop_counter == max_iter){
      break;
    }
  }
  cl=apply(tau_0, 1, which.max)
  info = c(e_step, m_step, list(Q=Q))
  info$tau_0 = NULL #tau_0 is not interesting, since cl is already computed
  res = new("Var_Bayes_output", x=X, cl=cl, info=info, Q=Q)
  res
}


#' @title The Variational EM-Algorithm with model selection
#'
#' @description Cycles between the Variational EM-Algorithm for different numbers
#' of classes until the class with the best ICL-value is found. 
#' It is assumed that the objective function of the ICL is monotonically decreasing
#' Still quite slow; Only good for small value of Q (<=6)
#' @param sbm An \code{adj_matrix_type}-object for X
#' @param Q number of clusters
#' @param epsilon=1e-3 convergence threshold for the algorithm
#' @param max_iter Maximal number of iterations for the algorithm
#' @param conv_threshold_m_step=1e-3 Convergence threshold for the fixed-point-solver of the m-step
#' @param max_iter_m_step=1e-3 Maximal number of iterations for the fixed-point-solver of the m-step
#' @return a \code{Var_Bayes_output} object all relevant informations + estimates
var_bayes_model_selection = function(sbm, Q=6, epsilon=1e-3, max_iter=5,
                                     conv_threshold_m_step=1e-3, 
                                     max_iter_m_step=5){
  if(class(sbm) == "SBM"){
    X = sbm$x
  } else{
    X = sbm
  }
  Q_dir = Q
  in_cluster = agnes(X, method = "ward")
  var_bayes_with_prefix = function(Q_dir){
    var_bayes(X, Q_dir, epsilon, in_cluster, max_iter,
              conv_threshold_m_step, max_iter_m_step)
  }
  var_bayes_res_0 = var_bayes_with_prefix(Q_dir)
  ICL_prot_1 = ICL(var_bayes_res_0)
  dir = 1
  var_bayes_res = var_bayes_with_prefix(Q_dir + dir)
  ICL_prot_2 = ICL(var_bayes_res)
  if (ICL_prot_2 <= ICL_prot_1){
    dir = -1
    var_bayes_res = var_bayes_with_prefix(Q_dir + dir)
    ICL_prot_2 = ICL(var_bayes_res) 
  }
  while(ICL_prot_2 > ICL_prot_1){
    var_bayes_res_0 = var_bayes_res
    Q_dir = Q_dir + dir
    ICL_prot_1 = ICL_prot_2
    var_bayes_res = var_bayes_with_prefix(Q_dir + dir)
    ICL_prot_2 = ICL(var_bayes_res)
  }
  var_bayes_res_0
}
