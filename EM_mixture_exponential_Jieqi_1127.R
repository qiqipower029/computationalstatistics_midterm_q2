library(magrittr)
############################################################################################### Generate Data 
############################################################################################### Scenario setting 
Scenario = function(n, k, pi.TRUE, lambda.TRUE, seed){
  set.seed(seed)
  z.TRUE = sample(1:k, size = n, replace = T, prob = pi.TRUE)
  y = sapply(z.TRUE, function(z) rexp(1, lambda.TRUE[z])) # specific for POI
  
  set.seed(1)
  cluster.kmeans = kmeans(y, centers = k)
  init = list(pi = as.numeric(table(cluster.kmeans$cluster)/n)[order(cluster.kmeans$centers[,1])],
              lambda = as.numeric(sort(cluster.kmeans$centers[,1]))) # specific for POI
  return(list(y=y,init=init,pi.TRUE=pi.TRUE,lambda.TRUE=lambda.TRUE))
}

Scenario_f = Scenario(n=10000, k=3, pi.TRUE=c(0.40, 0.30, 0.3), lambda.TRUE = c(0.5, 10, 60), seed = 3) # much much more sample size




mixture_exponential <- function(y, init, tol = 1e-06, max.iter = 10000000000000) {
  ptm = proc.time()
  n = length(y)
  k = length(init$pi)
  
  # E-step
  pi_after      = init$pi
  lambda_after  = init$lambda
  
  # M-step
  pi_before     = numeric(k)
  lambda_before = numeric(k)
  
  iter = 0
  while ((any(abs(lambda_before - lambda_after)/(lambda_before + 1e-04) > tol) | any(abs(pi_before - pi_after)/(pi_before + 1e-04) > tol)) & iter < max.iter) {
    # E-step
    pi_before = pi_after
    lambda_before = lambda_after
    w = t(replicate(n, pi_before)) * outer(y, lambda_before, dexp) / replicate(k, rowSums(t(replicate(n, pi_before)) * outer(y, lambda_before, dexp)))
    
    # M-step
    pi_after = colSums(w)/n
    lambda_after = colSums(w)/colSums(w*y)
    iter <- iter + 1
    
    cat('Iteration: ', iter, '\n')
  }
  
  # log likelihood
  w = t(replicate(n, pi_after)) * outer(y, lambda_after, dexp) / replicate(k, rowSums(t(replicate(n, pi_after)) * outer(y, lambda_after, dexp)))
  loglik = sum(w * log(t(replicate(n, pi_after)) * outer(y, lambda_after, dexp)), na.rm = T)
  
  ############################################################################################### Variance Estimate
  ############################################################################### Information Matrix for parameters
  w.matrix = w
  p.est = pi_after
  lambda.est = lambda_after
  ############################################################################### PSI matrix
  # -----psi 11-----#
  psi.11 = matrix(nrow = k-1, ncol = k-1)
  for (l in 1:k-1) {
    # off-diagonal entry
    b = l+1
    if(b < k){
      for (j in b:k-1) {
        psi.11[l,j] = sum( -w.matrix[,l]*w.matrix[,j]/(p.est[l]*p.est[j])
                           +w.matrix[,l]*w.matrix[,k]/(p.est[l]*p.est[k])
                           +w.matrix[,k]*w.matrix[,j]/(p.est[k]*p.est[j])
                           +w.matrix[,k]*(1-w.matrix[,k])/(p.est[k])^2   )
        psi.11[j,l] = psi.11[l,j]
      }
    } 
    
    # diagonal entry
    psi.11[l,l] = sum(   w.matrix[,l]*(1-w.matrix[,l])/(p.est[l])^2 
                         + w.matrix[,k]*(1-w.matrix[,k])/(p.est[k])^2 
                         + 2*w.matrix[,l]*w.matrix[,k]/(p.est[l]*p.est[k])  )
  }
  
  # -----psi 22-----#
  psi.22 = matrix(nrow = k, ncol = k)
  for (l in 1:k) {
    # off-diagonal entry
    b = l+1
    if(b <= k){
      for (j in b:k) {
        psi.22[l,j] = sum((1/lambda.est[l]-y)*(1/lambda.est[j]-y)*(-w.matrix[,l]*w.matrix[,j]))
        psi.22[j,l] = psi.22[l,j]
      }
    }
    # diagonal entry
    psi.22[l,l] = sum((1/lambda.est[l]-y)^2*w.matrix[,l]*(1-w.matrix[,l]))
  }
  
  # -----psi 12-----#
  psi.12 = matrix(nrow = k-1, ncol = k)
  for (l in 1:k-1) {
    # off-diagonal entry
    b = l+1
    for (j in b:k) {
      psi.12[l,j] = sum((1/lambda.est[j]-y)*(-w.matrix[,l]*w.matrix[,j]/p.est[l] + w.matrix[,k]*w.matrix[,j]/p.est[k]))
      if(j < k) {psi.12[j,l] = psi.12[l,j]}
    }
    
    # diagonal entry
    psi.12[l,l] = sum((1/lambda.est[l]-y)*(w.matrix[,l]*(1-w.matrix[,l])/p.est[l] + w.matrix[,l]*w.matrix[,k]/p.est[k]))
  }
  psi.21 = t(psi.12)
  
  # Generate psi matrix
  psi.upper = cbind(psi.11, psi.12);psi.lower = cbind(psi.21, psi.22); psi = rbind(psi.upper, psi.lower)
  
  ############################################################################### GAMMA matrix
  # -----gamma 11-----#
  gamma.11 = matrix(nrow = k-1, ncol = k-1)
  for (l in 1:k-1) {
    # off-diagonal entry
    b = l+1
    if(b <= k-1){
      for (j in b:k-1) {
        gamma.11[l,j] = sum(w.matrix[,k]/(p.est[k])^2)
        gamma.11[j,l] = gamma.11[l,j]
      }
    }
    # diagonal entry
    gamma.11[l,l] = sum(w.matrix[,l]/(p.est[l]^2)+w.matrix[,k]/(p.est[k])^2)
  }
  
  # -----gamma 22-----#
  gamma.22 = matrix(nrow = k, ncol = k)
  for (l in 1:k) {
    # off-diagonal entry
    b = l+1
    if( b<= k){
      for (j in b:k) {
        gamma.22[l,j] = 0
        gamma.22[j,l] = 0
      }
    }
    # diagonal entry
    gamma.22[l,l] = sum(1/(lambda.est[l])^2*w.matrix[,l])
  }
  
  # -----gamma 12-----#
  gamma.12 = matrix(0, nrow = k-1, ncol = k)
  gamma.21 = t(gamma.12)
  gamma.upper = cbind(gamma.11, gamma.12);gamma.lower = cbind(gamma.21, gamma.22); gamma = rbind(gamma.upper, gamma.lower)
  
  ############################################################################### Generate the final variance matrix
  variance.nopk = solve(-psi + gamma)
  
  # Add elements for p_k
  variance.matrix = matrix(0, ncol = 2*k, nrow = 2*k)
  variance.matrix[-k,-k] = variance.nopk
  vec1 = matrix(1, nrow = 1, ncol = k-1)
  var_pp = variance.nopk[1:k-1, 1:k-1]
  var_plambda = variance.nopk[1:k-1, k:(2*k-1)]
  variance.matrix[k,k] = vec1 %*% var_pp %*% t(vec1)
  
  for (j in 1:k-1) {
    variance.matrix[k, j] = sum(var_pp[,j])
    variance.matrix[j, k] = variance.matrix[k, j]
  }
  
  for (j in 1:k) {
    variance.matrix[k, k+j] = sum(var_plambda[,j])
    variance.matrix[k+j, k] = variance.matrix[k, k+j]
  }
  
  ############################################################################### Return final estimate
  return(list(comp_time = as.vector((proc.time() - ptm)[3]),
              sample_size = n,
              iter = iter, 
              loglik = loglik, 
              init_pi = round(init$pi,4), 
              esti_pi = round(pi_after,4),
              true_pi = round(Scenario_f$pi.TRUE,4),
              init_lambda = round(init$lambda,4), 
              esti_lambda = round(lambda_after,4),
              true_lambda = round(Scenario_f$lambda.TRUE,4), 
              esti_var = variance.matrix))
}

mixture_exponential(y = Scenario_f$y, init= Scenario_f$init)
