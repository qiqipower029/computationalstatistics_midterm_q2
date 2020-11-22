#-----------------Variance Estimation-----------------#
#----------------------Jieqi--------------------------#

# First, get the point estimates and the W matrix
library(magrittr)
k = 4
# Data Generation
n = 10000                    
pi.TRUE = c(0.2, 0.2, 0.3, 0.3) 
lambda.TRUE = c(1, 10, 20, 28)   
set.seed(1)
z.TRUE = sample(1:k, size = n, replace = T, prob = pi.TRUE) 
y = sapply(z.TRUE, function(z) rpois(1, lambda.TRUE[z])) 
set.seed(1)
cluster.kmeans = kmeans(y, centers = k)
init = list(pi = as.numeric(table(cluster.kmeans$cluster)/n)[order(cluster.kmeans$centers[,1])],
            lambda = as.numeric(sort(cluster.kmeans$centers[,1])))
mixture_poisson = function(y, init, tol = 1e-06, max.iter = 100000000) {
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

    w = t(replicate(n, pi_before)) * outer(y, lambda_before, dpois) / replicate(k, rowSums(   t(replicate(n, pi_before))*outer(y, lambda_before, dpois)  )  )
    
    # M-step
    pi_after = colSums(w)/n
    lambda_after = colSums(w*y)/colSums(w)
    
    iter = iter + 1
    cat('Iteration: ', iter, '\n')
  }
  
  # log likelihood
  w = t(replicate(n, pi_after)) * outer(y, lambda_after, dpois) / replicate(k, rowSums(t(replicate(n, pi_after)) * outer(y, lambda_after, dpois)))
  loglik = sum(w * log(t(replicate(n, pi_after)) * outer(y, lambda_after, dpois)))
  
  return(list(iter = iter, loglik = loglik, init = init, pi = pi_after, lambda = lambda_after, wp = w))
}
result = mixture_poisson(y, init)
w.matrix = result$wp
p.est = result$pi
lambda.est = result$lambda

#--------Variance function-------------#
variance.function.poisson = function(w.matrix, p.est, lambda.est, y) {
  n = nrow(w.matrix)
  k = ncol(w.matrix)
  # -----psi 11-----#
  psi.11 = matrix(nrow = k-1, ncol = k-1)
  for (l in 1:k-1) {
    
    b = l+1
    if(b<k){
    # off-diagonal entry
      for (j in b:k-1) {
        psi.11[l,j] = sum(-w.matrix[,l]*w.matrix[,j]/(p.est[l]*p.est[j])
                        +w.matrix[,l]*w.matrix[,k]/(p.est[l]*p.est[k])
                        +w.matrix[,k]*w.matrix[,j]/(p.est[k]*p.est[j])
                        +w.matrix[,k]*(1-w.matrix[,k])/(p.est[k])^2)
        psi.11[j,l] = psi.11[l,j]
      }
    } 
    # diagonal entry
    psi.11[l,l] = sum(w.matrix[,l]*(1-w.matrix[,l])/(p.est[l])^2 
                      + w.matrix[,k]*(1-w.matrix[,k])/(p.est[k])^2 
                      + 2*w.matrix[,l]*w.matrix[,k]/(p.est[l]*p.est[k]))
  }
  
  # -----psi 22-----#
  psi.22 = matrix(nrow = k, ncol = k)
  for (l in 1:k) {
    b = l+1
    if(b<=k){
      # off-diagonal entry
      for (j in b:k) {
        psi.22[l,j] = sum((y/lambda.est[l]-1)*(y/lambda.est[j]-1)*(-w.matrix[,l]*w.matrix[,j]))
        psi.22[j,l] = psi.22[l,j]
      }
    }
    # diagonal entry
    psi.22[l,l] = sum((y/lambda.est[l]-1)^2*w.matrix[,l]*(1-w.matrix[,l]))
  }
  
  # -----psi 12-----#
  psi.12 = matrix(nrow = k-1, ncol = k)
  for (l in 1:k-1) {
    b = l+1
  
      # off-diagonal entry
      for (j in b:k) {
        psi.12[l,j] = sum((y/lambda.est[j]-1)*(-w.matrix[,l]*w.matrix[,j]/p.est[l]+w.matrix[,k]*w.matrix[,j]/p.est[k]))
        if(j<k) {psi.12[j,l] = psi.12[l,j]}
      }
    
    # diagonal entry
    psi.12[l,l] = sum((y/lambda.est[l]-1)*(w.matrix[,l]*(1-w.matrix[,l])/p.est[l]+w.matrix[,l]*w.matrix[,k]/p.est[k]))
  }
  psi.21 = t(psi.12)
  
  # Generate psi matrix
  psi.upper = cbind(psi.11, psi.12);psi.lower = cbind(psi.21, psi.22); psi = rbind(psi.upper, psi.lower)
  
  # -----gamma 11-----#
  gamma.11 = matrix(nrow = k-1, ncol = k-1)
  for (l in 1:k-1) {
    b = l+1
    if(b<=k-1){
      # off-diagonal entry
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
    b = l+1
    if(b<=k){
      # off-diagonal entry
      for (j in b:k) {
        gamma.22[l,j] = 0
        gamma.22[j,l] = 0
      }
    }
    # diagonal entry
    gamma.22[l,l] = sum(y/(lambda.est[l])^2*w.matrix[,l])
  }
  
  # -----gamma 12-----#
  gamma.12 = matrix(0, nrow = k-1, ncol = k)
  gamma.21 = t(gamma.12)
  gamma.upper = cbind(gamma.11, gamma.12);gamma.lower = cbind(gamma.21, gamma.22); gamma = rbind(gamma.upper, gamma.lower)
  
  # Generate the final variance matrix
  variance.nopk = solve(-psi + gamma)
  
  # Add elements for p_k
  variance.matrix = matrix(0, ncol = 2*k, nrow = 2*k)
  variance.matrix[-k,-k] = variance.nopk
  vec1 = matrix(1, nrow = 1, ncol = k-1)
  var_pp = variance.nopk[1:k-1, 1:k-1]
  var_plambda = variance.nopk[1:k-1, k:(2*k-1)]
  variance.matrix[k,k] = vec1 %*% var_pk_1 %*% t(vec1)
  
  for (j in 1:k-1) {
    variance.matrix[k, j] = sum(var_pp[,j])
    variance.matrix[j, k] = variance.matrix[k, j]
  }
  
  for (j in 1:k) {
    variance.matrix[k, k+j] = sum(var_plambda[,j])
    variance.matrix[k+j, k] = variance.matrix[k, k+j]
  }
  
  return(variance.matrix)
}

