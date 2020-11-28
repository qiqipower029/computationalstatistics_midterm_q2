library(magrittr)
library(VGAM)
############################################################################################### Generate Data 
############################################################################################### Scenario setting 
Scenario = function(n, k, pi.TRUE, sigma.TRUE, seed){
  set.seed(seed)
  z.TRUE = sample(1:k, size = n, replace = T, prob = pi.TRUE)
  y = sapply(z.TRUE, function(z) rrayleigh(1, sigma.TRUE[z])) # specific for rayleigh
  
  set.seed(1)
  cluster.kmeans = kmeans(y, centers = k)
  init = list(pi = as.numeric(table(cluster.kmeans$cluster)/n)[order(cluster.kmeans$centers[,1])],
              sigma = as.numeric(sort(cluster.kmeans$centers[,1]))) # specific for rayleigh
  return(list(y=y,init=init,pi.TRUE=pi.TRUE,sigma.TRUE=sigma.TRUE))
}

Scenario_f = Scenario(n=10000, k=4, pi.TRUE=c(0.40, 0.30, 0.20, 0.10), sigma.TRUE = c(0.5, 10, 30, 60), seed = 1) # much much more sample size

mixture_rayleigh = function(y, init, tol = 1e-06, max.iter = 100000000000) {
  ptm = proc.time()
  n = length(y)
  k = length(init$pi)
  
  # E-step
  newpi    = init$pi
  newsigma = init$sigma
  
  # M-step
  pi = numeric(k)
  sigma = numeric(k)
  
  iter = 0
  
  while ((any(abs(sigma - newsigma)/abs(sigma) > tol) | any(abs(pi - newpi)/abs(pi) > tol)) & (iter<max.iter)) {
    # E-step
    pi = newpi
    sigma = newsigma
    w = t(replicate(n, pi)) * outer(y, sigma, drayleigh) / replicate(k, rowSums(t(replicate(n, pi)) * outer(y, sigma, drayleigh)))
    # M-step
    newpi = colSums(w)/n
    newsigma = sqrt(colSums(w*y^2)/2/colSums(w))
    
    iter = iter + 1
    
    cat('Iteration: ', iter, '\n')
  }
  
  # log likelihood
  w = t(replicate(n, newpi)) * outer(y, newsigma, drayleigh) / replicate(k, rowSums(t(replicate(n, newpi)) * outer(y, newsigma, drayleigh)))
  loglik = sum(w * log(t(replicate(n, newpi)) * outer(y, newsigma, drayleigh)), na.rm = T)
  
  ############################################################################################### Variance Estimate
  ############################################################################### Information Matrix for parameters
  w.matrix = w
  p.est = newpi
  sigma.est = newsigma
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
        psi.22[l,j] = sum((-2/sigma.est[l]+y^2/(sigma.est[l])^3)*(-2/sigma.est[j]+y^2/(sigma.est[j])^3)*(-w.matrix[,l]*w.matrix[,j]))
        psi.22[j,l] = psi.22[l,j]
      }
    }
    # diagonal entry
    psi.22[l,l] = sum((-2/sigma.est[l]+y^2/(sigma.est[l])^3)^2*w.matrix[,l]*(1-w.matrix[,l]))
  }
  
  # -----psi 12-----#
  psi.12 = matrix(nrow = k-1, ncol = k)
  for (l in 1:k-1) {
    # off-diagonal entry
    b = l+1
    for (j in b:k) {
      psi.12[l,j] = sum((-2/sigma.est[j]+y^2/(sigma.est[j])^3)*(-w.matrix[,l]*w.matrix[,j]/p.est[l]+w.matrix[,k]*w.matrix[,j]/p.est[k]))
      if(j < k) {psi.12[j,l] = psi.12[l,j]}
    }
    
    # diagonal entry
    psi.12[l,l] = sum((-2/sigma.est[l]+y^2/(sigma.est[l])^3)*(w.matrix[,l]*(1-w.matrix[,l])/p.est[l])+w.matrix[,l]*w.matrix[,k]/p.est[k])
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
    gamma.22[l,l] = sum((-2/(sigma.est[l])^2+3*y^2/(sigma.est[l])^4)*w.matrix[,l])
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
              esti_pi = round(p.est,4),
              true_pi = round(Scenario_f$pi.TRUE,4),
              init_sigma = round(init$sigma,4), 
              esti_sigma = round(sigma.est,4),
              true_sigma = round(Scenario_f$sigma.TRUE,4), 
              esti_var = variance.matrix))
}

mixture_rayleigh(y = Scenario_f$y, init= Scenario_f$init)
