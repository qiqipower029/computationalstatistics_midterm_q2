#---------------------------EM function for Poisson--------------------------#
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
