#---------------------EM function for Rayleigh --------------------#
mixture_rayleigh = function(y, init, tol = 1e-06, max.iter = 100000000000) {
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
  v = log(t(replicate(n, newpi)) * outer(y, newsigma, drayleigh))
  loglik = sum(w * log(t(replicate(n, newpi)) * outer(y, newsigma, drayleigh)), na.rm = T)
  
  return(list(n = n, k = k, iter = iter, loglikelihood = loglik, init = init, pi = newpi, sigma = newsigma, max.iter = max.iter))
}
