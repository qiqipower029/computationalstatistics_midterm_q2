#-------------------Scenario for Exponential Distribution Point Estimation  Jieqi-------------------#



#-------------------Define EM function for exponential distribution-----------------------#
mixture_exponential <- function(y, init, tol = 1e-06, max.iter = 10000000) {
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
  
  return(list(iter = iter, loglik = loglik, init = init, pi = pi_after, lambda = lambda_after, iter = iter, n = n, k = k, max.iter = max.iter))
}

#---------Scenario setting----------#
n = 10000
k = 6
pi.TRUE     = c(0.166667, 0.166667, 0.166667, 0.166667, 0.166667, 0.166667)
lambda.TRUE = c(0.3, 0.8, 3, 8, 15, 25)

set.seed(1)
z.TRUE = sample(1:k, size = n, replace = T, prob = pi.TRUE)
y = sapply(z.TRUE, function(z) rexp(1, lambda.TRUE[z]))
set.seed(1)
cluster.kmeans = kmeans(y, centers = k)
init = list(pi = as.numeric(table(cluster.kmeans$cluster)/n)[order(cluster.kmeans$centers[,1])],
            lambda = as.numeric(1/sort(cluster.kmeans$centers[,1])))

#--------Execute the function----------#
ptm = proc.time()
mixture_exponential(y, init)
proc.time() - ptm


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
    
    # Find matrix w
    # method 1
    # w = outer(1:n, 1:k, Vectorize(function(i, j) (pi_before[j] * dpois(y[i], lambda_before[j]))/sum(pi_before * dpois(y[i], lambda_before))))
    # method 2
    # w = matrix(nrow = n, ncol = k)
    # for (i in 1:n) {
    #   for (j in 1:k) {
    #     w[i,j] = (pi_before[j] * dpois(y[i], lambda_before[j]))/sum(pi_before * dpois(y[i], lambda_before))
    #   }
    # }
    # method 3
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
