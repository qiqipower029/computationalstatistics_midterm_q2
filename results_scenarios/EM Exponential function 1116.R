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

