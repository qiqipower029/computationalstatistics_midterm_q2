library(magrittr)


# k-Mixture Exponential Model
k = 5

# Data Generation
n = 10000
pi.TRUE     = c(0.2, 0.05, 0.35, 0.1, 0.3)
lambda.TRUE = c(1, 1/10, 1/30, 1/50, 1/70)

k = 4
pi.TRUE     = c(0.2, 0.4, 0.3, 0.1)
lambda.TRUE = c(0.02, 10, 200, 500)

set.seed(1)
z.TRUE = sample(1:k, size = n, replace = T, prob = pi.TRUE)
y = sapply(z.TRUE, function(z) rexp(1, lambda.TRUE[z]))

plot(density(y))


##################################################
# EM Algorithm
# (1) Set initial values (by k-means)
set.seed(1)
cluster.kmeans = kmeans(y, centers = k)
init = list(pi = as.numeric(table(cluster.kmeans$cluster)/n)[order(cluster.kmeans$centers[,1])],
             lambda = as.numeric(1/sort(cluster.kmeans$centers[,1])))

mixture_exponential(y, init,  max.iter = 100000)


# (2) Algorithm
mixture_exponential <- function(y, init, tol = 1e-05, max.iter = 5000) {
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
  loglik = sum(w * log(t(replicate(n, pi_after)) * outer(y, lambda_after, dexp)))
  
  return(list(iter = iter, loglik = loglik, init = init, pi = pi_after, lambda = lambda_after))
}

# (3) Execute iteration
mixture_exponential(y, init,  max.iter = 100000)

