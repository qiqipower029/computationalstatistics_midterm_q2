library(magrittr)
library(VGAM)


# k-Mixture Rayleigh Model
k = 5

# Data Generation
n = 10000
pi.TRUE = c(0.2, 0.05, 0.35, 0.1, 0.3)
sigma.TRUE = c(1, 5, 10, 7, 20)

set.seed(1)
z.TRUE = sample(1:k, size = n, replace = T, prob = pi.TRUE)
y = sapply(z.TRUE, function(z) rrayleigh(1, sigma.TRUE[z]))

plot(density(y))


##################################################


# EM Algorithm

# (1) Set initial values (by k-means)
set.seed(1)
cluster.kmeans = kmeans(y, centers = k)
init = list(pi = as.numeric(table(cluster.kmeans$cluster)/n)[order(cluster.kmeans$centers[,1])],
             sigma = as.numeric(sort(cluster.kmeans$centers[,1]*sqrt(2/pi))))

# (2) Algorithm
mixture_rayleigh = function(y, init, tol = 1e-05, max.iter = 5000) {
  n = length(y)
  k = length(init$pi)
  
  # E-step
  pi_after    = init$pi
  sigma_after = init$sigma
  
  # M-step
  pi_before    = numeric(k)
  sigma_before = numeric(k)
  
  iter = 0
  while ((any(abs(sigma_before - sigma_after)/(sigma_before + 1e-04) > tol) | any(abs(pi_before - pi_after)/(pi_before + 1e-04) > tol)) & iter < max.iter) {
    # E-step
    pi_before    = pi_after
    sigma_before = sigma_after
    w = t(replicate(n, pi_before)) * outer(y, sigma_before, drayleigh) / replicate(k, rowSums(t(replicate(n, pi_before)) * outer(y, sigma_before, drayleigh)))
    
    # M-step
    pi_after = colSums(w)/n
    sigma_after = sqrt(colSums(w*y^2)/2/colSums(w))
    
    iter = iter + 1
    cat('Iteration: ', iter, '\n')
  }
  
  # log likelihood
  w = t(replicate(n, pi_after)) * outer(y, sigma_after, drayleigh) / replicate(k, rowSums(t(replicate(n, pi_after)) * outer(y, sigma_after, drayleigh)))
  loglik = sum(w * log(t(replicate(n, pi_after)) * outer(y, sigma_after, drayleigh)))
  
  return(list(iter = iter, loglik = loglik, init = init, pi = pi_after, lambda = lambda_after))
}

# (3) Execute iteration
mixture_rayleigh(y, init)

