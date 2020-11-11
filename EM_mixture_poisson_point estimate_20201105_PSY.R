library(magrittr)

# k-Mixture Poisson Model
k = 5                                   # number of mixture components

# Data Generation
n = 10000                               # sample size.          ; i = 1,2,...,n
pi.TRUE = c(0.2, 0.05, 0.35, 0.1, 0.3)  # mixture weights, p_j, ; j = 1,2,...,k
lambda.TRUE = c(1, 10, 30, 50, 70)      # lambda_j,             ; j = 1,2,...,k

set.seed(1)
z.TRUE = sample(1:k, size = n, replace = T, prob = pi.TRUE) # Indicator variable Z (= delta_ij)
y = sapply(z.TRUE, function(z) rpois(1, lambda.TRUE[z]))    # observations y_i, i = 1,2,...,n

plot(density(y))
##################################################

# EM Algorithm

# (1) Set initial values (by k-means method)
set.seed(1)
cluster.kmeans = kmeans(y, centers = k)
init = list(pi = as.numeric(table(cluster.kmeans$cluster)/n)[order(cluster.kmeans$centers[,1])],
            lambda = as.numeric(sort(cluster.kmeans$centers[,1])))

#Cluster means:
#1 54.85071
#2 73.05784
#3  2.55703
#4 34.66868
#5 24.97744

# as.numeric(table(cluster.kmeans$cluster))
# 1400 2386 2411 1986 1817

# as.numeric(table(cluster.kmeans$cluster)/n)
# 0.1400 0.2386 0.2411 0.1986 0.1817

# order(cluster.kmeans$centers[,1])
# 3 5 4 1 2

# pi
# 0.2411 0.1817 0.1986 0.1400 0.2386

# lambda
# 2.55703 24.97744 34.66868 54.85071 73.05784


# (2) Algorithm
mixture_poisson = function(y, init, tol = 1e-05, max.iter = 5000) {
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

  return(list(iter = iter, loglik = loglik, init = init, pi = pi_after, lambda = lambda_after))
}

# (3) Execute iteration
mixture_poisson(y, init)

