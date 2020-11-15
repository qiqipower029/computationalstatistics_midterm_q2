library(magrittr)


# k-Mixture Exponential Model
k = 2

# Data Generation
n = 10000
pi.TRUE     = c(0.2, 0.8)
lambda.TRUE = c(0.1, 10)

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

# (2) Algorithm
mixture_exponential <- function(y, init, tol = 1e-06, max.iter = 1000000) {
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
  
  return(list(iter = iter, loglik = loglik, init = init, pi = pi_after, lambda = lambda_after, iter = iter, n = n, k = k, max.iter = max.iter))
}

# (3) Execute iteration

#----------------------k=2, n=10000-------------------------#
re1 = unlist(mixture_exponential(y, init,  max.iter = 1))
re2 = unlist(mixture_exponential(y, init,  max.iter = 3))
re3 = unlist(mixture_exponential(y, init,  max.iter = 100))

#----------------------k=2, n=1000-------------------------#
n = 1000
pi.TRUE     = c(0.2, 0.8)
lambda.TRUE = c(0.1, 10)

set.seed(1)
z.TRUE = sample(1:k, size = n, replace = T, prob = pi.TRUE)
y = sapply(z.TRUE, function(z) rexp(1, lambda.TRUE[z]))
set.seed(1)
cluster.kmeans = kmeans(y, centers = k)
init = list(pi = as.numeric(table(cluster.kmeans$cluster)/n)[order(cluster.kmeans$centers[,1])],
            lambda = as.numeric(1/sort(cluster.kmeans$centers[,1])))
re4 = unlist(mixture_exponential(y, init,  max.iter = 1))
re5 = unlist(mixture_exponential(y, init,  max.iter = 3))
re6 = unlist(mixture_exponential(y, init,  max.iter = 100))

k2 = rbind(re1, re2, re3, re4, re5, re6)


#----------------------k=3, n=10000-------------------------#
k = 3
n = 10000
pi.TRUE     = c(0.2, 0.4, 0.4)
lambda.TRUE = c(0.1, 5, 20)

set.seed(1)
z.TRUE = sample(1:k, size = n, replace = T, prob = pi.TRUE)
y = sapply(z.TRUE, function(z) rexp(1, lambda.TRUE[z]))
set.seed(1)
cluster.kmeans = kmeans(y, centers = k)
init = list(pi = as.numeric(table(cluster.kmeans$cluster)/n)[order(cluster.kmeans$centers[,1])],
            lambda = as.numeric(1/sort(cluster.kmeans$centers[,1])))
re7 = unlist(mixture_exponential(y, init,  max.iter = 5))
re8 = unlist(mixture_exponential(y, init,  max.iter = 30))
re9 = unlist(mixture_exponential(y, init,  max.iter = 1000))

#----------------------k=3, n=1000-------------------------#
k = 3
n = 1000
pi.TRUE     = c(0.2, 0.4, 0.4)
lambda.TRUE = c(0.1, 5, 20)

set.seed(1)
z.TRUE = sample(1:k, size = n, replace = T, prob = pi.TRUE)
y = sapply(z.TRUE, function(z) rexp(1, lambda.TRUE[z]))
set.seed(1)
cluster.kmeans = kmeans(y, centers = k)
init = list(pi = as.numeric(table(cluster.kmeans$cluster)/n)[order(cluster.kmeans$centers[,1])],
            lambda = as.numeric(1/sort(cluster.kmeans$centers[,1])))
re10 = unlist(mixture_exponential(y, init,  max.iter = 5))
re11 = unlist(mixture_exponential(y, init,  max.iter = 30))
re12 = unlist(mixture_exponential(y, init,  max.iter = 1000))

k3 = rbind(re7, re8, re9, re10, re11, re12)


#----------------------k=4, n=10000-------------------------#
k = 4
n = 10000
pi.TRUE     = c(0.2, 0.1, 0.3, 0.4)
lambda.TRUE = c(0.05, 0.5, 5, 50)

set.seed(1)
z.TRUE = sample(1:k, size = n, replace = T, prob = pi.TRUE)
y = sapply(z.TRUE, function(z) rexp(1, lambda.TRUE[z]))
set.seed(1)
cluster.kmeans = kmeans(y, centers = k)
init = list(pi = as.numeric(table(cluster.kmeans$cluster)/n)[order(cluster.kmeans$centers[,1])],
            lambda = as.numeric(1/sort(cluster.kmeans$centers[,1])))
re13 = unlist(mixture_exponential(y, init,  max.iter = 1))
re14 = unlist(mixture_exponential(y, init,  max.iter = 5))

#----------------------k=4, n=1000-------------------------#
k = 4
n = 1000
pi.TRUE     = c(0.2, 0.1, 0.3, 0.4)
lambda.TRUE = c(0.05, 0.5, 5, 50)

set.seed(1)
z.TRUE = sample(1:k, size = n, replace = T, prob = pi.TRUE)
y = sapply(z.TRUE, function(z) rexp(1, lambda.TRUE[z]))
set.seed(1)
cluster.kmeans = kmeans(y, centers = k)
init = list(pi = as.numeric(table(cluster.kmeans$cluster)/n)[order(cluster.kmeans$centers[,1])],
            lambda = as.numeric(1/sort(cluster.kmeans$centers[,1])))
re15 = unlist(mixture_exponential(y, init,  max.iter = 1))
re16 = unlist(mixture_exponential(y, init,  max.iter = 5))

k4 = rbind(re13, re14, re15, re16)


# Then we want to see what will happen if we remove NaN
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

#----------------------k=2, n=10000-------------------------#
n = 10000
k = 2
pi.TRUE     = c(0.2, 0.8)
lambda.TRUE = c(0.1, 10)

set.seed(1)
z.TRUE = sample(1:k, size = n, replace = T, prob = pi.TRUE)
y = sapply(z.TRUE, function(z) rexp(1, lambda.TRUE[z]))
set.seed(1)
cluster.kmeans = kmeans(y, centers = k)
init = list(pi = as.numeric(table(cluster.kmeans$cluster)/n)[order(cluster.kmeans$centers[,1])],
            lambda = as.numeric(1/sort(cluster.kmeans$centers[,1])))
re17 = unlist(mixture_exponential(y, init,  max.iter = 1))
re18 = unlist(mixture_exponential(y, init,  max.iter = 3))
re19 = unlist(mixture_exponential(y, init,  max.iter = 100))

#----------------------k=2, n=1000-------------------------#
n = 1000
pi.TRUE     = c(0.2, 0.8)
lambda.TRUE = c(0.1, 10)

set.seed(1)
z.TRUE = sample(1:k, size = n, replace = T, prob = pi.TRUE)
y = sapply(z.TRUE, function(z) rexp(1, lambda.TRUE[z]))
set.seed(1)
cluster.kmeans = kmeans(y, centers = k)
init = list(pi = as.numeric(table(cluster.kmeans$cluster)/n)[order(cluster.kmeans$centers[,1])],
            lambda = as.numeric(1/sort(cluster.kmeans$centers[,1])))
re20 = unlist(mixture_exponential(y, init,  max.iter = 1))
re21 = unlist(mixture_exponential(y, init,  max.iter = 3))
re22 = unlist(mixture_exponential(y, init,  max.iter = 100))

k2_removena = rbind(re17, re18, re19, re20, re21, re22)


#----------------------k=3, n=10000-------------------------#
k = 3
n = 10000
pi.TRUE     = c(0.2, 0.4, 0.4)
lambda.TRUE = c(0.1, 5, 20)

set.seed(1)
z.TRUE = sample(1:k, size = n, replace = T, prob = pi.TRUE)
y = sapply(z.TRUE, function(z) rexp(1, lambda.TRUE[z]))
set.seed(1)
cluster.kmeans = kmeans(y, centers = k)
init = list(pi = as.numeric(table(cluster.kmeans$cluster)/n)[order(cluster.kmeans$centers[,1])],
            lambda = as.numeric(1/sort(cluster.kmeans$centers[,1])))
re23 = unlist(mixture_exponential(y, init,  max.iter = 5))
re24 = unlist(mixture_exponential(y, init,  max.iter = 30))
re25 = unlist(mixture_exponential(y, init,  max.iter = 1000))

#----------------------k=3, n=1000-------------------------#
k = 3
n = 1000
pi.TRUE     = c(0.2, 0.4, 0.4)
lambda.TRUE = c(0.1, 5, 20)

set.seed(1)
z.TRUE = sample(1:k, size = n, replace = T, prob = pi.TRUE)
y = sapply(z.TRUE, function(z) rexp(1, lambda.TRUE[z]))
set.seed(1)
cluster.kmeans = kmeans(y, centers = k)
init = list(pi = as.numeric(table(cluster.kmeans$cluster)/n)[order(cluster.kmeans$centers[,1])],
            lambda = as.numeric(1/sort(cluster.kmeans$centers[,1])))
re26 = unlist(mixture_exponential(y, init,  max.iter = 5))
re27 = unlist(mixture_exponential(y, init,  max.iter = 30))
re28 = unlist(mixture_exponential(y, init,  max.iter = 1000))

k3_removena = rbind(re23, re24, re25, re26, re27, re28)


#----------------------k=4, n=10000-------------------------#
k = 4
n = 10000
pi.TRUE     = c(0.2, 0.1, 0.3, 0.4)
lambda.TRUE = c(0.05, 0.5, 5, 50)

set.seed(1)
z.TRUE = sample(1:k, size = n, replace = T, prob = pi.TRUE)
y = sapply(z.TRUE, function(z) rexp(1, lambda.TRUE[z]))
set.seed(1)
cluster.kmeans = kmeans(y, centers = k)
init = list(pi = as.numeric(table(cluster.kmeans$cluster)/n)[order(cluster.kmeans$centers[,1])],
            lambda = as.numeric(1/sort(cluster.kmeans$centers[,1])))
re29 = unlist(mixture_exponential(y, init,  max.iter = 1))
re30 = unlist(mixture_exponential(y, init,  max.iter = 5))

#----------------------k=4, n=1000-------------------------#
k = 4
n = 1000
pi.TRUE     = c(0.2, 0.1, 0.3, 0.4)
lambda.TRUE = c(0.05, 0.5, 5, 50)

set.seed(1)
z.TRUE = sample(1:k, size = n, replace = T, prob = pi.TRUE)
y = sapply(z.TRUE, function(z) rexp(1, lambda.TRUE[z]))
set.seed(1)
cluster.kmeans = kmeans(y, centers = k)
init = list(pi = as.numeric(table(cluster.kmeans$cluster)/n)[order(cluster.kmeans$centers[,1])],
            lambda = as.numeric(1/sort(cluster.kmeans$centers[,1])))
re31 = unlist(mixture_exponential(y, init,  max.iter = 1))
re32 = unlist(mixture_exponential(y, init,  max.iter = 5))

k4_removena = rbind(re29, re30, re31, re32)

library(openxlsx)
write.xlsx(k2, file = "expk2_1111_Jieqi.xlsx")
write.xlsx(k3, file = "expk3_1111_Jieqi.xlsx")
write.xlsx(k4, file = "expk4_1111_Jieqi.xlsx")
write.xlsx(k2_removena, file = "expk2rm_1111_Jieqi.xlsx")
write.xlsx(k3_removena, file = "expk3rm_1111_Jieqi.xlsx")
write.xlsx(k4_removena, file = "expk4rm_1111_Jieqi.xlsx")
