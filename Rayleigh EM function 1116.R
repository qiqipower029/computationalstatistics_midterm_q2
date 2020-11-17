library(magrittr)
library(VGAM)


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

#---------Scenario setting----------#
n = 10000
k = 6
pi.TRUE     = c(0.166667, 0.166667, 0.166667, 0.166667, 0.166667, 0.166667)
sigma.TRUE = c(0.3, 8, 30, 80, 150, 250)

set.seed(1)
z.TRUE = sample(1:k, size = n, replace = T, prob = pi.TRUE)
y = sapply(z.TRUE, function(z) rexp(1, lambda.TRUE[z]))
set.seed(1)
cluster.kmeans = kmeans(y, centers = k)
init = list(pi = as.numeric(table(cluster.kmeans$cluster)/n)[order(cluster.kmeans$centers[,1])],
            sigma = as.numeric(1/sort(cluster.kmeans$centers[,1])))

#--------Execute the function----------#
ptm = proc.time()
mixture_rayleigh(y, init)
proc.time() - ptm