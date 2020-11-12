library(magrittr)
library(VGAM)


# k-Mixture Rayleigh Model
k = 3

# Data Generation
n = 10000
pi.TRUE = c(0.3, 0.4, 0.3)
sigma.TRUE = c(1, 10, 30)

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
mixture_rayleigh = function(y, init, tol = 1e-05, max.iter = 1000) {
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
  loglik = sum(w * log(t(replicate(n, newpi)) * outer(y, newsigma, drayleigh)))
  
  return(list(n = n, k = k, iter = iter, loglikelihood = loglik, init = init, pi = newpi, sigma = newsigma))
}

# (3) Execute iteration
mixture_rayleigh(y, init, max.iter = 30)
result$loglikelihood
v = result$v
w = result$w
product = w * v
sum(product[,1], na.rm = T)
sum(product[,1], na.rm = T) + sum(product[,2], na.rm = T)
output = cbind(v, w, n, )
# Tried k = 2, 
x1 = c(rep(1, 900)) %>% as.numeric()
x2 = c(rep(2, 900)) %>% as.numeric()
x3 = data.frame(x1, x2)
sum(x3[1]) + sum(x3[2])


# ------------------------------ Generate output--------------------------------#
# k = 2, n = 10000
k = 2
n = 10000
pi.TRUE = c(0.3, 0.7)
sigma.TRUE = c(1, 10)
set.seed(1)
z.TRUE = sample(1:k, size = n, replace = T, prob = pi.TRUE)
y = sapply(z.TRUE, function(z) rrayleigh(1, sigma.TRUE[z]))
plot(density(y))
cluster.kmeans = kmeans(y, centers = k)
init = list(pi = as.numeric(table(cluster.kmeans$cluster)/n)[order(cluster.kmeans$centers[,1])],
            sigma = as.numeric(sort(cluster.kmeans$centers[,1]*sqrt(2/pi))))
result1 = unlist(mixture_rayleigh(y, init, max.iter = 30))

#-------------k = 2, n = 10000--------------#
k = 2
n = 10000
pi.TRUE = c(0.3, 0.7)
sigma.TRUE = c(1, 10)
set.seed(1)
z.TRUE = sample(1:k, size = n, replace = T, prob = pi.TRUE)
y = sapply(z.TRUE, function(z) rrayleigh(1, sigma.TRUE[z]))
plot(density(y))
cluster.kmeans = kmeans(y, centers = k)
init = list(pi = as.numeric(table(cluster.kmeans$cluster)/n)[order(cluster.kmeans$centers[,1])],
            sigma = as.numeric(sort(cluster.kmeans$centers[,1]*sqrt(2/pi))))
result1 = unlist(mixture_rayleigh(y, init, max.iter = 15))
result2 = unlist(mixture_rayleigh(y, init, max.iter = 50))

#-------------k = 2, n = 1000--------------#
k = 2
n = 1000
pi.TRUE = c(0.3, 0.7)
sigma.TRUE = c(1, 10)
set.seed(1)
z.TRUE = sample(1:k, size = n, replace = T, prob = pi.TRUE)
y = sapply(z.TRUE, function(z) rrayleigh(1, sigma.TRUE[z]))
plot(density(y))
cluster.kmeans = kmeans(y, centers = k)
init = list(pi = as.numeric(table(cluster.kmeans$cluster)/n)[order(cluster.kmeans$centers[,1])],
            sigma = as.numeric(sort(cluster.kmeans$centers[,1]*sqrt(2/pi))))
result3 = unlist(mixture_rayleigh(y, init, max.iter = 15))
result4 = unlist(mixture_rayleigh(y, init, max.iter = 50))
result5 = unlist(mixture_rayleigh(y, init, max.iter = 5000))

#-------------k = 3, n = 10000--------------#
k = 3
n = 10000
pi.TRUE = c(0.3, 0.4, 0.3)
sigma.TRUE = c(1, 10, 20)
set.seed(1)
z.TRUE = sample(1:k, size = n, replace = T, prob = pi.TRUE)
y = sapply(z.TRUE, function(z) rrayleigh(1, sigma.TRUE[z]))
plot(density(y))
cluster.kmeans = kmeans(y, centers = k)
init = list(pi = as.numeric(table(cluster.kmeans$cluster)/n)[order(cluster.kmeans$centers[,1])],
            sigma = as.numeric(sort(cluster.kmeans$centers[,1]*sqrt(2/pi))))
result6 = unlist(mixture_rayleigh(y, init, max.iter = 15))
result7 = unlist(mixture_rayleigh(y, init, max.iter = 50))
result8 = unlist(mixture_rayleigh(y, init, max.iter = 5000))

#-------------k = 3, n = 1000 --------------#
n = 1000
result9 = unlist(mixture_rayleigh(y, init, max.iter = 15))
result10 = unlist(mixture_rayleigh(y, init, max.iter = 50))
result11 = unlist(mixture_rayleigh(y, init, max.iter = 5000))

k2 = rbind(result1, result2, result3, result4, result5)
k3 = rbind(result6, result7, result8, result9, result10, result11)

#------------if we remove NaN---------------#
mixture_rayleigh = function(y, init, tol = 1e-05, max.iter = 1000) {
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
  
  return(list(n = n, k = k, iter = iter, loglikelihood = loglik, init = init, pi = newpi, sigma = newsigma))
}

# k = 2, n = 10000
k = 2
n = 10000
pi.TRUE = c(0.3, 0.7)
sigma.TRUE = c(1, 10)
set.seed(1)
z.TRUE = sample(1:k, size = n, replace = T, prob = pi.TRUE)
y = sapply(z.TRUE, function(z) rrayleigh(1, sigma.TRUE[z]))
plot(density(y))
cluster.kmeans = kmeans(y, centers = k)
init = list(pi = as.numeric(table(cluster.kmeans$cluster)/n)[order(cluster.kmeans$centers[,1])],
            sigma = as.numeric(sort(cluster.kmeans$centers[,1]*sqrt(2/pi))))
result12 = unlist(mixture_rayleigh(y, init, max.iter = 30))
result13 = unlist(mixture_rayleigh(y, init, max.iter = 50))

#-------------k = 2, n = 1000--------------#
k = 2
n = 1000
pi.TRUE = c(0.3, 0.7)
sigma.TRUE = c(1, 10)
set.seed(1)
z.TRUE = sample(1:k, size = n, replace = T, prob = pi.TRUE)
y = sapply(z.TRUE, function(z) rrayleigh(1, sigma.TRUE[z]))
plot(density(y))
cluster.kmeans = kmeans(y, centers = k)
init = list(pi = as.numeric(table(cluster.kmeans$cluster)/n)[order(cluster.kmeans$centers[,1])],
            sigma = as.numeric(sort(cluster.kmeans$centers[,1]*sqrt(2/pi))))
result14 = unlist(mixture_rayleigh(y, init, max.iter = 15))
result15 = unlist(mixture_rayleigh(y, init, max.iter = 50))
result16 = unlist(mixture_rayleigh(y, init, max.iter = 5000))

#-------------k = 3, n = 10000--------------#
k = 3
n = 10000
pi.TRUE = c(0.3, 0.4, 0.3)
sigma.TRUE = c(1, 10, 20)
set.seed(1)
z.TRUE = sample(1:k, size = n, replace = T, prob = pi.TRUE)
y = sapply(z.TRUE, function(z) rrayleigh(1, sigma.TRUE[z]))
plot(density(y))
cluster.kmeans = kmeans(y, centers = k)
init = list(pi = as.numeric(table(cluster.kmeans$cluster)/n)[order(cluster.kmeans$centers[,1])],
            sigma = as.numeric(sort(cluster.kmeans$centers[,1]*sqrt(2/pi))))
result17 = unlist(mixture_rayleigh(y, init, max.iter = 15))
result18 = unlist(mixture_rayleigh(y, init, max.iter = 50))
result19 = unlist(mixture_rayleigh(y, init, max.iter = 5000))

#-------------k = 3, n = 1000 --------------#
n = 1000
result20 = unlist(mixture_rayleigh(y, init, max.iter = 15))
result21 = unlist(mixture_rayleigh(y, init, max.iter = 50))
result22 = unlist(mixture_rayleigh(y, init, max.iter = 5000))

k2_removena = rbind(result12, result13, result14, result15, result16)
k3_removena = rbind(result17, result18, result19, result20, result21, result22)

library(openxlsx)
write.xlsx(k2, file = "rayleighk2.xlsx")
write.xlsx(k3, file = "rayleighk3.xlsx")
write.xlsx(k2_removena, file = "rayleighk2rm.xlsx")
write.xlsx(k3_removena, file = "rayleighk3rm.xlsx")
