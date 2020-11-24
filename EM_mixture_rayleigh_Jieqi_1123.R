library(magrittr)
library(VGAM)
############################################################################################### Generate Data 
############################################################################################### Scenario setting 
Scenario = function(n, k, pi.TRUE, sigma.TRUE, seed){
  set.seed(seed)
  z.TRUE = sample(1:k, size = n, replace = T, prob = pi.TRUE)
  y = sapply(z.TRUE, function(z) rrayleigh(1, lambda.TRUE[z])) # specific for rayleigh
  
  set.seed(1)
  cluster.kmeans = kmeans(y, centers = k)
  init = list(pi = as.numeric(table(cluster.kmeans$cluster)/n)[order(cluster.kmeans$centers[,1])],
              lambda = as.numeric(sort(cluster.kmeans$centers[,1]))) # specific for rayleigh
  return(list(y=y,init=init,pi.TRUE=pi.TRUE,lambda.TRUE=lambda.TRUE))
}

Scenario_f = Scenario(n=10000, k=4, pi.TRUE=c(0.40, 0.30, 0.20, 0.10), lambda.TRUE = c(0.5, 10, 30, 60), seed = 1) # much much more sample size