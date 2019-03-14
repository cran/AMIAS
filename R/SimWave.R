# simulate wave
# n determines length
# sigma is residual sd
# output is data (y),  input location (x), mean function (y0), location of knots in x scale (tau) and in sample scale (A)
SimuWave <- function (n, sigma = 0.1, seed=NA) 
{
  if (!is.na(seed)) set.seed(seed)
  x = seq(1/n, 1,length.out = n)
  tau = c(256, 512, 768, 1024, 1152, 1280, 1344, 1408)/1472
  nknot = length(tau)
  A = sapply(tau, function(z) which(x>=z)[1])
  tau=x[A]
  tau1 = c(0, tau, 1)
  h =  cumsum(c(1, (-1)^(1:nknot)*(1:nknot)))
  phi = rep(1, n)
  for (j in 1:(nknot+1)) phi = cbind(phi,pmin(pmax(x-tau1[j], 0), tau1[j+1]-tau1[j]))
  y0 = as.vector(phi%*%c(0, h))
  y <- y0 + sigma*rnorm(n)
  return(list(y = y, x = x, y0 = y0, tau = tau, SetA = A))  
}
