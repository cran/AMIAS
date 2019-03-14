
# simulate blocks
# n determines length
# sigma is residual sd
# output is data (y),  input location (x), mean function (y0), location of knots in x scale (tau) and in sample scale (A)

SimuBlocks <- function (n, sigma = 0.1, seed=NA) 
{
  if (!is.na(seed)) set.seed(seed)
  x = seq(1/n, 1,length.out = n)
  tau <- c(0.1, 0.13, 0.15, 0.23, 0.25, 0.4, 0.44, 0.65, 
           0.76, 0.78, 0.81)
  h <- c(0, 4, -5, 3, -4, 5, -4.2, 2.1, 4.3, -3.1, 2.1, -4.2)/5
  A = sapply(tau, function(z) which(x>=z)[1])
  tau1=c(0,tau,1)
  y0 = 0*x
  for (j in 1:(length(A)+1)) y0[x>tau1[j] & x<=tau1[j+1]] = h[j]
  y <- y0 + sigma*rnorm(n)
  return(list(y = y, x = x, y0 = y0, tau = tau, SetA = A))  
}