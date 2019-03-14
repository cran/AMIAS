# simulate doppler
# n determines length 
# sigma is residual sd
# output is data (y),  input location (x), mean function (y0), location of knots in x scale (tau) and in sample scale (A)


SimuDoppler <- function(n, sigma = 0.1, seed=NA){
  if (!is.na(seed)) set.seed(seed)
  x <- seq(1/n, 1,length.out = n)
  y0 <- sqrt(x*(1-x))*sin(2*pi*(1+0.05)/(x+0.05))
  y <- y0 + sigma*rnorm(n)
  return(list(y = y, y0 = y0, x=x))
}
