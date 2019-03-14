plot.amias <- function(x, add.knots = TRUE,...){
  
  n <- length(x$y)
  y0 <- seq(1/n, 1, length.out = n)
  yname = "y"
  if(x$smooth)yname <- "smooth-y"
  plot.default(y0, x$y, ylab = yname, col="lightgrey", pch = 20, xlab="",...)
  
  ## Plot the estimates
  type <- ifelse(x$q==0, 's', "l") # Specify the plot type for the estimate
  lines(y0, x$alpha, col = "red", type = type, lwd = 2, lty = 1)
  
  if(add.knots == TRUE){
    ## Add the locations for the detected knots
    knot <- y0[which(x$v!=0)]
    points(knot, rep(par("usr")[4], length(knot)), pch=3, cex=2, col="red")
  }
  
}