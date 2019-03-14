plot.samias <- function(x, type = c("coef", "vpath"), k, add.label = TRUE, add.knots = TRUE,...){
  
  type <- match.arg(type)
  if (!(type %in% c("coef", "vpath"))) {
    stop("Invalid type, must be \"coef\" or \"vpath\".")
  }
  
  n <- length(x$y)
  k.max <- dim(x$alpha.all)[2]
  if(type == "vpath"){
    vall <- x$v.all
    kopt <- x$kopt
    n <- length(x$y)
    
    k = dim(vall)[2]
    A1 = which(rowSums(abs(vall))>0)
    TmpX = data.frame(ID=A1,V=vall[A1,])
    newA1 = TmpX$ID[which(abs(TmpX[,2])>0)]  ## Estimated active set for each cardinality k
    for (j in 3:(k+1)){
      tmp = TmpX[!is.element(TmpX$ID,newA1),]
      newA1 = c(newA1, tmp$ID[which(abs(tmp[,j])>0)])
    }
    newV1 = rbind(0,t(vall[newA1,]))
    matplot(seq(0, k), newV1, type="l", pch=1, lty=1, col = rainbow(length(newA1)),
            xlim = c(0,k+1),  xlab = "Cardinality k", ylab = " Estimated v", ...)
    abline(v=kopt, lty=2, col=1) ## add a vertical line
    if (add.label) {
      d <- as.matrix(dist(newV1[k+1,]))
      idx <- which(d<0.1, arr.ind = TRUE)
      idx1 <- idx[idx[,1]<idx[,2],,drop=FALSE]
      xpos <- rep(k, length(newA1))
      xpos[idx1[,2]] <- xpos[idx1[,2]] + 0.5
      text(xpos,newV1[k+1,], round(newA1/n, 2), cex=0.6, pos=4, col = 1)
    }
    
    
  }
  if(type == "coef"){
    y0 <- seq(1/n, 1, length.out = n)
    yname = "y"
    if(x$smooth)yname <- "smooth-y"
    plot.default(y0, x$y, ylab = yname, col="lightgrey", pch = 20, xlab="",...)
    type.est <- ifelse(x$q==0, 's', "l") # Specify the plot type for the estimate
    
    if(missing(k)){
      lines(y0, x$alpha, col = "red", type=type.est, lwd=2, lty=1, cex = 0.5)
      
      if(add.knots == TRUE){
        ## Add the locations for the detected knots
        knot <- y0[which(x$v!=0)]
        points(knot, rep(par("usr")[4], length(knot)), pch=3, cex=2, col="red")
      }
    }else{
      lines(y0, x$alpha.all[,k], col = "red", type=type.est, lwd=2, lty=1, cex = 0.5)
      
      if(add.knots == TRUE){
        ## Add the locations for the detected knots
        knot <- y0[which(x$v.all[,k]!=0)]
        points(knot, rep(par("usr")[4], length(knot)), pch=3, cex=2, col="red")
      }
    }
   
  }
  
}
