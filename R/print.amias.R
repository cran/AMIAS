print.amias <- function(x, ...) {
  cat("\nCall:\n")
  dput(x$call)
  cat("\nOutput:\n")
  cat(paste("Generalized L0 coefficients with number of knots being ", x$k,".",
            "\n\n", sep=""))
}
