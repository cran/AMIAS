print.samias <- function (x, digits = max(3, getOption("digits") - 3), ...) 
{
  cat("\nCall: ", deparse(x$call), "\n\n")
  print(cbind(k = x$k, Df = x$df, `BIC` = signif(x$bic, digits)))
}