#------------------------------------------------------------------------------#
#                        S3 Methods for htln_fit                               #
#------------------------------------------------------------------------------#

#' @export
coef.htln_fit <- function(object, ...) {
  object$coef
}

#' @export
logLik.htln_fit <- function(object, ...) {
  ll <- object$logLik
  
  if (is.na(ll)) {
    return(NA_real_)
  }
  
  # phi, meanlog, sdlog are treated as free parameters when non-missing
  free_par <- c("phi", "meanlog", "sdlog")
  df <- sum(!is.na(object$coef[free_par]))
  
  attr(ll, "df")   <- df
  attr(ll, "nobs") <- object$nobs
  class(ll) <- "logLik"
  ll
}

#' @export
AIC.htln_fit <- function(object, ..., k = 2) {
  ll <- logLik(object)
  if (is.na(ll)) {
    return(NA_real_)
  }
  df <- attr(ll, "df")
  -2 * as.numeric(ll) + k * df
}

#' @export
BIC.htln_fit <- function(object, ...) {
  ll <- logLik(object)
  if (is.na(ll)) {
    return(NA_real_)
  }
  df   <- attr(ll, "df")
  nobs <- attr(ll, "nobs")
  -2 * as.numeric(ll) + log(nobs) * df
}

#' @export
nobs.htln_fit <- function(object, ...) {
  object$nobs
}

#' @export
print.htln_fit <- function(x, ...) {
  cat("HTLN model fit\n")
  cat("  distribution :", x$dist, "\n")
  
  # calcolo della "convergenza" dalla componente truncnorm
  conv_ok <- is.list(x$optim_trunc) && isTRUE(x$optim_trunc$convergence == 0L)
  cat("  convergence  :", if (conv_ok) "converged" else "not converged", "\n")
  cat("  nobs         :", x$nobs, "\n\n")
  
  cat("Coefficients:\n")
  print(x$coef)
  
  # Compute logLik / AIC / BIC on the fly usando i metodi S3
  ll  <- try(logLik(x), silent = TRUE)
  aic <- try(AIC(x),    silent = TRUE)
  bic <- try(BIC(x),    silent = TRUE)
  
  cat("\nlogLik :", if (inherits(ll,  "try-error")) NA_real_ else as.numeric(ll), "\n")
  cat("AIC    :", if (inherits(aic, "try-error")) NA_real_ else aic, "\n")
  cat("BIC    :", if (inherits(bic, "try-error")) NA_real_ else bic, "\n")
  
  invisible(x)
}


