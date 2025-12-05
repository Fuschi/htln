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
  cat("  success      :", x$success, "\n")
  cat("  nobs         :", x$nobs, "\n\n")
  
  cat("Coefficients:\n")
  print(x$coef)
  
  # Compute AIC/BIC on the fly
  aic_val <- try(AIC(x), silent = TRUE)
  bic_val <- try(BIC(x), silent = TRUE)
  
  cat("\nlogLik :", x$logLik, "\n")
  cat("AIC    :", if (inherits(aic_val, "try-error")) NA_real_ else aic_val, "\n")
  cat("BIC    :", if (inherits(bic_val, "try-error")) NA_real_ else bic_val, "\n")
  
  invisible(x)
}

