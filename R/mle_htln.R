#------------------------------------------------------------------------------#
#                           MLE for HTLN Model                                 #
#------------------------------------------------------------------------------#

#' @rdname htln
#' @export
mle_htln <- function(
    y,
    b_log = function(z) max(z) + 3 * stats::sd(z),
    warning.silent = TRUE
) {
  # ARGUMENT CHECKS -----------------------------------------------------------
  if (!is.numeric(y) || !is.vector(y)) {
    cli::cli_abort("{.arg y} must be a numeric vector.")
  }
  if (length(y) == 0L) {
    cli::cli_abort("{.arg y} must have length > 0.")
  }
  if (any(!is.finite(y))) {
    cli::cli_abort("{.arg y} must not contain NA, NaN or Inf values.")
  }
  if (any(y < 0)) {
    cli::cli_abort("{.arg y} must contain only values >= 0.")
  }
  if (!is.logical(warning.silent) || length(warning.silent) != 1L || is.na(warning.silent)) {
    cli::cli_abort("{.arg warning.silent} must be a single non-missing logical value.")
  }
  
  # log-transform: Z = log(Y + 1)
  y_log     <- log1p(y)
  y_log_pos <- y_log[y_log > 0]
  
  n  <- length(y_log)
  n1 <- length(y_log_pos)
  n0 <- n - n1
  
  # Hurdle probability: mass at zero
  phi <- n0 / n
  
  # UPPER BOUND ON THE LOG SCALE ----------------------------------------------
  # b_log can be:
  #  - numeric scalar (fixed bound, possibly Inf)
  #  - function: b_log(z) returning a scalar bound from positive values
  b_log_val <- NA_real_
  
  if (n1 > 0L) {
    if (is.function(b_log)) {
      # b_log is a function of the positive log-values
      b_log_val <- b_log(y_log_pos)
      if (!is.numeric(b_log_val) || length(b_log_val) != 1L || !is.finite(b_log_val)) {
        cli::cli_abort(
          "When {.arg b_log} is a function, it must return a single finite numeric value."
        )
      }
    } else if (is.numeric(b_log) && length(b_log) == 1L) {
      # fixed numeric bound (can be Inf)
      b_log_val <- b_log
    } else {
      cli::cli_abort(
        "{.arg b_log} must be either a numeric scalar or a function."
      )
    }
    
    # Basic sanity check: upper bound must be >= 0 and >= max positive log-value
    if (b_log_val < 0 || b_log_val < max(y_log_pos)) {
      cli::cli_abort(c(
        "{.arg b_log} defines an invalid upper truncation bound on the log scale.",
        "x" = "The computed bound is smaller than 0 or smaller than the maximum positive log-value.",
        "i" = "Please provide {.arg b_log} such that b_log >= 0 and b_log >= max(log1p(y[y > 0]))."
      ))
    }
  }
  
  # FIT TRUNCATED NORMAL ------------------------------------------------------
  param.trnorm.mle <- NULL
  
  # Fit the continuous component only if we have at least two positive values
  if (n1 > 1L) {
    sd_pos <- stats::sd(y_log_pos)
    
    # Fit only when sd is finite and strictly positive
    if (is.finite(sd_pos) && sd_pos > 0) {
      try(
        param.trnorm.mle <- MASS::fitdistr(
          x       = y_log_pos,
          densfun = truncnorm::dtruncnorm,
          start   = list(mean = mean(y_log_pos), sd = sd_pos),
          a       = 0,          # lower bound on log scale
          b       = b_log_val   # fixed upper bound (may be Inf)
        ),
        silent = warning.silent
      )
    }
  }
  # If n1 <= 1, or sd_pos is invalid (Inf/NA/0),
  # param.trnorm.mle remains NULL and the MLE is considered unsuccessful.
  
  # COMPUTE LOG-LIKELIHOOD ----------------------------------------------------
  if (!is.null(param.trnorm.mle) && phi != 0) {
    loglik <- n0 * log(phi) + n1 * log(1 - phi) + param.trnorm.mle$loglik
  } else if (!is.null(param.trnorm.mle) && phi == 0) {
    # No mass at zero, pure truncated log-normal on log(Y+1)
    loglik <- param.trnorm.mle$loglik
  } else {
    loglik <- NA_real_
  }
  
  # PARAMETER VECTOR ----------------------------------------------------------
  if (!is.null(param.trnorm.mle)) {
    
    # Extract mean and sd by name for safety
    est_mean <- as.numeric(param.trnorm.mle$estimate["mean"])
    est_sd   <- as.numeric(param.trnorm.mle$estimate["sd"])
    
    coef_vec <- c(
      "phi"      = phi,
      "meanlog"  = est_mean,
      "sdlog"    = est_sd,
      "b_log"    = b_log_val
    )
    success <- 1L
  } else {
    coef_vec <- c(
      "phi"      = phi,
      "meanlog"  = NA_real_,
      "sdlog"    = NA_real_,
      "b_log"    = NA_real_
    )
    success <- 0L
  }
  
  # BUILD RESULT OBJECT -------------------------------------------------------
  result <- list(
    dist    = "htln",
    coef    = coef_vec,
    logLik  = as.numeric(loglik),
    nobs    = n,
    success = success
  )
  
  class(result) <- "htln_fit"
  return(result)
}
