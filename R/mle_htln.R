#------------------------------------------------------------------------------#
#                           MLE for HTLN Model                                 #
#------------------------------------------------------------------------------#

#' @rdname htln
#' @export
mle_htln <- function(
    y,
    b_log = function(z) max(z) + 3 * stats::sd(z),
    optim_method = "BFGS",
    optim_control = list(),
    optim_hessian = FALSE,
    optim_lower = -Inf, optim_upper = Inf
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
  
  # LOG-TRANSFORM: Z = log(Y + 1) ---------------------------------------------
  y_log      <- log1p(y)
  y_log_pos  <- y_log[y_log > 0]
  sd_log_pos <- stats::sd(y_log_pos)
  
  n  <- length(y_log)
  n1 <- length(y_log_pos)
  n0 <- n - n1
  
  # Hurdle probability: mass at zero
  phi <- n0 / n
  
  # CHECK POSITIVE LOG-TRANSFORMED ELEMENT ------------------------------------
  if (n1 <= 1L)
    cli::cli_abort("At least two positive observations are required for truncated normal fitting (`n1 > 1`).")
  
  if (!is.finite(sd_log_pos) || sd_log_pos <= 0)
    cli::cli_abort("Standard deviation of positive values must be finite and > 0. Found {.val {sd_log_pos}}.")
  
  # UPPER BOUND ON THE LOG SCALE ----------------------------------------------
  b_log_val <- NA_real_
  
  if (n1 > 0L) {
    if (is.function(b_log)) {
      b_log_val <- b_log(y_log_pos)
      if (!is.numeric(b_log_val) || length(b_log_val) != 1L || !is.finite(b_log_val)) {
        cli::cli_abort(
          "When {.arg b_log} is a function, it must return a single finite numeric value."
        )
      }
    } else if (is.numeric(b_log) && length(b_log) == 1L) {
      b_log_val <- b_log
    } else {
      cli::cli_abort("{.arg b_log} must be either a numeric scalar or a function.")
    }
    
    if (b_log_val < 0 || b_log_val < max(y_log_pos)) {
      cli::cli_abort(c(
        "{.arg b_log} defines an invalid upper truncation bound on the log scale.",
        "x" = "The computed bound is smaller than 0 or smaller than the maximum positive log-value.",
        "i" = "Please provide {.arg b_log} such that b_log >= 0 and b_log >= max(log1p(y[y > 0]))."
      ))
    }
  }
  
  # FIT TRUNCATED NORMAL (CONTINUOUS COMPONENT) -------------------------------
  fit_trnorm <- stats::optim(
    par = c(mean = mean(y_log_pos), sd = sd_log_pos),
    fn  = function(par){
      -sum(log(truncnorm::dtruncnorm(
        x    = y_log_pos,
        a    = 0, 
        b    = b_log_val,
        mean = par[1], 
        sd   = par[2]
      )))
    },
    method = optim_method,
    control = optim_control,
    hessian = optim_hessian,
    lower = optim_lower, upper = optim_upper
  )
  
  # COMPUTE LOG-LIKELIHOOD ----------------------------------------------------
  loglik_trunc <- -fit_trnorm$value
  
  if (phi != 0) {
    loglik <- n0 * log(phi) + n1 * log(1 - phi) + loglik_trunc
  } else {
    loglik <- loglik_trunc
  }
  
  # BUILD RESULT OBJECT -------------------------------------------------------
  result <- list(
    dist       = "htln",
    coef       = c(
      phi     = phi, 
      meanlog = fit_trnorm$par[[1]],
      sdlog   = fit_trnorm$par[[2]], 
      b_log   = b_log_val
    ),
    logLik     = as.numeric(loglik),
    nobs       = n,
    optim_trunc = fit_trnorm
  )
  
  class(result) <- "htln_fit"
  return(result)
}


