#------------------------------------------------------------------------------#
#                            DHTLN                                             #
#------------------------------------------------------------------------------#

#' @rdname htln
#' @export
dhtln <- function(x, phi = 0, meanlog = 0, sdlog = 1, a_log = 0, b_log = Inf) {
  
  # ARGUMENT CHECKS -----------------------------------------------------------
  # x must be numeric, finite, non-empty
  if (!is.numeric(x) || !is.vector(x)) {
    cli::cli_abort("{.arg x} must be a numeric vector.")
  }
  if (length(x) == 0L) {
    cli::cli_abort("{.arg x} must have length > 0.")
  }
  if (any(!is.finite(x))) {
    cli::cli_abort("{.arg x} must not contain NA, NaN or Inf values.")
  }
  
  # phi checks
  if (!is.numeric(phi) || length(phi) != 1L || !is.finite(phi)) {
    cli::cli_abort("{.arg phi} must be a single finite numeric value.")
  }
  if (phi < 0 || phi > 1) {
    cli::cli_abort("{.arg phi} must be in the range [0, 1].")
  }
  
  # meanlog and sdlog checks
  if (!is.numeric(meanlog) || length(meanlog) != 1L || !is.finite(meanlog)) {
    cli::cli_abort("{.arg meanlog} must be a single finite numeric value.")
  }
  if (!is.numeric(sdlog) || length(sdlog) != 1L || !is.finite(sdlog)) {
    cli::cli_abort("{.arg sdlog} must be a single finite numeric value.")
  }
  if (sdlog <= 0) {
    cli::cli_abort("{.arg sdlog} must be greater than 0.")
  }
  
  # truncation bounds on the log scale
  if (!is.numeric(a_log) || length(a_log) != 1L || !is.finite(a_log)) {
    cli::cli_abort("{.arg a_log} must be a single finite numeric value.")
  }
  if (!is.numeric(b_log) || length(b_log) != 1L || is.na(b_log)) {
    cli::cli_abort("{.arg b_log} must be a single numeric value (finite or Inf).")
  }
  if (b_log < a_log) {
    cli::cli_abort("{.arg b_log} must be >= {.arg a_log}.")
  }
  
  # DENSITY -------------------------------------------------------------------
  # Continuous part (truncated normal on log scale)
  ans <- (1 - phi) * truncnorm::dtruncnorm(
    x    = x,
    a    = a_log,
    b    = b_log,
    mean = meanlog,
    sd   = sdlog
  )
  
  # Point mass at 0 on the log scale
  ans[x == 0] <- phi
  
  return(ans)
}


#------------------------------------------------------------------------------#
#                            QHTLN                                             #
#------------------------------------------------------------------------------#

#' @rdname htln
#' @export
qhtln <- function(p, phi = 0, meanlog = 0, sdlog = 1,
                  a_log = 0, b_log = Inf, integer = TRUE) {
  
  # ARGUMENT CHECKS -----------------------------------------------------------
  # p must be numeric, finite, non-empty
  if (!is.numeric(p) || !is.vector(p)) {
    cli::cli_abort("{.arg p} must be a numeric vector.")
  }
  if (length(p) == 0L) {
    cli::cli_abort("{.arg p} must have length > 0.")
  }
  if (any(!is.finite(p))) {
    cli::cli_abort("{.arg p} must not contain NA, NaN or Inf values.")
  }
  if (any(p < 0 | p > 1)) {
    cli::cli_abort("All probabilities in {.arg p} must be in the range [0, 1].")
  }
  
  # phi checks
  if (!is.numeric(phi) || length(phi) != 1L || !is.finite(phi)) {
    cli::cli_abort("{.arg phi} must be a single finite numeric value.")
  }
  if (phi < 0 || phi > 1) {
    cli::cli_abort("{.arg phi} must be in the range [0, 1].")
  }
  
  # meanlog and sdlog checks
  if (!is.numeric(meanlog) || length(meanlog) != 1L || !is.finite(meanlog)) {
    cli::cli_abort("{.arg meanlog} must be a single finite numeric value.")
  }
  if (!is.numeric(sdlog) || length(sdlog) != 1L || !is.finite(sdlog)) {
    cli::cli_abort("{.arg sdlog} must be a single finite numeric value.")
  }
  # Original logic: allow sdlog == 0
  if (sdlog < 0) {
    cli::cli_abort("{.arg sdlog} must be >= 0.")
  }
  
  # truncation bounds on the log scale
  if (!is.numeric(a_log) || length(a_log) != 1L || !is.finite(a_log)) {
    cli::cli_abort("{.arg a_log} must be a single finite numeric value.")
  }
  if (!is.numeric(b_log) || length(b_log) != 1L || is.na(b_log)) {
    cli::cli_abort("{.arg b_log} must be a single numeric value (finite or Inf).")
  }
  if (b_log < a_log) {
    cli::cli_abort("{.arg b_log} must be >= {.arg a_log}.")
  }
  
  # integer flag
  if (!is.logical(integer) || length(integer) != 1L || is.na(integer)) {
    cli::cli_abort("{.arg integer} must be a single non-missing logical value.")
  }
  
  # QUANTILE COMPUTATION ------------------------------------------------------
  if (phi == 1) {
    # Pure hurdle at zero
    ans <- rep(0, length(p))
  } else {
    ans <- rep(NA_real_, length(p))
    
    # Hurdle part (point mass at 0 on original scale)
    ans[p <= phi] <- 0
    
    # Continuous part (truncated normal on log scale)
    idx <- is.na(ans)
    if (any(idx)) {
      ans[idx] <- truncnorm::qtruncnorm(
        p    = (p[idx] - phi) / (1 - phi),
        a    = a_log,
        b    = b_log,
        mean = meanlog,
        sd   = sdlog
      )
    }
    
    # Transform back to original scale
    ans <- expm1(ans)
    
    if (integer) {
      ans[ans > 0 & ans < 1] <- 1
      ans <- round(ans)
    }
  }
  
  return(ans)
}


#------------------------------------------------------------------------------#
#                            RHTLN                                             #
#------------------------------------------------------------------------------#

#' @rdname htln
#' @export
rhtln <- function(n, phi = 0, meanlog = 0, sdlog = 1,
                  a_log = 0, b_log = Inf, integer = TRUE) {
  
  # ARGUMENT CHECKS -----------------------------------------------------------
  # n must be a positive integer
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n)) {
    cli::cli_abort("{.arg n} must be a single finite numeric value.")
  }
  if (round(n) != n) {
    cli::cli_abort("{.arg n} must be an integer.")
  }
  if (n < 1) {
    cli::cli_abort("{.arg n} must be a positive integer.")
  }
  
  # phi checks
  if (!is.numeric(phi) || length(phi) != 1L || !is.finite(phi)) {
    cli::cli_abort("{.arg phi} must be a single finite numeric value.")
  }
  if (phi < 0 || phi > 1) {
    cli::cli_abort("{.arg phi} must be in the range [0, 1].")
  }
  
  # meanlog and sdlog checks
  if (!is.numeric(meanlog) || length(meanlog) != 1L || !is.finite(meanlog)) {
    cli::cli_abort("{.arg meanlog} must be a single finite numeric value.")
  }
  if (!is.numeric(sdlog) || length(sdlog) != 1L || !is.finite(sdlog)) {
    cli::cli_abort("{.arg sdlog} must be a single finite numeric value.")
  }
  # Original logic: allow sdlog == 0
  if (sdlog < 0) {
    cli::cli_abort("{.arg sdlog} must be >= 0.")
  }
  
  # truncation bounds on the log scale
  if (!is.numeric(a_log) || length(a_log) != 1L || !is.finite(a_log)) {
    cli::cli_abort("{.arg a_log} must be a single finite numeric value.")
  }
  if (!is.numeric(b_log) || length(b_log) != 1L || is.na(b_log)) {
    cli::cli_abort("{.arg b_log} must be a single numeric value (finite or Inf).")
  }
  if (b_log < a_log) {
    cli::cli_abort("{.arg b_log} must be >= {.arg a_log}.")
  }
  
  # integer flag
  if (!is.logical(integer) || length(integer) != 1L || is.na(integer)) {
    cli::cli_abort("{.arg integer} must be a single non-missing logical value.")
  }
  
  # RANDOM GENERATION ---------------------------------------------------------
  # Bernoulli draw: 1 = positive, 0 = zero
  ans <- stats::rbinom(n = n, size = 1, prob = 1 - phi)
  m <- sum(ans > 0)
  
  if (m > 0) {
    ans[ans == 1] <- truncnorm::rtruncnorm(
      n    = m,
      a    = a_log,
      b    = b_log,
      mean = meanlog,
      sd   = sdlog
    )
  }
  
  # Back to original scale
  ans <- expm1(ans)
  
  if (integer) {
    ans[ans > 0 & ans < 1] <- 1
    ans <- round(ans)
  }
  
  return(ans)
}
