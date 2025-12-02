#' @name htln
#' @rdname htln
#'
#' @aliases mle_htln dhtln qhtln rhtln
#'
#' @title Hurdle Truncated Log-Normal Model
#'
#' @description
#' Functions to fit and work with a hurdle model whose continuous component
#' follows a truncated log-normal distribution.
#'
#' The model is tailored for non-negative data with many zeros and heavy tails,
#' such as NGS counts after a \code{log(x + 1)} transformation. Let
#' \eqn{Y} denote the original-scale data (e.g. counts) and
#' \eqn{Z = \log(Y + 1)}. The distribution is
#' \deqn{
#'   P(Y = 0) = \phi, \qquad
#'   Z \mid (Y > 0) \sim \mathrm{TruncNorm}(\mu, \sigma^2; [a, b]),
#' }
#' where \eqn{\phi} is the probability of an excess zero (hurdle component)
#' and \eqn{Z} is truncated to the interval \code{[a, b]} on the log scale.
#'
#' On the original scale, \code{a} and \code{b} correspond to truncation
#' bounds applied to \eqn{Z = \log(Y + 1)}, not directly to \eqn{Y}.
#'
#' @param x vector of quantiles on the log-transformed scale
#'   (i.e., values of \eqn{Z = \log(Y + 1)}).
#' @param p vector of probabilities.
#' @param n number of observations to generate.
#' @param phi probability mass at zero (hurdle component), \eqn{P(Y = 0)}.
#' @param meanlog mean of \eqn{Z = \log(Y + 1)} on the log scale.
#' @param sdlog standard deviation of \eqn{Z = \log(Y + 1)} on the log scale.
#' @param a lower bound for the truncation interval on the log scale.
#' @param b upper bound for the truncation interval on the log scale.
#' @param y numeric vector on the original scale (non-negative), used for
#'   parameter estimation. Internally transformed as \code{log1p(y)}.
#' @param integer logical (default \code{TRUE}). If \code{TRUE}, values
#'   on the original scale are coerced to integers: positive values in
#'   \code{(0, 1)} are set to 1 and all values are rounded.
#' @param warning.silent logical. If \code{TRUE}, suppress warning messages
#'   from \code{\link[MASS]{fitdistr}} during the MLE optimization.
#'
#' @return
#' \itemize{
#'   \item \code{mle_htln()} returns a list with components
#'     \code{estimate}, \code{loglik} and \code{success}.
#'   \item \code{dhtln()} returns the density evaluated at \code{x}.
#'   \item \code{qhtln()} returns the quantiles for probabilities \code{p}.
#'   \item \code{rhtln()} returns random draws from the hurdle truncated
#'     log-normal model, on the original scale (after inverse \code{log1p}).
#' }
#'
#' @seealso \code{\link[truncnorm]{dtruncnorm}},
#'   \code{\link[truncnorm]{qtruncnorm}}, \code{\link[truncnorm]{rtruncnorm}},
#'   \code{\link[MASS]{fitdistr}}
#'
#' @importFrom cli cli_abort
NULL


#------------------------------------------------------------------------------#
#                         MLE_HTLN                                             #
#------------------------------------------------------------------------------#
#' @rdname htln
#' @export
mle_htln <- function(y, warning.silent = TRUE) {
  
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
  
  phi <- n0 / n
  
  # FIT TRUNCATED NORMAL ------------------------------------------------------
  param.trnorm.mle <- NULL
  b_val <- NA_real_  

  # Fit the continuous component only if we have at least two positive values
  if (n1 > 1L) {
    sd_pos <- stats::sd(y_log_pos)
    b_val <- max(y_log_pos) + 3 * sd_pos
    
    # Fit only when sd is finite and strictly positive.
    # If sd is 0 (all positive values identical), the truncated normal cannot be estimated.
    if (is.finite(sd_pos) && sd_pos > 0) {
      
      # Use a finite upper bound instead of Inf to improve numerical stability.
      # max(y_log_pos) + 3 * sd_pos approximates the infinite support
      # without causing overflow or optimization failures.
      try(
        param.trnorm.mle <- MASS::fitdistr(
          x       = y_log_pos,
          densfun = truncnorm::dtruncnorm,
          start   = list(mean = mean(y_log_pos), sd = sd_pos),
          a       = 0,                         # lower bound (log(Y+1) >= 0)
          b       = b_val                      # stable approximation of +Inf
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
    loglik <- param.trnorm.mle$loglik
  } else {
    loglik <- NA_real_
  }
  
  # OUTPUT --------------------------------------------------------------------
  if (!is.null(param.trnorm.mle)) {
    estimate <- c(
      "phi"     = phi,
      "meanlog" = as.numeric(param.trnorm.mle$estimate[1]),
      "sdlog"   = as.numeric(param.trnorm.mle$estimate[2]),
      "b" = b_val
    )
    result <- list(
      estimate = estimate,
      loglik   = loglik,
      success  = 1L
    )
  } else {
    estimate <- c(
      "phi"     = phi,
      "meanlog" = NA_real_,
      "sdlog"   = NA_real_,
      "b" = NA_real_
    )
    result <- list(
      estimate = estimate,
      loglik   = loglik,
      success  = 0L
    )
  }
  
  return(result)
}


#------------------------------------------------------------------------------#
#                            DHTLN                                             #
#------------------------------------------------------------------------------#
#' @rdname htln
#' @importFrom truncnorm dtruncnorm
#' @export
dhtln <- function(x, phi = 0, meanlog = 0, sdlog = 1, a = 0, b = Inf) {
  
  # ARGUMENT CHECKS -----------------------------------------------------------
  if (!is.numeric(x) || !is.vector(x)) {
    cli::cli_abort("{.arg x} must be a numeric vector.")
  }
  if (any(!is.finite(x))) {
    cli::cli_abort("{.arg x} must not contain NA, NaN or Inf values.")
  }
  
  if (!is.numeric(phi) || length(phi) != 1L || !is.finite(phi)) {
    cli::cli_abort("{.arg phi} must be a single finite numeric value.")
  }
  if (phi < 0 || phi > 1) {
    cli::cli_abort("{.arg phi} must be in the range [0, 1].")
  }
  
  if (!is.numeric(meanlog) || length(meanlog) != 1L || !is.finite(meanlog)) {
    cli::cli_abort("{.arg meanlog} must be a single finite numeric value.")
  }
  if (!is.numeric(sdlog) || length(sdlog) != 1L || !is.finite(sdlog)) {
    cli::cli_abort("{.arg sdlog} must be a single finite numeric value.")
  }
  if (sdlog <= 0) {
    cli::cli_abort("{.arg sdlog} must be greater than 0.")
  }
  
  if (!is.numeric(a) || length(a) != 1L || !is.finite(a)) {
    cli::cli_abort("{.arg a} must be a single finite numeric value.")
  }
  if (!is.numeric(b) || length(b) != 1L || !is.finite(b)) {
    cli::cli_abort("{.arg b} must be a single finite numeric value.")
  }
  
  # DENSITY -------------------------------------------------------------------
  ans <- (1 - phi) * truncnorm::dtruncnorm(
    x    = x,
    a    = a,
    b    = b,
    mean = meanlog,
    sd   = sdlog
  )
  ans[x == 0] <- phi
  
  return(ans)
}


#------------------------------------------------------------------------------#
#                            QHTLN                                             #
#------------------------------------------------------------------------------#
#' @rdname htln
#' @export
qhtln <- function(p, phi = 0, meanlog = 0, sdlog = 1,
                  a = 0, b = Inf, integer = TRUE) {
  
  # ARGUMENT CHECKS -----------------------------------------------------------
  if (!is.numeric(p) || !is.vector(p)) {
    cli::cli_abort("{.arg p} must be a numeric vector.")
  }
  if (any(!is.finite(p))) {
    cli::cli_abort("{.arg p} must not contain NA, NaN or Inf values.")
  }
  if (any(p < 0 | p > 1)) {
    cli::cli_abort("All probabilities in {.arg p} must be in the range [0, 1].")
  }
  
  if (!is.numeric(phi) || length(phi) != 1L || !is.finite(phi)) {
    cli::cli_abort("{.arg phi} must be a single finite numeric value.")
  }
  if (phi < 0 || phi > 1) {
    cli::cli_abort("{.arg phi} must be in the range [0, 1].")
  }
  
  if (!is.numeric(meanlog) || length(meanlog) != 1L || !is.finite(meanlog)) {
    cli::cli_abort("{.arg meanlog} must be a single finite numeric value.")
  }
  if (!is.numeric(sdlog) || length(sdlog) != 1L || !is.finite(sdlog)) {
    cli::cli_abort("{.arg sdlog} must be a single finite numeric value.")
  }
  # mantengo la logica originale: sdlog < 0 (non <= 0)
  if (sdlog < 0) {
    cli::cli_abort("{.arg sdlog} must be greater than or equal to 0.")
  }
  
  if (!is.numeric(a) || length(a) != 1L || !is.finite(a)) {
    cli::cli_abort("{.arg a} must be a single finite numeric value.")
  }
  if (!is.numeric(b) || length(b) != 1L || !is.finite(b)) {
    cli::cli_abort("{.arg b} must be a single finite numeric value.")
  }
  
  if (!is.logical(integer) || length(integer) != 1L || is.na(integer)) {
    cli::cli_abort("{.arg integer} must be a single non-missing logical value.")
  }
  
  # QUANTILES -----------------------------------------------------------------
  if (phi == 1) {
    ans <- rep(0, length(p))
  } else {
    ans <- rep(NA_real_, length(p))
    ans[p <= phi] <- 0
    
    idx <- is.na(ans)
    ans[idx] <- truncnorm::qtruncnorm(
      p    = (p[idx] - phi) / (1 - phi),
      a    = a,
      b    = b,
      mean = meanlog,
      sd   = sdlog
    )
    
    # back to original scale
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
                  a = 0, b = Inf, integer = TRUE) {
  
  # ARGUMENT CHECKS -----------------------------------------------------------
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n)) {
    cli::cli_abort("{.arg n} must be a single finite numeric value.")
  }
  if (round(n) != n) {
    cli::cli_abort("{.arg n} must be an integer.")
  }
  if (n < 1) {
    cli::cli_abort("{.arg n} must be a positive integer.")
  }
  
  if (!is.numeric(phi) || length(phi) != 1L || !is.finite(phi)) {
    cli::cli_abort("{.arg phi} must be a single finite numeric value.")
  }
  if (phi < 0 || phi > 1) {
    cli::cli_abort("{.arg phi} must be in the range [0, 1].")
  }
  
  if (!is.numeric(meanlog) || length(meanlog) != 1L || !is.finite(meanlog)) {
    cli::cli_abort("{.arg meanlog} must be a single finite numeric value.")
  }
  if (!is.numeric(sdlog) || length(sdlog) != 1L || !is.finite(sdlog)) {
    cli::cli_abort("{.arg sdlog} must be a single finite numeric value.")
  }
  # mantengo la logica originale: sdlog < 0 (non <= 0)
  if (sdlog < 0) {
    cli::cli_abort("{.arg sdlog} must be a positive number (or 0).")
  }
  
  if (!is.numeric(a) || length(a) != 1L || !is.finite(a)) {
    cli::cli_abort("{.arg a} must be a single finite numeric value.")
  }
  if (!is.numeric(b) || length(b) != 1L) {
    cli::cli_abort("{.arg b} must be a single numeric value or Inf.")
  }
  
  if (!is.logical(integer) || length(integer) != 1L || is.na(integer)) {
    cli::cli_abort("{.arg integer} must be a single non-missing logical value.")
  }
  
  # RANDOM GENERATION ---------------------------------------------------------
  # Bernoulli draw for the hurdle part (1 = positive, 0 = zero)
  ans <- stats::rbinom(n = n, size = 1, prob = 1 - phi)
  m <- sum(ans > 0)
  
  if (m > 0) {
    ans[ans == 1] <- truncnorm::rtruncnorm(
      n    = m,
      a    = a,
      b    = b,
      mean = meanlog,
      sd   = sdlog
    )
  }
  
  # back to original scale
  ans <- expm1(ans)
  
  if (integer) {
    ans[ans > 0 & ans < 1] <- 1
    ans <- round(ans)
  }
  
  return(ans)
}
