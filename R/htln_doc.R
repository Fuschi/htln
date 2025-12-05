#------------------------------------------------------------------------------#
#                  Hurdle Truncated Log-Normal (HTLN) Model                    #
#------------------------------------------------------------------------------#

#' Hurdle Truncated Log-Normal (HTLN) Family
#'
#' Functions to fit, evaluate, and simulate from a hurdle model whose
#' continuous component follows a truncated log-normal distribution.
#'
#' This model is designed for non-negative data with excess zeros and heavy
#' tails—common in microbiome and NGS contexts—where the positive component is
#' modeled on the transformed scale \eqn{Z = \log(Y + 1)}.
#'
#' The model is:
#' \deqn{
#'   P(Y = 0) = \phi, \qquad
#'   Z \mid (Y > 0) \sim \mathrm{TruncNorm}(\mu, \sigma^2; [a_{\log}, b_{\log}]),
#' }
#'
#' where:
#' * \eqn{\phi} is the hurdle probability of observing a structural zero,
#' * \eqn{Z = \log(Y+1)} is truncated to the interval \eqn{[a_{\log}, b_{\log}]},
#' * \eqn{\mu} and \eqn{\sigma} are the mean and standard deviation on the
#'   log-transformed scale.
#'
#' Note that the truncation bounds \code{a_log} and \code{b_log} apply on the
#' \eqn{Z = \log(Y+1)} scale, not directly to \eqn{Y}.
#'
#' @section Fitting:
#' \code{mle_htln()} performs maximum likelihood estimation of the hurdle
#' probability and the parameters of the truncated normal component using
#' \code{\link[MASS]{fitdistr}}. The function returns an object of class
#' \code{"htln_fit"}, for which standard S3 generics \code{coef()},
#' \code{logLik()}, \code{AIC()}, \code{BIC()} and \code{nobs()} are available.
#'
#' The argument \code{b_log} controls the upper truncation bound on the
#' \eqn{Z = \log(Y+1)} scale and can be:
#' \itemize{
#'   \item a numeric scalar (fixed bound, possibly \code{Inf});
#'   \item a function of the positive log-values (default), e.g.
#'     \code{function(z) max(z) + 3 * sd(z)}.
#' }
#'
#' @param y Numeric vector of non-negative observations (used for fitting by
#'   \code{mle_htln()}). Values are internally transformed as \code{log1p(y)}.
#' @param x Numeric vector of points at which to evaluate the density
#'   \code{dhtln()} (on the \eqn{Z = \log(Y+1)} scale).
#' @param p Numeric vector of probabilities for \code{qhtln()}.
#' @param n Number of observations to generate in \code{rhtln()}.
#' @param phi Hurdle probability \eqn{P(Y = 0)}.
#' @param meanlog Mean parameter of the truncated log-normal component.
#' @param sdlog Standard deviation on the log-transformed scale.
#' @param a_log,b_log Lower and upper truncation bounds on the \eqn{Z = \log(Y+1)} scale.
#'   In \code{mle_htln()}, \code{b_log} can also be a function of the positive
#'   log-values controlling the upper truncation bound (see Details).
#' @param integer Logical. If \code{TRUE} (default), random variates generated
#'   by \code{rhtln()} are coerced to integers on the original scale.
#' @param warning.silent Logical. If \code{TRUE}, suppress warnings produced by
#'   \code{\link[MASS]{fitdistr}} during optimisation.
#'
#' @return
#' \itemize{
#'   \item \code{mle_htln()} returns an object of class \code{"htln_fit"} containing:
#'         a named coefficient vector (\code{phi}, \code{meanlog}, \code{sdlog},
#'         \code{b_log}), the maximised log-likelihood, the number of observations,
#'         and a convergence flag. Standard generics \code{coef()},
#'         \code{logLik()}, \code{AIC()}, \code{BIC()} and \code{nobs()}
#'         can be applied to this object.
#'
#'   \item \code{dhtln()} returns density values of the hurdle truncated
#'         log-normal model evaluated at \code{x} on the log scale
#'         (\eqn{Z = \log(Y+1)}), mixing a point mass \code{phi} at \eqn{Z = 0}
#'         and a truncated normal density for \eqn{Z > 0}.
#'
#'   \item \code{qhtln()} returns quantiles on the original scale for
#'         probabilities \code{p}.
#'
#'   \item \code{rhtln()} generates random samples from the hurdle truncated
#'         log-normal distribution on the original scale.
#' }
#'
#' @seealso
#'   \code{\link[MASS]{fitdistr}},
#'   \code{\link[truncnorm]{dtruncnorm}},
#'   \code{\link[truncnorm]{qtruncnorm}},
#'   \code{\link[truncnorm]{rtruncnorm}}
#'
#' @aliases
#'   mle_htln
#'   dhtln
#'   qhtln
#'   rhtln
#'
#' @importFrom cli cli_abort
#' @importFrom stats AIC BIC logLik nobs
#'
#' @name htln
NULL


