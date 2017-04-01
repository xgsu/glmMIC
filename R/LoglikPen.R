#' Compute the penalized log likelihood for GLM with MIC penalty (Self-Written)
#'
#' @param beta A p-dimensional vector containing the regression ceofficients.
#' @param I.preselect A 0-1 indicator vector of same length as \code{beta} showing preselected variables (value 0 for no penalty) that will not be penalized.
#' @param preselect.beta0 Indicator of whether or not the intercept is pre-selected. Default is \code{FALSE}.
#' @param X An \eqn{n} by \eqn{p} design matrix.
#' @param y The \eqn{n} by 1 response vector
#' @param lambda The penalty parameter euqals either 2 in AIC or ln(n) in BIC (by default).
#' It can be specified as any value of the user's own choice.
#' @param a The scale parameter in the hyperbolic tangent function of the MIC penalty. By default, \eqn{a = 50}.
#' @param family a description of the error distribution. To use this function, \code{family} needs to be one of \code{"gaussian"},
#'
#' @return The value of the penalized log likelihood function evaluated at beta.
#' @details
#' This function is much faster than \code{\link{LoglikPenGLM}}, but it is only applicable for Gaussian linear regression, logistic regression,
#' and loglinear or Poisson regression models. To take advantage, sepcify  \code{family} as \code{family="gaussian"}, \code{family="binomial"}, or \code{family="poisson"} only.
#' @seealso \code{\link[stats]{glm}}, \code{\link[stats]{family}}
#' @references
#'\itemize{
#' \item Su, X. (2015). Variable selection via subtle uprooting.
#' \emph{Journal of Computational and Graphical Statistics}, \bold{24}(4): 1092--1113.
#' URL \url{http://www.tandfonline.com/doi/pdf/10.1080/10618600.2014.955176}
#' \item Su, X., Fan, J., Levine, R. A., Nunn, M. E., and Tsai, C.-L. (2016+). Sparse estimation of generalized linear
#' models via approximated information criteria. Submitted, \emph{Statistica Sinica}.
#' }


LoglikPen <- function(beta, I.preselect=NULL, preselect.beta0=FALSE, X, y, lambda, a,
                      family = "gaussian")
{
  # PRE-SELECTED VARIABLES
  if (is.null(I.preselect))  {
	I.preselect <- rep(1, length(beta))
  	if (preselect.beta0) I.preselect[1] <- 0
  }
  # THE PANALTY PART
  w <- tanh(a*beta^2); w[as.logical(1-I.preselect)] <- 1; beta.prime <- beta*w
  pen <- lambda*sum(w*I.preselect)
  # pen <- ifelse(select.beta0, lambda*sum(w*I.preselect), lambda*sum((w*I.preselect)[-1]))

  # THE (-2)*LOGLIKELIHOOD PART
  eta <- X%*%beta.prime
  if (family=="gaussian")
    loglik <- NROW(X)* log(sum((y - eta)^2))
  else if (family=="binomial")
    loglik <- -2*sum(y*eta - log(1+exp(eta)))
  else if (family=="poisson")
    loglik <- -2*sum(y*eta - exp(eta))
  else stop("Hmmm. How did you get here? Wrong specification of family.")

  # OUTPUT THE PENALZIED LOGLIKELIHOOD OR SMOOTHED AIC/BIC
  loglik.pen <- loglik + pen;
  return(loglik.pen)
}
