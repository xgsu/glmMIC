#' Compute the penalized log likelihood for GLM with MIC penalty via R function \code{glm}
#'
#' @param beta A p-dimensional vector containing the regression ceofficients.
#' @param I.preselect A 0-1 indicator vector of same length as \code{beta} showing preselected variables (value 0 for no penalty) that will not be penalized.
#' @param preselect.beta0 Indicator of whether or not the intercept is pre-selected. Default is \code{FALSE}.
#' @param X An \eqn{n} by \eqn{p} design matrix.
#' @param y The \eqn{n} by 1 response vector
#' @param lambda The penalty parameter euqals either 2 in AIC or ln(n) in BIC (by default).
#' It can be specified as any value of the user's own choice.
#' @param a The scale parameter in the hyperbolic tangent function of the MIC penalty. By default, \eqn{a = 50}.
#' @param family a description of the error distribution and link function to be used in the model. It needs to be the result
#' of a call to a family function since \code{glm.fit} is used here. In other words, it can NOT be specified as a character
#' string naming a family function nor a family function. See \code{\link[stats]{glm}} and \code{\link[stats]{family}} for
#' details of family functions.
#' @return The value of the penalized log likelihood function evaluated at beta.
#' @seealso \code{\link[stats]{glm}}, \code{\link[stats]{family}}
#' @references
#'\itemize{
#' \item Su, X. (2015). Variable selection via subtle uprooting.
#' \emph{Journal of Computational and Graphical Statistics}, \bold{24}(4): 1092--1113.
#' URL \url{http://www.tandfonline.com/doi/pdf/10.1080/10618600.2014.955176}
#' \item Su, X., Fan, J., Levine, R. A., Nunn, M. E., and Tsai, C.-L. (2016+). Sparse estimation of generalized linear
#' models via approximated information criteria. Submitted, \emph{Statistica Sinica}.
#' }

LoglikPenGLM <- function(beta, I.preselect=NULL, preselect.beta0=FALSE,
                           X, y, lambda, a=20,
                           family = gaussian(link = "identity"))
{
  # PRE-SELECTED VARIABLES
  if (is.null(I.preselect))  {
	I.preselect <- rep(1, length(beta))
  	if (preselect.beta0) I.preselect[1] <- 0
  }
  # THE PANALTY PART
  w <- tanh(a*beta^2); w[as.logical(1-I.preselect)] <- 1; beta.prime <- beta*w
  pen <- lambda*sum(w*I.preselect)
  # pen <- ifelse(preselect.beta0, lambda*sum(w*I.preselect), lambda*sum((w*I.preselect)[-1]))

  # THE LOGLIKELIHOOD PART OBTAINED FROM R FUNCTION glm()
  eta <- X%*%beta.prime; n <- NROW(X);
  fit <- eval(call("glm.fit", x = matrix(, n, 0L), y = y, offset = eta,
                   family = family,
                   control = glm.control(epsilon = 1e+8, maxit = 0.5, trace = T),
                   intercept = F))
  dev <- fit$deviance

  # SIMILAR YET SLOWER WAY
  # dev <- glm(formula=y~offset(eta)- 1, family = family, etastart=eta,
  #	control = glm.control(epsilon = 1e+8, maxit = 1, trace = TRUE))$"deviance"

  # OUTPUT THE PENALZIED LOGLIKELIHOOD
  loglik.pen <- dev + pen;  # print(loglik.pen)
  return(loglik.pen)
}
