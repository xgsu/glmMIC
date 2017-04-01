#' The Generic \code{print} Function for Object of \code{glmMIC} Class
#'
#' @name print.glmMIC
#' @param fit an object of \code{glmMIC} class.
#' @param digits	the minimal number of significant digits. See \code{\link[base]{print.default}}.
#' @param \code{...} further arguments passed to or from other methods.
#' @details
#' The (generic) print method for an \code{glmMIC} object. The results include info on the estimated gamma and beta.
#' Depending on the options, significance testing and confidence intervals are also provided.
#' @return The table of estimated regression coefficients beta and the reparameterized gamma.
#' @references
#'\itemize{
#' \item Su, X. (2015). Variable selection via subtle uprooting.
#' \emph{Journal of Computational and Graphical Statistics}, \bold{24}(4): 1092--1113.
#' URL \url{http://www.tandfonline.com/doi/pdf/10.1080/10618600.2014.955176}
#' \item Su, X., Fan, J., Levine, R. A., Nunn, M. E., and Tsai, C.-L. (2016+). Sparse estimation of generalized linear
#' models via approximated information criteria. Submitted, \emph{Statistica Sinica}.
#' }
#' #' @seealso \code{\link[glmMIC]{glmMIC}}
#' @export

print.glmMIC <- function(fit, digits = max(3L, getOption("digits") - 3L), ...)
{
  cat("\nCall:\n", paste(deparse(fit$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  if (NROW(fit$result)>=1){
    cat("Table of Estimated Coefficients via MIC:\n\n")
    print(fit$result, digits=digits)
    if (!is.null(fit$conf.level)) cat("\n* Note that LB and UB are ",
        paste(fit$conf.level*100, sep="\n", collapse="\n"), "% confidence intervals for gamma. \n\n", sep = "")
  }
  else cat("No coefficients. Wasn't MIC successfully run? \n")
  cat("\n")
  invisible(fit)
}
