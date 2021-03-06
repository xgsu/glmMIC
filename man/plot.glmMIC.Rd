% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot-glmMIC.R
\name{plot.glmMIC}
\alias{plot.glmMIC}
\title{The Generic \code{plot} Function for Object of \code{glmMIC} Class}
\usage{
\method{plot}{glmMIC}(fit, conf.level = 0.95, horizontal = TRUE,
  mar = rep(4.5, 4), ...)
}
\arguments{
\item{fit}{an object of \code{glmMIC} class.}

\item{conf.level}{confidence level used for error bar plots. Default is 0.95.}

\item{horizontal}{Logical indicator of horizontal alignment. Default is \code{TRUE}.}

\item{mar}{margin in terms of the number of lines to be specified on the four sides of the plot}

\item{\code{...}}{further arguments passed to or from other methods.}
}
\value{
Error bar plots for estimated gamma and beta at a given confidenve level.
}
\description{
The Generic \code{plot} Function for Object of \code{glmMIC} Class
}
\details{
The (generic) plot method for an \code{glmMIC} object. It plots MIC estimates of
gamma and beta. For 0 beta estimates, their corresponding
SE are reset to 0 to make the plot.
}
\references{
\itemize{
\item Su, X. (2015). Variable selection via subtle uprooting.
\emph{Journal of Computational and Graphical Statistics}, \bold{24}(4): 1092--1113.
URL \url{http://www.tandfonline.com/doi/pdf/10.1080/10618600.2014.955176}
\item Su, X., Fan, J., Levine, R. A., Nunn, M. E., and Tsai, C.-L. (2016+). Sparse estimation of generalized linear
models via approximated information criteria. Submitted, \emph{Statistica Sinica}.
}
}
\seealso{
\code{\link[glmMIC]{glmMIC}}
}

