#' Sparse Estimation of a GLM via Minimum approximated Information Criterion
#'
#' @name glmMIC
#' @aliases glm.MIC
#' @aliases glmmic
#' @param formula An object of class \code{\link[stats]{formula}}, with the response on the left of a \code{~} operator, and the terms on the right.
#' @param preselection A formula of form, e.g., \code{~ x1 + x2} that gives the pre-selected variables. In this case, no penalty will be applied to their slope parameters. These variables must be contained in the \code{formula} argument; otherwise, an error message shows up.
#' @param family a description of the error distribution and link function to be used in the model. Preferably for computational speed,
#' this is a character string naming a family function among the following three choices: \code{"gaussian"},
#' \code{"binomial"}, or \code{"poisson"}.  Otherwise, it has to be a family function or the result of a call to a family function that
#' can be called for by \code{\link[stats]{glm.fit}}. See \code{\link[stats]{family}} for details of family functions.
#' @param data A data.frame in which to interpret the variables named in the \code{formula} argument.
#' @param beta0 User-supplied beta0 value, the starting point for optimization. If missing or \code{NULL} (by default), the maximum likelihood estimator (MLE) will be used.
#' @param preselect.intercept A logical value indicating whether the intercept term is pre-selected. By default, it is \code{FALSE}.
#' @param criterion Specifies the model selection criterion used. If \code{"AIC"}, the complexity penalty parameter (lambda)
#' equals 2; if \code{"BIC"}, lambda equals ln(n), where n is the sample size. You may specify
#' the penalty parameter of your choice by setting \code{lambda0}.
#' @param lambda0 User-supplied penalty parameter for model complexity. If \code{method="AIC"} or \code{"BIC"}, the value
#' of \code{lambda0} will be ignored.
#' @param a0 The scale (or sharpness) parameter used in the hyperbolic tangent penalty. By default, \code{a0=min(n, 100)} is used.
#' @param rounding.digits Number of digits after the decimal point for rounding-up estiamtes. Default value is 4.
#' @param use.GenSA Logical value indicating if the generalized simulated annealing \code{GenSA} is used. The default is \code{FALSE}.
#' @param lower The lower bounds for the search space in \code{GenSA}. The default is -10 (\eqn{p} by \eqn{1} vector).
#' @param upper The upper bounds for the search space in \code{GenSA}. The default is +10 (\eqn{p} by \eqn{1} vector).
#' @param maxit.global  Maximum number of iterations allowed for the global optimization algorithm \code{SANN}. Default value is 100.
#' @param maxit.local Maximum number of iterations allowed for the local optimizaiton algorithm \code{BFGS}. Default value is 100.
#' @param epsilon Tolerance level for convergence. Default is 1e-6.
#' @param se.gamma Logical indicator of whether the standard error for \code{gamma} is computed. Default is \code{TRUE}.
#' @param CI.gamma Logical indicator of whether the confidence inverval for \code{gamma} is outputed. Default is \code{TRUE}.
#' @param conf.level Specifies the confidence level for \code{CI.gamma}. Defaulted as 0.95.
#' @param se.beta Logical indicator of whether the (post-selection) standard error for \code{beta} is computed. Default is \code{TRUE}.
#' @param fit.ML  Logical indicator of whether we fit the best selected model with full iteration of maximum likelihood (ML). Default is \code{FALSE}.
#' @param details Logical value: if \code{TRUE}, detailed results will be printed out when running \code{coxphMIC}.
#' @details
#' The main idea of MIC involves approximation of the l0 norm with a continuous or smooth
#' unit dent function. This method bridges the best subset selection and regularization by
#' borrowing strength from both. It mimics the best subset selection using a penalized likelihood
#' approach yet with no need of a tuning parameter.
#'
#' The problem is further reformulated with a reparameterization step by relating \code{beta}
#' to \code{gamma}. There are two benefits of doing so:  first, it reduces the optimization to
#' one unconstrained nonconvex yet smooth programming problem, which can be solved efficiently
#' as in computing the maximum likelihood estimator (MLE); furthermore, the
#' reparameterization tactic yields an additional advantage in terms of circumventing post-selection inference.
#' Significance testing on \code{beta} can be done through \code{gamma}.
#'
#' To solve the smooth yet nonconvex optimization, two options are available. The first is a simulated annealing (\code{method="SANN"} option
#' in \code{\link[stats]{optim}}) global optimization algorithm is first applied. The resultant estimator is then used
#' as the starting point for another local optimization algorithm, where the quasi-Newton BFGS method (\code{method="BFGS"}
#' in \code{\link{optim}}) by default. Optionally, the generalized simulated annealing, implemented in \code{\link[GenSA]{GenSA}},
#' can be used instead. This latter approach tends to be slower. However, it does not need to be combined with another local optimization;
#' besides, it often yields the same final solution with different runs. Thus, when \code{use.GenSA=TRUE},
#' the output includes \code{opt.global} only, without \code{opt.local}.
#'
#' In its current version, some appropriate data preparation might be needed. Most important of all,  X variables in all scenarios need to be
#' standardized or scaled. In the case of Gaussian linear regression, the response variable needs to be centered or even standardized.
#' In addition, missing values would cause errors too and hence need prehanlding too.
#'
#' @return An object of class \code{glmMIC} is returned, which may contain the following components depending on the options.
#' \describe{
#' \item{opt.global}{Results from the preliminary run of a global optimization procedure (\code{SANN} as default).}
#' \item{opt.local}{Results from the second run of a local optimization procedure (\code{BFGS} as default).}
#' \item{min.Q}{Value of the minimized objective function.}
#' \item{gamma}{Estimated gamma (reparameterized);}
#' \item{beta}{Estimated beta;}
#' \item{VCOV.gamma}{The estimated variance-covariance matrix for the (reparameterized) gamma estimate;}
#' \item{se.gamma}{Standard errors for the gamma estimate;}
#' \item{VCOV.beta}{The estimated variance-covariance matrix for the beat estimate;}
#' \item{se.beta}{Standard errors for the beta estimate (post-selection);}
#' \item{BIC}{The BIC value for the \emph{selected} model;}
#' \item{result}{A summary table of the fitting results;}
#' \item{fit.ML}{The \code{glm} fitting results with the selected model with full ML iterations;}
#' \item{call}{the matched call.}
#' }
#' @examples
#'   # Note that glmMIC works with standardized data only. See below for examples.
#'   # GAUSSIAN LINEAR REGRESSION
#'   library(lars); data(diabetes);
#'   dat <- cbind(diabetes$x, y=diabetes$y)
#'   dat <- as.data.frame(scale(dat))
#'   fit.MIC <- glmMIC(formula=y~.-1, family="gaussian", data=dat)
#'   names(fit.MIC)
#'   print(fit.MIC)
#'   plot(fit.MIC)
#'   # WITH PRE-SELECTED VARIABLES AND A DIFFERENT a VALUE
#'   fit.MIC <- glmMIC(formula=y~.-1, preselection=~age+sex, family="gaussian", a0=20, data=dat)
#'   fit.MIC
#'
#'   # LOGISTIC REGRESSION
#'   library(ncvreg); data(heart)
#'   dat <- as.data.frame(cbind(scale(heart[, -10]), chd=heart$chd)); names(dat)
#'   dat <- dat[, -c(4, 6)]; head(dat)
#'   fit.MIC <- glmMIC(formula= chd~., data=dat, family = "binomial")
#'   fit.MIC
#'
#'  # LOGLINEAR REGRESSION
#'  fish <- read.csv("http://www.ats.ucla.edu/stat/data/fish.csv")
#'  form <- count ~ . -1 + xb:zg
#'  y <- fish[,names(fish)==as.character(form)[2]]
#'  X <- model.matrix(as.formula(form),fish)
#'  dat <- data.frame(scale(X), count=fish$count); head(dat)
#'  fit.MIC <- glmMIC(formula= count~ ., family = "poisson", data=dat)
#'  fit.MIC
#'
#' @seealso \code{\link[stats]{glm}}, \code{\link[glmMIC]{print.glmMIC}}, \code{\link[glmMIC]{plot.glmMIC}},
#' @references
#'\itemize{
#' \item Su, X. (2015). Variable selection via subtle uprooting.
#' \emph{Journal of Computational and Graphical Statistics}, \bold{24}(4): 1092--1113.
#' URL \url{http://www.tandfonline.com/doi/pdf/10.1080/10618600.2014.955176}
#' \item Su, X., Fan, J., Levine, R. A., Nunn, M. E., and Tsai, C.-L. (2016+). Sparse estimation of generalized linear
#' models via approximated information criteria. Submitted, \emph{Statistica Sinica}.
#' }
#' @import numDeriv
#' @import GenSA
#' @export

glmMIC <- function(formula, preselection=NULL, family = c("gaussian", "binomial", "poisson"), data,
                   beta0=NULL, preselect.intercept=FALSE,
                   criterion ="BIC", lambda0=0,
                   a0=NULL, rounding.digits = 4,		# ROUNDING DIGITS FOR FINAL RESULTS
                   use.GenSA=FALSE, lower=NULL, upper=NULL,
                   maxit.global=100, maxit.local=100, epsilon=1e-6, 									# CONVERGENCE TOLERANCE
                   se.gamma=TRUE, CI.gamma=TRUE, conf.level=0.95,
                   se.beta=TRUE, fit.ML=FALSE, details=FALSE)
{
  call <- match.call()
  # CHECK THE family= ARGUMENT
  family0 <- family
  if (is.character(family)) family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) family <- family()
  if (is.null(family$family)) {print(family); stop("'family0' not recognized")}

  # CHECK THE data= ARGUMENT
  if (missing(data)) data <- environment(formula)
  # OBTAIN THE DESIGN MATRIX X AND RESPONSE y
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  Y <- model.response(mf, "any")
  if (length(dim(Y)) == 1L) {
    nm <- rownames(Y); dim(Y) <- NULL
    if (!is.null(nm)) names(Y) <- nm
  }
  yname <- as.character(formula[[2]]); # EXTRACT RESPONSE NAME FROM THE FORMULA
  X <- if (!is.empty.model(mt)) model.matrix(mt, mf)
  else matrix(, NROW(Y), 0L)
  Xnames <- colnames(X)
  if (details) print(X[1:10, ])  ################

  # HANDLE PRE-SELECTED VARIABLES
  I.preselect <- rep(1, length(Xnames))
  if (!is.null(preselection)) {
	x.preselect <- all.vars(preselection)
	if (sum(!is.element(x.preselect, Xnames)) >0) stop("At least one pre-selected variables are NOT in the model. Please check.")
	I.preselect <- as.numeric(!is.element(Xnames, x.preselect))
   }
   if (preselect.intercept) I.preselect[1] <- 0

  # CENTER Y IN GAUSSIAN LINEAR MODELS
  out <- as.list(NULL)

  # DETERMINE lambda
  n <- NROW(Y);
  if (is.null(a0)) a0 <- min(n, 100)
  if (criterion =="BIC") {lambda <- log(n); if (details) print(n)}
  else if (criterion =="AIC") lambda <- 2
  else lambda <- lambda0

  # OBTAIN MLE AS STARTING VALUE
  if (is.null(beta0)) {
    fit <- eval(call("glm.fit", x =X, y = Y, family = family))
    beta0 <- fit$coef
    if (details) print(beta0)
  }
  p <- length(beta0)

  # THE OBJECTIVE FUNCTION
  if (is.element(family0, c("gaussian", "binomial", "poisson"))) fun <- LoglikPen
  else fun <- LoglikPenGLM;
  #
  # GRADIENT
  grad <- NULL

  # USING GENERALIZED SIMULATED ANNEALING {GenSA} PLUS BFGS
  if (use.GenSA) {
    # require(GenSA)
    # THE LOWER AND UPPER BOUNDS FOR SEARCH SPACE IN OPTIMIZATION
    if (is.null(lower)) {lower=rep(-10, p); upper <- rep(10, p)}
    else if (length(lower)==1 && p > 1) {lower <- rep(lower, p); upper <- rep(upper, p)}
    else if (length(lower)!= p) stop("Wrong specification for lower= and upper=. Be aware of its appropriate dimension!")
    # NOTE GenSA MINIMIZES, SAME AS optim
    opt.fit1 <- GenSA(par=beta0, fn=fun, lower = lower, upper = upper,
                      control=list(maxit=maxit.global, nb.stop.improvement=5),
                      I.preselect=I.preselect, preselect.beta0=preselect.intercept, X=X, y=Y, lambda=lambda, a=a0, family = family0)
    min.Q <- opt.fit1$value
    betavec1 <- betavec2 <-  opt.fit1$par;
  } else {
    # OPTIMIZATION USING SIMULATED ANNEALING, FOLLOWED BY BFGS
    opt.fit1 <- optim(par=beta0, fn=fun, gr = grad,
                      method = "SANN", control = list(maxit=maxit.global, trace=F, reltol=epsilon),
                      I.preselect=I.preselect, preselect.beta0=preselect.intercept, X=X, y=Y, lambda=lambda, a=a0, family = family0)
    betavec1 <- opt.fit1$par; #
    opt.fit2 <- optim(par=betavec1, fn=fun, gr = grad,
                      method = "BFGS", control = list(maxit=maxit.local, trace=F, reltol=epsilon),
                      I.preselect=I.preselect, preselect.beta0=preselect.intercept, X=X, y=Y, lambda=lambda, a=a0, family = family0)
    betavec2 <- opt.fit2$par
    out$opt.local <- opt.fit2
    if (details) print(betavec2)
    min.Q <- opt.fit2$value
  }
  out <- c(out, list(opt.global=opt.fit1, min.Q=min.Q, gamma=betavec2))

  # OBTAIN SE FOR GAMMA
  gamma <- betavec2; # print(formula)
  result <- data.frame(beta0=beta0, gamma=gamma)
  if (se.gamma) {
    fit0 <- glm(formula, data = data, family = family,
                start=gamma, control=glm.control(epsilon=1e20, maxit=1, trace=F))
    VCOV.gamma <- summary(fit0)$"cov.scaled"
    out$VCOV.gamma <- VCOV.gamma
    se.gamma <- sqrt(diag(VCOV.gamma))
    out$se.gamma <- se.gamma
    # INFERENCE BASED ON GAMMA, NOT BETA
    z.gamma <- gamma/se.gamma
    pvalue.gamma <- 2*pnorm(abs(z.gamma), lower.tail=F)
    result <-  cbind(result, se.gamma=se.gamma, z.gamma=z.gamma, pvalue.gamma=pvalue.gamma)
    # CONFIDENCE INTERVAL FOR GAMMA (OPTIONAL)
    if (CI.gamma) {
      out$conf.level <- conf.level
      if (conf.level <= 0 || conf.level >=1) stop("The argument conf.level needs to within 0 and 1.")
      sig.level <- 1-conf.level
      z0 <- qnorm(1-sig.level/2)
      LB <- gamma - z0*se.gamma
      UB <- gamma + z0*se.gamma
      result <- cbind(result, LB=LB, UB=UB)
    }
  }

  # PREPARE THE VARIABLE SELECTION RESULTS
  w.hat <- tanh(a0*gamma^2);
  w.hat[as.logical(1-I.preselect)] <- 1;  # PRE-SELECTED VARIABLES HANDLED
  out$w.hat <- w.hat
  beta <- gamma*(w.hat); out$beta <- beta
  result <- cbind(result, w.hat=w.hat, beta=round(beta, rounding.digits))

  beta00 <- beta[abs(beta) > epsilon];
  if (se.beta) {
    # OBTAIN THE STANDARD ERROR OF BETA USING GLM (POST-SELECTION INFERENCE)
    q <- length(beta00)
    if (q==0) se.beta <- rep(NA, p)
    else {
      xnames.selected <- Xnames[abs(beta) > epsilon];
      # print(beta); print(xnames.selected)
      form.bmodel <- ifelse(is.element("(Intercept)", xnames.selected),
                            paste(yname, " ~ ", paste(xnames.selected[-1], collapse= " + ")),
                            paste(yname, " ~ ", paste(xnames.selected, collapse= " + "), "- 1"))
      form.bmodel <- as.formula(form.bmodel); # print(form.bmodel)
      fit <- glm(form.bmodel, data = data, family = family,
                 start=beta00, control=glm.control(epsilon=1e20, maxit=1, trace=F))
      out$BIC <- BIC(fit);
      VCOV.beta <- summary(fit)$"cov.scaled"
      out$VCOV.beta <- VCOV.beta
      se.beta <- sqrt(diag(VCOV.beta));
      # ADD NA FOR UNSELECTED VARIABLES
      tmp <- rep(NA, p); names(tmp) <- Xnames;
      tmp[names(se.beta)] <- se.beta
      se.beta <- tmp
      out$se.beta <- se.beta

      # FITTING THE SELECTED MODEL VIA MLE WITH FULL ITERATIONS
      if (fit.ML) out$fit.ML <- glm(form.bmodel, data = data, family = family);
    }
    result <- cbind(result, se.beta=se.beta, z.beta=beta/se.beta,
                    pvalue.beta=2*pnorm(abs(beta/se.beta), lower.tail=F))
  }
  result <- as.data.frame(result)
  row.names(result) <- Xnames
  out$result <- result
  out$call <- call
  class(out) <- "glmMIC"
  return(out)
}
