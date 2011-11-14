#' Estimate between-group, within-group, and overall correlation in grouped data
#' using mixed-effect models.
#'
#' Given grouped data, \code{corme} will estimate between-group, within-group,
#' and overall correlation by fitting each pair of variables with a mixed-effect
#' regression model where:
#' \itemize{
#' \item{group means are drawn from a bivariate normal distribution}
#' \item{observations within a group also have a bivariate normal distribution}
#' }
#'
#' As with \code{cor}, you can specify either a single set of variables
#' (\code{x}) or a second set (\code{y}) as well.  Only a single grouping
#' factor is currently supported.
#'
#' Significance testing is by chi-square likelihood ratio tests, which may
#' be inaccurate at small sample sizes.  Note that the significance test
#' \code{p.total} may be confusing: it is a model comparison between a model
#' with both between-group and within-group correlation and a null model with
#' no correlation.  It does \emph{not} test a hypothesis that \code{r.total}
#' is zero.  If between-group and within-group correlations have different signs,
#' then it is possible that \code{r.total} could be near zero while \code{p.total}
#' is highly significant.
#'
#' @title corme: correlation with mixed-effect models
#' @rdname corme
#' @export corme
#' @param ... additional arguments, not currently used
#' @return a list, where each item is either a scalar value (if we are
#'   comparing two variables only) or a matrix:
#'   \describe{
#'   \item{r.group}{correlation between groups}
#'   \item{p.group}{significance of r.group}
#'   \item{r.obs}{correlation within groups (accounting for r.group)}
#'   \item{p.obs}{significance of r.obs}
#'   \item{r.total}{overall correlation from all sources}
#'   \item{p.total}{overall signifance of any correlation in the model}
#'   }
#' @references A. Hamlett, L. Ryan, P. Serrano-Trespalacios, R. Wolfinger,
#'   Journal of the Air & Waste Management Association (1995) 53, 442 (2003).
corme <- function(...) UseMethod('corme')

corme1 <- function(g, x, y, REML=T) {
  stopifnot(length(x)==length(y), length(x)==length(g))
  group <- factor(rep(g, 2))
  obs <- factor(rep(1:length(x), 2))
  var <- factor(rep(c('x','y'), each=length(x)))
  val <- c(x,y)
  model <- lmer(val ~ 0+var +
                (0+var | obs:group) +
                (0+var | group),
                REML=REML)
  varx <- as.numeric(var=='x')
  vary <- as.numeric(var=='y')
  model0.group <- lmer(val ~ 0+var +
                       (0+var | obs:group) +
                       (0+varx | group) + (0+vary | group),
                       REML=REML)
  model0.obs <- lmer(val ~ 0+var +
                     (0+varx | obs:group) + (0+vary | obs:group) +
                     (0+var | group),
                     REML=REML)
  model0.both <- lmer(val ~ 0+var +
                      (0+varx | obs:group) + (0+vary | obs:group) +
                      (0+varx | group) + (0+vary | group),
                      REML=REML)
  pval <- function(m0, m1) anova(m0, m1)$'Pr(>Chisq)'[2]
  vc <- VarCorr(model)
  var.epsilon <- attr(vc, 'sc')**2
  vc.group <- vc$group
  vc.obs <- vc$'obs:group' + diag(var.epsilon, 2, 2)
  vc.total <- vc.group + vc.obs
  list(r.group=cov2cor(vc.group)[2],
       p.group=pval(model0.group, model),
       r.obs=cov2cor(vc.obs)[2],
       p.obs=pval(model0.obs, model),
       r.total=cov2cor(vc.total)[2],
       p.total=pval(model0.both, model))
}

#' @rdname corme
#' @method corme default
#' @S3method corme default
#' @import lme4
#' @param g a factor with at least 2 levels identifying which group each
#'   observation is from.
#' @param x a numeric vector, matrix, or data frame
#' @param y an optional second vector, matrix, or data frame.  If absent, this
#'   is the same as specifying y=x (but more than twice as fast).
#' @param REML if true use restricted maximum likelihood estimation
corme.default <- function(g, x, y=NULL, REML=T, ...) {
  init.result <- function(nrow, ncol, rnames, cnames) {
    m0 <- matrix(0, nrow, ncol, dim=list(rnames, cnames))
    m1 <- matrix(0, nrow, ncol, dim=list(rnames, cnames))
    diag(m1) <- 1
    list(r.group=m1, p.group=m0, r.obs=m1,
         p.obs=m0, r.total=m1, p.total=m0)
  }
  if (is.null(y)) {
    x <- as.matrix(x)
    stopifnot(ncol(x) > 1)
    nc <- ncol(x)
    if (nc==2)
      corme1(g, x[,1], x[,2], REML=REML)
    else {
      result <- init.result(nc, nc, colnames(x), colnames(x))
      for (i in seq_len(nc-1))
        for (j in seq(i+1, nc)) {
          rij <- corme1(g, x[,i], x[,j], REML=REML)
          result$r.group[i, j] <- rij$r.group
          result$p.group[i,j] <- rij$p.group
          result$r.obs[i,j] <- rij$r.obs
          result$p.obs[i,j] <- rij$p.obs
          result$r.total[i,j] <- rij$r.total
          result$p.total[i,j] <- rij$p.total
        }
      lapply(result, function(m) m + t(m) - diag(diag(m)))
    }
  } else {
    x <- as.matrix(x)
    y <- as.matrix(y)
    stopifnot(ncol(x) > 0, ncol(y) > 0, nrow(x)==nrow(y))
    nr <- ncol(x)
    nc <- ncol(y)
    if ((nc==1) && (nr==1))
      corme1(g, x[,1], y[,1], REML=REML)
    else {
      result <- init.result(nr, nc, colnames(x), colnames(y))
      for (i in seq_len(nr))
        for (j in seq_len(nc)) {
          rij <- corme1(g, x[,i], y[,j], REML=REML)
          result$r.group[i, j] <- rij$r.group
          result$p.group[i,j] <- rij$p.group
          result$r.obs[i,j] <- rij$r.obs
          result$p.obs[i,j] <- rij$p.obs
          result$r.total[i,j] <- rij$r.total
          result$p.total[i,j] <- rij$p.total
        }
      result
    }
  }
}

#' @rdname corme
#' @method corme formula
#' @S3method corme formula
#' @import Formula
#' @param formula a formula of the form y1+...+yn ~ x1+...+xm | g.  The left-hand
#'   side is optional, and if missing is treated as if it were the same as the
#'   right.
#' @param data a data frame containing the variables in 'formula'
corme.formula <- function(formula, data, REML=T, ...) {
  f <- Formula(formula)
  stopifnot(length(f)[1] <= 1, length(f)[2] == 2)
  y <- if (length(f)[1]==0) NULL else model.part(f, data, lhs=1)
  x <- model.part(f, data, rhs=1)
  g <- model.part(f, data, rhs=2, drop=T)
  corme.default(g, x, y, REML=REML, ...)
}
