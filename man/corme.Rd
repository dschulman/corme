\name{corme}
\alias{corme}
\alias{corme.default}
\alias{corme.formula}
\title{corme: correlation with mixed-effect models}
\usage{
  corme(...)

  \method{corme}{default} (g, x, y = NULL, REML = T, ...)

  \method{corme}{formula} (formula, data, REML = T, ...)
}
\arguments{
  \item{...}{additional arguments, not currently used}

  \item{g}{a factor with at least 2 levels identifying
  which group each observation is from.}

  \item{x}{a numeric vector, matrix, or data frame}

  \item{y}{an optional second vector, matrix, or data
  frame.  If absent, this is the same as specifying y=x
  (but more than twice as fast).}

  \item{REML}{if true use restricted maximum likelihood
  estimation}

  \item{formula}{a formula of the form y1+...+yn ~
  x1+...+xm | g.  The left-hand side is optional, and if
  missing is treated as if it were the same as the right.}

  \item{data}{a data frame containing the variables in
  'formula'}
}
\value{
  a list, where each item is either a scalar value (if we
  are comparing two variables only) or a matrix: \describe{
  \item{r.group}{correlation between groups}
  \item{p.group}{significance of r.group}
  \item{r.obs}{correlation within groups (accounting for
  r.group)} \item{p.obs}{significance of r.obs}
  \item{r.total}{overall correlation from all sources}
  \item{p.total}{overall signifance of any correlation in
  the model} }
}
\description{
  Estimate between-group, within-group, and overall
  correlation in grouped data using mixed-effect models.
}
\details{
  Given grouped data, \code{corme} will estimate
  between-group, within-group, and overall correlation by
  fitting each pair of variables with a mixed-effect
  regression model where: \itemize{ \item{group means are
  drawn from a bivariate normal distribution}
  \item{observations within a group also have a bivariate
  normal distribution} }

  As with \code{cor}, you can specify either a single set
  of variables (\code{x}) or a second set (\code{y}) as
  well.  Only a single grouping factor is currently
  supported.

  Significance testing is by chi-square likelihood ratio
  tests, which may be inaccurate at small sample sizes.
  Note that the significance test \code{p.total} may be
  confusing: it is a model comparison between a model with
  both between-group and within-group correlation and a
  null model with no correlation.  It does \emph{not} test
  a hypothesis that \code{r.total} is zero.  If
  between-group and within-group correlations have
  different signs, then it is possible that \code{r.total}
  could be near zero while \code{p.total} is highly
  significant.
}
\references{
  A. Hamlett, L. Ryan, P. Serrano-Trespalacios, R.
  Wolfinger, Journal of the Air & Waste Management
  Association (1995) 53, 442 (2003).
}

