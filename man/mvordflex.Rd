% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mvordflex.R
\name{mvordflex}
\alias{mvordflex}
\title{Multivariate Ordinal Regression Models.}
\usage{
mvordflex(
  formula,
  data,
  error.structure = cor_MMO3(~1),
  link = mvprobit(),
  response.levels = NULL,
  coef.constraints = NULL,
  coef.values = NULL,
  threshold.constraints = NULL,
  threshold.values = NULL,
  weights.name = NULL,
  offset = NULL,
  PL.lag = NULL,
  contrasts = NULL,
  control = mvord.control(se = TRUE, solver.optimx.control = list(maxit = 2e+05, trace =
    1, kkt = FALSE))
)
}
\arguments{
\item{formula}{an object of class \code{\link{formula}} of the form \code{y ~ X1 + ... + Xp}.}

\item{data}{\code{\link{data.frame}} containing a subject index, an index for the multiple measurements,
an ordinal response \code{y} and covariates \code{X1, ..., Xp}.}

\item{error.structure}{different \code{error.structures}: structures for three dimensional
panel data are \code{cor_MMO3(~ 1)}, \code{cor_MMO3_ar1(~ 1)}, \code{cor_MMO3_cross(~1)}.
See \code{\link{error_struct}} or 'Details'.}

\item{link}{specifies the link function by \code{mvprobit()} (multivariate normally distributed errors - default)
or \code{mvlogit(df = 8)} (multivariate logistically distributed errors), where \code{df} specifies the degrees of freedom of the t copula.}

\item{response.levels}{(optional) \code{\link{list}} of length equal to the number of multiple measurements to specify the category labels
in case of varying categories across multiple measurements}

\item{coef.constraints}{(optional) \code{\link{vector}} or \code{\link{matrix}} of constraints on the regression coefficients. See 'Details'.}

\item{coef.values}{(optional) \code{\link{matrix}} setting fixed values on the regression coefficients. See 'Details'.}

\item{threshold.constraints}{(optional) \code{\link{vector}} of constraints on the threshold parameters. See 'Details'.}

\item{threshold.values}{(optional) \code{\link{list}} of (optional) fixed values for the threshold parameters. See 'Details'.}

\item{weights.name}{(optional) character string with the column name of subject-specific weights in \code{data} which need to be
constant across multiple measurements. Negative weights are not allowed.}

\item{offset}{(optional) this can be used to specify an a priori known component to be included in the linear predictor during fitting.
This should be NULL or a numeric vector of length equal to the number of cases. One or more offset terms can be included
in the formula instead or as well, and if more than one is specified their sum is used. See model.offset.}

\item{PL.lag}{(optional) specifies the time lag of the pairs in the pairwise likelihood approach to be optimized.}

\item{contrasts}{(optional) an optional list. See the \code{contrasts.arg} of \code{\link{model.matrix.default}}.}

\item{control}{(optional) a list of parameters for controlling the fitting process. See \code{\link{mvord.control}} for details.}
}
\value{
The function \code{mvordflex} returns an object of \code{\link{class}} \code{"mvordflex"}.

The functions \code{summary} and \code{print} are used to display the results.
The function \code{coef} extracts the regression coefficients, a function \code{thresholds} the threshold coefficients
and the function \cr
\code{error_structure} returns the estimated parameters of the corresponding error structure.

An object of \code{\link{class}} \code{"mvordflex"} is a list containing the following components:

\itemize{
 \item{\code{beta}}{

 a named \code{\link{matrix}} of regression coefficients}
 \item{\code{theta}}

 a named \code{\link{list}}{ of threshold parameters}
  \item{\code{error.struct}}{

  an object of class \code{\link{error_struct}} containing the parameters of the error
  structure}
  \item{\code{sebeta}}{

    a named \code{\link{matrix}} of the standard errors of the regression coefficients}
  \item{\code{setheta}}{

    a named \code{\link{list}} of the standard errors of the threshold parameters}
  \item{\code{seerror.struct}}{

  a \code{vector} of standard errors for the parameters of the error structure}
  \item{\code{rho}}{

    a \code{\link{list}} of all objects that are used in \code{mvordflex()}}
}
}
\description{
Multivariate ordinal regression models for three dimensional panel data in
the R package mvordflex can be fitted using the function \code{mvordflex()}.
The data can be passed on to \code{mvordflex()} through the use
of multiple measurement objects MMO3 in the left-hand side of the model
formula. \code{MMO3} requires the specification of the name of the response
variable, the subject ID, the time ID and the multiple
measurement ID. \code{MMO3} uses a long data format, which has the advantage
that it allows for varying covariates across multiple measurements across
different time points.
}
\details{
\describe{
\item{Implementation \code{MMO3}:}{
  \itemize{

    \item{\code{data}:}{

    In MMO3 we use a long format for the input of data, where each row
    contains a subject index (\code{i}), a time index (\code{t}), a multiple
    measurement index (\code{j}), an ordinal response (\code{Y}) and all
    covariates (\code{X1} to \code{Xp}). This long format data structure
    is internally transformed to a matrix of responses which contains NA
    in the case of missing entries and a list of covariate matrices.
    This is performed by the multiple measurement object \code{MMO3(Y, i, t, j)}
    specifying the column names of the subject, time and the multiple
    measurement index in data. The column containing the ordinal observations
    can contain integer or character values or can be of class (ordered)
    `factor'. When using the long data structure, this column is basically
    a concatenated vector of each of the multiple ordinal responses.
    Internally, this vector is then split according to the multiple measurement and
    time indices. Then the ordinal variable corresponding to each multiple measurement index
    and time index is transformed into an ordered
    factor. For an integer or a character vector the natural ordering is used (ascending,
    or alphabetical). If for character vectors the alphabetical order does not correspond to
    the ordering of the categories, the optional argument response.levels allows to specify the
    levels for each response explicitly. This is performed by a list of length q * T
    (number of multiple measurements multiplied to the number of time points), where each
    element contains the names of the levels of the ordered categories in ascending
    (or if desired descending) order. If all the multiple measurements use the same number
    of classes and same labelling of the classes, the column Y can be stored as an ordered `factor'.
    The order of the multiple measurements is needed when specifying constraints on the
    threshold or regression parameters. This order is based on the type of the multiple
    measurement index column in data. For `integer', `character' or `factor' the
    natural ordering is used (ascending, or alphabetical). If a different
    order of the multiple responses is desired, the multiple measurement index
    column should be an ordered factor with a corresponding ordering of the levels.
    If the categories differ across multiple measurements (either the number of categories
    or the category labels) one needs to specify the response.levels explicitly. This is
    performed by a list of length q * T (number of multiple measurements multiplied to
    the number of time points), where each element contains the names of the levels of the
    ordered categories in ascending or descending order.}
\preformatted{response.levels = list(c("G","F","E", "D", "C", "B", "A"),
                       c("G","F","E", "D", "C", "B", "A"),
                       c("O","N","M","L", "K", "J", "I", "H"))}

\item{\code{formula}}{
The ordinal responses (e.g., rating) are passed by a formula object. Intercepts
can be included or excluded in the model depending on the model parameterization:
\itemize{
\item {Model without intercept:} If the intercept should be removed the \code{formula} for a given response (\code{rating})
and covariates (\code{X1} to \code{Xp}) has the following form:

     \code{formula = MMO3(rating,firm_id,year_id,rater_id) ~ 0 + X1 + ... + Xp}.

\item {Model with intercept:} If one wants to include an intercept in the model, there are two equivalent possibilities
to set the model \code{formula}. Either one includes the intercept explicitly by:

    \code{formula = MMO3(rating,firm_id,year_id,rater_id) ~ 1 + X1 + ... + Xp},

or by

  \code{formula = MMO3(rating,firm_id,year_id,rater_id) ~ X1 + ... + Xp}.
}
}
}
}
  \item{\code{error.structure}}{
 We allow for different error structures depending on the model parameterization:
\itemize{
  \item {Correlation:}
  \itemize{
  \item \code{cor_MMO3}
     This parameterization is available for responses which exhibit correlation
     in the cross-section as well as longitudinally. The number of parameters
      to be estimated is q(q-1)/2 + q, where the first term corresponds to
      the cross-sectional correlation parameters and the second one to the
      number of auto-regressive parameters for the q outcomes.

       \code{error.structure = cor_MMO3(~ 1)}
}

  }
  }

  \item{\code{coef.constraints}}{
  The package supports
  constraints on the regression coefficients. Firstly, the
  user can specify whether the regression coefficients should be equal
  across some or all response dimensions. Secondly, the values of some
  of the regression coefficients can be fixed.

  As there is no unanimous way to specify such constraints, we offer
  two options. The first option is similar to the specification of constraints on the thresholds.
   The constraints can be specified in this case as a vector or matrix of integers,
    where coefficients getting same integer value are set equal.
  Values of the regression coefficients can be fixed through a matrix.
  Alternatively constraints on the regression coefficients can be specified
  by using the design employed by the \pkg{VGAM} package.
  The constraints in this setting are set through a named list,
  where each element of the list contains a matrix full-column rank.
  If the values of some regression coefficients should be fixed, offsets can be used.
  This design has the advantage that it supports
  constraints on outcome-specific as well as category-specific
  regression coefficients. While the first option has the advantage of requiring a more concise input,
   it does not support category-specific coefficients.
  The second option offers a more flexible design in this respect. For further information
  on the second option we refer to the vignette and to the documentation of \code{\link[VGAM]{vglm}}.

Using the first option, constraints can be specified by a vector or a matrix \cr
   \code{coef.constraints}.
    First, a simple and less flexible way by specifying a vector \cr
    \code{coef.constraints}
     of dimension \eqn{q*T}. See \code{\link[mvord]{mvord}}.
 }
  \item{\code{coef.values}}{
  In addition, specific values on regression coefficients can be set in the matrix \cr
  \code{coef.values}. See \code{\link[mvord]{mvord}}.
}


  \item{\code{threshold.constraints}}{
  Similarly, constraints on the threshold parameters can be imposed by a vector of positive integers,
   where dimensions with equal threshold parameters get the same integer. When restricting the thresholds of two
    outcome dimensions to be the same, one has to be careful that the number of categories in
     the two outcome dimensions must be the same. See \code{\link[mvord]{mvord}}.
}

  \item{\code{threshold.values}}{
  In addition, threshold parameter values can be specified by \code{threshold.values}
   in accordance with identifiability constraints. For this purpose we use a \code{list}
    with \eqn{q*T} elements, where each element specifies the constraints of the particular
     dimension by a vector of length of the number of threshold parameters (number of categories - 1).
     A number specifies a threshold parameter to a specific value and \code{NA} leaves the parameter flexible.
}
}
}
\examples{
library(mvordflex)
data(data_toy_mvordflex)
q <- 3 # number of multiple measurements
TT <- 5 # number of years
# Note that the number of responses is q*TT
res <- mvordflex(
  formula = MMO3(response, firm_id, year_id, outcome_id) ~ 0 + X1 + X2,
  data = data_toy_mvordflex,
  coef.constraints = rep(1, q * TT),
  threshold.constraints = rep(1:q, TT))
print(res)
summary(res)
thresholds(res)
coefficients(res)
head(error_structure(res))

}
\references{
Hirk R, Hornik K, Vana L (2020). “\pkg{mvord}: An \strong{R} Package for Fitting Multivariate Ordinal Regression Models.” \emph{Journal of Statistical Software}, \strong{93}(4), 1–41. doi: 10.18637/jss.v093.i04 (URL: \url{https://doi.org/10.18637/jss.v093.i04}).
}
\seealso{
\code{\link{print.mvord}}, \code{\link{summary.mvord}}, \code{\link{coef.mvord}},
 \code{\link{thresholds.mvord}}, \code{\link{error_structure.mvord}}, \cr
 \code{\link{mvord.control}}, \code{\link{data_toy_mvordflex}}
}
