\name{genJS}
\alias{genJS}
\alias{predict.genJS}
\alias{print.genJS}
\alias{summary.genJS}
\alias{coef.genJS}
\alias{vcov.genJS}
\alias{pvalue_A}
\title{Generalized James--Stein Estimator}
\usage{
# Generalized James-Stein estimator
genJS(yi, di, center = TRUE, maxit = 100, tol = 1e-5, verbose = TRUE)

# Empirical Bayes predictor
\method{predict}{genJS}(object, method = "EB", alpha = NULL, ...)

# p-value under the hypothesis H0: A = 0
pvalue_A(yi, di, center = TRUE)

# Utility functions
\method{print}{genJS}(x, digits = max(1L, getOption("digits") - 2L), ...)
\method{summary}{genJS}(object, ..., digits = max(1L,
    getOption("digits") - 2L))
\method{coef}{genJS}(object, ...)
\method{vcov}{genJS}(object, ...)
}
\arguments{
    \item{yi}{\code{[numeric vector]} direct estimator.}
    \item{di}{\code{[numeric vector]} variance of direct estimator.}
    \item{center}{\code{[logical]} toggle for the estimation of center,
        otherwise the center is assumed known to be zero.}
    \item{method}{\code{[character]} estimation method: \code{"EB"},
        \code{"orris"}, \code{"pretest"}.}
    \item{alpha}{\code{[numeric]} used in pretest estimator,
        \eqn{0 < \alpha 1}{0 < \alpha < 1}.}
    \item{maxit}{\code{integer} maximum number of iterative updates for
        Fisher scoring (default: \code{100}).}
    \item{tol}{\code{[numeric]} numerical tolerance criterion to stop the
        iterations (default: \code{1e-05})..}
    \item{verbose}{\code{[logical]} indicating whether additional information
        is printed to the console (default: \code{TRUE}).}
    \item{x}{an object of class \code{genJS}.}
    \item{object}{an object of class \code{genJS}.}
    \item{digits}{number of digits.}
    \item{...}{additional arguments.}
}
\description{
Generalized James--Stein maximum-likelihood (M) or restricted ML (REML)
estimator
}
\details{
Function \code{genJS()} implements the generalized James--Stein
maximum-likelihood (M) or restricted ML (REML) estimator; see Efron and Morris
(1973, 1975).
    \itemize{
        \item \code{center = FALSE}: Location is assumed known; shrinkage
            toward the origin. The variance is estimated by maximum-likelihood
            (ML)
        \item \code{center = TRUE}: Location is estimated together with
            variance using the restricted ML estimator (REML)
    }

The components estimators (empirical Bayes predictors) are computed using
\code{predict()}.
    \itemize{
        \item \code{method = "EB" } is the Empirical Best predictor
        \item \code{method = "morris"} is a bias-corrected EB predictor
            suggested in Morris (1983)
        \item \code{method = "pretest"} is a pretest predictor; see Datta
            et al. (2011) and Molina et al. (2015)
    }

The "classical" James--Stein estimator assumes that the variances (\code{di})
are constant.

In the \code{summary} method, the p-value for the estimated variance is
computed under the hypotheses \eqn{H_0: A = 0}{H0: A = 0} and
\eqn{H_A: A > 0}{HA: A > 0} (Datta et al., 2011); see also function
\code{pvalue_A}.
}
\references{
Datta, G. S., Hall, P. and Mandal, A. (2011).
    Model Selection by Testing for the Presence of
    Small-Area Effects, and Application to Area-Level.
     \emph{Journal of the American Statistical Association} \bold{106},
     pp. 362--374. \doi{10.1198/jasa.2011.tm10036}

Efron, B. and Morris, C. (1973). Stein's Estimation Rule and Its
    Competitors: An Empirical Bayes Approach. \emph{Journal of the
    American Statistical Association} \bold{68}, pp. 117--130.
    \doi{10.1080/01621459.1973.10481350}

Efron, B. and Morris, C. (1975). Data Analysis Using Stein's Estimator
    and its Generalizations. \emph{Journal of the American Statistical
    Association} \bold{70}, pp. 311--319. \doi{10.1080/01621459.1975.10479864}

Molina, I., Rao, J. N. K. and  Datta, G. S. (2015).
    Small area estimation under a Fay-Herriot model with preliminary
    testing for the presence of random area effects.
    \emph{Survey Methodology} \bold{41}, pp. 1--19.

Morris, C.N. (1983). Parametric Empirical Bayes Inference: Theory and
    Applications. \emph{Journal of the American Statistical Association}
    \bold{78}, pp. 47--55. \doi{10.2307/2287098}
}
\author{
Tobias Schoch
}
\seealso{
    \code{\link[=mse]{mse()}} for MSE estimation
}
\examples{
data(TURP)

# parameter estimation (center kept fix at zero)
m <- genJS(TURP$yi, TURP$di, center = FALSE)

# empirical Bayes predictor
predict(m)
}
