\name{genJS}
\alias{genJS}
\alias{predict.genJS}
\alias{print.genJS}
\alias{summary.genJS}
\alias{coef.genJS}
\alias{vcov.genJS}
\title{Generalized James--Stein Estimator}
\usage{
# Generalized James-Stein estimator
genJS(yi, di, center = TRUE, maxit = 100, tol = 1e-5, verbose = TRUE)

# Empirical Bayes predictor
\method{predict}{genJS}(object, morris = FALSE, ...)

# Utility functions
\method{print}{genJS}(x, ...)
\method{summary}{genJS}(object, ...)
\method{coef}{genJS}(object, ...)
\method{vcov}{genJS}(object, ...)
}
\arguments{
    \item{yi}{\code{[numeric vector]} direct estimator.}
    \item{di}{\code{[numeric vector]} variance of direct estimator.}
    \item{center}{\code{[logical]} toggle for the estimation of center,
        otherwise the center is assumed known to be zero.}
    \item{maxit}{\code{integer} maximum number of iterative updates for
        Fisher scoring (default: \code{100}).}
    \item{tol}{\code{[numeric]} numerical tolerance criterion to stop the
        iterations (default: \code{1e-05})..}
    \item{verbose}{\code{[logical]} indicating whether additional information
        is printed to the console (default: \code{TRUE}).}
    \item{morris}{\code{[logical]} toggle to use the bias-corrected predictor
        of Morris (1983).}
    \item{x}{an object of class \code{genJS}.}
    \item{object}{an object of class \code{genJS}.}
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
\code{predict()}.  Morris (1983) proposed a bias-corrected empirical Bayes
predictor (for small samples); this predictor obtains by using
\code{morris = TRUE} in the call.

The "classical" James--Stein estimator assumes that the variances (\code{di})
are constant.
}
\references{
Efron, B. und Morris, C. (1973). Stein's Estimation Rule and Its
    Competitors: An Empirical Bayes Approach. \emph{Journal of the
    American Statistical Association} \bold{68}, pp. 117--130.
    \doi{10.1080/01621459.1973.10481350}

Efron, B. und Morris, C. (1975). Data Analysis Using Stein's Estimator
    and its Generalizations. \emph{Journal of the American Statistical
    Association} \bold{70}, pp. 311--319. \doi{10.1080/01621459.1975.10479864}

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
