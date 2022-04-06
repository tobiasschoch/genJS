\name{jackknife_var}
\alias{jackknife_var}
\title{Jackknife Variance Estimator}
\usage{
jackknife_var(x)
}
\arguments{
    \item{x}{an object of class \code{genJS}.}
}
\description{
Jackknife estimator of variance \eqn{A}{A}
}
\references{
Jiang, J., Lahiri, P., Wan, S.-M. (2002). A unified Jackknife Theory
    for Empirical Best Prediction With M-Estimation.
    \emph{The Annals of Statistics} \bold{30}, pp. 1782--1810.
    \doi{10.1214/aos/1043351257}
}
\author{
Tobias Schoch
}
\seealso{
    \code{\link[=genJS]{genJS()}} generalized James--Stein estimator
}
\examples{
data(TURP)

# parameter estimation (center kept fix at zero)
m <- genJS(TURP$yi, TURP$di, center = FALSE)

# jackknife variance estiamtor
jackknife_var(m)
}
