\name{Mean Square Error Estimation}
\alias{mse}
\title{Mean Square Error Estimation}
\usage{
mse(x, method = "analytic")
}
\arguments{
    \item{x}{an object of class \code{genJS}.}
    \item{method}{\code{[character]} estimation method \code{"analytic"}
        or \code{"jackknife"}.}
}
\description{
Mean square error estimation: Analytic approximation (Prasad and Rao, 1990)
or jackknife estimator (Jiang et al., 2002).
}
\references{
Jiang, J., Lahiri, P., Wan, S.-M. (2002). A unified Jackknife Theory
    for Empirical Best Prediction With M-Estimation.
    \emph{The Annals of Statistics} \bold{30}, pp. 1782--1810.
    \doi{10.1214/aos/1043351257}
Prasad, N. G. N. and Rao, J. N. K. Rao (1990).
    The Estimation of the Mean Squared Error of Small-Area Estimators.
    \emph{Journal of the American Statistical Association} \bold{85},
    pp. 163--171. \doi{10.2307/2289539},
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

# analytical approxiation of MSE
mse(m, method = "analytic")
}
