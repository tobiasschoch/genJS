\name{EB}
\alias{EB}

\title{Empirical Bayes Estimator (Generalized James--Stein Estimator)}
\usage{
EB(oi, ei, center = FALSE, maxit = 100, tol = 1e-5, verbose = FALSE)
}
\arguments{
    \item{oi}{\code{[numeric vector]} observed counts.}
    \item{ei}{\code{[numeric vector]} expected counts.}
    \item{center}{\code{[logical]} toggle for the estimation of center,
        otherwise the center is assumed known to be zero.}
    \item{maxit}{\code{integer} maximum number of iterative updates for
        Fisher scoring (default: \code{100}).}
    \item{tol}{\code{[numeric]} numerical tolerance criterion to stop the
        iterations (default: \code{1e-05})..}
    \item{verbose}{\code{[logical]} indicating whether additional information
        is printed to the console (default: \code{TRUE}).}
}
\description{
The function \code{EB()} is a wrapper function for
\code{\link[=genJS]{genJS()}}, the generalized James--Stein estimator
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
m <- EB(TURP$oi, TURP$ei, center = FALSE)

# empirical Bayes predictor
predict(m)
}
