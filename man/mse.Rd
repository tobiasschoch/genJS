\name{Mean Square Error Estimation}
\alias{mse}
\alias{mse.genJS}
\alias{mse.EB}
\title{Mean Square Error Estimation}
\usage{
mse(x, ...)

\method{mse}{genJS}(x, method = "analytic", alpha = NULL, ...)
\method{mse}{EB}(x, method = "analytic", alpha = NULL, ...)
}
\arguments{
    \item{x}{an object of class \code{genJS}.}
    \item{method}{\code{[character]} estimation method \code{"analytic"},
        \code{"jackknife"}, or \code{"pretest"} (for pretest estimators).}
    \item{alpha}{\code{[numeric]} used in pretest estimator,
        \eqn{0 < \alpha < 1}{0 < \alpha < 1}.}
    \item{\dots}{additional arguments.}
}
\description{
Mean square error estimation:
    \itemize{
        \item Analytic approximation (Prasad and Rao, 1990; Datta and
            Lahiri, 2002)
        \item Jackknife estimator (Jiang et al., 2002)
        \item Pretest estimator (Datta et al., 2011; Molina et al., 2015)
    }
}
\references{
Datta, G. S., Hall, P. and Mandal, A. (2011).
    Model Selection by Testing for the Presence of
    Small-Area Effects, and Application to Area-Level.
     \emph{Journal of the American Statistical Association} \bold{106},
     pp. 362--374. \doi{10.1198/jasa.2011.tm10036}

Datta, G. S. and Lahiri, P. (2000).
    A unified measure of uncertainty of estimated best linear unbiased
    predictors in small area estimation problems.
    \emph{Statistica Sinica} \bold{10}, pp. 613--627.

Jiang, J., Lahiri, P., Wan, S.-M. (2002). A unified Jackknife Theory
    for Empirical Best Prediction With M-Estimation.
    \emph{The Annals of Statistics} \bold{30}, pp. 1782--1810.
    \doi{10.1214/aos/1043351257}

Molina, I., Rao, J. N. K. and  Datta, G. S. (2015).
    Small area estimation under a Fay-Herriot model with preliminary
    testing for the presence of random area effects.
    \emph{Survey Methodology} \bold{41}, pp. 1--19.

Prasad, N. G. N. and Rao, J. N. K. Rao (1990).
    The Estimation of the Mean Squared Error of Small-Area Estimators.
    \emph{Journal of the American Statistical Association} \bold{85},
    pp. 163--171. \doi{10.2307/2289539}
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
