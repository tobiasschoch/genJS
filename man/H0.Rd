\name{H0}
\alias{H0_bootstrap}
\title{Parametric Bootstrap of the Null Hypothesis}
\usage{
H0_bootstrap(object, replicates = 1000, p = c(0.025, 0.975), seed = 1)
}
\arguments{
    \item{object}{an SCV estimate.}
    \item{replicates}{number of bootstrap replicates.}
    \item{p}{percentiles.}
    \item{seed}{random seed}
}
\description{
Parametric bootstrap under the null hypothesis of no (excess) regional
variation (Ibanez et al., 2009).

The method is available for
    \itemize{
        \item \code{\link[=EB]{EB()}}
        \item \code{\link[=SCV]{SCV()}}
    }
}
\author{
Tobias Schoch
}
\references{
Ibanez, B., Librero, J., Bernal-Delgado, E., Peiro, S.,
    Gonzalez Lopez-Valcarcel, B., Martainez, N. and Aizpuru, F. (2009).
    Is there much variation in variation? Revisiting statistics of
    small area variation in health services research
    \emph{BMC Health Service Research} \bold{9}, pp. 1--12.
    \doi{10.1186/1472-6963-9-60}
}
\examples{
data(TURP)

# parameter estimation
m <- EB(TURP$oi, TURP$ei)

# empirical Bayes predictor
H0_bootstrap(m)
}
