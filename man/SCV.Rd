\name{SCV}
\alias{SCV}
\alias{print.scv}
\alias{summary.scv}
\alias{vcov.scv}
\title{Systematic Component of Variation}
\usage{
SCV(oi, ei, approx = TRUE)

# Utility functions
\method{print}{scv}(x, digits = max(1L, getOption("digits") - 2L), ...)
\method{summary}{scv}(object, ..., digits = max(1L,
    getOption("digits") - 2L))
\method{vcov}{scv}(object, ...)
}
\arguments{
    \item{oi}{\code{[numeric vector]} observed counts.}
    \item{ei}{\code{[numeric vector]} expected counts.}
    \item{approx}{\code{[logical]} if \code{TRUE}, then
        \eqn{log(o_i/e_i)}{log(o[i]/e[i])} is approximate by
        \eqn{(o_i - e_i) / e_i}{(o[i] - e[i]) / e[i]};
        default: \code{FALSE}.}
    \item{x}{an object of class \code{scv}.}
    \item{object}{an object of class \code{scv}.}
    \item{digits}{number of digits.}
    \item{...}{additional arguments.}
}
\description{
Estimator of Systematic Component of Variation (SCV)
}
\references{
    McPherson, K., Wennberg, J. E., Hovind, O. B., Clifford, P. (1982).
    Small-Area Variations in the Use of Common Surgical Procedures: An
    International Comparison of New England, England, and Norway
    \emph{New England Journal of Medicine} \bold{307}, pp. 1310--1314.
    \doi{10.1056/NEJM198211183072104}
}
