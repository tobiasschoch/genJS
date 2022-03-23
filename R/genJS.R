#===============================================================================
# SUBJECT R functions to quantify regional variation
# PROJECT Methodenmandat zum Versorgungsatlas im Auftrg des Obsan
# AUTHORS Tobias Schoch, tobias.schoch@fhnw.ch, January 18, 2022
# LICENSE GPL >= 2
#-------------------------------------------------------------------------------
# Copyright (C) 2022 Tobias Schoch (e-mail: tobias.schoch@fhnw.ch)
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Library General Public License for more details.
#
# You should have received a copy of the GNU Library General Public
# License along with this library; if not, a copy is available at
# https://www.gnu.org/licenses/

#-------------------------------------------------------------------------------
# genJS: Empirical Bayes estimators: (1) restricted maximum likelihood
#     estimator (REML) of location and variance; (2) maximum likelihood
#     estimator of variance (MLE), location set to zero
#-------------------------------------------------------------------------------
#   yi      Direct estimator
#   di      Variance of direct estimator
#   center  Logical; if TRUE center is estimated; otherwise center = 0
#   maxit   Maximum number of iterations
#   tol     Numeric tolerance (termination criterion)
#   verbose Indicating whether warnings should be printed to the console
genJS <- function(yi, di, center = TRUE, maxit = 100, tol = 1e-5,
    verbose = TRUE)
{
    stopifnot(all(di > 0), length(yi) == length(di))
    # initialization
    if (center) {
        # REML estimator
        scoring <- function(w, ri)
        {
            sw <- sum(w); sw2 <- sum(w^2)
            J <- sw2
            S <- -sw + sw2 / sw  + sum((w * ri)^2)
            S / J
        }
        mu <- mean(yi)
    } else {
        # ML estimator (location kept fixed at zero)
        scoring <- function(w, ri)
        {
            J <- sum(w^2)
            S <- -sum(w) + sum((w * ri)^2)
            S / J
        }
        mu <- 0
    }
    A <- 1; w <- 1 / (A + di)
    # Fisher scoring
    negflag <- FALSE; niter <- 0; converged <- FALSE
    while (niter <= maxit) {
        niter <- niter + 1
        Anew <- A + scoring(w, yi - mu)
        # if Anew is negative, we put Anew <- 0 (and set the negflag), the
        # second time we encounter a negative value, we terminate
        if (Anew < 0) {
            Anew <- 0
            if (negflag)
                break
            negflag <- TRUE
        }
        if (abs(Anew - A) < tol) {
            converged <- TRUE
            break
        } else {
            A <- Anew
            w <- 1 / (A + di)
            if (center)
                mu <- sum(w * yi) / sum(w)
        }
    }
    if (!converged && verbose)
        cat("\nAlgorithm DID NOT converge in", niter, "iterations\n\n")

    res <- list(mu = mu, A = A, center = center,
            model = list(yi = yi, di = di, n = length(yi)),
        method = ifelse(center, "REML", "MLE"), converged = converged,
        niter = niter, tol = tol, negflag = negflag)
    class(res) <- "genJS"
    res
}
print.genJS <- function(x, ...)
{
    cat(paste0("\nGeneralized James-Stein ", x$method, " estimator\n"))
    if (!x$converged) {
        cat("\nFisher scoring algorithm did not converge\n\n")
    } else {
        if (x$center) {
            cat("\nLocation\n")
            print(x$mu)
        } else {
            cat("\nLocation kept fix at zero\n")
        }
        cat("Variance:\n")
        print(x$A)
    }
}
summary.genJS <- function(object, ...)
{
    cat(paste0("\nGeneralized James-Stein ", object$method, " estimator\n"))
    if (!object$converged) {
        cat("\nFisher scoring algorithm did not converge\n\n")
    } else {
        cat(paste0("Number of iterations: ", object$niter,
            "\nNumerical tolerance for Fisher scoring: ", object$tol,"\n\n"))
        if (object$negflag)
            cat("Variance estimate was negative (negflag)\n")
        stats::pnorm(0.975)
    }
}
coef.genJS <- function(object, ...)
{
    object$mu
}
vcov.genJS <- function(object, ...)
{
    object$A
}
#-------------------------------------------------------------------------------
# EBP: Empirical best predictor
#-------------------------------------------------------------------------------
#   x       An estimated model computed with EB_est
#   Morris  Logical; bias correction of Morris (1983, JASA, Vol. 78)
predict.genJS <- function(object, morris = FALSE, ...)
{
    Bi <- object$model$di / (object$A + object$model$di)
    if (morris) {
        n <- object$model$n
        Bi <- (n - 3) / (n  - 2) * Bi
    }
    # empirical best predictor
    (1 - Bi) * object$model$yi + Bi * object$mu
}
#-------------------------------------------------------------------------------
# MSE estimation
#-------------------------------------------------------------------------------
#   x       model estimated genJS
#   method  mse estimation method
mse <- function(x, method = "analytic")
{
    if (class(x) != "genJS")
        stop("Argument 'x' must be an estimate computed by genJS()\n")
    switch(match.arg(method, c("analytic", "jackknife")),
        "analytic" = .mse_analytic(x),
        "jackknife" = .mse_jackknife(x))
}
