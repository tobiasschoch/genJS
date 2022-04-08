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
    at <- complete.cases(yi, di)
    if (sum(at) != length(at))
        stop("Some values are missing or not a number\n")
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

    structure(list(mu = mu, A = A, center = center,
        model = list(yi = yi, di = di, n = length(yi)),
        method = ifelse(center, "Generalized James-Stein REML estimator",
            "General James-Stein MLE estimator"), converged = converged,
        optim = list(niter = niter, tol = tol, negflag = negflag),
        call = match.call()), class = "genJS")
}
# print method
print.genJS <- function(x, digits = max(1L, getOption("digits") - 2L), ...)
{
    if (!x$converged) {
        cat("\nFisher scoring algorithm did not converge\n\n")
    } else {
        cat(paste0(x$method, "\n"))
        if (x$center) {
            cat(paste0("\nLocation: ", round(x$mu, digits),"\n"))
        } else {
            cat("\nLocation kept fix at zero\n")
        }
        cat(paste0("Variance: ", round(x$A, digits), "\n"))
    }
}
# summary method
summary.genJS <- function(object, ..., digits = max(1L,
    getOption("digits") - 2L))
{
    cat(paste0(object$method, "\n"))
    if (!object$converged) {
        cat("\nFisher scoring algorithm did not converge\n\n")
    } else {
        if (object$center) {
            cat(paste0("\nLocation: ", round(object$mu, digits), "\n"))
        } else {
            cat("\nLocation kept fix at zero\n")
        }
        pval <- pvalue_A(object$model$yi, object$model$di, object$center)
        cat(paste0("Variance: ", round(object$A, digits), " (p-value:",
            format.pval(pval), " for H0: A = 0)\n"))
        if (!is.null(object$optim)) {
            cat(paste0("\nNumber of iterations: ", object$optim$niter,
                "\nNumerical tolerance for Fisher scoring: ", object$optim$tol,
                    "\n"))
            if (object$optim$negflag)
                cat("NOTE: Variance estimate was negative (negflag)\n")
        }
    }
    object$pvalue <- pval
    invisible(object)
}
# extract estimate
coef.genJS <- function(object, ...)
{
    object$mu
}
# extract covariance matrix
vcov.genJS <- function(object, ...)
{
    object$A
}
#-------------------------------------------------------------------------------
# EBP: Empirical best predictor
#-------------------------------------------------------------------------------
#   x       An estimated model
#   method  "EB", "morris", "pretest"
predict.genJS <- function(object, method = "EB", alpha = NULL, ...)
{
    match.arg(method, c("EB", "morris", "pretest"))
    n <- object$model$n
    yi <- object$model$yi; di <- object$model$di
    Bi <- di / (object$A + di)
    mu <- object$mu
    switch(method,
        "morris" = {                    # Morris (1983, JASA)
            Bi <- (n - 3) / (n  - 2) * Bi
        },
        "pretest" = {                   # Datta, Hall & Mandal (2011, JASA)
            if (is.null(alpha))
                stop("Argument 'alpha' must be defined\n")
            stopifnot(alpha > 0, alpha < 1)
            pvalue <- pvalue_A(yi, di, object$center)
            if (pvalue > alpha) {
                Bi <- rep(1, n)
                mu <- ifelse(object$center, sum(yi / di) / sum(1 / di), 0)
            }
        })
    # empirical best predictor
    (1 - Bi) * object$model$yi + Bi * mu
}
#-------------------------------------------------------------------------------
# jackknife variance estimator of A
#-------------------------------------------------------------------------------
#FIXME: EB and genJS have different arguments: yi, di and oi, ei
jackknife_var <- function(x)
{
    if (!inherits(x, "genJS"))
        stop("Argument 'x' must be of class 'genJS'\n")
    yi <- x$model$yi
    di <- x$model$di
    n <- x$model$n
    call <- x$call
    call[[2]] <- substitute(yi[-i])
    call[[3]] <- substitute(di[-i])

    theta <- rep(0, n)
    for (i in 1:n) {
        tmp <- eval(call)
        theta[i] <- vcov(tmp)
    }
    # jackknife variance estimate
    var(theta) * (n - 1)^2 / n
}
#-------------------------------------------------------------------------------
# Test under the hypothesis H0: A = 0 and HA: A > 0; see Datta, Hall &
# Mandal (2011, JASA)
#-------------------------------------------------------------------------------
pvalue_A <- function(yi, di, center = TRUE){
    n <- length(yi)
    if (center) {
        # estimate location, mu, (under H0)
        w <- 1 / di
        mu <- sum(w * yi) / sum(w)
        df <- n - 1
    } else {
        mu <- 0
        df <- n
    }
    # test statistic under H0
    chi <- sum((yi - mu)^2 / di)
    # p-value
    stats::pchisq(chi, df = df, lower.tail = FALSE)
}
