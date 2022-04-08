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
# MSE estimation
mse <- function(x, ...)
{
    UseMethod("mse")
}
#-------------------------------------------------------------------------------
# mse method for EB
mse.EB <- function(x, method = "analytic", alpha = NULL, ...)
{
    # cast 'EB' object to a 'genJS' object
    x$call[[1]] <- substitute(genJS)
    names(x$call)[match(c("oi", "ei"), names(x$call))] <- c("yi", "di")
    if (is.null(x$call$center))
        x$call$center <- FALSE
    class(x) <- "genJS"
    # call mse
    mse(x, method, alpha, ...)
}
#-------------------------------------------------------------------------------
# mse method for genJS
mse.genJS <- function(x, method = "analytic", alpha = NULL, ...)
{
    switch(method,
        "analytic" = .mse_analytic(x),
        "jackknife" = .mse_jackknife(x),
        "pretest" = {
            if (is.null(alpha))
                stop("Argument 'alpha' must be defined\n")
            stopifnot(alpha > 0, alpha < 1)
            .mse_pretest(x, alpha)
        },
        stop(paste0("Argument method = '", method, "' is unknown\n")))
}
#-------------------------------------------------------------------------------
# mse analytic approximation (Prasad & Rao, 1990, JASA)
.mse_analytic <- function(x)
{
    di <- x$model$di; A <- x$A
    Bi <- di / (A + di)
    # first component
    g1 <- A * Bi
    # second component (if location is estimated, otherwise zero)
    g2 <- ifelse(x$center, Bi^2 / sum(1 / (A + di)), 0)
    # third component (the same for MLE and REML; see Datta & Lahiri, 2000,
    # Statistica Sinica)
    g3 <- 2 * Bi^2 / ((A + di) * sum(1 / (A + di)^2))
    # in total
    g1 + g2 + g3
}
#-------------------------------------------------------------------------------
# mse analytic approximation for pretest estimator; Molina, Rao & Datta (2015,
# Survey Methodology)
.mse_pretest <- function(x, alpha)
{
    pval <- pvalue_A(x$model$yi, x$model$di, x$center)
    if (pval > alpha) {
        di <- x$model$di; A <- x$A
        Bi <- di / (A + di)
        # second component (if location is estimated, otherwise zero)
        if (x$center)
            Bi^2 / sum(1 / (A + di))
        else
            rep(0, length(di))
    } else {
        .mse_analytic(x)
    }
}
#-------------------------------------------------------------------------------
# Jackknife mse estimator of Jiang, Lahiri & Wang  (2002, Annals)
.mse_jackknife <- function(x)
{
    yi <- x$model$yi
    di <- x$model$di; n <- x$model$n
    delta_base <- predict(x)                # EB predictor
    g1_base <- x$A * di / (x$A + di)        # g1 component of mse
    call <- x$call
    call[[2]] <- substitute(yi[-i])
    call[[3]] <- substitute(di[-i])

    m1 <- rep(0, n); m2 <- rep(0, n)
    for (i in 1:n) {
        tmp <- eval(call)
        tmp$model <- list(yi = yi, di = di)
        # jackknife estimator of rule delta
        m1 <- m1 + (tmp$A * di / (tmp$A + di) - g1_base)^2
        m2 <- m2 + (predict(tmp) - delta_base)^2
    }
    m1 <- m1 * (n - 1) / n
    m2 <- m2 * (n - 1) / n
    # mse estimator
    g1_base - m1 + m2
}
