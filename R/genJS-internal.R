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
# mse analytic approximation (Prasad & Rao, 1990)
.mse_analytic <- function(x)
{
    di <- x$model$di; A <- x$A
    bi <- di / (A + di)
    # first component
    g1 <- A * bi
    # second component (if location is estimated, otherwise zero)
    g2 <- ifelse(x$center, bi^2 / sum(1 / (A + di)), 0)
    # third component (the same for MLE and REML; see Datta & Lahiri, 2000)
    g3 <- 2 * bi^2 / ((A + di) * sum(1 / (A + di)^2))
    # in total
    g1 + g2 + g3
}
#-------------------------------------------------------------------------------
# Jackknife mse estimator of Jiang et al. (2002)
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
