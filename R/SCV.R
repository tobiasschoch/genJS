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
# SCV: systematic component of variation, McPerson (1982)
#-------------------------------------------------------------------------------
SCV <- function(oi, ei, approx = TRUE)
{
    stopifnot(all(oi > 0), all(ei > 0))
    yi <- if (approx)
        (oi - ei) / ei
    else
        log(oi) - log(ei)
    A <- mean(yi^2 - 1 / ei)
    structure(list(mu = 0, A = A, center = FALSE,
        model = list(oi = oi, ei = ei, n = length(oi)),
        method = "SCV (McPherson et al., 1982)", converged = TRUE,
        call = match.call()), class = "scv")
}
print.scv <- function(x, digits = max(1L, getOption("digits") - 2L), ...)
{
    cat(paste0(x$method, "\n"))
    cat(paste0("Variance: ", round(x$A, digits), "\n"))
}
# summary method
summary.scv <- function(object, ..., digits = max(1L,
    getOption("digits") - 2L))
{
    print(object)
}
# extract covariance matrix
vcov.scv <- function(object, ...)
{
    object$A
}
