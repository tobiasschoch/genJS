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
# Wrapper function for genJS
#-------------------------------------------------------------------------------
EB <- function(oi, ei, center = FALSE, maxit = 100, tol = 1e-5,
    verbose = FALSE)
{
    yi <- (oi - ei) / ei; di <- oi / ei^2
    tmp <- genJS(yi, di, center, maxit, tol, verbose)
    tmp$model$ei <- ei; tmp$model$oi <- oi
    tmp$call <- match.call()
    class(tmp) <- c("EB", "genJS")
    tmp
}
