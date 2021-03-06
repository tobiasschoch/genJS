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
# Parametric bootstrap (Ibanez et al. 2009, BMC Health Service Research)
#-------------------------------------------------------------------------------
H0_bootstrap <- function(object, replicates = 1000, p = c(0.025, 0.975),
    seed = 1)
{
    if (!(inherits(object, c("EB", "scv"))))
        stop("Bootstrap is not available for this object\n")
    call <- object$call
    call$oi <- substitute(oi)
    call$ei <- substitute(ei)
    ei <- object$model$ei
    n <- length(ei)
    set.seed(seed)
    res <- numeric(replicates)
    for (i in 1:replicates) {
        oi <- rpois(n, ei)          # observed counts, assumption rate = 1
		res[i] <- vcov(eval(call))  # compute variance estimate
	}
    quantile(res, probs = p)
}
