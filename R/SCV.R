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
# SCV: systematic component of variation, McPerson (1982)
#-------------------------------------------------------------------------------
#   yi  (Oi - Ei) / Ei, where Oi and Ei are, respectively, the i-th observed
#       and expected counts
SCV <- function(yi)
{
    mean(yi^2)
}
