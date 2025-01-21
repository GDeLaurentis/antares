#!/usr/bin/env python
# -*- coding: utf-8 -*-

#   ___ _           _     ___          _ _
#  / __(_)_ _  __ _| |___/ __| __ __ _| (_)_ _  __ _ ___
#  \__ \ | ' \/ _` | / -_)__ \/ _/ _` | | | ' \/ _` (_-<
#  |___/_|_||_\__, |_\___|___/\__\__,_|_|_|_||_\__, /__/
#             |___/                            |___/

# Author: Giuseppe


from __future__ import unicode_literals
from __future__ import print_function

from lips import Particles, myException

from antares.core.settings import settings
from antares.core.tools import mapThreads, log_linear_fit, retry, log_linear_fit_Exception


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def single_scalings(oUnknown, invariants, seed=0, verbose=False):
    if oUnknown.is_zero:
        return [0 for inv in invariants]
    else:
        return mapThreads(single_scaling, oUnknown, invariants, seed=seed, verbose=verbose,
                          UseParallelisation=settings.UseParallelisation, Cores=settings.Cores)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


@retry((myException, log_linear_fit_Exception, AssertionError, AttributeError, ), max_tries=2, silent=False)
def single_scaling(oUnknown, invariant, seed=0):

    oParticles = Particles(oUnknown.multiplicity, seed=seed, field=settings.field)

    if settings.field.name == "mpc":

        xaxis, yaxis = [], []

        for k in range(settings.ScalingsIterationsStart, settings.ScalingsIterationsStart + settings.ScalingsIterationsNumber):
            oParticles.variety((invariant, ), (10 ** -k,))
            xaxis += [abs(oParticles.compute(invariant))]
            yaxis += [abs(oUnknown(oParticles))]

        return log_linear_fit(xaxis, yaxis)

    elif settings.field.name == "padic":

        oParticles.variety((invariant, ), (1, ))
        res = oUnknown(oParticles)
        if res.k == 0:
            raise Exception("Lost all padic digits.")
        return res.n
