#!/usr/bin/env python
# -*- coding: utf-8 -*-

#   _  _                         _           ___ _ _
#  | \| |_  _ _ __  ___ _ _ __ _| |_ ___ _ _| __(_) |_
#  | .` | || | '  \/ -_) '_/ _` |  _/ _ \ '_| _|| |  _|
#  |_|\_|\_,_|_|_|_\___|_| \__,_|\__\___/_| |_| |_|\__|

# Author: Giuseppe


from __future__ import unicode_literals
from __future__ import print_function

import os

from sympy import pprint
from copy import deepcopy
from antares.core.bh_patch import accuracy
from antares.core.settings import settings
from lips.invariants import Invariants
from antares.scalings.single import single_scalings
from antares.scalings.pair import pair_scalings
# from antares.ansatze.interface import Ansatz
# from antares.ansatze.fitter import AnsatzeFit


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


class Terms_numerators_fit:

    def fit_numerators(self, llTrialNumInvs=[], accept_all_zeros=True):
        """Obtains numerators by fitting the coefficients of an ansatz. Updates self if succesful. Makes self an empty list otherwise."""#
        raise NotImplementedError  # coming soon

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

    def refine_fit(self, tempTerms, llTrialNumInvs):
        raise NotImplementedError  # coming soon

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

    def set_inversion_accuracy(self, silent=False):
        raise NotImplementedError  # coming soon