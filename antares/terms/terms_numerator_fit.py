#!/usr/bin/env python
# -*- coding: utf-8 -*-

#   _  _                         _           ___ _ _
#  | \| |_  _ _ __  ___ _ _ __ _| |_ ___ _ _| __(_) |_
#  | .` | || | '  \/ -_) '_/ _` |  _/ _ \ '_| _|| |  _|
#  |_|\_|\_,_|_|_|_\___|_| \__,_|\__\___/_| |_| |_|\__|

# Author: Giuseppe


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


class Terms_numerators_fit:

    def fit_numerators(self, llTrialNumInvs=[], accept_all_zeros=True):
        """Obtains numerators by fitting the coefficients of an ansatz. Updates self if succesful. Makes self an empty list otherwise."""
        raise NotImplementedError  # coming soon

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

    def refine_fit(self, tempTerms, llTrialNumInvs):
        raise NotImplementedError  # coming soon

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

    def set_inversion_accuracy(self, silent=False):
        raise NotImplementedError  # coming soon
