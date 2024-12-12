#!/usr/bin/env python
# -*- coding: utf-8 -*-

#                            _   _   _
#   __ ___ _ _ ___   ___ ___| |_| |_(_)_ _  __ _ ___
#  / _/ _ \ '_/ -_)_(_-</ -_)  _|  _| | ' \/ _` (_-<
#  \__\___/_| \___(_)__/\___|\__|\__|_|_||_\__, /__/
#                                          |___/

# Author: Giuseppe

from __future__ import print_function
from __future__ import unicode_literals

import sys
import os
import logging

from syngular import Field

from antares.core.bh_patch import BH_found

if BH_found:
    import BH

MainPythonDirectory = os.path.dirname(os.path.abspath(__file__))[:-5]


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


class Settings(object):

    def __init__(self):

        # Miscellaneous
        self.BHsettings = "USE_KNOWN_FORMULAE no \n SET_ALL_RAT_TO_ZERO no"
        self.UseParallelisation = True
        self.Cores = 6
        self.gmp_precision = 1024
        self.to_int_prec = "1e-6"
        self.LoggingLevel = "logging.WARNING"
        self.UsePartialResultsIfPossible = True

        # Collinear Limits
        self.DoScalings = True
        self.SingleScalingsUseCache = False
        self.ScalingsIterationsStart = 28
        self.ScalingsIterationsNumber = 2
        self.ScalingsMaxNumberOfFailes = 5
        self.SingleScalingsUse4Brackets = True
        self.Restrict3Brackets = True  # Restrict 3Brackets to neighbouring ones
        self.Restrict4Brackets = True  # Restrict 4Brackets to neighbouring ones
        self.FurtherRestrict4Brackets = True  # [|( ) & ( )|] independerly neighbouring

        # Ansatze Generation
        self.NumeratorAnsatzProvider = ["CP_SAT", "CP_SAT_HIGGS", "CP_SAT_W_BOSON", ][0]  # "SPINOR_SOLVE" - Deprecated
        self.AutomaticPartialFractioning = True
        self.ProceedByGuesses = True  # Specifies the type of automatic partial fractioning used
        self.UseAdHocNumeratorAnsatze = False
        self.ExploreDoubleScalings = False  # Requires UseAdHocNumeratorAnsatze
        self.SplitHigherOrderPoles = True
        self.AddNumeratorsFromPreviousRelevantTerms = True

        # Ansatze Fit
        self.PerformNumeratorFit = True
        self.RefineFit = True
        self.ForceUseOfJustOneDenominator = False
        self.MaximumMatrixSizeToInvert = 6000
        self.ObtainOnlyOneSolution = True
        self.InversionUsesGMPs = False
        self.UseGpu = False
        self.CacheMatrixInfo = False
        self.LoadMatrixInfo = False

        # Number field
        self.field = Field('mpc', 0, 300)   # Field('padic', prime, 6)

    def read_from_file(self, file_full_path):
        # Needs reworking
        return None
        import importlib
        if os.path.isfile(file_full_path) is False:
            raise Exception("Settings file with path {} not found.".format(file_full_path))
        base_settings = importlib.machinery.SourceFileLoader('', file_full_path).load_module()
        base_settings_dictionary = base_settings.__dict__
        [base_settings_dictionary.pop(key) for key in list(base_settings_dictionary.keys()) if key in ['os', 'sys', 'unicode_literals']]
        self.run_file_name = file_full_path.split("/")[-1]
        for key, value in base_settings_dictionary.items():
            self.__setattr__(key, value)

    @property
    def run_file_name(self):
        return self._run_file_name

    @run_file_name.setter
    def run_file_name(self, value):
        self._run_file_name = value
        self.logfile = self.run_file_name + ".log"

    @property
    def base_cache_path(self):
        return MainPythonDirectory + "/../.cache/"

    @property
    def base_res_path(self):
        return MainPythonDirectory + "/../Results/"

    @property
    def PWTestingCachePath(self):
        pwtestingcache = MainPythonDirectory + "/testing/cache/"
        if not os.path.exists(pwtestingcache):
            os.makedirs(pwtestingcache)
        return pwtestingcache

    @property
    def logfile(self):
        return self._logfile

    @logfile.setter
    def logfile(self, value):
        self._logfile = value
        self.logpath = MainPythonDirectory + "/../Logging"
        if not os.path.exists(self.logpath):
            os.makedirs(self.logpath)
        self.logpath = MainPythonDirectory + "/../Logging/" + self.logfile
        logging.basicConfig(filename=self.logpath, filemode='w', level=logging.WARNING)
        if sys.stdout.isatty() is False and hasattr(sys.stdout, "__str__") and "ipykernel.iostream.OutStream" not in sys.stdout.__str__():
            f = open(self.logpath, mode="w+")
            sys.stdout = f
            sys.stderr = f
            print.refs.append(sys.stdout.tell())

    @property
    def Cores(self):
        return self._Cores

    @Cores.setter
    def Cores(self, value):
        self._Cores = value
        os.environ["NUM_OF_CORES"] = str(self.Cores)

    @property
    def BHsettings(self):
        return self._BHsettings

    @BHsettings.setter
    def BHsettings(self, value):
        self._BHsettings = value
        if BH_found:
            BH.use_setting(str(self.BHsettings))

    @property
    def gmp_precision(self):
        return self._gmp_precision

    @gmp_precision.setter
    def gmp_precision(self, value):
        self._gmp_precision = value
        if BH_found:
            BH.RGMP.set_precision(self.gmp_precision)

    @property
    def to_int_prec(self):
        return self._to_int_prec

    @to_int_prec.setter
    def to_int_prec(self, value):
        self._to_int_prec = value
        if BH_found:
            BH.to_int_prec = BH.RGMP(str(self.to_int_prec))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


settings = Settings()
