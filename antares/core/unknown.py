#!/usr/bin/env python
# -*- coding: utf-8 -*-

#   _   _      _
#  | | | |_ _ | |___ _  _____ __ ___ _
#  | |_| | ' \| / / ' \/ _ \ V  V / ' \
#   \___/|_||_|_\_\_||_\___/\_/\_/|_||_|
#

# Author: Giuseppe

import os
import numpy
import mpmath

from sympy import pprint
from copy import deepcopy

from lips import Particles

from antares.core.bh_unknown import BHUnknown
from antares.core.settings import settings
from antares.core.tools import Generate_LaTeX_and_PDF, forbidden_ordering, flatten
from lips.invariants import Invariants
from antares.core.numerical_methods import Numerical_Methods

local_directory = os.path.dirname(os.path.abspath(__file__))
mpmath.mp.dps = 300

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def caching_decorator(func):

    def decorated(self):
        if (hasattr(self, "_{}_cached".format(func.__name__)) and hasattr(self, "_{}_hash_when_cached".format(func.__name__)) and
           getattr(self, "_{}_hash_when_cached".format(func.__name__)) == self.__hash__()):
            return getattr(self, "_{}_cached".format(func.__name__))
        else:
            result = func(self)
            setattr(self, "_{}_cached".format(func.__name__), result)
            setattr(self, "_{}_hash_when_cached".format(func.__name__), self.__hash__())
            return result

    return decorated


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


class Unknown(Numerical_Methods, object):

    def __init__(self, original_unknown, load_partial_results=False, silent=True):
        if callable(original_unknown):
            self.original_unknown = original_unknown
            if not hasattr(self, "__name__") and hasattr(original_unknown, "__name__"):
                self.__name__ = original_unknown.__name__
            self.partial_pieces = []
        else:
            raise Exception("Unknown initialisation argument must be callable, got {} of type {}.".format(original_unknown, type(original_unknown)))
        # Load partial results if asked to do so
        if load_partial_results is True:
            from antares.terms.terms import LoadResults
            PartialResults, partial = LoadResults(self.res_path, load_partial_results_only=True, silent=silent)
            if PartialResults is not None:
                self.add_partial_piece(PartialResults[0])
                self.original_unknown = deepcopy(self)
                self.partial_pieces = []

    def __hash__(self):
        return hash(self.recursively_extract_terms())

    def add_partial_piece(self, partial_piece):
        partial_piece.multiplicity = self.multiplicity
        self.partial_pieces += [partial_piece]

    def remove_last_partial_piece(self):
        self.partial_pieces = self.partial_pieces[:-1]

    def reset(self):
        self.partial_pieces = []

    @property
    def helconf(self):
        return self.recursively_extract_original_unknown().helconf

    @property
    def helconf_and_loopid(self):
        return self.recursively_extract_original_unknown().helconf_and_loopid

    @property
    def amppart(self):
        return self.recursively_extract_original_unknown().amppart

    @property
    def ampindex(self):
        return self.recursively_extract_original_unknown().ampindex

    @property
    def loopid(self):
        return self.recursively_extract_original_unknown().loopid

    @property
    def amppart_and_ampindex(self):
        return self.recursively_extract_original_unknown().amppart_and_ampindex

    @property
    def multiplicity(self):
        return self.original_unknown.multiplicity

    @property
    def internal_masses(self):
        return self.original_unknown.internal_masses

    @property
    def basis_functions(self):
        if hasattr(self.recursively_extract_original_unknown(), "basis_functions"):
            return self.recursively_extract_original_unknown().basis_functions
        else:
            return []

    @property
    def basis_functions_invs(self):
        return list(set(flatten(self.basis_functions)))

    @property
    def spurious_poles(self):
        if hasattr(self.recursively_extract_original_unknown(), "spurious_poles"):
            return self.recursively_extract_original_unknown().spurious_poles
        else:
            return []

    @property
    @caching_decorator
    def poles_to_be_eliminated(self):
        if settings.AutomaticPartialFractioning is False or settings.ProceedByGuesses is False:
            return [""]
        poles_to_be_eliminated = self._ordering_for_inversion_and_printing()
        return poles_to_be_eliminated

    def _ordering_for_inversion_and_printing(self, lInvariants=None, dExponents=None):
        if lInvariants is None:
            lInvariants = self.den_invs
        if dExponents is None:
            dExponents = self.den_exps
        oInvariants = Invariants(self.multiplicity, Restrict3Brackets=settings.Restrict3Brackets,
                                 Restrict4Brackets=settings.Restrict4Brackets, FurtherRestrict4Brackets=settings.FurtherRestrict4Brackets)
        spurious_poles = self.spurious_poles
        basis_functions_invs = self.basis_functions_invs
        return sorted(lInvariants, key=lambda inv:
                      (
                          inv not in [_inv for _inv in oInvariants.invs_3 if _inv[1] == _inv[-2]],  # e.g. ⟨1|(2+3)|1] before ⟨1|(2+3)|6], ⟨3|(1+2)|3] before ⟨3|(1+2)|4]
                          inv not in [_inv for _inv in oInvariants.invs_2 if _inv not in basis_functions_invs] and basis_functions_invs != [],
                          inv not in oInvariants.invs_s,
                          # poles forbidden with Delta if present <-- not necessary given the other criterea
                          # inv not in forbidden_invariants([_inv for _inv in lInvariants if _inv in oInvariants.invs_D][0], self) if any([
                          #     _inv in oInvariants.invs_D for _inv in lInvariants]) else True,
                          inv not in [_inv for _inv in oInvariants.invs_3 if _inv not in spurious_poles],
                          -dExponents[inv],
                          inv not in basis_functions_invs,
                          inv not in self.spurious_poles,
                          oInvariants.full.index(inv) if inv in oInvariants.full else -1
                      ))

    @property
    def easy_boxes(self):
        if hasattr(self.recursively_extract_original_unknown(), "easy_boxes"):
            return self.recursively_extract_original_unknown().easy_boxes

    @property
    @caching_decorator
    def mass_dimension(self):
        if self.what_am_I == "Numerical" and hasattr(self, "helconf") and hasattr(self, "amppart") and "external" not in self.amppart:
            if self.amppart == "box":
                return 8 - len(self.helconf)
            elif self.amppart == "triangle":
                return 6 - len(self.helconf)
            elif self.amppart == "bubble" or self.amppart == "tree" or self.amppart == "rational":
                return 4 - len(self.helconf)
        return super(Unknown, self).mass_dimension

    @property
    @caching_decorator
    def phase_weights(self):
        if self.what_am_I == "Numerical" and hasattr(self, "helconf") and hasattr(self, "amppart") and "external" not in self.amppart:
            phase_weights = []
            for entry in self.helconf:
                if entry == "p" or entry == "yp":
                    phase_weights += [-2]
                elif entry == "m" or entry == "ym":
                    phase_weights += [+2]
                elif entry == "qp" or entry == "qbp":
                    phase_weights += [-1]
                elif entry == "qm" or entry == "qbm":
                    phase_weights += [+1]
            return phase_weights
        return super(Unknown, self).phase_weights

    @property
    @caching_decorator
    def soft_weights(self):
        oParticles = Particles(self.multiplicity, seed=0)
        soft_weights = [[], []]
        for i, oParticle in enumerate(oParticles):
            a, b = self.multiplicity - 1, self.multiplicity
            if i in [a, b]:
                a, b = 1, 2
            # Angles
            oParticles.randomise_all()
            oParticle.r_sp_u = numpy.array([mpmath.mpc('1e-50', '2e-50'), mpmath.mpc('3e-50', '2e-50')])
            oParticles.fix_mom_cons(a, b)
            before = self(oParticles)
            oParticle.r_sp_u = oParticle.r_sp_u / 10
            oParticles.fix_mom_cons(a, b)
            after = self(oParticles)
            temp_p_w = mpmath.log((abs(after / before))) / mpmath.log(10)
            soft_weights[0] += [int(round(float(temp_p_w)))]
            # Squares
            oParticles.randomise_all()
            oParticle.l_sp_u = numpy.array([mpmath.mpc('1e-50', '2e-50'), mpmath.mpc('3e-50', '2e-50')])
            oParticles.fix_mom_cons(a, b)
            before = self(oParticles)
            oParticle.l_sp_u = oParticle.l_sp_u / 10
            oParticles.fix_mom_cons(a, b)
            after = self(oParticles)
            temp_p_w = mpmath.log((abs(after / before))) / mpmath.log(10)
            soft_weights[1] += [int(round(float(temp_p_w)))]
        return soft_weights

    @property
    def what_am_I(self):
        from antares.terms.terms import Terms
        if isinstance(self.original_unknown, Unknown):
            return self.original_unknown.what_am_I
        elif isinstance(self.original_unknown, BHUnknown):
            return "Numerical"
        elif isinstance(self.original_unknown, Terms):
            return "Analytical"

    def __call__(self, oParticles):
        result = self.original_unknown(oParticles)
        for partial_piece in self.partial_pieces:
            result -= partial_piece(oParticles)
        return result

    def print_partial_result(self, partial=True, compile_tex_to_pdf=False):
        oTerms = self.recursively_extract_terms()
        oTerms.rearrange_and_finalise()
        Generate_LaTeX_and_PDF(oTerms.Write_LaTex(), self.res_path, partial=partial, compile_tex_to_pdf=compile_tex_to_pdf)

    def recursively_extract_original_unknown(self):
        if hasattr(self.original_unknown, "original_unknown") and hasattr(self.original_unknown, "recursively_extract_original_unknown"):
            return self.original_unknown.recursively_extract_original_unknown()
        else:
            return self.original_unknown

    def save_original_unknown_call_cache(self):
        if hasattr(self.recursively_extract_original_unknown(), "save_call_cache"):
            self.recursively_extract_original_unknown().save_call_cache()

    def recursively_extract_terms(self):
        from antares.terms.terms import Terms
        if isinstance(self.original_unknown, Unknown):
            oTerms = self.original_unknown.recursively_extract_terms()
        else:
            oTerms = Terms([])
        for oPartialTerms in self.partial_pieces:
            oTerms += oPartialTerms
        if oTerms != []:
            oTerms.oUnknown = self
        return oTerms

    def reorder_invariants_for_partial_fractioning(self):
        # Optimal Ordering for Partial Fractioning
        self.den_invs.reverse()
        forbidden_ordering(self.den_invs, self.den_invs, self)
        unsorted_den_invs = [_entry for _entry in self.den_invs]
        self.den_invs = []
        while len(unsorted_den_invs) > 0:
            self.den_invs += [unsorted_den_invs[0]]
            unsorted_den_invs.remove(unsorted_den_invs[0])
            forbidden_ordering(unsorted_den_invs, self.den_invs, self)
        print("The reorganized denominator invariants are:")
        pprint(self.den_invs)
        print("")

    def does_not_require_partial_fractioning(self):
        if ((self.oTermsFromSingleScalings.ansatze_mass_dimensions[0], self.oTermsFromSingleScalings.ansatze_phase_weights[0]) == (0, [0 for i in range(self.multiplicity)]) or
           self.oTermFromSingleScalings.oDen.lInvs == [] or settings.ForceUseOfJustOneDenominator is True):
            if settings.ForceUseOfJustOneDenominator is True:
                print("Attempting to use a single denominator as required.")
            else:
                print("Both mass dimension and phase weights are zero. Attempting to finish here.")
            return True
        else:
            return False

    def fit_single_scalings_result(self, split_pole_orders=True):
        from antares.terms.terms import Terms
        self._pair_invs, self._pair_exps, self._pair_friends = None, None, None
        oTerms = Terms([self.oTermFromSingleScalings])
        oTerms.oUnknown = self
        if any([self.den_exps[inv] > 1 for inv in self.den_invs]) and split_pole_orders is True:
            higher_poles, higher_poles_exps = map(list, zip(*[(inv, self.den_exps[inv]) for inv in self.den_invs if self.den_exps[inv] > 1]))
            oTerms.oFittingSettings.lSmallInvs = [higher_poles[0]]
            oTerms.oFittingSettings.lSmallInvsExps = [higher_poles_exps[0]]
        else:
            oTerms.oFittingSettings.lSmallInvs = []
            oTerms.oFittingSettings.lSmallInvsExps = []
        oTerms.fit_numerators()
        if len(oTerms) >= 1:
            return [oTerms]
        else:
            raise Exception("Failed to fit single scalings result.")

    def get_partial_fractioned_terms(self, invariant, max_nbr_terms=1):
        from antares.terms.terms import fix_old_automatic_partial_fractioning_format
        from antares.partial_fractioning.automatic import AutomaticPartialFraction
        # from partial_fractioning.by_guesses import make_guesses
        from antares.partial_fractioning.v3 import partial_partial_fractioning
        lTerms = []
        if settings.AutomaticPartialFractioning is True:
            if settings.ProceedByGuesses is False:
                self.reorder_invariants_for_partial_fractioning()
                lTerms = list(AutomaticPartialFraction(self))
                lTerms.sort(key=lambda s: (s.ansatze_mass_dimensions[0], sum(s.ansatze_mass_dimensions)))
                fix_old_automatic_partial_fractioning_format(lTerms)
            else:
                # lTerms = make_guesses(self, invariant)
                lTerms = partial_partial_fractioning(self, invariant, max_nbr_terms=max_nbr_terms)
        elif settings.AutomaticPartialFractioning is False:
            lTerms = settings.ManualPartialFraction(self)
        return lTerms

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

    @property
    def res_path(self):
        if not hasattr(self, "_res_path"):
            res_path = settings.base_res_path + self.original_unknown.__name__
            if "/" in self.original_unknown.__name__:
                upper_path = settings.base_res_path + "/".join(self.original_unknown.__name__.split("/")[:-1])
                if not os.path.exists(upper_path):
                    os.makedirs(upper_path)
                if not os.path.exists(upper_path + "/cpp"):
                    os.makedirs(upper_path + "/cpp")
                if not os.path.exists(upper_path + "/math"):
                    os.makedirs(upper_path + "/math")
            return res_path
        else:
            return self._res_path

    @res_path.setter
    def res_path(self, temp_res_path):
        self._res_path = temp_res_path


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
