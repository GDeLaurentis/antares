#!/usr/bin/env python
# -*- coding: utf-8 -*-

#    ___             _ _
#   | _ \___ ____  _| | |_
#   |   / -_|_-< || | |  _|
#   |_|_\___/__/\_,_|_|\__|


# Author: Giuseppe

from __future__ import unicode_literals
from __future__ import print_function

import sys
import math
import re
import copy
import numpy
import sympy
import functools

from lips import Particles
from lips.invariants import Invariants
from lips.particles_eval import non_unicode_powers, pA2, pS2, pSijk, pDijk_non_adjacent
from lips.tools import pNB as pNB_internal
from lips.particles_eval import pNB as pNB_overall
from lips.symmetries import inverse

from antares.core.tools import Write, LaTeXToPython, flatten, get_common_Q_factor, get_max_abs_numerator, get_max_denominator
from antares.core.settings import settings
from antares.core.numerical_methods import Numerical_Methods
from antares.core.bh_patch import accuracy
from antares.core.unknown import Unknown
from antares.scalings.pair import pair_scalings
from antares.terms.terms_numerator_fit import Terms_numerators_fit
from antares.terms.term import Term, Numerator, Denominator


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


class Terms(Numerical_Methods, Terms_numerators_fit, list):

    def __init__(self, lllNumInvs_or_Terms, lllNumExps=None, llNumCoefs_or_Sym=None, llDenInvs=None, llDenExps=None):
        list.__init__(self)
        if isinstance(lllNumInvs_or_Terms, str):
            self.__rstr__(lllNumInvs_or_Terms)
        elif all([isinstance(oTerm, Term) for oTerm in lllNumInvs_or_Terms]):
            for oTerm in lllNumInvs_or_Terms:
                self.append(oTerm)
            for oTerm in lllNumInvs_or_Terms:
                if hasattr(oTerm, "multiplicity"):
                    self.multiplicity = oTerm.multiplicity
        else:
            for i in range(len(lllNumInvs_or_Terms)):
                self.append(Term(Numerator(lllNumInvs_or_Terms[i], lllNumExps[i], llNumCoefs_or_Sym[i]), Denominator(llDenInvs[i], llDenExps[i])) if (
                    lllNumInvs_or_Terms[i] != [[]] or llDenInvs[i] != [] or (lllNumInvs_or_Terms[i] == [[]] and lllNumExps[i] == [[]] and llNumCoefs_or_Sym[i] == [] and
                                                                             llDenInvs[i] == [] and llDenExps[i] == [])) or
                            len(llNumCoefs_or_Sym[i][0]) == 2  # i.e. just a number
                            else Term(llNumCoefs_or_Sym[i][0]))
        self.oFittingSettings = FittingSettings()
        self.useful_symmetries = []
        self.symmetries = []
        self.update_invs_set()

    def __and__(self, other):
        assert isinstance(other, Terms)
        appendedTerms = copy.copy(self)
        for entry in other:
            appendedTerms.append(entry)
        return appendedTerms

    def __add__(self, other):
        if other is None:
            return self
        if any([oTerm.am_I_a_symmetry for oTerm in self]):
            if any([oTerm.am_I_a_symmetry for oTerm in other]):
                last_symmetry_1 = self.index_of_last_symmetry()
                last_symmetry_2 = other.index_of_last_symmetry()
                oSumTerms = Terms(list(self[0:last_symmetry_1 + 1]))
                for oTerm in other[0:last_symmetry_2 + 1]:
                    oSumTerms.append(oTerm)
                oSumTerms.compactify_symmetries()
                for oTerm in self[last_symmetry_1 + 1:]:
                    oSumTerms.append(oTerm)
                for oTerm in other[last_symmetry_2 + 1:]:
                    oSumTerms.append(oTerm)
            else:
                oSumTerms = Terms(list(self))
                for oTerm in other:
                    oSumTerms.append(oTerm)
        else:
            if any([oTerm.am_I_a_symmetry for oTerm in other]):
                last_symmetry_2 = other.index_of_last_symmetry()
                oSumTerms = Terms(list(other[0:last_symmetry_2 + 1]))
                for oTerm in self:
                    oSumTerms.append(oTerm)
                for oTerm in other[last_symmetry_2 + 1:]:
                    oSumTerms.append(oTerm)
            else:
                oSumTerms = Terms(list(self))
                for oTerm in other:
                    oSumTerms.append(oTerm)
        oSumTerms.update_invs_set()
        return oSumTerms

    def __sub__(self, other):
        return self + (-1 * other)

    def __radd__(self, other):
        return self.__add__(other)

    def __iadd__(self, other):
        oSumTerms = self.__add__(other)
        del self[:]
        for oTerm in oSumTerms:
            self.append(oTerm)
        self.update_invs_set()
        return self

    def __mul__(self, other):
        return Terms([oTerm * other for oTerm in self])

    def __rmul__(self, other):
        return self * other

    def __neg__(self):
        return -1 * self

    def __truediv__(self, other):
        return self.__div__(other)

    def __div__(self, other):
        return Terms([oTerm / other for oTerm in self])

        # print("\n\n\n\n WTF!! \n\n\n\n\n")
        # if type(other) is int:
        #     new = copy.deepcopy(self)
        #     for lCoefs in new[-1:].llCoefs:
        #         for i, coef in enumerate(lCoefs):
        #             lCoefs[i] = tuple([coef[0] / other, coef[1] / other])
        #     return new
        # else:
        #     raise Exception("Terms division implemented only for integers.")

    def __getitem__(self, item):
        if isinstance(list.__getitem__(self, item), Term):
            oTerm = list.__getitem__(self, item)
            if hasattr(self, "multiplicity"):
                oTerm.multiplicity = self.multiplicity
            return oTerm
        else:
            oTerms = Terms(list.__getitem__(self, item))
            if hasattr(self, "oFittingSettings"):
                oTerms.oFittingSettings = self.oFittingSettings
            if hasattr(self, "oUnknown"):
                oTerms.oUnknown = self.oUnknown
            if hasattr(self, "multiplicity"):
                oTerms.multiplicity = self.multiplicity
            return oTerms

    def __getslice__(self, i, j):
        NewTerms = Terms(list.__getslice__(self, i, j))
        if hasattr(self, "oFittingSettings"):
            NewTerms.oFittingSettings = self.oFittingSettings
        if hasattr(self, "oUnknown"):
            NewTerms.oUnknown = self.oUnknown
        if hasattr(self, "multiplicity"):
            NewTerms.multiplicity = self.multiplicity
        return NewTerms

    def __contains__(self, other):  # is other in self?
        if isinstance(other, Term):
            if any([other in oTerm for oTerm in self]):
                return True
            else:
                return False
        if isinstance(other, Terms):
            if not all([other_term in self for other_term in other]):
                return False
            else:
                return True

    def __str__(self):
        string = "\n".join([str(oTerm) if hasattr(oTerm, 'oNum') else "+" + str(oTerm) for oTerm in self])
        if sys.version_info.major > 2:
            return string
        else:
            return string.encode('utf-8')

    def __rstr__(self, string):
        self.__init__([Term(entry) for entry in string.replace(" ", "").split("\n") if entry != ''])

    def __repr__(self):
        return f"Terms(\"\"\"{self}\"\"\")"

    def __hash__(self):
        return hash(tuple(self))

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

    def index_of_last_symmetry(self):
        for index, oTerm in list(enumerate(self))[::-1]:
            if oTerm.am_I_a_symmetry is True:
                return index
        else:
            return None

    def compactify_symmetries(self):
        oTerms = copy.copy(self)
        something_was_rearranged = True
        while something_was_rearranged is True:
            something_was_rearranged = False
            indices_of_symmetries = [i for i, iTerm in enumerate(oTerms) if hasattr(iTerm, "tSym")]
            if indices_of_symmetries != []:
                grouped_indices_of_symmetries = [[indices_of_symmetries[0]]]
                for ind in indices_of_symmetries[1:]:
                    if ind - 1 in grouped_indices_of_symmetries[-1]:
                        grouped_indices_of_symmetries[-1] += [ind]
                    else:
                        grouped_indices_of_symmetries += [[ind]]
                sets_grouped_symmetries = list(map(set, [[oTerms[index_of_symmetry].tSym for index_of_symmetry in group_of_indices]
                                                         for group_of_indices in grouped_indices_of_symmetries]))
                for symmetry_group in sets_grouped_symmetries:
                    if sets_grouped_symmetries.count(symmetry_group) > 1:
                        something_was_rearranged = True
                        first_symmetry_to_move_to = grouped_indices_of_symmetries[sets_grouped_symmetries.index(symmetry_group)][0]
                        first_symmetry_to_remove = grouped_indices_of_symmetries[::-1][sets_grouped_symmetries[::-1].index(symmetry_group)][0]
                        last_symmetry_to_remove = grouped_indices_of_symmetries[::-1][sets_grouped_symmetries[::-1].index(symmetry_group)][-1]
                        first_term_to_move = grouped_indices_of_symmetries[::-1][sets_grouped_symmetries[::-1].index(symmetry_group) + 1][-1] + 1
                        oTerms = Terms(list(oTerms[0:first_symmetry_to_move_to]) + list(oTerms[first_term_to_move:first_symmetry_to_remove]) + list(oTerms[
                            first_symmetry_to_move_to:first_term_to_move]) + list(oTerms[last_symmetry_to_remove + 1:]))
                        break
        del self[:]
        for oTerm in oTerms:
            self.append(oTerm)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

    @property
    def multiplicity(self):
        if hasattr(self, "_multiplicity"):
            return self._multiplicity
        elif hasattr(self, "oUnknown"):
            return self.oUnknown.multiplicity
        else:
            raise AttributeError("Terms object has no attribute multiplicity")

    @multiplicity.setter
    def multiplicity(self, value):
        self._multiplicity = value
        for oTerm in self:
            oTerm.multiplicity = self.multiplicity

    @property
    def internal_masses(self):
        if hasattr(self, "_internal_masses"):
            return self._internal_masses
        elif hasattr(self, "oUnknown"):
            return self.oUnknown.internal_masses
        else:
            return set()

    @internal_masses.setter
    def internal_masses(self, value):
        self._internal_masses = value
        for oTerm in self:
            oTerm.internal_masses = self.internal_masses

    @property
    def are_fully_reduced(self):
        are_fully_reduced = True
        for oTerm in self:
            if oTerm.is_fully_reduced is False:
                are_fully_reduced = False
                break
        return are_fully_reduced

    @property
    @caching_decorator
    def mass_dimensions(self):
        lM = []
        for oTerm in self:
            lM += [oTerm.mass_dimension]
        return lM

    @property
    @caching_decorator
    def lphase_weights(self):
        lPW = []
        for oTerm in self:
            lPW += [oTerm.phase_weights]
        return lPW

    @property
    @caching_decorator
    def ansatze_mass_dimensions(self):
        # if self.oUnknown.den_invs == [] and self.oUnknown.num_invs == []:
        #     return [self.oUnknown.mass_dimension]
        lM = []
        oUnknown_mass_dimension = self.oUnknown.mass_dimension
        for i in range(len(self)):
            if self[i].am_I_a_symmetry is True:
                continue
            mdPartial = self[i].mass_dimension
            lM += [oUnknown_mass_dimension - mdPartial]
        return lM

    @property
    @caching_decorator
    def ansatze_phase_weights(self):
        # if self.oUnknown.den_invs == [] and self.oUnknown.num_invs == []:
        #    return [self.oUnknown.phase_weights]
        lPW = []
        oUnknown_phase_weights = self.oUnknown.phase_weights
        for i in range(len(self)):
            if self[i].am_I_a_symmetry is True:
                continue
            pwPartial = self[i].phase_weights
            lPW += [[_i - _j for _i, _j in zip(oUnknown_phase_weights, pwPartial)]]
        return lPW

    @property
    def ansatze_angle_degrees(self):
        lmd, lpws = self.ansatze_mass_dimensions, self.ansatze_phase_weights
        return [int(md + sum(pws) / 2) // 2 for md, pws in zip(lmd, lpws)]

    @property
    def ansatze_square_degrees(self):
        lmd, lpws = self.ansatze_mass_dimensions, self.ansatze_phase_weights
        return [int(md - sum(pws) / 2) // 2 for md, pws in zip(lmd, lpws)]

    @property
    @caching_decorator
    def lllNumInvs(self):
        return [oTerm.oNum.llInvs for oTerm in self if oTerm.am_I_a_symmetry is False]

    @property
    @caching_decorator
    def lllNumExps(self):
        return [oTerm.oNum.llExps for oTerm in self if oTerm.am_I_a_symmetry is False]

    @property
    @caching_decorator
    def llCoefs(self):
        return [oTerm.oNum.lCoefs for oTerm in self if oTerm.am_I_a_symmetry is False]

    @property
    @caching_decorator
    def llDenInvs(self):
        return [oTerm.oDen.lInvs for oTerm in self if oTerm.am_I_a_symmetry is False]

    @property
    @caching_decorator
    def llDenExps(self):
        return [oTerm.oDen.lExps for oTerm in self if oTerm.am_I_a_symmetry is False]

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

    @functools.cached_property
    def common_coeff_factor(self):
        arr = numpy.array([entry[0] for entry in flatten(self.llCoefs)])
        sign = -1 if (arr < 0).tolist().count(True) > (len(arr) / 2) else 1
        # Reason for removing outliers: if a numerator accidentaly has a prime factor which cancels with the denominator
        # then it may spoil the gcd computation for the common factor. Since this is unlikely, removing outliers fixes it.
        outlier_fractions = [0, 0.01, 0.05, 0.1]
        common_Q_factors = [get_common_Q_factor(arr, outlier_fraction) for outlier_fraction in outlier_fractions]
        max_denominators = [get_max_denominator(arr / common_Q_factor) for common_Q_factor in common_Q_factors]
        min_max_denom = min(max_denominators)
        return sign * common_Q_factors[max_denominators.index(min_max_denom)]

    @functools.cached_property
    def with_normalized_coeffs(self):
        return self / self.common_coeff_factor

    @property
    def max_abs_numerator(self):
        array = numpy.array([entry[0] for entry in flatten(self.llCoefs)])
        return get_max_abs_numerator(array)

    @property
    def max_denominator(self):
        array = numpy.array([entry[0] for entry in flatten(self.llCoefs)])
        return get_max_denominator(array)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

    def __call__(self, oParticles):
        InvsDict = {}
        for inv in self.invs_set:
            InvsDict[inv] = oParticles.compute(inv)
        SymInvsDict = {}
        for oTerm in self:
            if oTerm.am_I_a_symmetry is True and oTerm.tSym not in SymInvsDict.keys():
                NewDict = {}
                NewParticles = oParticles.image(inverse(oTerm.tSym[0]))
                if oTerm.tSym[1] is True:
                    NewParticles.angles_for_squares()
                for inv in self.invs_set:
                    NewDict[inv] = NewParticles.compute(inv)
                SymInvsDict[oTerm.tSym] = NewDict
        NumericalResult = 0
        for i, iTerm in enumerate(self):
            if iTerm.am_I_a_symmetry is True:
                something_was_not_a_symmetry = False
                for j, jTerm in enumerate(self[:i][::-1]):
                    if jTerm.am_I_a_symmetry and something_was_not_a_symmetry:
                        break
                    elif jTerm.am_I_a_symmetry is False:
                        something_was_not_a_symmetry = True
                        if iTerm.tSym[2] == "-":
                            NumericalResult -= jTerm(SymInvsDict[iTerm.tSym])
                        else:
                            NumericalResult += jTerm(SymInvsDict[iTerm.tSym])
            else:
                NumericalResult += iTerm(InvsDict)
        return NumericalResult

    def Image(self, Rule):
        if not hasattr(self, "multiplicity"):  # deduce multiplicity from length of permutation string.
            self.multiplicity = len(Rule[0])
        return Terms([oTerm.Image(Rule) for oTerm in self])

    def cluster(self, rule):
        oClusteredTerms = Terms([oTerm.cluster(rule) for oTerm in self])
        oClusteredTerms.multiplicity = len(rule)
        return oClusteredTerms

    def explicit_representation(self, raw=False):
        oTermsExplicit = copy.copy(self)
        del oTermsExplicit[:]
        for i, iTerm in enumerate(self):
            if iTerm.am_I_a_symmetry is True:
                something_was_not_a_symmetry = False
                to_be_appended = []
                for j, jTerm in enumerate(self[:i][::-1]):
                    if jTerm.am_I_a_symmetry and something_was_not_a_symmetry:
                        break
                    elif jTerm.am_I_a_symmetry is False:
                        something_was_not_a_symmetry = True
                        if raw is False:
                            to_be_appended.append(copy.copy(jTerm.Image(iTerm.tSym)))
                        else:
                            to_be_appended.append(copy.copy(jTerm.rawImage(iTerm.tSym)))
                if iTerm.tSym[-1] == "-":
                    oTermsExplicit += [-1 * entry for entry in to_be_appended[::-1]]
                else:
                    oTermsExplicit += to_be_appended[::-1]
            else:
                oTermsExplicit.append(copy.copy(iTerm))
        oTermsExplicit.update_invs_set()
        return oTermsExplicit

    def canonical_ordering(self):
        for oTerm in self:
            oTerm.canonical_ordering()

    def update_invs_set(self):
        invs_set = []
        for oTerm in self:
            if oTerm.am_I_a_symmetry is True:
                continue
            invs_set += flatten(oTerm.oNum.llInvs)
            invs_set += flatten(oTerm.oNum.lCommonInvs)
            invs_set += flatten(oTerm.oDen.lInvs)
        if hasattr(self, "multiplicity"):
            for oTerm in self:
                oTerm.multiplicity = self.multiplicity
        self.invs_set = set(invs_set)

    def relevant_old_Terms(self, i):
        oOldTerms = self.oUnknown.recursively_extract_terms()
        return [_oTerm for _oTerm in oOldTerms if not _oTerm.am_I_a_symmetry and _oTerm.is_fully_reduced and
                # all numerator invariants have changed by at most 1 power
                all([True if inv in _oTerm.oNum.llInvs[0] and abs(_oTerm.oNum.llExps[0][_oTerm.oNum.llInvs[0].index(inv)] - self[i].oNum.llExps[0][self[i].oNum.llInvs[0].index(inv)]) <= 1 or
                     inv not in _oTerm.oNum.llInvs[0] and self[i].oNum.llExps[0][self[i].oNum.llInvs[0].index(inv)] == 1 else False for inv in self[i].oNum.llInvs[0]]) and
                # all denomiantor invariants have changed by at most 1 power
                all([True if inv in _oTerm.oDen.lInvs and abs(_oTerm.oDen.lExps[_oTerm.oDen.lInvs.index(inv)] - self[i].oDen.lExps[self[i].oDen.lInvs.index(inv)]) <= 1 or
                     inv not in _oTerm.oDen.lInvs and self[i].oDen.lExps[self[i].oDen.lInvs.index(inv)] == 1 else False for inv in self[i].oDen.lInvs])]

    def remove_duplicates(self):
        NewTerms = Terms(list(set(self)))
        if hasattr(self, "oFittingSettings"):
            NewTerms.oFittingSettings = self.oFittingSettings
        if hasattr(self, "oUnknown"):
            NewTerms.oUnknown = self.oUnknown
        return NewTerms

    def do_exploratory_double_collinear_limits(self):
        oInvariants = Invariants(self.multiplicity, Restrict3Brackets=settings.Restrict3Brackets,
                                 Restrict4Brackets=settings.Restrict4Brackets, FurtherRestrict4Brackets=settings.FurtherRestrict4Brackets)
        other_invariants = [inv for inv in oInvariants.invs_2 + oInvariants.invs_3 + oInvariants.invs_s if all(
            [inv not in oTerm.oDen.lInvs + flatten(oTerm.oNum.llInvs) for oTerm in self if not oTerm.am_I_a_symmetry])]
        # check all double scalings of the invariants only in iDen with the ones not appearing anywhere
        _pair_invs, _pair_exps, _pair_friends = pair_scalings(self.oUnknown, self[0].oDen.lInvs, other_invariants, oInvariants.full)
        if hasattr(self.oUnknown, "save_original_unknown_call_cache"):
            self.oUnknown.save_original_unknown_call_cache()
        sys.stdout.write("\r Finished calculating pair scalings.",)
        sys.stdout.flush()
        __pair_invs = [entry for entry in _pair_invs]
        __pair_exps = [entry for entry in _pair_exps]
        for j, jpair in enumerate(__pair_invs):
            if (__pair_exps[j] == "F" or __pair_exps[j] <= -self[0].oDen.lExps[self[0].oDen.lInvs.index(jpair[0])]):
                index = _pair_invs.index(jpair)
                _pair_invs.pop(index)
                _pair_exps.pop(index)
                _pair_friends.pop(index)
        print(" The interesting ones are:                        ")
        for j, _pair_exp in enumerate(_pair_exps):
            if isinstance(_pair_exp, int):
                _pair_exps[j] = abs(_pair_exps[j])

        if _pair_exps != []:
            for i, _pair_exp in enumerate(_pair_exps):
                if _pair_exp != "F" and (_pair_exp == int(_pair_exp) or _pair_exp * 2 == int(_pair_exp * 2)):
                    _pair_exps[i] = abs(_pair_exps[i])
            col_width = max([len(_pair_inv[0]) + len(_pair_inv[1]) + 6 for _pair_inv in _pair_invs])
            for i in range(len(_pair_invs)):
                print(("[" + _pair_invs[i][0] + ", " + _pair_invs[i][1] + "]:").ljust(col_width) + str(_pair_exps[i]) + ", " + str(len(_pair_friends[i])))
            print("")
        else:
            print("Non-existent")

        # col_width = max([len(unicode(_pair_inv)) for _pair_inv in _pair_invs])
        # for j in range(len(_pair_invs)):
        #     print "{} ---> {} --- {}".format(unicode(_pair_invs[j]).ljust(col_width), _pair_exps[j], len(_pair_friends[j]))

        # try weights
        # _flat_friends = [item for entry in _pair_friends for item in entry if item not in self[0].oDen.lInvs]
        # _all_friends = list(set(_flat_friends))
        # _weights = []
        # for _friend in _all_friends:
        #     _weights += [_flat_friends.count(_friend)]
        # _switch = False
        # while _switch is False:
        #     _switch = True
        #     for j in range(len(_all_friends) - 1):
        #         if _weights[j + 1] > _weights[j]:
        #             _weights[j + 1], _weights[j] = (_weights[j], _weights[j + 1])
        #             (_all_friends[j + 1], _all_friends[j]) = (_all_friends[j], _all_friends[j + 1])
        #             _switch = False
        # col_width = max(len(friend) for friend in _all_friends)
        # for j in range(len(_all_friends)):
        #     print(unicode(_all_friends[j]).ljust(col_width) + " ---> " + unicode(_weights[j]))

    def check_md_pw_consistency(self):
        raise NotImplementedError

    def order_terms_by_complexity(self):
        lM = self.ansatze_mass_dimensions
        lPW = self.ansatze_phase_weights
        labsPW = [sum(map(abs, entry)) for entry in lPW]
        [self, lM, labsPW] = map(list, zip(*sorted(zip(self, lM, labsPW), key=lambda s: (s[1], -s[2]))))
        self = Terms(self)

    def rearrange_and_finalise(self):
        if len(self) < 2:
            return
        invariants_for_ordering = []
        powers_for_ordering = {}
        for oTerm in self:
            if oTerm.am_I_a_symmetry is True:
                continue
            for inv in oTerm.oDen.lInvs:
                if inv not in invariants_for_ordering:
                    invariants_for_ordering += [inv]
                    powers_for_ordering[inv] = oTerm.oDen.lExps[oTerm.oDen.lInvs.index(inv)]
                elif powers_for_ordering[inv] < oTerm.oDen.lExps[oTerm.oDen.lInvs.index(inv)]:
                    powers_for_ordering[inv] = oTerm.oDen.lExps[oTerm.oDen.lInvs.index(inv)]
        invariants_for_ordering = self.oUnknown._ordering_for_inversion_and_printing(invariants_for_ordering, powers_for_ordering)

        def switch(oDen1, oDen2):
            need_to_switch = False
            know_what_to_do = False
            i = 0
            zipped_invs_weights = list(zip(sorted([invariants_for_ordering.index(inv) for inv in oDen1.lInvs]), sorted([invariants_for_ordering.index(inv) for inv in oDen2.lInvs])))
            while know_what_to_do is False and i < len(zipped_invs_weights):
                # first compare the invariants (that is, their weight according to the ordering)
                if zipped_invs_weights[i][0] < zipped_invs_weights[i][1]:
                    know_what_to_do = True
                    need_to_switch = False
                elif zipped_invs_weights[i][0] > zipped_invs_weights[i][1]:
                    know_what_to_do = True
                    need_to_switch = True
                else:
                    know_what_to_do = False
                # if that is not conclusive compare exponents
                if know_what_to_do is False:
                    exp1 = oDen1.lExps[oDen1.lInvs.index(invariants_for_ordering[zipped_invs_weights[i][0]])]
                    exp2 = oDen2.lExps[oDen2.lInvs.index(invariants_for_ordering[zipped_invs_weights[i][1]])]
                    if exp1 > exp2:
                        know_what_to_do = True
                        need_to_switch = False
                    elif exp1 < exp2:
                        know_what_to_do = True
                        need_to_switch = True
                    else:
                        know_what_to_do = False
                # if even the exponents where the same, go to the next invariant
                if know_what_to_do is True:
                    break
                else:
                    i += 1
            return need_to_switch

        switched = True
        while switched is True:
            switched = False
            i = 0
            while i < len(self) - 1:
                if self[i].am_I_a_symmetry is True or self[i + 1].am_I_a_symmetry is True:
                    pass
                elif switch(self[i].oDen, self[i + 1].oDen) is True:
                    self[i], self[i + 1] = self[i + 1], self[i]
                    switched = True
                i += 1

    def summarise(self, print_ansatze_info=True):
        if print_ansatze_info is True:
            lM = self.ansatze_mass_dimensions
            lPW = self.ansatze_phase_weights
            for i, iTerm in enumerate(self):
                print(iTerm,)
                if iTerm.am_I_a_symmetry is False:
                    print(" ------- [{}, {}]".format(lM[i], lPW[i]))
                else:
                    print("")
        else:
            for i, iTerm in enumerate(self):
                print(iTerm)

    def collapse(self):
        if len(self) > 2:
            collapsed_den = self[1].oDen
            for oTerm in self[2:]:
                for i, inv in enumerate(oTerm.oDen.lInvs):
                    if inv in collapsed_den.lInvs and collapsed_den.lExps[collapsed_den.lInvs.index(inv)] < oTerm.oDen.lExps[i]:
                        collapsed_den.lExps[collapsed_den.lInvs.index(inv)] = oTerm.oDen.lExps[i]
                    elif inv not in collapsed_den.lInvs:
                        collapsed_den.lInvs += [inv]
                        collapsed_den.lExps += [oTerm.oDen.lExps[i]]
            self = self[0:2]

    def invs_only_in_iDen(self, i):
        oInvariants = Invariants(self.multiplicity, Restrict3Brackets=settings.Restrict3Brackets,
                                 Restrict4Brackets=settings.Restrict4Brackets, FurtherRestrict4Brackets=settings.FurtherRestrict4Brackets)
        iDen = self[i].oDen
        invs_only_in_iDen = [inv for inv in iDen.lInvs if all([j == i or inv not in jTerm.oDen.lInvs or len(self) == 1 for j, jTerm in enumerate(self)])]
        exps_only_in_iDen = [jExp for (j, jExp) in enumerate(self[i].oDen.lExps) if iDen.lInvs[j] in invs_only_in_iDen]
        if len(invs_only_in_iDen) > 1:
            # order them from most complicated to less complicated
            [invs_only_in_iDen, exps_only_in_iDen] = map(
                list, zip(*sorted(zip(invs_only_in_iDen, exps_only_in_iDen),
                                  key=lambda s: (s[0] in oInvariants.invs_2, s[1] == 1,
                                                 # s[0] not in oInvariants.invs_P, s[0] not in oInvariants.invs_O, s[0] not in oInvariants.invs_D,
                                                 s[0] not in oInvariants.invs_s, s[0] not in oInvariants.invs_5, s[0] not in oInvariants.invs_4,
                                                 s[0] not in oInvariants.invs_3, -s[1]))))
        return invs_only_in_iDen, exps_only_in_iDen

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

    def check_against_unknown(self):
        oParticles = Particles(self.multiplicity)
        for i in range(5):
            sys.stdout.write("\rChecking the analytical result against the numerical one: {}/5.".format(i))
            sys.stdout.flush()
            oParticles.randomise_all()
            oParticles.fix_mom_cons()
            # compare to numerical result
            if abs(self(oParticles) - self.oUnknown(oParticles)) > 10 ** -(0.9 * accuracy()):
                raise Exception("\nThis result can't be trusted. Something went wrong.")
        else:
            print("\rChecked the analytical result 5 times against the numerical one.")
            print("Never did the difference exceed 10^-{}.\n".format(int(math.floor(0.9 * accuracy()))))

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

    def Write_LaTex(self):
        # generate latex message and print the full result
        latex_result = ""
        for i, iTerm in enumerate(self):
            if hasattr(iTerm, "oNum"):
                latex_result += "("
                # build numerator, if there is more then one term in this numerator write the common invs at the beginning
                if len(iTerm.oNum.lCommonInvs) >= 1:
                    latex_result += Write(iTerm.oNum.lCommonInvs, iTerm.oNum.lCommonExps)[1:-3]
                    latex_result += "("
                for j in range(len(iTerm.oNum.llInvs)):
                    if len(iTerm.oNum.lCoefs) > 0:
                        real = iTerm.oNum.lCoefs[j][0]
                        imag = iTerm.oNum.lCoefs[j][1]
                        if real != 0 or imag != 0:
                            if real == 0:
                                latex_result += "{}i".format(imag)
                            elif imag == 0:
                                latex_result += "{}".format(real)
                            elif imag > 0:
                                latex_result += "{}+{}i".format(real, imag)
                            elif imag < 0:
                                latex_result += "{}{}i".format(real, imag)
                    else:
                        assert len(iTerm.oNum.llInvs) == 1
                        latex_result += "1"
                    jInvs = copy.copy(iTerm.oNum.llInvs[j])
                    jExps = copy.copy(iTerm.oNum.llExps[j])
                    latex_result += Write(jInvs, jExps)[1:-3]
                    latex_result += " + "
                latex_result = latex_result[:-3]
                if len(iTerm.oNum.lCommonInvs) >= 1:
                    latex_result += ")"
                if iTerm.oDen.lInvs != []:
                    latex_result += ")/"
                    # add denominator
                    latex_result += Write(iTerm.oDen.lInvs, iTerm.oDen.lExps)[:-2]
                else:
                    latex_result += ")"
            elif hasattr(iTerm, "tSym"):
                if iTerm.tSym[2] == "-":
                    latex_result += r"\scriptscriptstyle({},\;\text{{{}}},\;{})".format(iTerm.tSym[0], iTerm.tSym[1], iTerm.tSym[2])
                else:
                    latex_result += r"\scriptscriptstyle({},\;\text{{{}}})".format(iTerm.tSym[0], iTerm.tSym[1])
            else:
                raise Exception("Counld't understand a term while writing latex expression.")
            latex_result += " +\\\\\n"
        latex_result = latex_result[:-4]
        if len(self) > 1:
            latex_result += "\\phantom{+}"
        latex_result = latex_result.replace(".0", "")
        latex_result = latex_result.replace(" ", "")
        latex_result = latex_result.replace("|(", "|")
        latex_result = latex_result.replace(")|", "|")
        latex_result = re.sub(r"\[(?P<a>\d)\|(?P<b>\d)\]", r"[\g<a>\g<b>]", latex_result)
        latex_result = re.sub(r"⟨(?P<a>\d)\|(?P<b>\d)⟩", r"⟨\g<a>\g<b>⟩", latex_result)
        latex_result = latex_result.replace("+-", "-")
        latex_frac = re.compile(r"\((?P<num>[/i\-\+\d⟨⟩\[\]\|a-zA-ZΔΩΠ_\^\(\)(?!/).]*)\)/\((?P<den>[/i\-\+\d⟨⟩\[\]\|a-zA-ZΔΩΠ_\^\(\)(?!/)]*)\)")
        latex_result = re.sub(latex_frac, r"\\scriptscriptstyle\\frac{\g<num>}{\g<den>}", latex_result)
        latex_frac = re.compile(r"\((?P<number>[/i\-\+\d⟨⟩\[\]\|a-zA-ZΔΩΠ_\^\(\)(?!/).]*)\)(?!/|\}\{)")
        latex_result = re.sub(latex_frac, r"\\scriptscriptstyle(\g<number>)", latex_result)
        latex_s_ijk = re.compile(r"s_(?P<numbers>[\d]*)")
        latex_result = re.sub(latex_s_ijk, r"s_{\g<numbers>}", latex_result)
        latex_D_ijk_new = re.compile(r"Δ_(?P<numbers1>[\d]*)\|(?P<numbers2>[\d]*)\|(?P<numbers3>[\d]*)")
        latex_result = re.sub(latex_D_ijk_new, r"Δ_{\g<numbers1>|\g<numbers2>|\g<numbers3>}", latex_result)
        latex_D_ijk = re.compile(r"Δ_(?P<numbers>[\d]+)")
        latex_result = re.sub(latex_D_ijk, r"Δ_{\g<numbers>}", latex_result)
        latex_O_ijk = re.compile(r"Ω_(?P<numbers>[\d]*)")
        latex_result = re.sub(latex_O_ijk, r"Ω_{\g<numbers>}", latex_result)
        latex_P_ijk = re.compile(r"Π_(?P<numbers>[\d]*)")
        latex_result = re.sub(latex_P_ijk, r"Π_{\g<numbers>}", latex_result)
        latex_tr5_ijk = re.compile(r"tr5_(?P<numbers>[\d]*)")
        latex_result = re.sub(latex_tr5_ijk, r"tr5_{\g<numbers>}", latex_result)
        latex_result = latex_result.replace("+\n", "+")
        latex_result = re.sub(r'(?<!\n)(?<!^)\\scriptscriptstyle', '', latex_result)
        latex_message = "\\begin{my}\n$\\begin{gathered}\n" + latex_result + "\n\\end{gathered}$\n\\end{my}\n"
        return latex_message

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

    @staticmethod
    def _parse_term(iTerm, dsubs=None):
        """Returns a functional equivalent represerentation of the imputted string."""
        iTerm = non_unicode_powers(iTerm)
        iTerm = re.sub(r"⟨(\d)\|(\d)⟩", r"za(p\1,p\2)", iTerm)
        iTerm = re.sub(r"\[(\d)\|(\d)\]", r"zb(p\1,p\2)", iTerm)
        iTerm = re.sub(r"s_{0,1}(\d)(\d)(\d)", r"s(p\1,p\2,p\3)", iTerm)
        iTerm = re.sub(r"s_{0,1}(\d)(\d)", r"s(p\1,p\2)", iTerm)
        iTerm = pNB_overall.sub(lambda match: parse_pNBs_to_functions(match), iTerm)

        iTerm = iTerm.replace("Δ_", "delta").replace("Π_", "pi")
        iTerm = iTerm.replace("^", "**")
        iTerm = iTerm.replace("tr5_1234", "tr5").replace("tr5(1|2|3|4)", "tr5")
        iTerm = re.sub(r"\)([a-zA-Z])", r")*\1", iTerm)
        iTerm = iTerm.replace(")(", ")*(")
        iTerm = re.sub(r"(\d)([a-z])", r"\1*\2", iTerm)
        iTerm = re.sub(r"(?<!b)(?<!b\d)(?<!a\d)(?<!b\d\d)(?<!a\d\d)(?<!b\d\d\d)(?<!a\d\d\d)(?<!tr)(\d)\(", r"\1*(", iTerm)
        iTerm = iTerm.replace(r"pi135", r"(s561-s562)")
        iTerm = iTerm.replace(r"pi351", r"(s123-s124)")
        iTerm = iTerm.replace(r"pi513", r"(s345-s346)")
        iTerm = iTerm.replace(r"/z", r")/(z").replace("_", "").replace(")(", ")*(").replace("|", "c")

        if dsubs is not None:
            from .term import subs_dict
            iTerm = subs_dict(iTerm, dsubs)

        return iTerm

    def toFORM(self):
        form_expr = []
        for oTerm in self:
            num_expr = self._parse_term(str(oTerm.oNum))
            if oTerm.oNum.lCommonInvs == [] and len(oTerm.oNum.llInvs) > 1:
                num_expr = "(" + num_expr + ")"
            den_expr = self._parse_term(str(oTerm.oDen))
            den_expr = den_expr.replace('s', 'is').replace('z', 'iz').replace('delta', 'idelta')
            form_expr += [num_expr + '*' + den_expr]
        return "+".join(form_expr).replace('++', '+').replace('+-', '-')

    def toFortran(self, dsubs=None):
        fortran_exprs = []
        for oTerm in self:
            if oTerm.am_I_a_symmetry:
                fortran_exprs += [str(oTerm)]
            else:
                num_expr = self._parse_term(str(oTerm.oNum), dsubs=dsubs)
                if oTerm.oNum.lCommonInvs == [] and len(oTerm.oNum.llInvs) > 1:
                    num_expr = "(" + num_expr + ")"
                if oTerm.oDen != Denominator():
                    den_expr = "(" + self._parse_term(str(oTerm.oDen), dsubs=dsubs) + ")"
                    fortran_exprs += [num_expr + '/' + den_expr]
                else:
                    fortran_exprs += [num_expr]
        fortran_expr = "\n+".join(fortran_exprs).replace('++', '+').replace('+-', '-')
        fortran_expr = re.sub(r"(?<![a-z])(?<!\*\*)([\+|\-]{0,1}\d+)/(\d+)\*", r"\1.0_dp/\2.0_dp*", fortran_expr)   # parse rational numbers
        # fortran_expr = re.sub(r"(?<![a-z])(?<!\*\*)([\+|\-]{0,1}\d+)\*", r"\1.0_dp*", fortran_expr)   # parse integers
        return fortran_expr

    def toSympy(self):
        txt = self.toFortran()
        txt = re.sub(r'za\(p(\d),p(\d)\)', r'za\1\2', txt)
        txt = re.sub(r'zb\(p(\d),p(\d)\)', r'zb\1\2', txt)
        all_zas = set(re.findall(r'za\d\d', txt))
        all_zbs = set(re.findall(r'zb\d\d', txt))
        print(all_zas, all_zbs)
        sympy.var(' '.join(all_zas))
        sympy.var(' '.join(all_zbs))
        expr = eval(txt)
        return expr

    def toSaM(self):
        explicit_repr = self.explicit_representation()
        return string_toSaM(str(explicit_repr))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

def parse_pNBs_to_functions(match):
    """Converts a pNB to a function with the zab12 notation."""
    match = match.group(0)
    if "-" in match:
        return match
    abcd = pNB_internal.search(match)
    a = abcd.group('start')
    bc = abcd.group('middle').replace("(", "").replace(")", "").split("|")
    bcs = [entry.split("+") for entry in bc]
    d = abcd.group('end')
    if match[0] == "⟨":
        res = "za"
        res += "".join(['b' if i % 2 == 0 else 'a' for i in range(len(bcs))])
    else:
        res = "zb"
        res += "".join(['a' if i % 2 == 0 else 'b' for i in range(len(bcs))])
    res += "".join(map(str, [len(bc) for bc in bcs]))
    res += "(" + ",".join(flatten([a, ] + bcs + [d, ])) + ")"
    return res


def string_toSaM(string):
    string = non_unicode_powers(string)
    # expand differences
    string = re.sub(r"⟨(\d)\|(?:\(){0,1}(\d)\-(\d)(?:\)){0,1}\|(\d)\]", r"(⟨\1|\2⟩[\2|\4]-⟨\1|\3⟩[\3|\4])", string)
    string = re.sub(r"\[(\d)\|(?:\(){0,1}(\d)\-(\d)(?:\)){0,1}\|(\d)⟩", r"([\1|\2]⟨\2|\4⟩-[\1|\3]⟨\3|\4⟩)", string)
    # parse invariants
    string = pA2.sub(r"Spaa[Sp[\1],Sp[\2]]", string)
    string = pS2.sub(r"Spbb[Sp[\1],Sp[\2]]", string)
    string = pSijk.sub(lambda match: f"S[{','.join(f'Sp[{num}]' for num in match.group(1))}]", string)
    string = pDijk_non_adjacent.sub(lambda match: f"Delta[{','.join(['{' + ','.join(f'Sp[{num}]' for num in group) + '}' for group in match.groups() ])}]", string)

    def parse_pNB(match):
        abcd = pNB_internal.search(match.group())
        a = int(abcd.group('start'))
        bcs = abcd.group('middle').replace("(", "").replace(")", "").split("|")
        d = int(abcd.group('end'))
        for i, bc in enumerate(bcs):
            if "-" in bc:
                raise Exception(f"Can't export to S@M a spinor string containing a difference. Failed to convert {match.group()}.")
            bcs[i] = [f"Sp[{i}]" for i in bc.split("+")]
        bcs = ",".join(["{" + ",".join(bc) + "}" for bc in bcs])
        begin = match.group()[0]
        end = match.group()[-1]
        begin = "a" if begin == "⟨" else "b"
        end = "a" if end == "⟩" else "b"
        return f"Sp{begin}{end}[Sp[{a}],{bcs},Sp[{d}]]"

    string = pNB_overall.sub(lambda match: parse_pNB(match), string)

    string = string.replace("**", "^").replace("++", "+").replace("+-", "-")
    string = re.sub(r"tr5_(\d)(\d)(\d)(\d)",
                    r"(Spbb[Sp[\1],Sp[\2]]Spaa[Sp[\2],Sp[\3]]Spbb[Sp[\3],Sp[\4]]Spaa[Sp[\4],Sp[\1]]-Spaa[Sp[\1],Sp[\2]]Spbb[Sp[\2],Sp[\3]]Spaa[Sp[\3],Sp[\4]]Spbb[Sp[\4],Sp[\1]])",
                    string)

    if any([char in "\n".join(string) for char in "⟨⟩"]):
        raise Exception(f"Failed to S@M export for:\n{string}\nDetected an left-over angle bracket in parsed expression.")

    return string

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def LoadResults(res_path, load_partial_results_only=False, silent=True, callable_to_check_with=None):

    lResults, loaded_result_is_partial = LaTeXToPython(res_path, load_partial_results_only)

    if lResults is None:
        return None, None
    elif lResults == 0:
        if callable_to_check_with is not None:
            oUnknown = Unknown(callable_to_check_with)
            if not oUnknown.is_zero:
                raise Exception("The result at {} doesn't match the function {}.".format(res_path, callable_to_check_with))
            elif silent is False:
                print("The loaded result is correct.")
        return lResults, loaded_result_is_partial
    for i, oRes in enumerate(lResults):
        lResults[i] = oRes
        if callable_to_check_with is not None:
            lResults[i].multiplicity = callable_to_check_with.multiplicity
    if silent is False:
        print("Loaded:\n", lResults)
    if loaded_result_is_partial is False and callable_to_check_with is not None:
        oUnknown = Unknown(callable_to_check_with)
        for i, oRes in enumerate(lResults):
            oUnknown.add_partial_piece(oRes)
            if not oUnknown.is_zero:
                if len(lResults) == 1:
                    raise Exception("The result at {} doesn't match the function {}.".format(res_path, callable_to_check_with))
                else:
                    raise Exception("The {} result at {} doesn't match the function {}.".format(i, res_path, callable_to_check_with))
            oUnknown.reset()
        if silent is False and len(lResults) == 1:
            print("The loaded result is correct.")
        elif silent is False:
            print("The loaded results are correct.")
    return lResults, loaded_result_is_partial


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


class FittingSettings(object):
    def __init__(self):
        self.lSmallInvs = []
        self.lSmallInvsExps = []
        self.lSymmetries = []


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def fix_old_automatic_partial_fractioning_format(lTerms):
    for i, iTerms in enumerate(lTerms):
        iTerms.collapse()
        oFittingSettings = FittingSettings()
        oFittingSettings.lSmallInvs, oFittingSettings.lSmallInvsExps = iTerms.invs_only_in_iDen(0)
        if len(oFittingSettings.lSmallInvs) >= 1 and (oFittingSettings.lSmallInvsExps[0] != 1 or len(iTerms) != 1):
            oFittingSettings.lSmallInvs, oFittingSettings.lSmallInvsExps = [oFittingSettings.lSmallInvs[0]], [oFittingSettings.lSmallInvsExps[0]]
        else:
            oFittingSettings.lSmallInvs, oFittingSettings.lSmallInvsExps = [], []
        iTerms.oFittingSettings = oFittingSettings
        lTerms[i] = iTerms[:1]
