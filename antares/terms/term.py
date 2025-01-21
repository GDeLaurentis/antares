#!/usr/bin/env python
# -*- coding: utf-8 -*-

#   _  _                         _           __     ___                     _           _
#  | \| |_  _ _ __  ___ _ _ __ _| |_ ___ _ _/ _|___|   \ ___ _ _  ___ _ __ (_)_ _  __ _| |_ ___ _ _
#  | .` | || | '  \/ -_) '_/ _` |  _/ _ \ '_> _|_ _| |) / -_) ' \/ _ \ '  \| | ' \/ _` |  _/ _ \ '_|
#  |_|\_|\_,_|_|_|_\___|_| \__,_|\__\___/_| \_____||___/\___|_||_\___/_|_|_|_|_||_\__,_|\__\___/_|


# Author: Giuseppe

from __future__ import unicode_literals
from __future__ import print_function

import lips
import mpmath
import pyadic
import re
import sys

from functools import reduce
from fractions import Fraction
from operator import mul
from copy import copy, deepcopy

from antares.core.settings import settings
from antares.core.numerical_methods import Numerical_Methods
from antares.scalings.single import single_scalings

from lips.tools import subs_dict
from lips.invariants import Invariants
from lips.particles_eval import unicode_powers, non_unicode_powers

if sys.version_info.major > 2:
    unicode = str

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


class Term(Numerical_Methods, object):

    def __init__(self, object1, object2=None):
        if isinstance(object1, str):
            self.__rstr__(object1)
        elif isinstance(object1, Numerator):
            if not isinstance(object2, Denominator):
                raise Exception("Term object created with Numerator but no Denominator.")
            else:
                self.oNum = object1
                self.oDen = object2
        elif isinstance(object1, tuple):
            if object2 is not None:
                raise Exception("Term created with symmetry and denominator.")
            elif not isinstance(object1[0], unicode) or not isinstance(object1[1], bool):
                raise Exception("Term created with non vailid symmetry.")
            else:
                self.tSym = object1
        else:
            raise Exception("Bad constructor.")
        self.simplify_factored_monomials()

    @classmethod
    def from_single_scalings(cls, oUnknown, invariants=None, verbose=False):

        # Choose the variables
        if invariants is None:
            oInvariants = Invariants(oUnknown.multiplicity, Restrict3Brackets=settings.Restrict3Brackets,
                                     Restrict4Brackets=settings.Restrict4Brackets, FurtherRestrict4Brackets=settings.FurtherRestrict4Brackets)
            if settings.SingleScalingsUse4Brackets is True:
                invariants = oInvariants.full
            else:
                invariants = oInvariants.full_minus_4_brackets
        else:
            invariants = copy(invariants)

        # Calculate the exponents
        exponents = single_scalings(oUnknown, invariants, verbose)

        # Clean invariants and exponents of zeros and failed scalings
        for i in range(len(invariants)):
            if exponents[i] == 0 or exponents[i] is None or exponents[i] == "F":
                invariants[i] = None
                exponents[i] = None
        invariants = list(filter(None, invariants))
        exponents = list(filter(None, exponents))

        # Split invariants in numerator and denominator
        num_invs, num_exps, den_invs, den_exps = [], [], [], []
        for i in range(len(invariants)):
            if exponents[i] > 0:
                num_invs += [invariants[i]]
                num_exps += [exponents[i]]
            else:
                den_invs += [invariants[i]]
                den_exps += [abs(exponents[i])]

        # results at four-point phase space need tweaking because of non-uniqueness of singular limits
        # if oUnknown.multiplicity == 4:
        #     num_invs, num_exps = [], []

        oTerm = cls(Numerator([num_invs], [num_exps]), Denominator(den_invs, den_exps))
        oTerm.multiplicity = oUnknown.multiplicity

        if verbose:
            print("The partial result is:   ")
            print(oTerm)

        return oTerm

    @property
    def am_I_a_symmetry(self):
        if hasattr(self, "tSym"):
            return True
        else:
            return False

    @property
    def multiplicity(self):
        if hasattr(self, "_multiplicity"):
            return self._multiplicity
        elif hasattr(self, "oUnknown"):
            return self.oUnknown.multiplicity
        else:
            raise AttributeError("Term object has no attribute multiplicity")

    @multiplicity.setter
    def multiplicity(self, value):
        self._multiplicity = value
        if not self.am_I_a_symmetry:
            self.oNum.multiplicity = value
            self.oDen.multiplicity = value

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
        if not self.am_I_a_symmetry:
            self.oNum.internal_masses = value
            self.oDen.internal_masses = value

    def __call__(self, InvsDict_or_Particles=None):
        from antares.terms.terms import Terms
        if isinstance(InvsDict_or_Particles, dict):
            if len(list(InvsDict_or_Particles.values())) > 0:
                dtype = type(list(InvsDict_or_Particles.values())[0])
            else:
                dtype = mpmath.mpc
            NumericalDenominator = 1
            for inv, exp in zip(self.oDen.lInvs, self.oDen.lExps):
                NumericalDenominator *= InvsDict_or_Particles[inv]**exp
            NumericalCommonNumerator = 1
            for inv, exp in zip(self.oNum.lCommonInvs, self.oNum.lCommonExps):
                NumericalCommonNumerator *= InvsDict_or_Particles[inv]**exp
            NumericalNumerator = 0
            for lInvs, lExps, coef in zip(self.oNum.llInvs, self.oNum.llExps, self.oNum.lCoefs if self.oNum.lCoefs != [] else [(Fraction(1, 1), Fraction(0, 1))]):
                if dtype == mpmath.mpc or dtype == mpmath.mpf:
                    coef = [make_proper(coef[0]), make_proper(coef[1])]
                    NumericalTerm = mpmath.mpc(mpmath.mpf(coef[0][0]) + mpmath.mpf(coef[0][1]) / mpmath.mpf(coef[0][2]),
                                               mpmath.mpf(coef[1][0]) + mpmath.mpf(coef[1][1]) / mpmath.mpf(coef[1][2]))
                else:
                    assert coef[1] == 0
                    NumericalTerm = coef[0]
                for inv, exp in zip(lInvs, lExps):
                    NumericalTerm *= InvsDict_or_Particles[inv]**exp
                NumericalNumerator += NumericalTerm
            return NumericalNumerator * NumericalCommonNumerator / NumericalDenominator
        elif isinstance(InvsDict_or_Particles, lips.Particles):
            return Terms([self])(InvsDict_or_Particles)

    def Image(self, Rule):
        from antares.ansatze.eigenbasis import Image
        if self.am_I_a_symmetry:
            newSym = Image(self.tSym, Rule)
            oSymTerm = Term(newSym)
            if hasattr(self, "multiplicity"):
                oSymTerm.multiplicity = self.multiplicity
        else:
            den_sign = int(reduce(mul, [Image(inv, Rule)[1] ** exp for inv, exp in zip(self.oDen.lInvs, self.oDen.lExps)]))
            num_common_sign = int(reduce(mul, [Image(inv, Rule)[1] ** exp for inv, exp in zip(self.oNum.lCommonInvs, self.oNum.lCommonExps)]) if self.oNum.lCommonInvs != [] else 1)
            num_signs = map(int, ([num_common_sign * reduce(mul, [Image(inv, Rule)[1] ** exp for inv, exp in zip(lInvs, lExps)]) for lInvs, lExps in zip(self.oNum.llInvs, self.oNum.llExps)]
                                  if self.oNum.llInvs != [[]] else []))
            signs = [den_sign * num_sign for num_sign in num_signs]
            oSymNum = Numerator([[Image(inv, Rule)[0] for inv in lInvs] for lInvs in self.oNum.llInvs], self.oNum.llExps, [
                (sign * coef[0], sign * coef[1]) for sign, coef in zip(signs, self.oNum.lCoefs)], [Image(inv, Rule)[0] for inv in self.oNum.lCommonInvs], self.oNum.lCommonExps)
            oSymDen = Denominator([Image(inv, Rule)[0] for inv in self.oDen.lInvs], self.oDen.lExps)
            oSymTerm = Term(oSymNum, oSymDen)
            if hasattr(self, "multiplicity"):
                oSymTerm.multiplicity = self.multiplicity
            oSymTerm.canonical_ordering()
        return oSymTerm

    def rawImage(self, Rule):
        from antares.topologies.topology import convert_invariant
        if self.am_I_a_symmetry:
            raise Exception("rawImage of symmetry not implemented")
        else:
            oSymNum = Numerator([[convert_invariant(inv, Rule) for inv in lInvs] for lInvs in self.oNum.llInvs], self.oNum.llExps, self.oNum.lCoefs,
                                [convert_invariant(inv, Rule) for inv in self.oNum.lCommonInvs], self.oNum.lCommonExps)
            oSymDen = Denominator([convert_invariant(inv, Rule) for inv in self.oDen.lInvs], self.oDen.lExps)
            oSymTerm = Term(oSymNum, oSymDen)
        return oSymTerm

    def cluster(self, rule):
        if self.am_I_a_symmetry:
            return Term(cluster_symmetry(self.tSym, rule))
        else:
            oNum = Numerator([[cluster_invariant(inv, rule) for inv in lInvs] for lInvs in self.oNum.llInvs], self.oNum.llExps, self.oNum.lCoefs,
                             [cluster_invariant(inv, rule) for inv in self.oNum.lCommonInvs], self.oNum.lCommonExps)
            oDen = Denominator([cluster_invariant(inv, rule) for inv in self.oDen.lInvs], self.oDen.lExps)
            return Term(oNum, oDen)

    @property
    def is_fully_reduced(self):
        if not self.am_I_a_symmetry and len(self.oNum.lCoefs) == 1:
            return True
        else:
            return False

    def simplify_factored_monomials(self):
        """Cancels powers of manifestly common factors between numerator and denominator.
        Lighter than rerunning single scaling study, but less powerful, since single scalings
        can also handle cancellations involving non-trivial rewritings (e.g. shoutens).
        """
        if self.am_I_a_symmetry:
            return
        if self.oNum.lCommonInvs == [] and len(self.oNum.llInvs) > 1:
            return
        if self.oNum.lCommonInvs != []:
            common_num_invs = copy(self.oNum.lCommonInvs)
            common_num_exps = copy(self.oNum.lCommonExps)
        elif len(self.oNum.llInvs) == 1:
            common_num_invs = copy(self.oNum.llInvs[0])
            common_num_exps = copy(self.oNum.llExps[0])
        shared_factors = list(set(self.oDen.lInvs) & set(common_num_invs))
        if shared_factors == []:
            return
        for shared_factor in shared_factors:
            if common_num_exps[common_num_invs.index(shared_factor)] == self.oDen.lExps[self.oDen.lInvs.index(shared_factor)]:
                common_num_exps.pop(common_num_invs.index(shared_factor))
                common_num_invs.pop(common_num_invs.index(shared_factor))
                self.oDen.lExps.pop(self.oDen.lInvs.index(shared_factor))
                self.oDen.lInvs.pop(self.oDen.lInvs.index(shared_factor))
            elif common_num_exps[common_num_invs.index(shared_factor)] > self.oDen.lExps[self.oDen.lInvs.index(shared_factor)]:
                common_num_exps[common_num_invs.index(shared_factor)] -= self.oDen.lExps[self.oDen.lInvs.index(shared_factor)]
                self.oDen.lExps.pop(self.oDen.lInvs.index(shared_factor))
                self.oDen.lInvs.pop(self.oDen.lInvs.index(shared_factor))
            elif common_num_exps[common_num_invs.index(shared_factor)] < self.oDen.lExps[self.oDen.lInvs.index(shared_factor)]:
                self.oDen.lExps[self.oDen.lInvs.index(shared_factor)] -= common_num_exps[common_num_invs.index(shared_factor)]
                common_num_exps.pop(common_num_invs.index(shared_factor))
                common_num_invs.pop(common_num_invs.index(shared_factor))
        if self.oNum.lCommonInvs != []:
            self.oNum.lCommonInvs = common_num_invs
            self.oNum.lCommonExps = common_num_exps
        else:
            self.oNum.llInvs[0] = common_num_invs
            self.oNum.llExps[0] = common_num_exps

    def canonical_ordering(self):
        oInvariants = Invariants(self.multiplicity, Restrict3Brackets=settings.Restrict3Brackets,
                                 Restrict4Brackets=settings.Restrict4Brackets, FurtherRestrict4Brackets=settings.FurtherRestrict4Brackets)
        if self.am_I_a_symmetry is False:
            if len(self.oDen.lInvs) >= 1:
                self.oDen.lInvs, self.oDen.lExps = map(
                    list, zip(*sorted(zip(self.oDen.lInvs, self.oDen.lExps), key=lambda x: oInvariants.full.index(x[0]) if x[0] in oInvariants.full else 999)))
            for i, (lInvs, lExps) in enumerate(zip(self.oNum.llInvs, self.oNum.llExps)):
                if len(lInvs) >= 1:
                    self.oNum.llInvs[i], self.oNum.llExps[i] = map(
                        list, zip(*sorted(zip(lInvs, lExps), key=lambda x: oInvariants.full.index(x[0]) if x[0] in oInvariants.full else 999)))

    def __mul_or_div__(self, other, operation):
        if self.am_I_a_symmetry is True or other.am_I_a_symmetry is True:
            raise Exception("Division not defined for term containing symmetry.")
        if len(self.oNum.llInvs) > 1 or len(other.oNum.llInvs) > 1:
            raise Exception("Division not defined for terms with more than one set of invariants in numerator.")
        result = Term(Numerator(), Denominator())
        for inv in self.oNum.llInvs[0] + self.oDen.lInvs + other.oNum.llInvs[0] + other.oDen.lInvs:
            if inv in result.oNum.llInvs[0] or inv in result.oDen.lInvs:
                continue
            exp1 = self.oNum.llExps[0][self.oNum.llInvs[0].index(inv)] if inv in self.oNum.llInvs[0] else 0
            exp2 = self.oDen.lExps[self.oDen.lInvs.index(inv)] if inv in self.oDen.lInvs else 0
            exp3 = other.oNum.llExps[0][other.oNum.llInvs[0].index(inv)] if inv in other.oNum.llInvs[0] else 0
            exp4 = other.oDen.lExps[other.oDen.lInvs.index(inv)] if inv in other.oDen.lInvs else 0
            if operation == "mul":
                exp = exp1 - exp2 + exp3 - exp4
            elif operation == "div":
                exp = exp1 - exp2 - exp3 + exp4
            if exp > 0:
                result.oNum.llInvs[0] += [inv]
                result.oNum.llExps[0] += [exp]
            elif exp < 0:
                result.oDen.lInvs += [inv]
                result.oDen.lExps += [-exp]
        return result

    def __neg__(self):
        return -1 * self

    def __mul__(self, other):
        if (isinstance(other, int) or isinstance(other, Fraction) or (
                isinstance(other, complex) and other.real.is_integer() and other.imag.is_integer())):
            oResTerm = deepcopy(self)
            if oResTerm.am_I_a_symmetry is True:
                return oResTerm
            for i, coef in enumerate(oResTerm.oNum.lCoefs):
                new_coef = pyadic.GaussianRational(coef[0], coef[1]) * other
                oResTerm.oNum.lCoefs[i] = (new_coef.real, new_coef.imag)
            return oResTerm
        return self.__mul_or_div__(other, "mul")

    def __rmul__(self, other):
        if type(other) is int:
            return self * other

    def __truediv__(self, other):
        return self.__div__(other)

    def __div__(self, other):
        if (isinstance(other, int) or isinstance(other, Fraction) or (
                isinstance(other, complex) and other.real.is_integer() and other.imag.is_integer())):
            oResTerm = deepcopy(self)
            if oResTerm.am_I_a_symmetry is True:
                return oResTerm
            for i, coef in enumerate(oResTerm.oNum.lCoefs):
                oResTerm.oNum.lCoefs[i] = (coef[0] / other, coef[1] / other)
            return oResTerm
        return self.__mul_or_div__(other, "div")

    def __contains__(self, other):  # is other in self?
        if self.am_I_a_symmetry is True and other.am_I_a_symmetry is True:
            return True if self == other else False
        elif self.am_I_a_symmetry is True or other.am_I_a_symmetry is True:
            return False
        if other.oNum in self.oNum and other.oDen in self.oDen:
            return True
        else:
            return False

    def __unicode__(self):
        return str(self)

    def __str__(self):
        if hasattr(self, "tSym"):
            string = str(self.tSym)
        elif str(self.oDen) != "":
            string = str(self.oNum) + "/(" + str(self.oDen) + ")"
        else:
            string = str(self.oNum)
        if sys.version_info.major > 2:
            return string
        else:
            return string.encode('utf-8')

    def __rstr__(self, string):
        if "True" in string or "False" in string:
            # this is a symmetry
            string = string.replace("+(", "(")
            symmetry = eval(string)
            self.__init__(symmetry)
        else:
            # this is kinematic expression
            if ")/(" in string:
                numerator, denominator = string.split(")/(")
                numerator = numerator + ")"
                denominator = "(" + denominator
            else:
                numerator = string
                denominator = ""
            self.__init__(Numerator(numerator), Denominator(denominator))

    def __repr__(self):
        return str(self)

    def __hash__(self):
        return hash(str(self))

    def __eq__(self, other):
        return self.__hash__() == other.__hash__() and isinstance(other, Term)

    def __ne__(self, other):
        return self.__hash__() != other.__hash__() and isinstance(other, Term)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def make_proper(fraction):
    numerator = abs(fraction.numerator)
    denominator = fraction.denominator
    integer_part = numerator // denominator
    proper_numerator = numerator % denominator
    if fraction > 0:
        return (integer_part, proper_numerator, denominator)
    else:
        return (-integer_part, -proper_numerator, denominator)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


class Numerator(object):

    def __init__(self, llInvs=[[]], llExps=[[]], lCoefs=[], lCommonInvs=[], lCommonExps=[]):
        if isinstance(llInvs, str):
            self.__rstr__(llInvs)
        else:
            self.llInvs = deepcopy(llInvs)
            self.llExps = deepcopy(llExps)
            self.lCoefs = deepcopy(lCoefs)
            if lCommonInvs == []:
                self.isolate_common_invariants_and_exponents()
            else:
                self.lCommonInvs = deepcopy(lCommonInvs)
                self.lCommonExps = deepcopy(lCommonExps)

    def isolate_common_invariants_and_exponents(self):
        if len(self.llInvs) > 1:
            self.lCommonInvs = [inv for inv in self.llInvs[0] if all([inv in lInvs for lInvs in self.llInvs])]
            self.lCommonExps = [min([self.llExps[j][self.llInvs[j].index(inv)] for j in range(len(self.llInvs))]) for inv in self.lCommonInvs]
            self.llExps = [[Exp if self.llInvs[i][j] not in self.lCommonInvs else Exp - self.lCommonExps[self.lCommonInvs.index(self.llInvs[i][j])] for (j, Exp) in enumerate(lExps)]
                           for (i, lExps) in enumerate(self.llExps)]
            self.llInvs = [[inv for (j, inv) in enumerate(lInv) if self.llExps[i][j] != 0] for (i, lInv) in enumerate(self.llInvs)]
            self.llExps = [[exp for exp in lExp if exp != 0] for lExp in self.llExps]
        else:
            self.lCommonInvs = []
            self.lCommonExps = []

    def __unicode__(self):
        return str(self)

    def __repr__(self):
        return str(self)

    def __str__(self):
        lCoefsString = self.lCoefsString
        lInvsString = self.lInvsString
        if len(lInvsString) == 1 and len(lCoefsString) == 0:
            lCoefsString = [""]
        if self.lCommonInvs == []:
            string = "".join(lCoefsString[i] + lInvsString[i] for i in range(len(lCoefsString)))
            if string == "":
                string = "+(1)"
            else:
                string = "+(" + string + ")"
        else:
            string = "+" + self.CommonInvsString + "(" + "".join(lCoefsString[i] + lInvsString[i] for i in range(len(lCoefsString))) + ")"
        if sys.version_info.major > 2:
            return string
        else:
            return string.encode('utf-8')

    def __rstr__(self, string):
        if string == "+(1)" or string == "(1)" or string == "1":
            self.__init__()
            return
        string = non_unicode_powers(string)
        string = string.replace("|(", "|").replace(")|", "|")
        if string[0] == "+":
            string = string[1:]
        split_numerator = re.split(r"(?<!tr)(\()(?=[\+\-]{0,1}\d)", string)
        # print(split_numerator)
        if len(split_numerator) == 1:  # single monomial without (gaussian) rational coefficient
            lInvs, lExps = parse_monomial(split_numerator[0][1:-1])
            llInvs, llExps = [lInvs], [lExps]
            lCommonInvs, lCommonExps, lCoefs = [], [], []
        elif len(split_numerator) == 3:   # factored monomial times polynomial - factored monomial may be empty
            # print(split_numerator)
            common_numerator = split_numerator[0]
            if common_numerator.count("(") > common_numerator.count(")") and common_numerator[0] == "(":
                common_numerator = common_numerator[1:]
            lCommonInvs, lCommonExps = parse_monomial(common_numerator)
            lCommonExps = list(map(int, lCommonExps))
            rest_numerator = split_numerator[1] + split_numerator[2]
            if rest_numerator.count("(") < rest_numerator.count(")") and rest_numerator[-1] == ")":
                rest_numerator = rest_numerator[:-1]
            rest_numerator = rest_numerator[1:-1]
            split_rest_numerator = re.split(r"(?<!\||_)(?<![\(\|\+\-]\d)((?:^|[\+\-])\d+[/\d+]*[IiJj]{0,1})(?!\||_)", rest_numerator)
            # print(split_rest_numerator)
            split_rest_numerator_fixed = []  # rejoin split pieces until parentheses are balanced
            comulative_string = ""
            for entry in split_rest_numerator:
                comulative_string += entry
                if comulative_string.count("(") == comulative_string.count(")"):
                    split_rest_numerator_fixed += [comulative_string]
                    comulative_string = ""
            split_rest_numerator = split_rest_numerator_fixed
            split_rest_numerator = [entry for entry in split_rest_numerator if entry != '']
            coeffs, monomials = split_rest_numerator[::2], split_rest_numerator[1::2]
            llInvs, llExps = map(list, zip(*list(map(parse_monomial, monomials))))
            lCoefs = list(map(lambda x: (Fraction(x), 0) if not any([imag_unit_char in x for imag_unit_char in "IiJj"])
                              else (0, Fraction(x.replace("I", "").replace("i", "").replace("J", "").replace("j", ""))), coeffs))
        else:
            raise Exception("Numerator string not understood (split).")
        llInvsCompactified = [list(set(lInvs)) for lInvs in llInvs]
        llExpsCompactified = [[sum([exp for j, exp in enumerate(lExps) if llInvs[i][j] == inv]) for inv in llInvsCompactified[i]]
                              for i, lExps in enumerate(llExps)]
        self.__init__(llInvsCompactified, llExpsCompactified, lCoefs, lCommonInvs, lCommonExps)

    def __contains__(self, other):
        if not all([inv in self.llInvs[0] for inv in other.llInvs[0]]):
            return False  # all zeros of other are also zeros of self
        if not all([other.llExps[0][other.llInvs[0].index(inv)] <= self.llExps[0][self.llInvs[0].index(inv)] for inv in other.llInvs[0]]):
            return False  # for each zero in other, that same zero in self has at least the degree of other
        if other.lCommonInvs != [] or self.lCommonInvs != []:
            print("detected common invs!")
            return False
        return True

    @property
    def lCoefsString(self):
        lCoefsString = []
        for i, complex_number in enumerate(self.lCoefs):
            complex_fraction_as_string = ""
            if complex_number[0] != 0:
                if (i != 0 or len(self.lCoefs) == 1) and complex_number[0] > 0:
                    complex_fraction_as_string += "+" + list(map(str, complex_number))[0]
                else:
                    complex_fraction_as_string += list(map(str, complex_number))[0]
            if complex_number[1] > 0:
                complex_fraction_as_string += "+" + list(map(str, complex_number))[1] + "I"
            elif complex_number[1] < 0:
                complex_fraction_as_string += list(map(str, complex_number))[1] + "I"
            lCoefsString += [complex_fraction_as_string]
        return lCoefsString

    @property
    def CommonInvsString(self):
        return unicode_powers(re.sub(r"\^1(?!\.|\d)", "", "".join(["^".join(map(unicode, entry)) for entry in zip(self.lCommonInvs, self.lCommonExps)]).replace(".0", "")))

    @property
    def lInvsString(self):
        return [unicode_powers(re.sub(r"\^1(?!\.|\d)", "", "".join(["^".join(map(unicode, entry)) for entry in zip(self.llInvs[i], self.llExps[i])]).replace(".0", "")))
                for i in range(len(self.llInvs))]


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


class Denominator(object):

    def __init__(self, lInvs=[], lExps=[]):
        if isinstance(lInvs, str):
            self = self.__rstr__(lInvs)
        else:
            self.lInvs = deepcopy(lInvs)
            self.lExps = deepcopy(lExps)

    def __unicode__(self):
        return str(self)

    def __repr__(self):
        return str(self)

    def __str__(self):
        string = unicode_powers(re.sub(r"\^1(?!\.)", "", "".join(["^".join(map(str, entry)) for entry in zip(self.lInvs, self.lExps)]).replace(".0", "")))
        if sys.version_info.major > 2:
            return string
        else:
            return string.encode('utf-8')

    def __rstr__(self, string):
        string = non_unicode_powers(string)[1:-1]
        string = string.replace("|(", "|").replace(")|", "|")
        lInvs, lExps = parse_monomial(string)
        self.__init__(lInvs, lExps)

    def __hash__(self):
        return hash(str(self))

    def __eq__(self, other):
        return self.__hash__() == other.__hash__() and isinstance(other, Denominator)

    def __ne__(self, other):
        return self.__hash__() != other.__hash__() and isinstance(other, Denominator)

    def __contains__(self, other):  # is other in self?
        if not all([inv in self.lInvs for inv in other.lInvs]):
            return False  # all poles of other are also poles of self
        if not all([other.lExps[other.lInvs.index(inv)] <= self.lExps[self.lInvs.index(inv)] for inv in other.lInvs]):
            return False  # for each pole in other, that same pole in self has at least the degree of other
        return True

    def as_dict(self):
        return dict(zip(self.lInvs, self.lExps))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def parse_monomial(string):
    if string == '':
        return [], []
    splitted_string = [entry for entry in re.split(r"(?<![\+\-\(])(?<![a-zA-Z])(?<!tr5)([⟨\[\(a-zA-ZΔΩΠ])", string) if entry != '']
    splitted_string = [splitted_string[i] + splitted_string[i + 1] for i in range(len(splitted_string))[::2]]
    # sqeuentially remerge strings until parenthesis are (minimally) balanced
    splitted_string_partially_remerged = [splitted_string[0]]
    for entry in splitted_string[1:]:
        if splitted_string_partially_remerged[-1].count("(") != splitted_string_partially_remerged[-1].count(")"):
            splitted_string_partially_remerged[-1] += entry
        else:
            splitted_string_partially_remerged += [entry]
    invs, exps = list(map(list, zip(*[re.split(r"\^(?=\d+$)", entry) if re.findall(r"\^(?=\d+$)", entry) != []
                                      else [entry, '1'] for entry in splitted_string_partially_remerged])))
    return invs, list(map(int, exps))


def cluster_symmetry(symmetry, rule):
    drule = dict(zip(["".join(map(str, entry)) for entry in rule], map(str, range(1, len(rule) + 1))))
    return (subs_dict(symmetry[0], drule), ) + symmetry[1:]


def cluster_invariant(invariant, rule):
    drule1 = dict(zip(["".join(map(str, entry)) for entry in rule], map(str, range(1, len(rule) + 1))))
    drule2 = dict(zip(["+".join(map(str, entry)) for entry in rule], map(str, range(1, len(rule) + 1))))
    drule3 = dict(zip(["-" + "-".join(map(str, entry)) for entry in rule], [f"-{i}" for i in range(1, len(rule) + 1)]))
    drule = drule1 | drule2 | drule3
    return subs_dict(invariant, drule)
