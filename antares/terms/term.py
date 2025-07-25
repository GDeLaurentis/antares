#   _  _                         _           __     ___                     _           _
#  | \| |_  _ _ __  ___ _ _ __ _| |_ ___ _ _/ _|___|   \ ___ _ _  ___ _ __ (_)_ _  __ _| |_ ___ _ _
#  | .` | || | '  \/ -_) '_/ _` |  _/ _ \ '_> _|_ _| |) / -_) ' \/ _ \ '  \| | ' \/ _` |  _/ _ \ '_|
#  |_|\_|\_,_|_|_|_\___|_| \__,_|\__\___/_| \_____||___/\___|_||_\___/_|_|_|_|_||_\__,_|\__\___/_|

# Author: Giuseppe

import numpy
import mpmath
import re
import functools
import operator
import warnings

from collections.abc import Sequence
from functools import reduce
from fractions import Fraction
from operator import mul
from copy import copy, deepcopy

from lips.tools import subs_dict
from lips.invariants import Invariants

from syngular import Field, Monomial, Polynomial

from ..core.tools import flatten
from ..core.settings import settings
from ..core.numerical_methods import Numerical_Methods
from ..scalings.single import single_scalings


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
            elif not isinstance(object1[0], str) or not isinstance(object1[1], bool):
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
        exponents = single_scalings(oUnknown, invariants, verbose=verbose)

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

        oTerm = cls(Numerator([num_invs], [num_exps], [1]), Denominator(den_invs, den_exps))
        oTerm.multiplicity = oUnknown.multiplicity

        if verbose:
            print("The partial result is:   ")
            print(oTerm)

        return oTerm

    @property
    def is_symmetry(self):
        return hasattr(self, "tSym")

    @property
    def am_I_a_symmetry(self):
        warnings.warn(
            "am_I_a_symmetry is deprecated, use is_symmetry instead.",
            category=DeprecationWarning,
            stacklevel=2,
        )
        return self.is_symmetry

    @property
    def ansatz(self):
        return self.oNum.polynomial.monomials

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
        if not self.is_symmetry:
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
        if not self.is_symmetry:
            self.oNum.internal_masses = value
            self.oDen.internal_masses = value

    def __getitem__(self, item):

        # NumPy boolean mask
        if (
            (isinstance(item, numpy.ndarray) and item.dtype == bool) or
            (isinstance(item, Sequence) and all(isinstance(i, bool) for i in item))
        ):
            result = deepcopy(self)
            result.oNum.polynomial = self.oNum.polynomial[item]

        else:
            raise NotImplementedError(f"Unsupported indexing type: {type(item)}")

        return result

    def __call__(self, InvsDict_or_Particles):
        from antares.terms.terms import Terms
        if type(InvsDict_or_Particles) is dict:
            field = InvsDict_or_Particles['field']
            # denominator monomial & common (factored) numerator monomial
            NumericalDenominator = self.oDen(InvsDict_or_Particles)
            NumericalCommonNumerator = self.oNum.monomial(InvsDict_or_Particles)
            # numerator polynomial (with or without contraction with (gaussian) rational coefficients)
            if all(entry is None for entry in self.oNum.lCoefs):
                rat_coefs = None
            else:
                rat_coefs = []
                for coef in self.oNum.lCoefs:
                    if field.characteristic == 0:
                        coef = [make_proper(coef[0]), make_proper(coef[1])]   # is this make proper story really needed ?!
                        coef = mpmath.mpc(mpmath.mpf(coef[0][0]) + mpmath.mpf(coef[0][1]) / mpmath.mpf(coef[0][2]),
                                          mpmath.mpf(coef[1][0]) + mpmath.mpf(coef[1][1]) / mpmath.mpf(coef[1][2]))
                    else:
                        assert coef[1] == 0
                        coef = field(coef[0])
                    rat_coefs += [coef]
                rat_coefs = numpy.array(rat_coefs)

            # print(InvsDict_or_Particles)

            run_on_gpu = settings.UseGpu and (len(self.oNum.llInvs) > 1000 or rat_coefs is None)

            if run_on_gpu:
                from linac import load_matrices
                monomials = [flatten([[inv, ] * exp for inv, exp in zip(lInvs, lExps)]) for lInvs, lExps in zip(self.oNum.llInvs, self.oNum.llExps)]
                prefactor = NumericalCommonNumerator / NumericalDenominator
                numerical_monomials = load_matrices(
                    ['prefactor'], [monomials], {'prefactor': prefactor} | InvsDict_or_Particles, use_cuda=True)[0].T
            else:
                numerical_monomials = []
                for lInvs, lExps in zip(self.oNum.llInvs, self.oNum.llExps):
                    numerical_monomial = 1
                    for inv, exp in zip(lInvs, lExps):
                        numerical_monomial = numerical_monomial * InvsDict_or_Particles[inv] ** exp
                    numerical_monomials += [numerical_monomial]
                numerical_monomials = numpy.array(numerical_monomials)

            if rat_coefs is None:
                # print(numpy.atleast_2d(numerical_monomials).shape)
                if not run_on_gpu:
                    return ((NumericalCommonNumerator / NumericalDenominator) * numpy.atleast_2d(numerical_monomials)).T
                else:
                    return numpy.atleast_2d(numerical_monomials).T
            else:
                numerical_poly = numpy.einsum("i,i...->...", rat_coefs, numerical_monomials)  # @ dot product fails with tensor monomials
                if not run_on_gpu:
                    return NumericalCommonNumerator / NumericalDenominator * numerical_poly
                else:
                    return numerical_poly

        elif callable(InvsDict_or_Particles):
            return Terms([self])(InvsDict_or_Particles)

    def Image(self, Rule):
        from antares.ansatze.eigenbasis import Image
        if self.is_symmetry:
            newSym = Image(self.tSym, Rule)
            oSymTerm = Term(newSym)
            if hasattr(self, "multiplicity"):
                oSymTerm.multiplicity = self.multiplicity
        else:
            den_sign = int(reduce(mul, [Image(inv, Rule)[1] ** exp for inv, exp in zip(self.oDen.lInvs, self.oDen.lExps)], 1))
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
        if self.is_symmetry:
            raise Exception("rawImage of symmetry not implemented")
        else:
            oSymNum = Numerator([[convert_invariant(inv, Rule) for inv in lInvs] for lInvs in self.oNum.llInvs], self.oNum.llExps, self.oNum.lCoefs,
                                [convert_invariant(inv, Rule) for inv in self.oNum.lCommonInvs], self.oNum.lCommonExps)
            oSymDen = Denominator([convert_invariant(inv, Rule) for inv in self.oDen.lInvs], self.oDen.lExps)
            oSymTerm = Term(oSymNum, oSymDen)
        return oSymTerm

    def cluster(self, rule):
        if self.is_symmetry:
            return Term(cluster_symmetry(self.tSym, rule))
        else:
            oNum = Numerator([[cluster_invariant(inv, rule) for inv in lInvs] for lInvs in self.oNum.llInvs], self.oNum.llExps, self.oNum.lCoefs,
                             [cluster_invariant(inv, rule) for inv in self.oNum.lCommonInvs], self.oNum.lCommonExps)
            oDen = Denominator([cluster_invariant(inv, rule) for inv in self.oDen.lInvs], self.oDen.lExps)
            return Term(oNum, oDen)

    @property
    def is_fully_reduced(self):
        if not self.is_symmetry and len(self.oNum.lCoefs) == 1:
            return True
        else:
            return False

    def simplify_factored_monomials(self):
        """Cancels powers of manifestly common factors between numerator and denominator.
        Lighter than rerunning single scaling study, but less powerful, since single scalings
        can also handle cancellations involving non-trivial rewritings (e.g. shoutens).
        """
        if self.is_symmetry:
            return
        elif len(self.oNum.monomial) == 0 and len(self.oNum.polynomial) > 1:
            return
        else:
            if len(self.oNum.monomial) != 0:
                monomial_to_cancel = self.oNum.monomial & self.oDen
            else:
                monomial_to_cancel = self.oNum.polynomial.monomials[0] & self.oDen
            self.oNum /= monomial_to_cancel
            self.oDen /= monomial_to_cancel

    def canonical_ordering(self):
        oInvariants = Invariants(self.multiplicity, Restrict3Brackets=settings.Restrict3Brackets,
                                 Restrict4Brackets=settings.Restrict4Brackets, FurtherRestrict4Brackets=settings.FurtherRestrict4Brackets)
        if self.is_symmetry is False:
            if len(self.oDen.lInvs) >= 1:
                self.oDen = Denominator(sorted(zip(self.oDen.lInvs, self.oDen.lExps),
                                               key=lambda x: oInvariants.full.index(x[0]) if x[0] in oInvariants.full else 999))
            for i, (lInvs, lExps) in enumerate(zip(self.oNum.llInvs, self.oNum.llExps)):
                if len(lInvs) >= 1:
                    self.oNum.llInvs[i], self.oNum.llExps[i] = map(
                        list, zip(*sorted(zip(lInvs, lExps), key=lambda x: oInvariants.full.index(x[0]) if x[0] in oInvariants.full else 999)))

    def __mul_or_div__(self, other, operation):
        from .terms import Terms
        if isinstance(other, Terms) and len(other) == 1:
            other = other[0]
        if not isinstance(other, Term):
            raise Exception(f"Attempted to {operation} Term object by {other} of type {type(other)}.")
        if self.is_symmetry is True or other.is_symmetry is True:
            raise Exception("Division not defined for term containing symmetry.")
        if len(self.oNum.llInvs) > 1 or len(other.oNum.llInvs) > 1:
            raise Exception("Division not defined for terms with more than one set of invariants in numerator.")
        if operation == "divide":
            return Term(Numerator(Monomial(), Polynomial([
                (self.oNum.polynomial.coeffs[0] / other.oNum.polynomial.coeffs[0],
                 self.oNum.polynomial.monomials[0] * other.oDen)],
                self.oNum.polynomial.field)),
                Denominator(self.oDen * other.oNum.polynomial.monomials[0]))
        elif operation == "multiply":
            return Term(Numerator(Monomial(), Polynomial([
                (self.oNum.polynomial.coeffs[0] * other.oNum.polynomial.coeffs[0],
                 self.oNum.polynomial.monomials[0] * other.oNum.polynomial.monomials[0])],
                self.oNum.polynomial.field)),
                Denominator(self.oDen * other.oDen))
        else:
            raise Exception(f"Operation {operation} between Terms not understood.")

    def __neg__(self):
        return -1 * self

    def __mul__(self, other):
        if (isinstance(other, int) or isinstance(other, Fraction) or (
                isinstance(other, complex) and other.real.is_integer() and other.imag.is_integer())):
            oResTerm = deepcopy(self)
            if oResTerm.is_symmetry is True:
                return oResTerm
            oResTerm.oNum.polynomial = Polynomial(
                [(coeff * other, monomial) for coeff, monomial in oResTerm.oNum.polynomial.coeffs_and_monomials],
                oResTerm.oNum.polynomial.field
            )
            return oResTerm
        return self.__mul_or_div__(other, "multiply")

    def __rmul__(self, other):
        if type(other) is int:
            return self * other

    def __truediv__(self, other):
        return self.__div__(other)

    def __div__(self, other):
        if (isinstance(other, int) or isinstance(other, Fraction) or (
                isinstance(other, complex) and other.real.is_integer() and other.imag.is_integer())):
            oResTerm = deepcopy(self)
            if oResTerm.is_symmetry is True:
                return oResTerm
            oResTerm.oNum.polynomial = Polynomial(
                [(coeff / other, monomial) for coeff, monomial in oResTerm.oNum.polynomial.coeffs_and_monomials],
                oResTerm.oNum.polynomial.field
            )
            return oResTerm
        return self.__mul_or_div__(other, "divide")

    def __contains__(self, other):  # is other in self?
        if self.is_symmetry is True and other.is_symmetry is True:
            return True if self == other else False
        elif self.is_symmetry is True or other.is_symmetry is True:
            return False
        if other.oNum in self.oNum and other.oDen in self.oDen:
            return True
        else:
            return False

    def __str__(self):
        if hasattr(self, "tSym"):
            return str(self.tSym)
        elif str(self.oDen) != "":
            return str(self.oNum) + "/(" + str(self.oDen) + ")"
        else:
            return str(self.oNum)

    def __rstr__(self, string):
        if "True" in string or "False" in string:  # this is a symmetry
            string = string.replace("+(", "(")
            symmetry = eval(string)
            self.__init__(symmetry)
        else:  # this is kinematic expression
            if ")/(" in string:
                numerator, denominator = string.split(")/(")
                numerator = numerator + ")"
                denominator = "(" + denominator
            else:
                numerator = string
                denominator = ""
            # print(numerator, denominator)
            self.__init__(Numerator(numerator), Denominator(denominator))

    def __repr__(self):
        return f"Term(\"{self}\")"

    def __hash__(self):
        return hash(str(self))

    def __eq__(self, other):
        return self.__hash__() == other.__hash__() and isinstance(other, Term)

    def __ne__(self, other):
        return self.__hash__() != other.__hash__() and isinstance(other, Term)

    @property
    def variables(self):
        if self.is_symmetry:
            return set()
        return self.oNum.variables | self.oDen.variables


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

    def __init__(self, linvs=[[]], lexps=[[]], coeffs=[], common_invs=[], common_exps=[], field=Field("Qi", 0, 0)):
        # print("in init", linvs, lexps, coeffs, common_invs, common_exps)
        if isinstance(linvs, str) and isinstance(lexps, list):
            self.__rstr__(linvs)
        elif isinstance(linvs, Monomial) and isinstance(lexps, Polynomial):
            self.monomial = linvs
            self.polynomial = lexps
        elif isinstance(linvs, str) and isinstance(lexps, str):
            self.monomial = Monomial(linvs)
            self.polynomial = Polynomial(lexps, field)
        else:
            self.monomial = Monomial(common_invs, common_exps)
            self.polynomial = Polynomial(list(zip(coeffs, [Monomial(invs, exps) for invs, exps in zip(linvs, lexps)])), field=field)
        if len(self.polynomial) == 1 and len(self.monomial) != 0:
            self.polynomial *= self.monomial
            self.monomial /= self.monomial
        if len(self.polynomial) > 1:
            extra_common_monomial = functools.reduce(operator.and_, self.polynomial.monomials)
            self.polynomial = self.polynomial / extra_common_monomial
            self.monomial = self.monomial * extra_common_monomial
        else:
            assert self.monomial == Monomial("")

    def __eq__(self, other):
        return isinstance(other, Numerator) and self.monomial == other.monomial and self.polynomial == other.polynomial

    def __repr__(self):
        return f"Numerator(\"{self}\")"

    def __str__(self):
        return f"{self.monomial}({self.polynomial})"

    def __rstr__(self, string):
        string = " ".join(string.split())
        string = string.replace("|(", "|").replace(")|", "|")
        if string[0] == "+":
            string = string[1:]
        split_numerator = re.split(r"(?<!tr)(?<!tr5)(\()(?=[\+\-]{0,1}[?\d])", string)
        # print(split_numerator)
        if len(split_numerator) == 1:  # single monomial without gaussian rational coefficient
            split_numerator = split_numerator[0][1:-1]  # remove parenthesis
            self.__init__('', split_numerator)
        elif len(split_numerator) == 3:   # factored monomial times polynomial - factored monomial may be empty
            common_numerator = split_numerator[0]
            if common_numerator.count("(") > common_numerator.count(")") and common_numerator[0] == "(":
                common_numerator = common_numerator[1:]
            rest_numerator = split_numerator[1] + split_numerator[2]
            if rest_numerator.count("(") < rest_numerator.count(")") and rest_numerator[-1] == ")":
                rest_numerator = rest_numerator[:-1]
            rest_numerator = rest_numerator[1:-1]
            self.__init__(common_numerator, rest_numerator)
        else:
            raise Exception("Numerator string not understood (split).")

    def __contains__(self, other):
        if not all([inv in self.llInvs[0] for inv in other.llInvs[0]]):
            return False  # all zeros of other are also zeros of self
        if not all([other.llExps[0][other.llInvs[0].index(inv)] <= self.llExps[0][self.llInvs[0].index(inv)] for inv in other.llInvs[0]]):
            return False  # for each zero in other, that same zero in self has at least the degree of other
        if other.lCommonInvs != [] or self.lCommonInvs != []:
            print("detected common invs!")
            return False
        return True

    def __truediv__(self, other):
        if isinstance(other, Monomial):
            if len(self.monomial) != 0 and other.issubset(self.monomial):
                return Numerator(self.monomial / other, self.polynomial)
            elif len(self.polynomial) == 1 and other.issubset(self.polynomial.monomials[0]):
                return Numerator(self.monomial, self.polynomial / other)
            else:
                return NotImplemented
        else:
            return NotImplemented

    @property
    def lCommonInvs(self):
        return self.monomial.invs

    @lCommonInvs.setter
    def lCommonInvs(self, value):
        raise AttributeError("lCommonInvs is a legacy read-only property.")

    @property
    def lCommonExps(self):
        return self.monomial.exps

    @lCommonExps.setter
    def lCommonExps(self, value):
        raise AttributeError("lCommonExps is a legacy read-only property.")

    @property
    def llInvs(self):
        return self.polynomial.linvs

    @llInvs.setter
    def llInvs(self, value):
        raise AttributeError("llInvs is a legacy read-only property.")

    @property
    def llExps(self):
        return self.polynomial.lexps

    @llExps.setter
    def llExps(self, value):
        raise AttributeError("llExps is a legacy read-only property.")

    @property
    def lCoefs(self):
        return self.polynomial.coeffs

    @lCoefs.setter
    def lCoefs(self, value):
        raise AttributeError("lCoefs is a legacy read-only property.")

    @property
    def coeffs(self):
        return self.polynomial.coeffs

    @coeffs.setter
    def coeffs(self, value):
        self.polynomial.coeffs = value

    @property
    def variables(self):
        return self.polynomial.variables | self.monomial.variables

    def as_term(self):
        return Term(str(self))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


class Denominator(Monomial):

    @property
    def lInvs(self):
        return self.invs

    @lInvs.setter
    def lInvs(self, value):
        raise AttributeError("lInvs is a legacy read-only property.")

    @property
    def lExps(self):
        return self.exps

    @lExps.setter
    def lExps(self, value):
        raise AttributeError("lExps is a legacy read-only property.")

    def __rstr__(self, string):
        string = string[1:-1]
        string = string.replace("|(", "|").replace(")|", "|")  # fixes spinor chain notation
        return super(Denominator, self).__rstr__(string)

    def __repr__(self):
        return f"Denominator(\"{self}\")"

    def as_term(self):
        return Term('(1' + str(self) + ')')


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def cluster_symmetry(symmetry, rule):
    drule = dict(zip(["".join(map(str, entry)) for entry in rule], map(str, range(1, len(rule) + 1))))
    return (subs_dict(symmetry[0], drule), ) + symmetry[1:]


def cluster_invariant(invariant, rule):
    drule1 = dict(zip(["".join(map(str, entry)) for entry in rule], map(str, range(1, len(rule) + 1))))
    drule2 = dict(zip(["+".join(map(str, entry)) for entry in rule], map(str, range(1, len(rule) + 1))))
    drule3 = dict(zip(["-" + "-".join(map(str, entry)) for entry in rule], [f"-{i}" for i in range(1, len(rule) + 1)]))
    drule = drule1 | drule2 | drule3
    return subs_dict(invariant, drule)
