#!/usr/bin/env python
# -*- coding: utf-8 -*-

#     _                 _          ___     _            __
#    /_\  _ _  ___ __ _| |_ ______|_ _|_ _| |_ ___ _ _ / _|__ _ __ ___
#   / _ \| ' \(_-</ _` |  _|_ / -_)| || ' \  _/ -_) '_|  _/ _` / _/ -_)
#  /_/ \_\_||_/__/\__,_|\__/__\___|___|_||_\__\___|_| |_| \__,_\__\___|

# Author: Giuseppe
# Created: 10/07/2018

from math import factorial
from copy import deepcopy

from pycoretools import flatten
from ..core.settings import settings


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


class Ansatz(list):
    """Ansatz for numerator spinorial expression: lllInvariants. Obtained through Google's CP SAT + Singular."""  # or Daniel's Spinor Solve

    @staticmethod
    def minimal_mass_dimension(phase_weights):
        half_abs_sum = sum([abs(pw) for pw in phase_weights]) / 2.0
        total_positive = sum([pw for pw in phase_weights if pw > 0])
        total_negative = sum([-pw for pw in phase_weights if pw < 0])
        highest_positive = max(phase_weights) if max(phase_weights) > 0 else 0
        lowest_negative = min(phase_weights) if min(phase_weights) < 0 else 0
        abs_lowest_negative = abs(lowest_negative)
        minimum_mass_dimension = half_abs_sum
        if highest_positive > (total_positive - highest_positive):
            unbalanced_positive = (2 * highest_positive - total_positive)
        else:
            unbalanced_positive = total_positive % 2
        if abs_lowest_negative > (abs_lowest_negative - total_negative):
            unbalanced_negative = (2 * abs_lowest_negative - total_negative)
        else:
            unbalanced_negative = total_negative % 2
        minimum_mass_dimension += max([unbalanced_positive, unbalanced_negative])
        return minimum_mass_dimension

    @classmethod
    def exists(cls, mass_dimension, phase_weights):
        """Checks if an anstaz can be built, without actually trying to build it."""
        half_abs_sum = sum([abs(pw) for pw in phase_weights]) / 2.0

        # odd number of contractions means it is impossible
        if not half_abs_sum.is_integer():
            return False

        # naive mass dimension check
        if not mass_dimension >= half_abs_sum:
            return False

        # impossible combination
        if (sum(phase_weights) / 2) % 2 != mass_dimension % 2:
            return False

        # advanced mass dimension check
        minimum_mass_dimension = cls.minimal_mass_dimension(phase_weights)
        if not mass_dimension >= minimum_mass_dimension:
            return False

        return True

    @staticmethod
    def estimate_size(mass_dimension, phase_weights):
        """See Eq.s 3.1 to 3.3 of arXiv:2010.14525"""
        minimal_mass_dimension = Ansatz.minimal_mass_dimension(phase_weights)
        m = len(phase_weights)
        d = int(mass_dimension - minimal_mass_dimension)
        return sum([CwR(n_s(m), i // 2) if i == d else n_t(m) if i == 0 else CwR(n_s(m), i // 2) * n_t(m) for i in range(d, -1, -4)[0:2]])

    # @staticmethod
    # def from_spinor_solve(md, pw):
    #     ansatz_result = ansatz(md, pw)
    #     ansatz_terms = ansatz_result.terms()
    #     ansatz_terms = list(map(str, ansatz_terms))
    #     ansatz_terms = [ansatz_term.replace("<", "⟨").replace(">", "⟩").replace("'", "").split(" ") for ansatz_term in ansatz_terms]
    #     return ansatz_terms

    @staticmethod
    def from_cp_sat_higgs(md, pw, verbose=True):
        from antares.ansatze.Numerator_Ansatz_Generator_CP_SAT_Higgs import Ansatz as ansatze_higgs
        return ansatze_higgs([md], [pw], verbose=verbose)[0]

    @staticmethod
    def from_cp_sat_HHj(md, pw, verbose=True):
        from antares.ansatze.Numerator_Ansatz_Generator_CP_SAT_2Higgses import Ansatz as ansatze_HHj
        return ansatze_HHj([md], [pw], verbose=verbose)[0]

    @staticmethod
    def from_cp_sat_HHH(md, pw, verbose=True):
        from antares.ansatze.Numerator_Ansatz_Generator_CP_SAT_3Higgses import Ansatz as ansatze_HHH
        return ansatze_HHH([md], [pw], verbose=verbose)[0]

    @staticmethod
    def from_cp_sat_HHHj(md, pw, verbose=True):
        from antares.ansatze.Numerator_Ansatz_Generator_CP_SAT_3Hj import Ansatz as ansatze_HHHj
        return ansatze_HHHj([md], [pw], verbose=verbose)[0]

    @staticmethod
    def from_cp_sat_ttH(md, pw, verbose=True):
        from antares.ansatze.Numerator_Ansatz_Generator_CP_SAT_ttH import Ansatz as ansatze_ttH
        return ansatze_ttH([md], [pw], verbose=verbose)[0]

    @staticmethod
    def from_cp_sat_w_boson(md, pw, verbose=True):
        from antares.ansatze.Numerator_Ansatz_Generator_CP_SAT_W_Boson import AnsatzFlipped as ansatze_w_boson
        return ansatze_w_boson([md], [pw], verbose=verbose)[0]

    @staticmethod
    def from_cp_sat(md, pw, verbose=True):
        from antares.ansatze.Numerator_Ansatz_Generator_CP_SAT import Ansatz as ansatze
        return ansatze([md], [pw], verbose=verbose)[0]

    def get_ansatz(self, md, pw, verbose=True):
        # if settings.NumeratorAnsatzProvider == "SPINOR_SOLVE":
        #     return self.from_spinor_solve(md, pw)
        if settings.NumeratorAnsatzProvider == "CP_SAT":
            return self.from_cp_sat(md, pw, verbose=verbose)
        elif settings.NumeratorAnsatzProvider == "CP_SAT_HIGGS":
            return self.from_cp_sat_higgs(md, pw, verbose=verbose)
        elif settings.NumeratorAnsatzProvider == "CP_SAT_HHJ":
            return self.from_cp_sat_HHj(md, pw, verbose=verbose)
        elif settings.NumeratorAnsatzProvider == "CP_SAT_HHH":
            return self.from_cp_sat_HHH(md, pw, verbose=verbose)
        elif settings.NumeratorAnsatzProvider == "CP_SAT_3Hj":
            return self.from_cp_sat_HHHj(md, pw, verbose=verbose)
        elif settings.NumeratorAnsatzProvider == "CP_SAT_ttH":
            return self.from_cp_sat_ttH(md, pw, verbose=verbose)
        elif settings.NumeratorAnsatzProvider == "CP_SAT_W_BOSON":
            return self.from_cp_sat_w_boson(md, pw, verbose=verbose)
        else:
            raise Exception("Ansatz provider not understood.")

    def __init__(self, lM, lPW, MaximumAnsatzeSize=10 ** 6, verbose=True):
        """Initialisation. Requires list of mass dimensions (lM), list of phase weights (lPW)."""
        if verbose:
            print(f"Obtaining ansatz from {settings.NumeratorAnsatzProvider} with lM, lPW: {lM}, {lPW}.", end=" ")

        if False in map(self.exists, lM, lPW):
            self._nInput = 0
            if verbose:
                print("\nImpossible ansatze! {}           ".format(list(map(self.exists, lM, lPW))))
            # return
        elif sum(map(self.estimate_size, lM, lPW)) > MaximumAnsatzeSize:
            self._nInput = sum(map(self.estimate_size, lM, lPW))
            if verbose:
                print("\nEstimated ansatz size is {}. It exceeds the set maximum                                       ".format(self._nInput))
            return

        list.__init__(self)
        ansatz_terms = [(self.get_ansatz(md, pw, verbose=verbose) if md != 0 else [[]])
                        if self.exists(md, pw) else [] for md, pw in zip(lM, lPW)]
        for ansatz_term in ansatz_terms:
            self.append(ansatz_term)

        self._set_nInput()
        self._set_lInvariants()
        self._set_as_indices()

        if verbose:
            print(f"\rObtained ansatz from {settings.NumeratorAnsatzProvider} with lM, lPW: {lM}, {lPW}. Size: {self.nInput}.",
                  end="                                                       \n")

    def _set_nInput(self):
        self._nInput = len(flatten(self, max_recursion=1))

    def _set_lInvariants(self):
        self._lInvariants = list(set(flatten(self)) if '1' in flatten(self) else set(flatten(self) + ['1']))
        self._lInvariants.sort()

    def _set_as_indices(self):
        self._as_indices = [[[self.lInvariants.index(entry) for entry in term] for term in a_term] for a_term in self]

    @property
    def nInput(self):
        """Length of flattened Ansatz object."""
        return self._nInput

    @nInput.setter
    def nInput(self, other):
        raise Exception("Illegal access: nInput in Ansatz object is read only.")

    @property
    def lInvariants(self):
        """List of invariants appearing in Ansatz object."""
        return self._lInvariants

    @lInvariants.setter
    def lInvariants(self, other):
        raise Exception("Illegal access: lInvariants in Ansatz object is read only.")

    @property
    def as_indices(self):
        """as_indices on sInvariants representation Ansatz object."""
        return self._as_indices

    @as_indices.setter
    def as_indices(self, other):
        raise Exception("Illegal access: as_indices in Ansatz object is read only.")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def SymmetriseAnsatze(oAnsatz, oTerms):
    from antares.topologies.topology import convert_invariant

    oSymmetrisedAnsatz = deepcopy(oAnsatz)
    del oSymmetrisedAnsatz[:]
    for i, iTerm in enumerate(oTerms):
        if iTerm.is_symmetry is True:
            oSymmetrisedAnsatz.append(iTerm.tSym)
        else:
            oSymmetrisedAnsatz.append(oAnsatz[i])

    oSymmetrisedAnsatzExplicit = deepcopy(oAnsatz)
    del oSymmetrisedAnsatzExplicit[:]
    for i, iAnsatz in enumerate(oSymmetrisedAnsatz):
        if type(iAnsatz) is tuple:
            something_was_not_a_symmetry = False
            to_be_appended = []
            for j, jAnsatz in enumerate(oSymmetrisedAnsatz[:i][::-1]):
                if type(jAnsatz) is tuple and something_was_not_a_symmetry:
                    break
                elif type(jAnsatz) is not tuple:
                    something_was_not_a_symmetry = True
                    to_be_appended.append([[convert_invariant(inv, iAnsatz) for inv in term] for term in jAnsatz])
            oSymmetrisedAnsatz += to_be_appended[::-1]
        else:
            oSymmetrisedAnsatzExplicit.append(iAnsatz)

    oSymmetrisedAnsatzExplicit._set_nInput()
    oSymmetrisedAnsatzExplicit._set_lInvariants()
    oSymmetrisedAnsatzExplicit._set_as_indices()

    return oSymmetrisedAnsatzExplicit


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def C(n, r):
    """Combinations."""
    if n >= 0 and r >= 0 and (n - r) >= 0:
        return factorial(n) // factorial(r) // factorial(n - r)
    else:
        return


def CwR(n, r):
    """Combinations with replacement."""
    return C(r + n - 1, r)


def n_s(m):
    """Number of independents Mandelstams."""
    return m * (m - 3) // 2


def n_t(m):
    """Number of trace_5's."""
    if m > 4:
        return C(m - 1, 4)
    else:
        return 0


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# Symmetric Ansatze Example   ~~~~~   for future work

# # Symmetrised ansatze for 6, [-3, 0, 0, 0, 0, 3] w/ sym: (654321, True)
# invs = {"⟨2|6⟩", "⟨3|6⟩", "⟨4|6⟩", "⟨5|6⟩", "[1|2]", "[1|5]", "[1|4]", "[1|3]"}
# terms = ["⟨2|6⟩⟨2|6⟩⟨5|6⟩[1|2][1|2][1|5]+⟨2|6⟩⟨5|6⟩⟨5|6⟩[1|2][1|5][1|5]", "⟨2|6⟩⟨2|6⟩⟨2|6⟩[1|2][1|2][1|2]+⟨5|6⟩⟨5|6⟩⟨5|6⟩[1|5][1|5][1|5]",
#          "⟨2|6⟩⟨2|6⟩⟨4|6⟩[1|2][1|2][1|4]+⟨3|6⟩⟨5|6⟩⟨5|6⟩[1|3][1|5][1|5]", "⟨2|6⟩⟨4|6⟩⟨4|6⟩[1|2][1|4][1|4]+⟨3|6⟩⟨3|6⟩⟨5|6⟩[1|3][1|3][1|5]"]
# nInput = 4


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
