# Author: Giuseppe

import sys
import multiset
import functools
import itertools

from sympy import pprint
from copy import deepcopy
from pycoretools import flatten, mapThreads, filterThreads, all_non_empty_subsets
from antares.core.settings import settings
from antares.core.tools import p3B, pOijk, pDijk
from lips.invariants import Invariants
from antares.terms.term import Term, Numerator, Denominator
from antares.terms.terms import Terms
from antares.topologies.topology import internal_symmetry
from antares.ansatze.eigenbasis import PermutationCycles
from antares.partial_fractioning.by_guesses import forced_invariants, forbidden_invariants, optional_tuples, improve_from_past_terms


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def partial_partial_fractioning(oUnknown, invariant, max_nbr_terms=1, silent=False):
    if silent is False:
        print("\rForced:                                             ")
        pprint(forced_invariants(invariant, oUnknown))
        print("Forbidden:                                            ")
        pprint(forbidden_invariants(invariant, oUnknown))
        print("Optional:                                             ")
        pprint(optional_tuples(invariant, oUnknown))
        print("")

    for nbr_terms in range(1, max_nbr_terms + 1):
        # seed the oTerms with the correct number of denominators
        lTerms = [Terms([Term(Numerator(), Denominator([invariant], [oUnknown.den_exps[invariant]])) for i in range(nbr_terms)])]
        for oTerms in lTerms:
            oTerms.multiplicity = oUnknown.multiplicity

        # add all forced invariants in all possible combinations
        for forced_invariant in forced_invariants(invariant, oUnknown):
            print("\rAdding forced invariant: {}.                                             ".format(forced_invariant), end="\r")
            sys.stdout.flush()
            add_forced_invariant_partial = functools.partial(add_forced_invariant, oUnknown, invariant, forced_invariant)
            lTerms = list(set(map(Terms, map(multiset.FrozenMultiset, flatten(mapThreads(add_forced_invariant_partial, lTerms), treat_list_subclasses_as_list=False)))))
            if lTerms == []:
                break
        # add all optional tuples in all possible combinations
        else:
            for optional_tuple in optional_tuples(invariant, oUnknown):
                print("\rUnwrapping optional tuple: {}.                                       ".format("(" + ", ".join(map(str, optional_tuple)) + ")"), end="\r")
                sys.stdout.flush()
                add_optional_tuple_partial = functools.partial(add_optional_tuple, oUnknown, invariant, optional_tuple)
                lTerms = list(set(map(Terms, map(multiset.FrozenMultiset, flatten(mapThreads(add_optional_tuple_partial, lTerms), treat_list_subclasses_as_list=False)))))

        # if any two high order poles can be splitted by a spurious pole allow it
        print("\rSplitting higher order poles by means of additional spurious poles.                                  ", end="\r")
        sys.stdout.flush()
        for oTerms in lTerms:
            if len(oTerms) == 1:
                high_order_poles = [inv for inv in oTerms[0].oDen.lInvs if inv in oUnknown.den_invs and oUnknown.den_exps[inv] > 1]
                for comb in itertools.combinations(high_order_poles, 2):
                    _optional_spurious_poles = [pole for pole in oUnknown.true_friends[comb[0], comb[1]] if pole in oUnknown.spurious_poles]
                    if len(_optional_spurious_poles) == 0:
                        continue
                    else:
                        for spurious_pole in _optional_spurious_poles:
                            _oTerms = deepcopy(oTerms) + deepcopy(oTerms)
                            # pop comb[1] from _oTerms[0]
                            _oTerms[0].oDen.lExps.pop(_oTerms[0].oDen.lInvs.index(comb[1]))
                            _oTerms[0].oDen.lInvs.pop(_oTerms[0].oDen.lInvs.index(comb[1]))
                            if spurious_pole not in _oTerms[0].oDen.lInvs:
                                _oTerms[0].oDen.lInvs.append(spurious_pole)
                                _oTerms[0].oDen.lExps.append(1)
                            # pop comb[0] from _oTerms[1]
                            _oTerms[1].oDen.lExps.pop(_oTerms[1].oDen.lInvs.index(comb[0]))
                            _oTerms[1].oDen.lInvs.pop(_oTerms[1].oDen.lInvs.index(comb[0]))
                            if spurious_pole not in _oTerms[1].oDen.lInvs:
                                _oTerms[1].oDen.lInvs.append(spurious_pole)
                                _oTerms[1].oDen.lExps.append(1)
                            _oTerms.canonical_ordering()
                            lTerms.append(_oTerms)
                            # !try trimming the optional tuples of these higher order poles!
                            _oTerms = deepcopy(oTerms) + deepcopy(oTerms) + deepcopy(oTerms)
                            # pop both comb[0] and comb[1] from _oTerms[0]
                            _oTerms[0].oDen.lExps.pop(_oTerms[0].oDen.lInvs.index(comb[0]))
                            _oTerms[0].oDen.lInvs.pop(_oTerms[0].oDen.lInvs.index(comb[0]))
                            _oTerms[0].oDen.lExps.pop(_oTerms[0].oDen.lInvs.index(comb[1]))
                            _oTerms[0].oDen.lInvs.pop(_oTerms[0].oDen.lInvs.index(comb[1]))
                            if spurious_pole not in _oTerms[0].oDen.lInvs:
                                _oTerms[0].oDen.lInvs.append(spurious_pole)
                                _oTerms[0].oDen.lExps.append(1)
                            # pop comb[1] from _oTerms[1]
                            _oTerms[1].oDen.lExps.pop(_oTerms[1].oDen.lInvs.index(comb[1]))
                            _oTerms[1].oDen.lInvs.pop(_oTerms[1].oDen.lInvs.index(comb[1]))
                            if spurious_pole not in _oTerms[1].oDen.lInvs:
                                _oTerms[1].oDen.lInvs.append(spurious_pole)
                                _oTerms[1].oDen.lExps.append(1)
                            for den_inv in _oTerms[1].oDen.lInvs:  # trimming
                                if den_inv not in oUnknown.spurious_poles and (den_inv, None) in optional_tuples(comb[0], oUnknown):
                                    del _oTerms[1].oDen.lExps[_oTerms[1].oDen.lInvs.index(den_inv)]
                                    del _oTerms[1].oDen.lInvs[_oTerms[1].oDen.lInvs.index(den_inv)]
                            # pop comb[0] from _oTerms[2]
                            _oTerms[2].oDen.lExps.pop(_oTerms[2].oDen.lInvs.index(comb[0]))
                            _oTerms[2].oDen.lInvs.pop(_oTerms[2].oDen.lInvs.index(comb[0]))
                            if spurious_pole not in _oTerms[2].oDen.lInvs:
                                _oTerms[2].oDen.lInvs.append(spurious_pole)
                                _oTerms[2].oDen.lExps.append(1)
                            for den_inv in _oTerms[2].oDen.lInvs:  # trimming
                                if den_inv not in oUnknown.spurious_poles and (den_inv, None) in optional_tuples(comb[1], oUnknown):
                                    del _oTerms[2].oDen.lExps[_oTerms[2].oDen.lInvs.index(den_inv)]
                                    del _oTerms[2].oDen.lInvs[_oTerms[2].oDen.lInvs.index(den_inv)]
                            _oTerms.canonical_ordering()
                            # keep the first term only if an invariant has been removed from both [1] and [2]
                            if all([inv in _oTerms[1].oDen.lInvs or inv in _oTerms[2].oDen.lInvs for inv in oTerms[0].oDen.lInvs]):
                                del _oTerms[0]
                            lTerms.append(_oTerms)

        # check for 0.5 exponents and add the correct Omega or Pi in numerator
        print("\rChecking for .5 scalings in double collinear limits.                                                     ", end="\r")
        sys.stdout.flush()
        oInvariants = Invariants(oUnknown.multiplicity,
                                 Restrict3Brackets=settings.Restrict3Brackets, Restrict4Brackets=settings.Restrict4Brackets, FurtherRestrict4Brackets=settings.FurtherRestrict4Brackets)
        for oTerms in lTerms:
            for oTerm in oTerms:
                # delta (or other?) .5 in denominator
                for i, iInv in enumerate(oTerm.oDen.lInvs):
                    if round(oTerm.oDen.lExps[i]) != oTerm.oDen.lExps[i]:
                        delta = [entry for entry in oInvariants.invs_D if entry in oTerm.oDen.lInvs][0]
                        if p3B.findall(invariant) == []:
                            raise Exception("Please eliminate .5 behviour first.")
                        for inv in oInvariants.invs_O:
                            if set(pDijk.findall(delta)[0]) == set(pOijk.findall(inv)[0]) and (pOijk.findall(inv)[0][0] == p3B.findall(invariant)[0][0] or
                                                                                               pOijk.findall(inv)[0][0] == p3B.findall(invariant)[0][-1]):
                                break
                        oTerm.oDen.lExps[i] = round(oTerm.oDen.lExps[i])
                        oTerm.oNum.llInvs[0] += [inv]
                        oTerm.oNum.llExps[0] += [1]
                # delta 0.5 in numerator
                for i, iInv in enumerate(oTerm.oNum.llInvs[0]):
                    if round(oTerm.oNum.llExps[0][i]) != oTerm.oNum.llExps[0][i]:
                        delta = [entry for entry in oInvariants.invs_D if entry in oTerm.oDen.lInvs or entry in oTerm.oNum.llInvs[0]][0]
                        if p3B.findall(invariant) == []:
                            raise Exception("Please eliminate .5 behviour first.")
                        for inv in oInvariants.invs_O:
                            if set(pDijk.findall(delta)[0]) == set(pOijk.findall(inv)[0]) and (pOijk.findall(inv)[0][0] == p3B.findall(invariant)[0][0] or
                                                                                               pOijk.findall(inv)[0][0] == p3B.findall(invariant)[0][-1]):
                                break
                        oTerm.oNum.llExps[0][i] = round(oTerm.oNum.llExps[0][i])
                        oTerm.oNum.llInvs[0][i] = inv

        # add known numerators to everyone and add inversion information
        print("\rAdding numerators from collinear limits.                                                                     ", end="\r")
        sys.stdout.flush()
        for i, iTerms in enumerate(lTerms):
            for oTerm in iTerms:
                oTerm.oNum.llInvs[0] += oUnknown.num_invs
                oTerm.oNum.llExps[0] += [oUnknown.num_exps[inv] for inv in oUnknown.num_invs]
            iTerms.oFittingSettings.lSmallInvs = [invariant]
            iTerms.oFittingSettings.lSmallInvsExps = [oUnknown.den_exps[invariant]]
            iTerms.oUnknown = oUnknown

        # check mass dimension and phase weights consistency
        print("\rChecking lM, lPW consistency                                                                                   ", end="\r")
        sys.stdout.flush()
        lTerms = filterThreads(lambda oTerms: oTerms.check_md_pw_consistency() is not False, lTerms)

        # add numerator guesses from previously obtained pieces of the one loop amplitude
        print("\rAdding numerators from previously obtained pieces of the one loop amplitude.                                                    ", end="\r")
        sys.stdout.flush()
        for i in range(len(lTerms)):
            iTerms = lTerms[i]
            if len(iTerms) > 1:   # atm can only do it for one term ansatz (otherwise requires considering all possible combinations...)
                continue
            for easy_box in oUnknown.easy_boxes:
                if easy_box[0].oDen == iTerms[0].oDen:   # require perfect match in denominators
                    new_oTerms = Terms([deepcopy(_Term) for _Term in iTerms])
                    new_oTerms.oUnknown = oUnknown
                    new_oTerms.oFittingSettings.lSmallInvs = [invariant]
                    new_oTerms.oFittingSettings.lSmallInvsExps = [oUnknown.den_exps[invariant]]
                    for j, jInv in enumerate(easy_box[0].oNum.llInvs[0]):
                        if easy_box[0].oNum.llExps[0][j] > 1 and jInv not in new_oTerms[0].oNum.llInvs[0]:
                            new_oTerms[0].oNum.llInvs[0].append(jInv)
                            new_oTerms[0].oNum.llExps[0].append(easy_box[0].oNum.llExps[0][j] - 1)
                    if new_oTerms.check_md_pw_consistency() is True and new_oTerms != iTerms:
                        lTerms.append(new_oTerms)

        # add numerators from old relevant terms
        if settings.AddNumeratorsFromPreviousRelevantTerms is True:
            print("\rLooking for previous relevant terms...                                                                                                   ", end="\r")
            sys.stdout.flush()
            lOldTerms = oUnknown.recursively_extract_terms()
            lRelevantOldTerms = [deepcopy(oTerm) for oTerm in lOldTerms if hasattr(oTerm, "oDen") and invariant in oTerm.oDen.lInvs and len(oTerm.oNum.lCoefs) < 10 and
                                 oUnknown.den_exps[invariant] < oTerm.oDen.lExps[oTerm.oDen.lInvs.index(invariant)]]
            if lRelevantOldTerms != []:
                print("\r                                                                                                                                        ", end="\n")
                print("Previous relevant terms:                                                                                                                  ", end="\n")
                # keep only common invariants from past terms
                for i, oTerm in enumerate(lRelevantOldTerms):
                    if len(oTerm.oNum.llInvs) > 1:
                        if oTerm.oNum.lCommonInvs == []:
                            lRelevantOldTerms[i] = None
                        else:
                            oTerm.oNum.llInvs = [deepcopy(oTerm.oNum.lCommonInvs)]
                            oTerm.oNum.llExps = [deepcopy(oTerm.oNum.lCommonExps)]
                            oTerm.oNum.lCommonInvs = []
                            oTerm.oNum.lCommonExps = []
                            oTerm.oNum.lCoefs = []
                lRelevantOldTerms = filterThreads(lambda x: x is not None, lRelevantOldTerms)
                for relevant_old_term in lRelevantOldTerms:
                    print(relevant_old_term)
            if len(lRelevantOldTerms) > 0:
                print("\rAdding numerators from previous relevant terms.                                                                                    ", end="\r")
                for i in range(len(lTerms)):
                    new_oTerms = improve_from_past_terms(Terms([deepcopy(_Term) for _Term in lTerms[i]]), invariant, lRelevantOldTerms)
                    if new_oTerms is not None:
                        new_oTerms.oUnknown = oUnknown
                        new_oTerms.oFittingSettings.lSmallInvs = [invariant]
                        new_oTerms.oFittingSettings.lSmallInvsExps = [oUnknown.den_exps[invariant]]
                        if new_oTerms.check_md_pw_consistency() is True and new_oTerms != iTerms:
                            lTerms.append(new_oTerms)

        # a^2 b^2 c^2 in c-tower becomes a^2 b c^2, a b^2 c^2
        if settings.SplitHigherOrderPoles is True:
            print("\rSplitting higher poles...                                                                                                                          ", end="\r")
            sys.stdout.flush()
            for oTerms in lTerms:
                if len(oTerms) < 1:
                    continue
                invs_to_be_split = [iInv for i, iInv in enumerate(oTerms[0].oDen.lInvs) if oTerms[0].oDen.lExps[i] > 1 and iInv != invariant]
                if len(oTerms) == 1 and len(invs_to_be_split) > 1:
                    for inv_to_be_split in invs_to_be_split:
                        oTerms += [deepcopy(oTerms[0])]
                        oTerms[-1].oDen.lExps = [exp if exp == 1 or inv == inv_to_be_split or inv == invariant else exp - 1 for exp, inv in zip(oTerms[-1].oDen.lExps, oTerms[-1].oDen.lInvs)]
                    del oTerms[0]

        # sort + canonical ordering + print
        print("\rSorting results...                                                                                                                      ", end="\r")
        sys.stdout.flush()
        lTerms = list(set(lTerms))
        lTerms = filterThreads(lambda _oTerms: (_oTerms is not None and _oTerms != [] and all([entry in oUnknown.den_invs + oUnknown.spurious_poles
                                                                                               for _oTerm in _oTerms for entry in _oTerm.oDen.lInvs])), lTerms)
        lTerms.sort(key=lambda oTerms: (
            # priority to spurious poles when higher poles are involved
            max([sum([1 if entry in oUnknown.spurious_poles else (
                max([oUnknown.den_exps[entry], max([0] + [oTerm.oDen.lExps[oTerm.oDen.lInvs.index(entry)] for oTerm in oUnknown.recursively_extract_terms()
                                                          if oTerm.am_I_a_symmetry is False and entry in oTerm.oDen.lInvs])]) - 1) for entry in oTerm.oDen.lInvs]) for oTerm in oTerms]),
            # by worst mass dimension
            -min(oTerms.mass_dimensions),
            # by overall mass dimension
            -sum(oTerms.mass_dimensions))
        )

        # limit the number of guesses to the best
        lTerms = lTerms[:]

        # make sure they are internally ordered
        for oTerms in lTerms:
            oTerms.canonical_ordering()

        # symmetries handling
        if ((hasattr(oTerms, "oUnknown") and hasattr(oTerms.oUnknown.recursively_extract_original_unknown(), "amppart") and
           oTerms.oUnknown.recursively_extract_original_unknown().amppart not in [None, "external"]) or
           (hasattr(oUnknown.recursively_extract_original_unknown(), "symmetries"))):
            specify_symmetries(invariant, lTerms, oUnknown)

        # print them
        for oTerms in lTerms:
            print(oTerms, end="\n")

        # try with one more denominator or exit
        if lTerms != []:
            print("\rFound at least a result.                                                       ")
            break
        else:
            print("\rNo result found with {} terms.                                                 ".format(nbr_terms))

    return lTerms


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def add_forced_invariant(oUnknown, invariant, forced_invariant, oTerms):
    lTerms = []
    for indices in all_non_empty_subsets(range(len(oTerms))):
        _oTerms = deepcopy(oTerms)
        for index in indices:
            # check that forced_invariant is not forbidden with other invariants already in denominator
            if any([forced_invariant in forbidden_invariants(_inv, oUnknown, exponent=1) for _inv in _oTerms[index].oDen.lInvs if
                    _inv in oUnknown.den_invs and forced_invariant not in oUnknown.spurious_poles and _inv not in oUnknown.spurious_poles]):
                break
            # append the invariant with an appropriate exponent
            else:
                _oTerms[index].oDen.lInvs.append(forced_invariant)
                _oTerms[index].oDen.lExps.append(
                    # forced_invariant exponent is the difference between that of the pair and that of the invariant
                    oUnknown.pair_exps[(invariant, forced_invariant)] - oUnknown.den_exps[invariant] if
                    # check that the subtraction makes sense
                    type(oUnknown.pair_exps[(invariant, forced_invariant)]) in [float, int] and
                    # if forced_invariant is in the denominator
                    forced_invariant in oUnknown.den_invs and
                    # and if it does not exceed the exponent of forced_invariant by itself
                    oUnknown.pair_exps[(invariant, forced_invariant)] - oUnknown.den_exps[invariant] <= oUnknown.den_exps[forced_invariant] and
                    # and if it is at least 0.5 (generally positive)
                    oUnknown.pair_exps[(invariant, forced_invariant)] - oUnknown.den_exps[invariant] >= 0.5 else
                    # otherwise forced_invariant exponent is 1
                    1)
        else:
            _oTerms.canonical_ordering()
            lTerms += [_oTerms]
    return lTerms


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def add_optional_tuple(oUnknown, invariant, optional_tuple, oTerms):
    lTerms = []
    for optional_invariant in optional_tuple:
        if optional_invariant is None:
            lTerms += [deepcopy(oTerms)]
            continue
        for indices in all_non_empty_subsets(range(len(oTerms))):
            _oTerms = deepcopy(oTerms)
            _oTerms2 = deepcopy(oTerms)  # Naively some invariants may want to go in the numerator, but in that case allow them to be on the bottom as well (other friends could be in num)
            for index in indices:
                # if this invariant has already been added than skip it
                if optional_invariant in oTerms[index].oDen.lInvs or optional_invariant in oTerms[index].oNum.llInvs[0]:
                    continue
                # if it is forbidden with any other invariant already added than skip it (unless one of them a spurious pole)
                if any([optional_invariant in forbidden_invariants(_inv, oUnknown, exponent=1) for _inv in oTerms[index].oDen.lInvs if
                        _inv in oUnknown.den_invs and optional_invariant not in oUnknown.spurious_poles and _inv not in oUnknown.spurious_poles]):
                    break
                # if the naive exponenet is zero allow it in the denominator with power of 1
                if (optional_invariant in oUnknown.den_invs and oUnknown.pair_exps[(invariant, optional_invariant)] != "F" and
                   oUnknown.pair_exps[(invariant, optional_invariant)] - oUnknown.den_exps[invariant] == 0):
                    _oTerms2[index].oDen.lInvs.append(optional_invariant)
                    _oTerms2[index].oDen.lExps.append(1)
                    _oTerms2[index].canonical_ordering()
                # if the naive exponent is less than zero then it mean it goes in the numerator.. consider both cases in numerator with naive exponent or in denominator
                elif (optional_invariant in oUnknown.den_invs and oUnknown.pair_exps[(invariant, optional_invariant)] != "F" and
                      oUnknown.pair_exps[(invariant, optional_invariant)] - oUnknown.den_exps[invariant] < 0):
                    _oTerms[index].oNum.llInvs[0].append(optional_invariant)
                    _oTerms[index].oNum.llExps[0].append(abs(oUnknown.pair_exps[(invariant, optional_invariant)] - oUnknown.den_exps[invariant]))
                    _oTerms[index].canonical_ordering()
                    _oTerms2[index].oDen.lInvs.append(optional_invariant)
                    _oTerms2[index].oDen.lExps.append(1)
                    _oTerms2[index].canonical_ordering()
                else:
                    _oTerms[index].oDen.lInvs.append(optional_invariant)
                    _oTerms[index].oDen.lExps.append(oUnknown.pair_exps[(invariant, optional_invariant)] - oUnknown.den_exps[invariant] if
                                                     optional_invariant in oUnknown.den_invs and oUnknown.pair_exps[(invariant, optional_invariant)] != "F" and
                                                     oUnknown.pair_exps[(invariant, optional_invariant)] - oUnknown.den_exps[invariant] <= oUnknown.den_exps[optional_invariant] else
                                                     oUnknown.den_exps[optional_invariant] if optional_invariant in oUnknown.den_invs else 1)
                    _oTerms[index].canonical_ordering()
                    _oTerms2[index].oDen.lInvs.append(optional_invariant)
                    _oTerms2[index].oDen.lExps.append(oUnknown.pair_exps[(invariant, optional_invariant)] - oUnknown.den_exps[invariant] if
                                                      optional_invariant in oUnknown.den_invs and oUnknown.pair_exps[(invariant, optional_invariant)] != "F" and
                                                      oUnknown.pair_exps[(invariant, optional_invariant)] - oUnknown.den_exps[invariant] <= oUnknown.den_exps[optional_invariant] else
                                                      oUnknown.den_exps[optional_invariant] if optional_invariant in oUnknown.den_invs else 1)
                    _oTerms2[index].canonical_ordering()
            else:
                if _oTerms == _oTerms2:
                    lTerms += [_oTerms2]
                elif _oTerms != _oTerms2:
                    lTerms += [_oTerms, _oTerms2]
    return lTerms


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def specify_symmetries(invariant, lTerms, oUnknown):
    from antares.terms.terms import Term
    if hasattr(oUnknown.recursively_extract_original_unknown(), "symmetries"):
        symmetries = oUnknown.recursively_extract_original_unknown().symmetries
    else:
        symmetries = internal_symmetry(oUnknown.recursively_extract_original_unknown())
    pprint("\rThere are {} symmetries:                                                         ".format(len(symmetries)))
    pprint(symmetries)
    print("")
    for symmetry in symmetries:
        cycles = PermutationCycles(symmetry, oUnknown.den_invs, ordering=False)
        for cycle in cycles:
            if invariant not in cycle:
                continue
            for oTerms in lTerms:
                oImmagedTerms = oTerms.Image(symmetry)
                # These symmetries are respected by the terms
                if oImmagedTerms == oTerms:
                    oTerms.symmetries += [symmetry]
                # These symmetries will be appended to the result in the end
                elif len(cycle) > 1 and not any([oTerms.Image(sym) == oImmagedTerms for sym in oTerms.useful_symmetries]):
                    oTerms.useful_symmetries += [symmetry]
                # These symmetries will be used in the inversion
                elif len(cycle) == 1 and not any([oTerms.Image(sym) == oImmagedTerms for sym in oTerms.oFittingSettings.lSymmetries]):
                    oTerms.append(Term(symmetry))
