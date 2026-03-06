
#   ___          _   _      _ ___            _   _          _           ___       ___
#  | _ \__ _ _ _| |_(_)__ _| | __| _ __ _ __| |_(_)___ _ _ (_)_ _  __ _| _ )_  _ / __|_  _ ___ ______ ___ ___
#  |  _/ _` | '_|  _| / _` | | _| '_/ _` / _|  _| / _ \ ' \| | ' \/ _` | _ \ || | (_ | || / -_|_-<_-</ -_|_-<
#  |_| \__,_|_|  \__|_\__,_|_|_||_| \__,_\__|\__|_\___/_||_|_|_||_\__, |___/\_, |\___|\_,_\___/__/__/\___/__/
#                                                                 |___/     |__/

# Author: Giuseppe

from sympy import pprint
from copy import deepcopy
from pycoretools import flatten
from antares.core.settings import settings
from lips.invariants import Invariants
from antares.terms.terms import Terms, FittingSettings
from antares.topologies.topology import internal_symmetry, get_label
from antares.ansatze.eigenbasis import PermutationCycles
from numbers import Number


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def forced_invariants(invariant, oUnknown):
    return list(set(flatten([
        # loop over all invarants in denomiantor except the invariant itself. If conditions are met save the true friends cleaned of the invariant itself
        list(filter(lambda x: x != invariant, oUnknown.true_friends[(invariant, other_invariant)])) for other_invariant in filter(lambda x: x != invariant, oUnknown.den_invs) if
        # criteron 1: sum of all true friends exps adds up to the pair_exp from double collinear limits (eg: 3 invariants all appearing with power 1 and pair_exp = 3)
        sum([oUnknown.den_exps[friend_invariant] if friend_invariant in oUnknown.den_invs else 1
             for friend_invariant in oUnknown.true_friends[(invariant, other_invariant)]]) == oUnknown.pair_exps[(invariant, other_invariant)] or
        # critereon 2: no additional true friends and pair_exp exceeds the exp from single collinear limits (eg: 2 invariants both with power 2 but pair_exp is 3)
        isinstance(oUnknown.pair_exps[(invariant, other_invariant)], Number) and
        (oUnknown.pair_exps[(invariant, other_invariant)] > oUnknown.den_exps[invariant] and oUnknown.pair_exps[(invariant, other_invariant)] > oUnknown.den_exps[other_invariant] and
         len(oUnknown.true_friends[(invariant, other_invariant)]) == 2)
    ])))
    # old forced without knowing spurious poles
    # return [inv1 if inv1 != invariant else inv2 for i, (inv1, inv2) in enumerate(oUnknown.pair_invs)
    #         if invariant in [inv1, inv2] and oUnknown.pair_exps[i] > oUnknown.den_exps[oUnknown.den_invs.index(inv1)] and
    #         oUnknown.pair_exps[i] > oUnknown.den_exps[oUnknown.den_invs.index(inv2)] and
    #         len(oUnknown.pair_friends[i]) == 2]


def optional_tuples(invariant, oUnknown):
    higher_order_poles = [inv for inv in oUnknown.den_invs if oUnknown.den_exps[inv] > 1]
    optional_invs = [tuple(filter(lambda x: x != invariant, sorted(oUnknown.true_friends[(invariant, other_invariant)], key=lambda y: y != other_invariant))) if
                     isinstance(oUnknown.pair_exps[(invariant, other_invariant)], Number) and
                     # pair_exps exceeds exp from single collinear limits
                     oUnknown.pair_exps[(invariant, other_invariant)] > oUnknown.den_exps[invariant] and
                     # no higher order pole in true friends (e.g. a, b have pair_exp 4 but c is a true friend and has exp 4, than c is enough to explain it without a or b)
                     not any([inv in higher_order_poles and oUnknown.den_exps[inv] >= oUnknown.pair_exps[(invariant, other_invariant)]
                              for inv in [_inv for _inv in oUnknown.true_friends[(invariant, other_invariant)] if _inv not in (invariant, other_invariant)]]) and
                     # not critereon 1 in forced invs
                     not sum([oUnknown.den_exps[friend_invariant] if friend_invariant in oUnknown.den_invs else 1
                              for friend_invariant in oUnknown.true_friends[(invariant, other_invariant)]]) == oUnknown.pair_exps[(invariant, other_invariant)] and
                     # not critereon 2 in forced invs
                     not len(oUnknown.true_friends[(invariant, other_invariant)]) == 2 and
                     # not failed
                     tuple(filter(lambda x: x != invariant, oUnknown.true_friends[(invariant, other_invariant)])) != ("F",) else
                     # allow it if failed
                     (other_invariant, ) if tuple(filter(lambda x: x != invariant, oUnknown.true_friends[(invariant, other_invariant)])) == ("F",)
                     else None for other_invariant in oUnknown.den_invs if other_invariant != invariant]
    optional_invs = list(filter(lambda x: x is not None, optional_invs))
    # if doesn't appear at the start of a tuple then allow it here
    optional_invs += [(inv, None) for inv in oUnknown.den_invs if inv not in [invariant] + forced_invariants(invariant, oUnknown) +
                      forbidden_invariants(invariant, oUnknown) + [optional_tuple[0] for optional_tuple in optional_invs]]
    # allow spurious poles between forced and forbidden invariants
    for forced_inv in forced_invariants(invariant, oUnknown):
        for forbidden_inv in forbidden_invariants(invariant, oUnknown):
            for spurious_pole in oUnknown.spurious_poles:
                if spurious_pole in oUnknown.true_friends[(forced_inv, forbidden_inv)] and not any([spurious_pole in optional_tuple for optional_tuple in optional_invs]):
                    optional_invs += [(spurious_pole, None)]
    # clean up + sort
    oInvariants = Invariants(oUnknown.multiplicity,
                             Restrict3Brackets=settings.Restrict3Brackets, Restrict4Brackets=settings.Restrict4Brackets, FurtherRestrict4Brackets=settings.FurtherRestrict4Brackets)
    optional_invs = [tuple(sorted(list(optional_tuple), key=lambda inv: (oInvariants.full.index(inv) if inv in oInvariants.full else 9999))) for optional_tuple in optional_invs]
    optional_invs = [tuple([inv for inv in optional_tuple if inv not in forbidden_invariants(invariant, oUnknown)]) for optional_tuple in optional_invs]
    optional_invs = [optional_tuple if len(optional_tuple) >= 2 else optional_tuple + (None,) for optional_tuple in optional_invs]
    return list(set(optional_invs))


def suggested_invariants(invariant, oUnknown):
    return [inv1 if inv1 != invariant else inv2 for i, (inv1, inv2) in enumerate(oUnknown.pair_invs)
            if invariant in [inv1, inv2] and oUnknown.pair_exps[i] > oUnknown.den_exps[oUnknown.den_invs.index(inv1)] and
            oUnknown.pair_exps[i] > oUnknown.den_exps[oUnknown.den_invs.index(inv2)] and
            len(oUnknown.pair_friends[i]) > 2]


def forbidden_invariants(invariant, oUnknown, exponent="Maximal"):
    if exponent == "Maximal":
        exponent = oUnknown.den_exps[invariant]
    return [other_invariant for other_invariant in oUnknown.den_invs if other_invariant != invariant and
            isinstance(oUnknown.pair_exps[(invariant, other_invariant)], Number) and
            oUnknown.pair_exps[(invariant, other_invariant)] <= exponent and
            len(oUnknown.pair_friends[(invariant, other_invariant)]) == 2 and other_invariant not in oUnknown.spurious_poles]


def discouraged_invariants(invariant, oUnknown):
    return [inv1 if inv1 != invariant else inv2 for i, (inv1, inv2) in enumerate(oUnknown.pair_invs)
            if invariant in [inv1, inv2] and oUnknown.pair_exps[i] <= [oUnknown.den_exps[oUnknown.den_invs.index(inv1)]
                                                                       if inv1 == invariant else oUnknown.den_exps[oUnknown.den_invs.index(inv2)]][0] and
            len(oUnknown.pair_friends[i]) > 2]


def allowed_basis_functions_with_clean_up(invariant, oUnknown):
    to_be_returned = [basis_function for basis_function in oUnknown.basis_functions
                      if all([inv in basis_function for inv in forced_invariants(invariant, oUnknown)])]
    if to_be_returned == []:
        to_be_returned = deepcopy(oUnknown.basis_functions)
    # clean up from forbidden
    to_be_returned = [[inv for inv in basis_function if inv not in forbidden_invariants(invariant, oUnknown) and (inv in oUnknown.den_invs or inv in oUnknown.spurious_poles)]
                      for basis_function in to_be_returned]
    if invariant in oUnknown.basis_functions_invs:
        to_be_returned = [cleaned_basis_function for cleaned_basis_function in to_be_returned if invariant in cleaned_basis_function]
    return sorted(to_be_returned, key=lambda entry: 0 -
                  sum([1 for inv in entry if inv in suggested_invariants(invariant, oUnknown)]) +
                  sum([1 for inv in entry if inv in discouraged_invariants(invariant, oUnknown)]) +
                  sum([1 for inv in entry if inv not in oUnknown.den_invs])
                  )


def dress(oTerms, invariant, oUnknown):
    # remove discouraged
    oTerms[0].oDen.lInvs = [inv for inv in oTerms[0].oDen.lInvs if inv not in discouraged_invariants(invariant, oUnknown)]
    # add suggested
    oTerms[0].oDen.lInvs = oTerms[0].oDen.lInvs + [inv for inv in suggested_invariants(invariant, oUnknown) if inv not in oTerms[0].oDen.lInvs and
                                                   all([inv2 not in oTerms[0].oDen.lInvs for inv2 in forbidden_invariants(inv, oUnknown)])]
    # if some suggested could not be added due to forbidden with other invariants already present generate a second denominator
    if any(True for inv in suggested_invariants(invariant, oUnknown) if inv not in oTerms[0].oDen.lInvs):
        missing_suggested = [inv for inv in suggested_invariants(invariant, oUnknown) + forced_invariants(invariant, oUnknown) if inv not in oTerms[0].oDen.lInvs]
        print("missing:", end="")
        pprint(missing_suggested)
        oTerms += [deepcopy(oTerms[0])]
        oTerms[1].oDen.lInvs = filter(lambda x: all([x not in forbidden_invariants(inv, oUnknown) for inv in missing_suggested]), oTerms[1].oDen.lInvs)
        oTerms[1].oDen.lInvs += missing_suggested
    # generate first ansatze without spurious poles that haven't appeared yet

    for oTerm in oTerms:
        if invariant not in oTerm.oDen.lInvs:
            oTerm.oDen.lInvs += [invariant]
        # add exponents
        oTerm.oDen.lExps = [oUnknown.pair_exps[oUnknown.pair_invs.index([inv, invariant])] - oUnknown.den_exps[oUnknown.den_invs.index(invariant)]
                            if [inv, invariant] in oUnknown.pair_invs and oUnknown.pair_exps[oUnknown.pair_invs.index([inv, invariant])] != "F" and
                            oUnknown.pair_exps[oUnknown.pair_invs.index([inv, invariant])] - oUnknown.den_exps[oUnknown.den_invs.index(invariant)] <
                            oUnknown.den_exps[oUnknown.den_invs.index(inv)] else oUnknown.den_exps[oUnknown.den_invs.index(inv)] if [inv, invariant] in oUnknown.pair_invs else
                            oUnknown.pair_exps[oUnknown.pair_invs.index([invariant, inv])] - oUnknown.den_exps[oUnknown.den_invs.index(invariant)]
                            if [invariant, inv] in oUnknown.pair_invs and oUnknown.pair_exps[oUnknown.pair_invs.index([invariant, inv])] != "F" and
                            oUnknown.pair_exps[oUnknown.pair_invs.index([invariant, inv])] - oUnknown.den_exps[oUnknown.den_invs.index(invariant)] <
                            oUnknown.den_exps[oUnknown.den_invs.index(inv)] else oUnknown.den_exps[oUnknown.den_invs.index(inv)] if [invariant, inv] in oUnknown.pair_invs else
                            oUnknown.den_exps[oUnknown.den_invs.index(invariant)]
                            if inv == invariant else 1
                            for inv in oTerm.oDen.lInvs]
        # if Delta appears with x.5 power, then take ceiling and add a Omega on top
        if any(round(exp) != exp for exp in oTerm.oDen.lExps):
            oTerm.oDen.lExps = map(round, oTerm.oDen.lExps)
        oTerm.canonical_ordering()
    # add inversion info
    oTerms.oUnknown = oUnknown
    oFittingSettings = FittingSettings()
    oFittingSettings.lSmallInvs, oFittingSettings.lSmallInvsExps = [invariant], [oUnknown.den_exps[oUnknown.den_invs.index(invariant)]]
    oTerms.oFittingSettings = oFittingSettings


def improve_from_past_terms(oTerms, invariant, lRelevantOldTerms):
    invariant_exps_in_old_terms = list(set([relevant_old_term.oDen.lExps[relevant_old_term.oDen.lInvs.index(invariant)] for relevant_old_term in lRelevantOldTerms]))
    if len(invariant_exps_in_old_terms) == 1 and len(lRelevantOldTerms) == len(oTerms):
        # add numerator if the structures match
        for i, oTerm in enumerate(oTerms):
            for j, oRelTerm in enumerate(lRelevantOldTerms):
                if all([inv in oRelTerm.oDen.lInvs and oRelTerm.oDen.lExps[oRelTerm.oDen.lInvs.index(inv)] in range(
                        int(oTerm.oDen.lExps[oTerm.oDen.lInvs.index(inv)]) - 1, int(oTerm.oDen.lExps[oTerm.oDen.lInvs.index(inv)]) + 1 + 1) or
                        0 in range(int(oTerm.oDen.lExps[oTerm.oDen.lInvs.index(inv)]) - 1, int(oTerm.oDen.lExps[oTerm.oDen.lInvs.index(inv)]) + 1 + 1)
                        for inv in oTerm.oDen.lInvs]):
                    # add numerator guess
                    numerator_guesses = [(inv, exp - 1) for inv, exp in zip(lRelevantOldTerms[j].oNum.llInvs[0], lRelevantOldTerms[j].oNum.llExps[0]) if exp > 1]
                    if numerator_guesses != []:
                        oTerms[i].oNum.llInvs[0], oTerms[i].oNum.llExps[0] = map(list, zip(*numerator_guesses))
        # # split denominator --- although it is done elsewhere, it could be beneficial to allow it here as well? maybe
        # _oTerms = deepcopy(oTerms)
        # for i, iTerm in enumerate(oTerms):
        #     invs_to_be_split = [jInv for j, jInv in enumerate(iTerm.oDen.lInvs) if iTerm.oDen.lExps[j] > 1 and jInv != invariant]
        #     for inv_to_be_split in invs_to_be_split:
        #         _oTerms = Terms([deepcopy(iTerm)] + list(_oTerms))
        #         _oTerms[0].oDen.lExps = [exp - 1 if inv == inv_to_be_split else exp for exp, inv in zip(_oTerms[0].oDen.lExps, _oTerms[0].oDen.lInvs)]
        # oTerms = _oTerms
    elif len(invariant_exps_in_old_terms) > 1:
        invariant_powers = sorted(list(set([oTerm.oDen.lExps[i] for oTerm in lRelevantOldTerms for i, inv in enumerate(oTerm.oDen.lInvs) if inv == invariant])), key=lambda x: -x)
        llRelevantOldTermsByInvariantPowers = [[oTerm for oTerm in lRelevantOldTerms if oTerm.oDen.lExps[oTerm.oDen.lInvs.index(invariant)] == invariant_power]
                                               for invariant_power in invariant_powers]
        steps = []
        if len(llRelevantOldTermsByInvariantPowers) < 2:
            return
        for relevant_old_term in llRelevantOldTermsByInvariantPowers[1]:
            steps += [relevant_old_term / llRelevantOldTermsByInvariantPowers[0][0]]
        # print "Steps:"
        # for step in steps:
        #     print step
        # print ""
        oTerms = oTerms[0:0]
        for relevant_old_term in llRelevantOldTermsByInvariantPowers[1]:
            for step in steps:
                oTerms += Terms([relevant_old_term * step])
                oTerms[-1].canonical_ordering()
        oTerms = oTerms.remove_duplicates()
    return oTerms


def make_guesses(oUnknown, invariant):
    print("\rThere are {} basis functions.".format(len(oUnknown.basis_functions)))
    for basis_function in oUnknown.basis_functions:
        pprint(basis_function)
    print("")

    forced_invs = forced_invariants(invariant, oUnknown)
    suggested_invs = suggested_invariants(invariant, oUnknown)
    forbidden_invs = forbidden_invariants(invariant, oUnknown)
    discouraged_invs = discouraged_invariants(invariant, oUnknown)

    print("Forced: ")
    pprint(forced_invs)
    print("Suggested: ")
    pprint(suggested_invs)
    print("Forbidden: ")
    pprint(forbidden_invs)
    print("Discouraged: ")
    pprint(discouraged_invs)
    print("")

    lOldTerms = oUnknown.recursively_extract_terms()
    print("Old relevant terms:")
    lRelevantOldTerms = [oTerm for oTerm in lOldTerms if hasattr(oTerm, "oDen") and invariant in oTerm.oDen.lInvs and len(oTerm.oNum.lCoefs) == 1 and
                         oUnknown.den_exps[oUnknown.den_invs.index(invariant)] < oTerm.oDen.lExps[oTerm.oDen.lInvs.index(invariant)] and
                         not any([True for _oTerm in lOldTerms if hasattr(_oTerm, "oDen") and invariant in _oTerm.oDen.lInvs and len(_oTerm.oNum.lCoefs) > 1 and
                                  oUnknown.den_exps[oUnknown.den_invs.index(invariant)] < oTerm.oDen.lExps[_oTerm.oDen.lInvs.index(invariant)]])]
    for relevant_old_term in lRelevantOldTerms:
        print(relevant_old_term)
    print("")

    if len(forced_invs) + len(suggested_invs) + len(forbidden_invs) + len(discouraged_invs) != len(oUnknown.den_invs) - 1 and lRelevantOldTerms == []:
        print("Warning: some invariants were not found in any of the four lists.")
        pprint([inv for inv in oUnknown.den_invs if inv not in forced_invs + suggested_invs + forbidden_invs + discouraged_invs and inv != invariant])

    allowed_basis_fncs = allowed_basis_functions_with_clean_up(invariant, oUnknown)
    print("Allowed basis functions with clean up: ")
    for allowed_basis_fnc in allowed_basis_fncs:
        pprint(allowed_basis_fnc)
    print("")

    pprint([[inv for inv in allowed_basis_fnc if inv not in oUnknown.spurious_poles] for allowed_basis_fnc in allowed_basis_fncs])
    pprint([allowed_basis_fnc for allowed_basis_fnc in allowed_basis_fncs])
    lTerms = [Terms([[[]]], [[[]]], [[]], [[inv for inv in allowed_basis_fnc if inv not in oUnknown.spurious_poles]], [[]]) for allowed_basis_fnc in allowed_basis_fncs if
              any([True for inv in allowed_basis_fnc if inv not in oUnknown.spurious_poles])]
    lTerms += [Terms([[[]]], [[[]]], [[]], [allowed_basis_fnc], [[]]) for allowed_basis_fnc in allowed_basis_fncs]
    for oTerms in lTerms:
        dress(oTerms, invariant, oUnknown)

    lTerms.sort(key=lambda oTerms: -sum(oTerms.mass_dimensions))

    print("Dressed terms:")
    for oTerms in lTerms:
        print(oTerms)

    if lRelevantOldTerms != []:
        for i, iTerms in enumerate(lTerms):
            lTerms[i] = improve_from_past_terms(oTerms, invariant, lRelevantOldTerms)

    symmetries = internal_symmetry(get_label(oUnknown.helconf, oUnknown.amppart)) if oUnknown.amppart in ["tree", "rational"] else internal_symmetry(
        get_label(oUnknown.helconf, oUnknown.amppart, oUnknown.ampindex))
    pprint("There are {} symmetries:".format(len(symmetries)))
    pprint(symmetries)
    for symmetry in symmetries:
        cycles = PermutationCycles(symmetry, oUnknown.den_invs, ordering=False)
        for cycle in cycles:
            if invariant in cycle and len(cycle) > 1:
                for oTerms in lTerms:
                    oTerms.useful_symmetries += [symmetry]
            elif invariant in cycle and len(cycle) == 1:
                for oTerms in lTerms:
                    oTerms.oFittingSettings.lSymmetries += [symmetry]
    print("")

    return lTerms
