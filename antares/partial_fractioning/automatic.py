#     _       _                  _   _    ___          _   _      _ ___            _   _          _
#    /_\ _  _| |_ ___ _ __  __ _| |_(_)__| _ \__ _ _ _| |_(_)__ _| | __| _ __ _ __| |_(_)___ _ _ (_)_ _  __ _
#   / _ \ || |  _/ _ \ '  \/ _` |  _| / _|  _/ _` | '_|  _| / _` | | _| '_/ _` / _|  _| / _ \ ' \| | ' \/ _` |
#  /_/ \_\_,_|\__\___/_|_|_\__,_|\__|_\__|_| \__,_|_|  \__|_\__,_|_|_||_| \__,_\__|\__|_\___/_||_|_|_||_\__, |
#                                                                                                       |___/

# Author: Giuseppe

import sys
import itertools
import logging

from antares.core.my_warnings import oWarning as w


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def AutomaticPartialFraction(oUnknown, starting_terms=1):
    Dens_Result = []

    for m_counter in range(starting_terms, 10):
        # seed the List with the correct number of denominators
        DenList = [[]]
        for i in range(m_counter - 1):
            DenList += [[]]
        # Recursive partial-fractioning routine
        yield from Evolve(0, DenList, Dens_Result, oUnknown)
        # Choose if to continue with 1 more denominator or to stop
        if Dens_Result != []:
            print("\r                                                                         ")
            print("Found at least a result for m = {}.\n".format(m_counter))
            return
        else:
            print("\r                                                                         ", end="")
            print("\rNo result found for m = {}.".format(m_counter))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def Evolve(i, DenList, Dens_Result, oUnknown):

    if i == len(oUnknown.den_invs):    # end of recursion condition
        yield from EndOfRecursion(i, DenList, Dens_Result, oUnknown)
        return

    yield from AddAnotherInvariant(i, DenList, Dens_Result, oUnknown)   # recursive bit


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def AddAnotherInvariant(i, DenList, Dens_Result, oUnknown):
    inv = oUnknown.den_invs[i]
    _DenListSet = []

    # print(i, den_invs[i], DenList)
    sys.stdout.write("\rDepth: {}/{}, adding invariant: {}.                       ".format(i, len(oUnknown.den_invs), oUnknown.den_invs[i]))
    sys.stdout.flush()

    # loop over all denominators in the given list and decide which ones are tagged, forced or forbidden
    # the invariant inv must be inserted in a forced denominator, must not be inserted in a forbidden denominator and may be insterted in a tagged denominator
    tagged, forced, forbidden = [], [], []
    reason = [[]]
    for __i in range(len(DenList) - 1):
        reason += [[]]
    for k, Den in enumerate(DenList):
        tagged += [k]
        # loop over all entries of that denominator to choose if it need to be evolved (recall: enumerate([]) doens't loop even once)
        for l, l_inv in enumerate(Den):
            # if the pair scales as a single one then check if they have any "friends"
            if oUnknown.pair_exps[(inv, l_inv)] == max([oUnknown.den_exps[inv], oUnknown.den_exps[l_inv]]):
                # if they don't have any friends then this pair is forbidden (from being in the same denominator)
                if len(oUnknown.pair_friends[(inv, l_inv)]) == 2:
                    if k not in forbidden:
                        forbidden += [k]
            # if the pair scales as the product of the two then check if they have any "friends"
            elif oUnknown.pair_exps[(inv, l_inv)] == (oUnknown.den_exps[inv] + oUnknown.den_exps[l_inv]):
                # if they don't have any friends then this is forced
                if len(oUnknown.pair_friends[(inv, l_inv)]) == 2:
                    if k not in forced:
                        forced += [k]
                    # save the reason why the two are forced
                    if l_inv not in reason[k]:
                        reason[k] += [l_inv]
        # if it is both forbidden and forced then this set of denominators is no good, don't evolve it at all
        if (k in forbidden and k in tagged):
            tagged.remove(k)
        elif (k in forced and k in tagged):
            tagged.remove(k)
    # print(forbidden, forced, tagged)

    # evolve it according to the decision made
    if forced != [] or tagged != []:
        # some forced...
        if forced != []:
            for iforced in forced:
                # local variable for checking reasons
                _loopf = [_entry for _entry in forced]
                _loopf.remove(iforced)
                # refresh evolution variable
                _DenList = [denominator[:] for denominator in DenList]
                if iforced in forbidden:
                    # bad starting denominator to evolve
                    continue
                else:
                    _DenList[iforced] += [inv]
                # new default for permutations
                _DenListPermutations = [denominator[:] for denominator
                                        in _DenList]
                # check if others are still forced (different reasons)
                # need to consider permutations here ...
                for permu in itertools.permutations(_loopf, len(_loopf)):
                    _DenList = [denominator[:] for denominator
                                in _DenListPermutations]
                    _further_reason = []
                    _forced = [_entry for _entry in forced]
                    _forced.remove(iforced)
                    _ireason = [_entry for _entry in reason[iforced]]
                    bad_permutation = False
                    for jforced in permu:
                        _jreason = [_entry for _entry in reason[jforced]]
                        _ireason += _further_reason
                        for _ielem in _ireason:
                            if _ielem in _jreason:
                                _jreason.remove(_ielem)
                        if _jreason != []:
                            if jforced in forbidden:
                                bad_permutation = True
                                break
                            else:
                                _DenList[jforced] += [inv]
                            _further_reason += _jreason
                            _forced.remove(jforced)
                    # print(_forced, forbidden)
                    if bad_permutation is True:
                        continue
                    # forced one & the ones without an overlapping reason
                    if _DenListSet == []:
                        yield from Evolve(i + 1, _DenList, Dens_Result, oUnknown)
                        _DenListSet = [[_Den[:] for _Den in _DenList]]
                    elif all(sorted(_DenList) != sorted(__DenList)
                             for __DenList in _DenListSet):
                        yield from Evolve(i + 1, _DenList, Dens_Result, oUnknown)
                        _DenListSet += [_DenList]
                    # new default List of Denominators for this set
                    _DenListForced = [denominator[:]
                                      for denominator in _DenList]
                    # add the overlapping reason forced to tagged
                    # once stripped of forbidden ones
                    __forced = [_entry for _entry in _forced]
                    for _entry in __forced:
                        if _entry in forbidden:
                            _forced.remove(_entry)
                    _tagged = tagged + _forced
                    # then all combinations of the tagged ones
                    for rlength in range(1, len(_tagged) + 1):
                        for comb in itertools.combinations(_tagged, rlength):
                            # refresh it to the forced one
                            _DenList = [denominator[:]
                                        for denominator in _DenListForced]
                            for comb_item in comb:
                                _DenList[comb_item] += [inv]
                            # ... do it ...
                            if _DenListSet == []:
                                yield from Evolve(i + 1, _DenList, Dens_Result, oUnknown)
                                _DenListSet = [[_Den[:] for _Den in _DenList]]
                            # check that this is not overcounting
                            elif all(sorted(_DenList) != sorted(__DenList)
                                     for __DenList in _DenListSet):
                                yield from Evolve(i + 1, _DenList, Dens_Result, oUnknown)
                                _DenListSet += [_DenList]
        # no forced...
        else:
            # all combinations of the tagged ones
            for rlength in range(1, len(tagged) + 1):
                for comb in itertools.combinations(tagged, rlength):
                    # refresh evolution variable
                    _DenList = [denominator[:] for denominator in DenList]
                    for comb_item in comb:
                        _DenList[comb_item] += [inv]
                    # ... do it ...
                    if _DenListSet == []:
                        yield from Evolve(i + 1, _DenList, Dens_Result, oUnknown)
                        _DenListSet = [[_Den[:] for _Den in _DenList]]
                    # check that this is not overcounting
                    elif all(sorted(_DenList) != sorted(__DenList)
                             for __DenList in _DenListSet):
                        yield from Evolve(i + 1, _DenList, Dens_Result, oUnknown)
                        _DenListSet += [_DenList]


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def EndOfRecursion(i, DenList, Dens_Result, oUnknown):
    from antares.terms.terms import Terms

    sys.stdout.write("\rDepth: {}/{}, looking for spurious poles.         ".format(i, len(oUnknown.den_invs)))
    sys.stdout.flush()

    if [] in DenList:                                                       # empty list in denominator might mean we should have found a result by now
        msg = "Empty list in ansatze for denominator!"
        w.warn(msg, msg)
        return

    # if oUnknown.amppart != 'rational' or False:
    logging.debug("\nLooking for spurious spoles in:")                      # spurious poles search
    logging.debug(DenList)
    DenList = Add_Spurious_Poles(DenList, oUnknown)

    if DenList == []:                                                       # if not empty list then this might be a solution
        return

    for _DenList in Dens_Result:                                            # additional check to avoid overcounting results:
        is_subset = True                                                    # if all denominators are a subset of a previously found result
        for j in range(len(DenList)):                                       # then this would lead to cancellation between num and den
            bigSet = set(DenList[j])
            smallSet = set(_DenList[j])
            if not smallSet.issubset(bigSet):
                is_subset = False
                break
        if is_subset is True:
            return

    sys.stdout.write("\rDepth: {}/{}, loading result and checking mass dimensions.".format(i, len(oUnknown.den_invs)))
    sys.stdout.flush()

    DensExps = []                                                           # RECONSTRUCTION OF DENOMINATOR EXPONENTS
    for j, jDen in enumerate(DenList):                                      # loop over all partial fractioned denominators
        DenExps = []
        for k, kInv in enumerate(jDen):                                     # loop over all invariants in that denominator
            if kInv in oUnknown.den_invs:                                       # --- real poles --- they could have their maximum power (from single scalings)
                free_powers = 100                                           # or a lower power due to constraints from double scalings
                for l, lInv in enumerate(jDen):
                    if l == k:                                              # look at terms which have already been inserted only
                        break
                    if (lInv, kInv) not in oUnknown.pair_exps.keys():
                        continue
                    if len(oUnknown.pair_friends[(lInv, kInv)]) != 2:
                        continue
                    combined_exponent = oUnknown.pair_exps[(lInv, kInv)]               # power from double scalings
                    if combined_exponent == 'F':
                        continue
                    lInvInsertedExp = DenExps[l]                            # inserted power of lInv
                    _free_powers = round(combined_exponent - lInvInsertedExp)  # according to this pair, kInvPower cannot exceed combined_exponent-lInvInsertedExp
                    if _free_powers < 0:
                        _free_powers = 0
                    if _free_powers < free_powers:
                        free_powers = _free_powers
                kInvMaxExp = oUnknown.den_exps[kInv]
                DenExps += [min([kInvMaxExp, free_powers])]                 # insert minimum between free_power and max_power
            else:                                                           # --- spurious poles --- they should only have power 1
                DenExps += [1]
        DensExps += [DenExps]

    lPRN_Invs = [[[entry for entry in oUnknown.num_invs]] for den in DenList]   # load the oRes object with the current partial result
    lPRN_Exps = [[[oUnknown.num_exps[entry] for entry in oUnknown.num_invs]] for den in DenList]
    lPRN_Coefs = [[1] for den in DenList]
    lPRD_Invs = DenList
    lPRD_Exps = DensExps
    oRes = Terms(lPRN_Invs, lPRN_Exps, lPRN_Coefs, lPRD_Invs, lPRD_Exps)
    oRes.oUnknown = oUnknown

    if oRes.check_md_pw_consistency() is False:                             # check the mass dimension and phase weights consistency
        return
    sys.stdout.write("\r                                                                   ")
    sys.stdout.flush()
    print("\nFound possible partial fractioning (#{}),".format(len(Dens_Result) + 1), end="")
    print("with combined mass dimension of {}:".format(sum(oRes.ansatze_mass_dimensions)))
    Dens_Result += [DenList]                                                # <--- needs changing?

    oRes.order_terms_by_complexity()                                        # order partial fractioned terms from most to less complicated
    oRes.summarise()                                                        # print a summary of the current partial result

    yield oRes


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def Add_Spurious_Poles(_DenList, oUnknown):
    # print("")
    DenList = [den[:] for den in _DenList]
    linked = []
    for i, iDen in enumerate(DenList):
        for j, jDen in enumerate(DenList):
            if (not j > i):
                continue
            linked += [False]
    switch = True
    warning = ""
    # try again if a spurious pole is found
    while switch is True:
        switch = False
        counter = 0
        # loop over all combinations of denominators in the each list
        for i, iDen in enumerate(DenList):
            for j, jDen in enumerate(DenList):
                # make sure no overcounting happens
                if (not j > i):
                    continue
                logging.debug("Checking link {}-{}".format(i, j))
                if linked[counter] is True:
                    logging.debug("Skipped because already linked.")
                    counter += 1
                    continue
                # build the fake/linked denominators
                counter2 = 0
                linked_dens = []
                for ii, iiDen in enumerate(DenList):
                    for jj, jjDen in enumerate(DenList):
                        if (not jj > ii):
                            continue
                        if linked[counter2] is True:
                            linked_dens += [iiDen + jjDen]
                        counter2 += 1
                spurious_pole = []
                # loop over all combinations of terms in the denominators
                for k, kInv in enumerate(iDen):
                    for l, lInv in enumerate(jDen):
                        if (kInv, lInv) not in oUnknown.pair_exps.keys():
                            continue
                        need_spurious_pole = False
                        if oUnknown.pair_exps[(kInv, lInv)] == 2:  # why not? if exp == exp1+exp2:
                            need_spurious_pole = True
                        # if the pair already appears somewhere then not
                        for r, rDen in enumerate(DenList + linked_dens):
                            if (kInv in rDen and lInv in rDen):
                                need_spurious_pole = False
                                break
                        for friend in oUnknown.pair_friends[(kInv, lInv)]:
                            if (friend in iDen and friend in jDen):
                                need_spurious_pole = False
                                break
                        # add these friends to the spurious pole list
                        if need_spurious_pole is True:
                            logging.debug("{} and {} should scale as 2".format(kInv, lInv))
                            # logging the pair_friends leads to massive logfile
                            # logging.debug(pair_friends[index])
                            if spurious_pole == []:
                                _spurious_pole = oUnknown.pair_friends[(kInv, lInv)]
                                if kInv in _spurious_pole:
                                    _spurious_pole.remove(kInv)
                                if lInv in _spurious_pole:
                                    _spurious_pole.remove(lInv)
                                spurious_pole = _spurious_pole
                            else:
                                # look for intersection
                                _spurious_pole = [sp for sp in spurious_pole if sp in oUnknown.pair_friends[(kInv, lInv)]]
                                if kInv in _spurious_pole:
                                    _spurious_pole.remove(kInv)
                                if lInv in _spurious_pole:
                                    _spurious_pole.remove(lInv)
                                spurious_pole = _spurious_pole
                                # if get back to empty list then break
                                if spurious_pole == []:
                                    logging.debug("No spurious pole found")
                                    warning = "No spurious pole"
                                    break
                    if warning == "No spurious pole":
                        break
                if warning == "No spurious pole":
                    pass
                elif len(spurious_pole) > 1:
                    logging.debug("Error: found multiple spurious poles.")
                    pass
                elif linked[counter] is False and spurious_pole != []:
                    linked[counter] = True
                    switch = True
                    iDen += spurious_pole
                    jDen += spurious_pole
                    logging.debug(spurious_pole)
                else:
                    logging.debug("No spurious pole needed for this link")
                    linked[counter] = True
                warning = ""
                counter += 1
    if False in linked:
        # remove this denominator list since I coundn't
        # reconstruct the spurious poles
        logging.debug("No valid solution has been found for this set.")
        logging.debug("Now removing it.")
        return []
    else:
        # print("")
        logging.debug("Possible valid solution has been found for this set.")
        logging.debug(DenList)
        return DenList


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
