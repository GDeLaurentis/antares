#!/usr/bin/env python
# -*- coding: utf-8 -*-

#    _      _ _  _           _                 _
#   /_\  __| | || |___  __  /_\  _ _  ___ __ _| |_ ______
#  / _ \/ _` | __ / _ \/ _|/ _ \| ' \(_-</ _` |  _|_ / -_)
# /_/ \_\__,_|_||_\___/\__/_/ \_\_||_/__/\__,_|\__/__\___|

# Author: Giuseppe

import sys
import itertools
import antares.core.tools as SPT

from lips.invariants import Invariants
from ..core.settings import settings
from ..scalings.pair import pair_scalings
from ..ansatze.interface import Ansatz


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


class Result_AdHocNumeratorAnsatze:

    def NumeratorAnsatze(self):
        n = self.multiplicity
        lM = self.mass_dimensions
        lPW = self.mass_dimensions
        oInvariants = Invariants(self.multiplicity, Restrict3Brackets=settings.Restrict3Brackets,
                                 Restrict4Brackets=settings.Restrict4Brackets, FurtherRestrict4Brackets=settings.FurtherRestrict4Brackets)
        invs_2, invs_3, invs_s = oInvariants.invs_2, oInvariants.invs_3, oInvariants.invs_s
        invariants_full = oInvariants.full
        pair_invs, pair_exps = self.oUnknown.pair_invs, self.oUnknown.pair_exps
        pair_friends = self.oUnknown.pair_friends
        den_invs, den_exps = self.oUnknown.den_invs, self.oUnknown.den_exps

        # look for invariants not appearing anywhere (dont put 4 brackets here)
        other_invariants = invs_2 + invs_3 + invs_s
        for i, iInv in enumerate(invs_2 + invs_3 + invs_s):
            nowhere = True
            for j, jDen in enumerate(self.dens):
                if iInv in jDen.lInvs:
                    nowhere = False
            if nowhere is False:
                other_invariants.remove(iInv)

        # loop over all denominators to try to guess the numerator
        for i, iDen in enumerate(self.dens):

            # here I check which pairs scale worse than they should
            print("In denominator number {}:".format(i + 1))
            ViciniVicini = []
            ViciniViciniPower = []
            for j, jInv in enumerate(iDen.lInvs):
                for k, kInv in enumerate(iDen.lInvs):
                    if not k > j:
                        continue
                    key = "{}&{}".format(jInv, kInv)
                    if [jInv, kInv] in pair_invs:
                        index = pair_invs.index([jInv, kInv])
                    elif [kInv, jInv] in pair_invs:
                        index = pair_invs.index([kInv, jInv])
                    else:
                        continue
                    counter = iDen.lExps[j]
                    counter += iDen.lExps[k]
                    for l, lInv in enumerate(iDen.lInvs):
                        if (lInv in pair_friends[index] and
                           lInv != jInv and lInv != kInv):
                            counter += iDen.lExps[l]
                    if counter > pair_exps[index]:
                        msg = " {}".format(key)
                        msg += " should scale as "
                        msg += "{} instead they scale as {}.".format(
                            pair_exps[index], counter)
                        print(msg)
                        ViciniViciniPower += [counter - pair_exps[index]]
                        _vicini = []
                        # both are ⟨⟩ or []
                        if ((SPT.pA2.findall(jInv) != [] or
                             SPT.pS2.findall(jInv) != []) and
                            (SPT.pA2.findall(kInv) != [] or
                             SPT.pS2.findall(kInv) != [])):
                            for l in range(1, n + 1):
                                if str(l) in key:
                                    _vicini += [l]
                            ViciniVicini += [_vicini]
                        else:
                            if SPT.pSijk.findall(jInv) != []:
                                for l in range(1, n + 1):
                                    if str(l) in jInv:
                                        _vicini += [l]
                                ViciniVicini += [_vicini]
                                _vicini = []
                            if SPT.pSijk.findall(kInv) != []:
                                for l in range(1, n + 1):
                                    if str(l) in kInv:
                                        _vicini += [l]
                                ViciniVicini += [_vicini]
                                _vicini = []
                            if SPT.p3B.findall(jInv) != []:
                                # print jInv
                                pass
                            if SPT.p3B.findall(kInv) != []:
                                # print kInv
                                pass
            ViciniVicini = map(list, list(set(map(frozenset, ViciniVicini))))
            print(ViciniVicini, ViciniViciniPower)

            # here I look for invariants in this denominator only
            invs_only_in_iDen, _ = self.invs_only_in_iDen(i)

            # First is a case by case very quick analysis,
            # then is a slow numerical double scalings analysis

            # If half absolute sum is the same as minimal mass dimension,
            # then this is a constrained bunch of ⟨⟩ and []
            if SPT.Sum_Abs(lPW[i]) / 2.0 == lM[i]:
                # This is a bunch of ⟨a|b⟩ and [c|d]
                pass
            # If the minimal mass dimension is twice the half absolute sum,
            # and there is a positive and a negative entry,
            # then this is ⟨ | | ] to some power
            elif (SPT.Sum_Abs(lPW[i]) == lM[i] and
                  SPT.Positive_Entries(lPW[i]) == 1 and
                  SPT.Negative_Entries(lPW[i]) == 1):
                if len(ViciniVicini) >= 1:
                    num_ansatz = ""
                    end = lPW[i].index(min(lPW[i])) + 1
                    start = lPW[i].index(max(lPW[i])) + 1
                    num_ansatz += "⟨{}|(".format(start)
                    for vicini in ViciniVicini:
                        if ((start in vicini or end in vicini) and
                           len(vicini) == 3):
                            for k in vicini:
                                if k != start and k != end:
                                    num_ansatz += "{}+".format(k)
                            num_ansatz = num_ansatz[:-1]
                            break
                    num_ansatz += ")|{}]".format(end)
                    self.nums[i].llInvs[0] += [num_ansatz]
                    power = max(lPW[i])
                    self.nums[i].llExps[0] += [power]
            # If the minimal mass dimension is twice the half absolute sum,
            # and there are two positive and two negative entries,
            # then this might still be ⟨ | | ] to some power
            elif (SPT.Sum_Abs(lPW[i]) == lM[i] and
                  SPT.Positive_Entries(lPW[i]) == 2 and
                  SPT.Negative_Entries(lPW[i]) == 2):
                if (len(ViciniViciniPower) == 2 and
                   ViciniViciniPower[0] in lPW[i] and
                   ViciniViciniPower[1] in lPW[i]):
                    for ij in range(2):
                        num_ansatz = ""
                        for j, jPw in enumerate(lPW[i]):
                            if jPw == ViciniViciniPower[ij]:
                                start = j + 1
                                break
                        for j, jPw in enumerate(lPW[i]):
                            if jPw == -ViciniViciniPower[ij]:
                                end = j + 1
                                break
                        middle = ""
                        for k in ViciniVicini[ij]:
                            if k != start and k != end:
                                middle += "{}+".format(k)
                        middle = middle[:-1]
                        num_ansatz = "⟨{}|({})|{}]".format(start, middle, end)
                        if len(num_ansatz) == 11:
                            self.nums[i].invs[0] += [num_ansatz]
                            self.nums[i].exps[0] += [ViciniViciniPower[ij]]
            # 1 spaa or spbb with a bunch of spab
            elif (SPT.Sum_Abs(lPW[i]) + 1 == lM[i] and
                  ((SPT.Positive_Entries(lPW[i]) == 2 and
                   SPT.Negative_Entries(lPW[i]) == 1) or
                  (SPT.Positive_Entries(lPW[i]) == 1 and
                   SPT.Negative_Entries(lPW[i]) == 2))):
                num_ansatz = ""
                end = lPW[i].index(min(lPW[i])) + 1
                start = lPW[i].index(max(lPW[i])) + 1
                num_ansatz += "⟨{}|(".format(start)
                for vicini in ViciniVicini:
                    if start in vicini or end in vicini:
                        for k in vicini:
                            if k != start and k != end:
                                num_ansatz += "{}+".format(k)
                        num_ansatz = num_ansatz[:-1]
                        break
                num_ansatz += ")|{}]".format(end)
                if len(num_ansatz) == 11:
                    self.nums[i].llInvs[0] += [num_ansatz]
                    power = max([max(lPW[i]), abs(min(lPW[i]))])
                    self.nums[i].llExps[0] += [power]
            elif (SPT.Sum_Abs(lPW[i]) * 1.5 == lM[i] and
                  ((SPT.Positive_Entries(lPW[i]) == 0 and SPT.Negative_Entries(lPW[i]) == 2) or
                   (SPT.Positive_Entries(lPW[i]) == 2 and SPT.Negative_Entries(lPW[i]) == 0)) and len(ViciniVicini) > 1):
                if SPT.Positive_Entries(lPW[i]) == 2:
                    [start, end] = [j + 1 for j, x in enumerate(lPW[i]) if x == max(lPW[i])]
                elif SPT.Negative_Entries(lPW[i]) == 2:
                    [start, end] = [j + 1 for j, x in enumerate(lPW[i]) if x == min(lPW[i])]
                if (len([entry for entry in ViciniVicini if (start in entry and len(entry) == 3)]) == 1 and
                   len([entry for entry in ViciniVicini if (end in entry and len(entry) == 3)]) == 1):
                    for vicini in ViciniVicini:
                        if len(vicini) != 3:
                            continue
                        if start in vicini:
                            middle1 = "("
                            for entry in vicini:
                                if entry == start:
                                    continue
                                middle1 += str(entry) + "+"
                            middle1 = middle1[:-1]
                            middle1 += ")"
                        if end in vicini:
                            middle2 = "("
                            for entry in vicini:
                                if entry == end:
                                    continue
                                middle2 += str(entry) + "+"
                            middle2 = middle2[:-1]
                            middle2 += ")"
                    if SPT.Positive_Entries(lPW[i]) == 0:
                        self.nums[i].llInvs[0] += ["[{}|{}|{}|{}]".format(start, middle1, middle2, end)]
                    else:
                        self.nums[i].llInvs[0] += ["⟨{}|{}|{}|{}⟩".format(start, middle1, middle2, end)]
                    power = max([max(lPW[i]), abs(min(lPW[i]))])
                    self.nums[i].llExps[0] += [power]

            elif settings.ExploreDoubleScalings is True:
                # check all double scalings of the invariants only in iDen with the ones not appearing anywhere
                _pair_invs, _pair_exps, _pair_friends = pair_scalings(self.oUnknown, invs_only_in_iDen, other_invariants, invariants_full)
                sys.stdout.write("\r Finished calculating pair scalings.",)
                sys.stdout.flush()
                __pair_invs = [entry for entry in _pair_invs]
                __pair_exps = [entry for entry in _pair_exps]
                for j, jpair in enumerate(__pair_invs):
                    if (__pair_exps[j] <= -den_exps[den_invs.index(jpair[0])] or __pair_exps[j] == "F"):
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

                # col_width = max([len(str(_pair_inv)) for _pair_inv in _pair_invs])
                # for j in range(len(_pair_invs)):
                #     print "{} ---> {} --- {}".format(str(_pair_invs[j]).ljust(col_width), _pair_exps[j], len(_pair_friends[j]))

                # try weights
                _flat_friends = [item for entry in _pair_friends for item in entry if item not in invs_only_in_iDen]
                _all_friends = list(set(_flat_friends))
                _weights = []
                for _friend in _all_friends:
                    _weights += [_flat_friends.count(_friend)]
                _switch = False
                while _switch is False:
                    _switch = True
                    for j in range(len(_all_friends) - 1):
                        if _weights[j + 1] > _weights[j]:
                            _weights[j + 1], _weights[j] = (_weights[j], _weights[j + 1])
                            (_all_friends[j + 1], _all_friends[j]) = (_all_friends[j], _all_friends[j + 1])
                            _switch = False
                col_width = max(len(friend) for friend in _all_friends)
                for j in range(len(_all_friends)):
                    print(str(_all_friends[j]).ljust(col_width) + " ---> " + str(_weights[j]))
                # get the minimal intersection set!
                num_ansatze = []
                for j in range(1, len(_all_friends)):
                    for comb in itertools.combinations(_all_friends, j):
                        if (sum([_weights[_all_friends.index(entry)] for entry in comb]) < len(_pair_invs)):
                            continue
                        elif all(list(set(comb).intersection(_pair_friend)) != [] for _pair_friend in _pair_friends):
                            num_ansatze = list(comb)
                            break
                    if num_ansatze != []:
                        break
                print(num_ansatze)
                # THIS PART IS NOT FIXED YET
                # for j, jnsatz in enumerate(num_ansatze):
                #     for k, knv in enumerate(invs_only_in_iDen):
                #         multiplier = 1
                #         while True:
                #             print "Checking scaling of {}&{} w/ multiplier {}.".format(invs_only_in_iDen[0], jnsatz, multiplier),
                #             (__pair_invs, __pair_exps, __pair_friends) = PairScalings(self.oEnv.oUnknown, [knv], [jnsatz], invariants_full, multiplier)
                #             print __pair_exps
                #             if __pair_exps[0] == "F":
                #                 multiplier -= 1
                #                 break
                #             elif __pair_exps[0] >= 0:
                #                 multiplier += 1
                #             else:
                #                 multiplier -= 1
                #                 break
                #         print "{} should be at power {}".format(
                #             jnsatz, multiplier)
                #         # if multiplier != 0 and jnsatz not in lPRN[i]:
                #         #     lPRN[i] += [jnsatz]
                #         #     lPRN_Exps[i] += [multiplier]

            # if the numerator has been partially or entirely guessed
            if self.nums[i].llInvs[0] != []:
                # Update the unknown mass dimension and phase weights
                lM = self.Mass_Dimension_List
                lPW = self.Phase_Weights_List
            else:
                # print "This numerator needs to be fitted entirely."
                pass

            # if the length of the ansatz is 1 then consider it here
            if SPT.Sum_Abs(lPW[i]) / 2.0 == lM[i]:
                if lM[i] != 0:
                    Ansatze = Ansatz([lM[i]], [lPW[i]])
                    len_ansatze = len(Ansatze)
                    if len_ansatze == 1:
                        terms = []
                        for j in range(len(Ansatze.terms())):
                            terms += [str(Ansatze.terms()[j]).split(' ')]
                        compact_terms = list(set(terms[j]))
                        compact_terms_exps = [terms[j].count(entry) for entry in compact_terms]
                        self.nums[i].invs[0] += compact_terms
                        self.nums[i].exps[0] += compact_terms_exps

            # if the numerator has been partially or entirely guessed
            if self.nums[i].llInvs[0] != []:
                print("\r The new partial result for this numerator is:                         ")
                print(SPT.Write(self.nums[i].llInvs[0], self.nums[i].llExps[0])[:-2] + " ------- ", end="")
                # Update the unknown mass dimension and phase weights
                lM = self.Mass_Dimension_List
                lPW = self.Phase_Weights_List
                print("[{}, {}]".format(lM[i], lPW[i]))
            else:
                print("\r This numerator needs to be fitted entirely.")
