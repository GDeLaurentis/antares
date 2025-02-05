#!/usr/bin/env python
# -*- coding: utf-8 -*-

#     _                 _          ___ _ _   _
#    /_\  _ _  ___ __ _| |_ ______| __(_) |_| |_ ___ _ _
#   / _ \| ' \(_-</ _` |  _|_ / -_) _|| |  _|  _/ -_) '_|
#  /_/ \_\_||_/__/\__,_|\__/__\___|_| |_|\__|\__\___|_|

# Author: Giuseppe


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def SolutionToResult(oTerms, oAnsatz, Solution):
    at_least_a_non_zero_solution = False
    lllNumInvs, lllNumExps, llNumCoefs = [], [], []
    for i in range(len(Solution)):
        if len(oTerms.lllNumInvs[i]) > 1:
            print("Error: Solution to result does not expect longer than 1 common invariants for numerators.")
        common_invs = [entry for entry in oTerms.lllNumInvs[i][0]]
        common_exps = [entry for entry in oTerms.lllNumExps[i][0]]
        lllNumInvs += [[]]
        lllNumExps += [[]]
        llNumCoefs += [[]]
        for j, coefficient in enumerate(Solution[i]):
            if not isinstance(coefficient, tuple):
                coefficient = (coefficient, 0)
            if coefficient[0] == 0 and coefficient[1] == 0:
                continue
            else:
                at_least_a_non_zero_solution = True
            if oAnsatz[i][j] == ['1']:
                # this is a single coefficient term
                lllNumInvs[i] = [common_invs]
                lllNumExps[i] = [common_exps]
                llNumCoefs[i] = [coefficient]
                continue
            # else there are multiple coefficients in the numerator
            compact_terms = list(set(oAnsatz[i][j]))
            compact_exps = [oAnsatz[i][j].count(entry) for entry in compact_terms]
            # merge compact and common
            if common_invs != []:
                for _i, common_inv in enumerate(common_invs):
                    if common_inv not in compact_terms:
                        compact_terms += [common_inv]
                        compact_exps += [common_exps[_i]]
                    else:
                        index = compact_terms.index(common_inv)
                        compact_exps[index] += common_exps[_i]
            lllNumInvs[i] += [compact_terms]
            lllNumExps[i] += [compact_exps]
            llNumCoefs[i] += [coefficient]
    if at_least_a_non_zero_solution:
        return lllNumInvs, lllNumExps, llNumCoefs
    else:
        return [[]], [[]], []


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
