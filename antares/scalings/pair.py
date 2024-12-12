#!/usr/bin/env python
# -*- coding: utf-8 -*-

#  ___      _     ___          _ _
# | _ \__ _(_)_ _/ __| __ __ _| (_)_ _  __ _ ___
# |  _/ _` | | '_\__ \/ _/ _` | | | ' \/ _` (_-<
# |_| \__,_|_|_| |___/\__\__,_|_|_|_||_\__, /__/
#                                      |___/

# Author: Giuseppe

from __future__ import unicode_literals
from __future__ import print_function

import sys

from mpmath.libmp.libhyper import NoConvergence
from lips import Particles, myException

from antares.core.settings import settings
from antares.core.tools import mapThreads, log_linear_fit, MyShelf, retry, log_linear_fit_Exception


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def pair_scalings(oUnknown, some_invs, other_invs, all_invariants, relative=1, seed=0, silent=True):

    # This concept of pair 'Friends' is actually that of ideal members.
    pw_friends = settings.base_cache_path + "n={}".format(oUnknown.multiplicity) + "/Friends"
    with MyShelf(pw_friends, "c"):
        pass   # just to make sure it exists

    with MyShelf(pw_friends, 'r') as _sFriends:
        sFriends = dict(_sFriends)

    if hasattr(oUnknown, "new"):
        raise myException("new encountered in oUnknown in set_pair... deprecated?")
        del oUnknown.new  # new is not picklable

    # build the tuples list
    invs_tuples = []
    for i, some_inv in enumerate(some_invs):
        if some_invs == other_invs:
            j_start = i + 1
        else:
            j_start = 0
        for j in range(j_start, len(other_invs)):
            other_inv = other_invs[j]
            invs_tuples += [(some_inv, other_inv)]

    data = mapThreads(pair_scaling, oUnknown, all_invariants, relative, sFriends, invs_tuples, seed=seed, UseParallelisation=settings.UseParallelisation, Cores=settings.Cores)

    # un pack the data and write it to cache if it is not there yet
    pair_invs, pair_exps, pair_friends = [], [], []
    with MyShelf(pw_friends, 'c') as sFriends:
        for i, entry in enumerate(data):
            if entry is None:
                pair_invs += [list(invs_tuples[i])]
                pair_exps += ["F"]
                pair_friends += [["F"]]
            else:
                pair_invs += [entry[0]]
                pair_exps += [entry[1]]
                pair_friends += [entry[2]]
                key1 = "{}&{}".format(entry[0][0], entry[0][1])
                key2 = "{}&{}".format(entry[0][1], entry[0][0])
                if sys.version_info.major < 2:
                    key1 = key1.encode("utf-8")
                    key2 = key2.encode("utf-8")
                if key1 not in sFriends:
                    sFriends[key1] = entry[2]
                    sFriends[key2] = entry[2]

    # Clean using known numerator information
    if oUnknown.num_invs != [] and oUnknown.num_invs is not None:
        clean_pair_scalings_from_numerator(oUnknown, pair_invs, pair_exps, pair_friends)

    return pair_invs, pair_exps, pair_friends


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


@retry((myException, log_linear_fit_Exception, AssertionError, TimeoutError, AttributeError, NoConvergence), max_tries=2, silent=False)
def pair_scaling(oUnknown, all_invariants, relative, sFriends, invs_tuple, seed=0):

    some_inv, other_inv = invs_tuple[0], invs_tuple[1]

    oParticles = Particles(oUnknown.multiplicity, seed=seed, field=settings.field)

    if settings.field.name == "mpc":

        xaxis, yaxis = [], []

        for k in range(settings.ScalingsIterationsStart, settings.ScalingsIterationsStart + settings.ScalingsIterationsNumber):
            oParticles.variety((some_inv, other_inv), (10 ** -(relative * k), 2 * 10 ** -k))
            xaxis += [abs(oParticles(some_inv))]
            yaxis += [abs(oUnknown(oParticles))]

        pair_exp = log_linear_fit(xaxis, yaxis)

    elif settings.field.name == "padic":

        oParticles.variety((some_inv, other_inv), (1, 1,))
        pair_exp = oUnknown(oParticles).n

    pair_invs = [some_inv, other_inv]
    key1 = "{}&{}".format(some_inv, other_inv)
    if sys.version_info.major < 2:
        key1 = key1.encode("utf-8")

    if key1 not in sFriends:                        # time consuming part, since it involves recalculating all invariants
        _pair_friends = oParticles.phasespace_consistency_check(all_invariants)[3]
    else:                                           # read it from the cache
        _pair_friends = sFriends[key1]
    pair_friends = _pair_friends

    return pair_invs, pair_exp, pair_friends


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def clean_pair_scalings_from_numerator(oUnknown, pair_invs, pair_exps, pair_friends):
    print("\rCleaning pair scalings results from known numerator information.                                                             ")
    for i, pair_inv in enumerate(pair_invs):
        for num_inv in oUnknown.num_invs:
            if num_inv in pair_friends[i]:
                pair_exps[i] -= oUnknown.num_exps[num_inv]
