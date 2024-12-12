#!/usr/bin/env python
# -*- coding: utf-8 -*-

#     _                 _          ___ _               _             _
#    /_\  _ _  ___ __ _| |_ ______| __(_)__ _ ___ _ _ | |__  __ _ __(_)___
#   / _ \| ' \(_-</ _` |  _|_ / -_) _|| / _` / -_) ' \| '_ \/ _` (_-< (_-<
#  /_/ \_\_||_/__/\__,_|\__/__\___|___|_\__, \___|_||_|_.__/\__,_/__/_/__/
#                                       |___/

# Author: Giuseppe

from lips import Particles
from lips.symmetries import inverse

from antares.core.settings import settings
from antares.core.tools import pSijk, pA2, pS2, pNB, pDijk
from lips.invariants import Invariants


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def CanonicalOrdering(ProductOfSpinorsAsString, DimOneBasis):

    ProductOfSpinorsAsList = []
    for j in range(len(ProductOfSpinorsAsString) / 5):
        ProductOfSpinorsAsList += [ProductOfSpinorsAsString[j * 5:j * 5 + 5]]

    swap = True
    while swap is True:
        swap = False
        for j in range(len(ProductOfSpinorsAsList) - 1):
            if DimOneBasis.index(ProductOfSpinorsAsList[j]) > DimOneBasis.index(ProductOfSpinorsAsList[j + 1]):
                swap = True
                ProductOfSpinorsAsList[j], ProductOfSpinorsAsList[j + 1] = ProductOfSpinorsAsList[j + 1], ProductOfSpinorsAsList[j]

    return "".join(ProductOfSpinorsAsList)


def Image(Spinor, Rule, verbose=False):
    assert type(Rule) is tuple
    if type(Spinor) is tuple:
        # given R (Rule) and S (Symmetry), the Image of S under R is given by: R.S.R-1
        permutation = Spinor[0]
        inverse_rule_permutation = inverse(Rule[0])
        new_permutation = "".join([Rule[0][int(permutation[int(inverse_rule_permutation[entry - 1]) - 1]) - 1] for entry in range(1, len(permutation) + 1)])
        new_conjugation = bool((Rule[1] + Spinor[1] + Rule[1]) % 2)
        sign = "+" if len(Spinor) == 2 else Spinor[2]
        return (new_permutation, new_conjugation, sign)
    else:
        from antares.topologies.topology import convert_invariant
        SpinorImage = convert_invariant(Spinor, Rule)

        n = len(Rule[0])
        oInvariants = Invariants(n, Restrict3Brackets=settings.Restrict3Brackets, Restrict4Brackets=settings.Restrict4Brackets, FurtherRestrict4Brackets=settings.FurtherRestrict4Brackets)

        NbrFlips = 0
        if SpinorImage not in oInvariants.full:

            if pA2.findall(SpinorImage) != [] or pS2.findall(SpinorImage) != []:
                SpinorImage = SpinorImage[::-1].replace("[", "A").replace("]", "[").replace("A", "]").replace("⟨", "A").replace("⟩", "⟨").replace("A", "⟩")
                NbrFlips += 1

            elif pSijk.findall(SpinorImage) != []:
                ijk = sorted(map(int, list(pSijk.findall(SpinorImage)[0])))
                SpinorImage = "s_" + "".join(map(str, ijk))
                if SpinorImage not in oInvariants.full:
                    oParticles = Particles(n)
                    SpinorImageAlternative = "s_" + "".join(oParticles._complementary(list(SpinorImage[2:])))
                    if len(SpinorImageAlternative) < len(SpinorImage):
                        SpinorImage = SpinorImageAlternative

            elif pNB.findall(SpinorImage) != []:
                start, middles, end = pNB.findall(SpinorImage)[0]
                middles = middles.split("|")
                middles = [sorted(middle.replace("(", "").replace(")", "").split("+")) for middle in middles]
                middles = ["(" + "+".join(middle) + ")" for middle in middles]
                SpinorImage = SpinorImage[0] + start + "|" + "|".join(middles) + "|" + end + SpinorImage[-1]
                if SpinorImage not in oInvariants.full and len(middles) == 1:
                    oParticles = Particles(n)
                    if "+" in middles[0] and "-" in middles[0]:
                        raise NotImplementedError("Mixed plus and minuses in sandwiched p-slashes not implemented.")
                    sign_operator = "+" if "+" in middles[0] else "-"
                    complementary_middle = oParticles._complementary(middles[0].replace("(", "").replace(")", "").split(sign_operator))
                    complementary_middle = [entry for entry in complementary_middle if entry not in [start, end]]
                    _SpinorImage = SpinorImage[0] + start + "|(" + sign_operator.join(complementary_middle) + ")|" + end + SpinorImage[-1]
                    if _SpinorImage in oInvariants.full:
                        NbrFlips += 1
                        SpinorImage = _SpinorImage
                if SpinorImage not in oInvariants.full:
                    SpinorImage = SpinorImage[-1].replace("]", "[").replace("⟩", "⟨") + end + "|" + "|".join(middles[::-1]) + "|" + start + SpinorImage[0].replace("[", "]").replace("⟨", "⟩")
                    if len(middles) % 2 == 0:
                        NbrFlips += 1

            elif pDijk.findall(SpinorImage) != [] and SpinorImage.count("|") == 1:
                for counter in range(10):
                    ijk = pDijk.findall(SpinorImage)[0][0]
                    SpinorImage = "Δ_{}".format(ijk[-1] + ijk[:-1])
                    if SpinorImage in oInvariants.full:
                        break
            elif SpinorImage == "Δ_26|15|43":
                SpinorImage = "Δ_15|26|34"

        if verbose and SpinorImage not in oInvariants.full:
            print("Warning: could not find image in oInvariants. {} -> {}.".format(Spinor, SpinorImage))

        if NbrFlips % 2 == 0:
            Sign = 1
        else:
            Sign = -1

        # make sure all spinors are in the correct order
        # ProductOfSpinorsAsStringImage = CanonicalOrdering(ProductOfSpinorsAsStringImage, D1Basis)

        return SpinorImage, Sign
