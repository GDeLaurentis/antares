#!/usr/bin/env python
# -*- coding: utf-8 -*-

#   _  _                    _         _ _____      _             _      _   _         _
#  | \| |_  _ _ __  ___ _ _(_)__ __ _| |_   _|__  /_\  _ _  __ _| |_  _| |_(_)__ __ _| |
#  | .` | || | '  \/ -_) '_| / _/ _` | | | |/ _ \/ _ \| ' \/ _` | | || |  _| / _/ _` | |
#  |_|\_|\_,_|_|_|_\___|_| |_\__\__,_|_| |_|\___/_/ \_\_||_\__,_|_|\_, |\__|_\__\__,_|_|
#                                                                  |__/

# Author: Giuseppe

import sys

from sympy import pprint
from copy import deepcopy

from .core.settings import settings
from .core.tools import generate_latex_and_pdf, configuration_unpacker
from .core.unknown import BHUnknown, Unknown

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def numerical_to_analytical(oUnknown):

    if not hasattr(oUnknown, 'depth'):
        oUnknown.depth = 0

    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print("{}ToAnalytical called at depth {}.".format(oUnknown.what_am_I, oUnknown.depth), end="\n\n")

    if oUnknown.is_zero:

        print("The requested quantity is identically zero.\n\n{}ToAnalytical exiting depth {}.".format(oUnknown.what_am_I, oUnknown.depth))
        print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        return "The requested quantity is identically zero."

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

    if settings.DoScalings is True:

        oUnknown.do_single_collinear_limits()

        if oUnknown.does_not_require_partial_fractioning():
            iTerms = oUnknown.fit_single_scalings_result()[0]
            if len(iTerms) >= 1:  # if the inversion was done in a collinear limit it might be necessary to remove subleading singularities.
                oNewUnknown = Unknown(oUnknown)
                oUnknown.print_partial_result()
                oNewUnknown.depth = oUnknown.depth + 1
                NewResults = numerical_to_analytical(oNewUnknown)
                if NewResults != "The requested quantity is identically zero.":
                    print("")
                    for oTerm in NewResults[0]:
                        iTerms.append(oTerm)
                    iTerms.update_variables()
                oUnknown.reset()
            else:
                raise Exception("Single denominator fit is not supposed to fail - unless the matrix is too big.")
            iTerms.check_against_unknown()
            print("{}ToAnalytical exiting depth {}.".format(oUnknown.what_am_I, oUnknown.depth))
            print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
            return [iTerms]

        oUnknown.do_double_collinear_limits()

    elif settings.DoScalings is False:

        print("\033[F", end="\r")
        oUnknown.num_invs, oUnknown.num_exps, oUnknown.den_invs, oUnknown.den_exps = None, None, None, None
        oUnknown.pair_invs, oUnknown.pair_exps, oUnknown.pair_friends = None, None, None

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

    print("Poles to be eliminated: ", end="")
    pprint(oUnknown.poles_to_be_eliminated)
    print("Spurious poles: ", end="")
    pprint(oUnknown.spurious_poles)
    print("")

    for i, invariant in enumerate(oUnknown.poles_to_be_eliminated):  # this loop is essentially used only by ProceedByGuesses (otherwise it runs only once)

        if settings.AutomaticPartialFractioning is True and settings.ProceedByGuesses is False:
            print("\rGenerating full partial fractioning.\n")
        elif settings.AutomaticPartialFractioning is True and settings.ProceedByGuesses is True:
            print("\rGenerating partial fractioning for {}.\n".format(invariant))
        else:
            print("\r\nReading partial fractioning from run file.\n")

        lTerms = oUnknown.get_partial_fractioned_terms(invariant)

        for i, iTerms in enumerate(lTerms):
            print("Considering partial fractioning #{}/{}:".format(i + 1, len(lTerms)), end="")
            iTerms.summarise()
            print("")

            if settings.UseAdHocNumeratorAnsatze is True:
                iTerms.NumeratorAnsatze()
                print("")
                iTerms.summarise()
                print("")
            if settings.PerformNumeratorFit is False:
                print("END OF PROGRAM - Without numerator fit.")
                sys.exit(0)

            iTerms.fit_numerators()

            if len(iTerms) > 0:
                oNewUnknown = Unknown(oUnknown)
                oUnknown.print_partial_result()
                oNewUnknown.depth = oUnknown.depth + 1
                NewResults = numerical_to_analytical(oNewUnknown)
                if NewResults != "The requested quantity is identically zero.":
                    print("")
                    iTerms += NewResults[0]
                    iTerms.update_variables()
                lTerms[i] = iTerms
                oUnknown.reset()
                if settings.ObtainOnlyOneSolution is True or oUnknown.depth > 1:
                    lTerms = lTerms[i:i + 1]
                    break
            else:
                lTerms[i] = None
        lTerms = list(filter(None, lTerms))
        if lTerms != []:
            break

    if lTerms == [] and settings.AutomaticPartialFractioning is True:
        print("As last resort I'll try to use a single denominator although double collinear limits imply the need for more then one.")
        settings.RefineFit = False
        settings.MaximumMatrixSizeToInvert = 4000
        iTerms = oUnknown.fit_single_scalings_result(split_pole_orders=False)[0]
        if len(iTerms) >= 1:  # if the inversion was done in a collinear limit it might be necessary to remove subleading singularities.
            oNewUnknown = Unknown(oUnknown)
            oUnknown.print_partial_result()
            oNewUnknown.depth = oUnknown.depth + 1
            NewResults = numerical_to_analytical(oNewUnknown)
            if NewResults != "The requested quantity is identically zero.":
                print("")
                for oTerm in NewResults[0]:
                    iTerms.append(oTerm)
                iTerms.update_variables()
            oUnknown.reset()
        else:
            raise Exception("Single denominator fit is not supposed to fail - unless the matrix is too big.")
        iTerms.check_against_unknown()
        print("{}ToAnalytical exiting depth {}.".format(oUnknown.what_am_I, oUnknown.depth))
        print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        return [iTerms]
        # return oUnknown.fit_single_scalings_result()

    for oTerms in lTerms:
        oTerms.check_against_unknown()

    print("{}ToAnalytical exiting depth {}.".format(oUnknown.what_am_I, oUnknown.depth))
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    return lTerms


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


if __name__ == "__main__":

    settings.read_from_file(sys.argv[1])

    oBHUnknown = BHUnknown(**configuration_unpacker(sys.argv))
    oBHUnknown.do_single_collinear_limits()
    oBHUnknown.do_double_collinear_limits()

    oUnknown = Unknown(oBHUnknown, load_partial_results=settings.UsePartialResultsIfPossible)

    if oUnknown.recursively_extract_terms() != []:
        oUnknown = Unknown(deepcopy(oUnknown))

    if oUnknown.is_zero:
        if oUnknown.recursively_extract_terms() == []:
            print("The requested quantity is identically zero.")
            generate_latex_and_pdf("\\begin{my}\n$\\begin{gathered}\\scriptscriptstyle\\text{The requested quantity is identically zero.}\\end{gathered}$\n\\end{my}\n",
                                   oUnknown.res_path)
            print("END OF PROGRAM")
            sys.exit(0)
        else:
            print("The partial result is really the full result.")
            lResults = [oUnknown.recursively_extract_terms()]
    else:
        oUnknown.depth = 0
        lResults = numerical_to_analytical(oUnknown)
        for oTerms in lResults:
            oTerms += oUnknown.recursively_extract_terms()
    latex_result = ""
    for oRes in lResults:
        oRes.compactify_symmetries()
        oRes.rearrange_and_finalise()
        latex_result += oRes.Write_LaTex()
    if len(lResults) >= 1:
        generate_latex_and_pdf(latex_result, oUnknown.res_path)

    print("END OF PROGRAM")
