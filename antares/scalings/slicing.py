import functools
import numpy
import sympy

from copy import copy
from multiset import Multiset
from pyadic import ModP, PAdic
from pyadic.interpolation import Newton_polynomial_interpolation, Thiele_rational_interpolation
from syngular import flatten

from ..terms.terms import Terms, Term, Numerator, Denominator


@functools.lru_cache(maxsize=1024)
def get_invariant_dict(possible_denominators, evaluator, field=None):
    if field is None:
        field = evaluator.field
    candidate_denoms_dict = {}
    for possible_denominator in possible_denominators:
        possible_denominator_of_t = evaluator(possible_denominator)
        if isinstance(possible_denominator_of_t, (ModP, PAdic)):
            candidate_denoms_dict[possible_denominator] = possible_denominator_of_t
        elif isinstance(possible_denominator_of_t, numpy.ndarray):
            pass
        else:
            candidate_denoms_dict[possible_denominator] = sympy.poly(
                4 * possible_denominator_of_t.expand(), modulus=field.characteristic).as_expr().factor(modulus=field.characteristic)
    if not set(flatten([list(val.as_powers_dict().values()) for key, val in candidate_denoms_dict.items()])) == {1}:
        raise NotImplementedError("Expecting only factors to power 1 on a generic slice.")
    # check that no factors are shared
    candidate_denom_factors_dict = {key: [entry for entry in val.as_ordered_factors() if 't' in str(entry)]
                                    for key, val in candidate_denoms_dict.items()}
    candidate_unique_factors = [key for key, val in Multiset(flatten([val for key, val in candidate_denom_factors_dict.items()])).items() if val == 1]
    unique_candidate_denom_factors_dict = {key: [entry for entry in val if entry in candidate_unique_factors] for key, val in candidate_denom_factors_dict.items()}
    unique_candidate_denom_factors_dict = {key: val for key, val in unique_candidate_denom_factors_dict.items() if val != []}
    non_unique_candidate_denom_factors_dict = {key: val for key, val in candidate_denom_factors_dict.items()
                                               if key not in unique_candidate_denom_factors_dict.keys()}
    if not non_unique_candidate_denom_factors_dict == {}:
        print("""Warning: non unique candidate denominator factors found.
              Are you sure they all generate distinct ideals at codimension 1?""", non_unique_candidate_denom_factors_dict)
    return candidate_denom_factors_dict


def match_factors(polynomial, candidate_factors_dict, field, assert_full_match=False, assert_factors=True, verbose=False):
    if isinstance(polynomial, sympy.Number):
        return {}
    polynomial_degree = polynomial.as_poly(modulus=field.characteristic).degree()
    # ensure the polynomial is normalized to have leading coefficient of 1
    t = sympy.symbols('t')
    polynomial = sympy.poly(polynomial * int(1 / ModP(int(polynomial.coeff(t ** polynomial_degree)), field.characteristic)),
                            modulus=field.characteristic).as_expr()
    matched_irreds, matched_degree = {}, 0
    polynomial_power_dict = polynomial.factor(modulus=field.characteristic).as_powers_dict()
    unmatched_polynomial_power_dict = copy(polynomial_power_dict)
    for candidate, all_factors in candidate_factors_dict.items():
        if any([factor in polynomial_power_dict.keys() for factor in all_factors]):  # check if any matches
            all_factors_appear = all([factor in polynomial_power_dict.keys() for factor in all_factors])
            # check that all factors actually do appear
            if assert_factors:
                assert all_factors_appear
            else:
                if not all_factors_appear:
                    print("Warning: would have failed assertion, some factors are missing:",
                          [(factor, factor in polynomial_power_dict.keys()) for factor in all_factors])
            # and check that they appear with the same power
            all_factors_powers = [polynomial_power_dict[factor] for factor in all_factors]
            if assert_factors:
                assert len(set(all_factors_powers)) == 1
            else:
                if not len(set(all_factors_powers)) == 1:
                    print("Warning: would have failed assertion, multiple powers appear:",
                          set([polynomial_power_dict[factor] for factor in all_factors]))
            # Following unmatched poly factors check could be improved when one of the two Warnings above are raised
            unmatched_polynomial_power_dict = {key: val for key, val in unmatched_polynomial_power_dict.items() if key not in all_factors}
            matched_irreds[candidate] = all_factors_powers[-1]
            matched_degree += sum([factor.as_poly(modulus=field.characteristic).degree() for factor in all_factors]) * matched_irreds[candidate]
    if verbose:
        print(f"Polynomial {polynomial_power_dict}")
        print(f"Matched {matched_degree} / {polynomial_degree}: {matched_irreds}")
    if assert_full_match:
        if polynomial_degree != matched_degree:
            print(f"Unmatched factors: {unmatched_polynomial_power_dict}")
            raise AssertionError
    return matched_irreds


def univariate_Newton_on_slice(oFunc, oSlice, verbose=False):
    field = oSlice.field

    def f(tval):
        oPoint = oSlice.copy()
        oPoint.subs({'t': tval})
        return ModP(int(oFunc(oPoint)), oPoint.field.characteristic)

    rat_func_t = Newton_polynomial_interpolation(f, field.characteristic, verbose=verbose)
    return rat_func_t


def univariate_Thiele_on_slice(oFunc, oSlice, verbose=False):
    field = oSlice.field

    def f(tval):
        oPoint = oSlice.copy()
        oPoint.subs({'t': tval})
        return ModP(int(oFunc(oPoint)), oPoint.field.characteristic)

    rat_func_t = Thiele_rational_interpolation(f, field.characteristic, verbose=verbose)
    return rat_func_t


def univariate_Thiele_on_slice_given_LCD(oFunc, oTerms, oSlice, verbose=False):
    """univariate rational functions of t via Newton polynomial interpolation of BlackBox function, with common rational factor pulled out"""
    # this is guaranteed to be a polynomial, if full set of denominator factors is known
    tnum = univariate_Newton_on_slice(lambda oPs: oFunc(oPs) / oTerms(oPs), oSlice, verbose=verbose)
    # this may be a rational function if common numerator factor is found
    tdenom = univariate_Thiele_on_slice(oTerms, oSlice, verbose=verbose)
    # univariate field of fraction of galois field
    FFGF = sympy.GF(oSlice.field.characteristic).frac_field(sympy.symbols('t'))
    return (FFGF(tnum) * FFGF(tdenom)).as_expr()
    # return FFGF(tnum) * FFGF(tdenom)


def do_codimension_one_study(oFunc, oSlice, denominator_candidates, assert_factors=True, verbose=False):
    field = oSlice.field
    rat_func_t = univariate_Thiele_on_slice(oFunc, oSlice, verbose=verbose)
    if verbose:
        print("\n", rat_func_t)
    numerator = rat_func_t.as_numer_denom()[0]
    denominator = rat_func_t.as_numer_denom()[1]
    invariant_dict = get_invariant_dict(tuple(denominator_candidates), oSlice)
    numerator_dict = match_factors(numerator, invariant_dict, field, assert_full_match=False, assert_factors=assert_factors, verbose=verbose)
    denominator_dict = match_factors(denominator, invariant_dict, field, assert_full_match=True, assert_factors=assert_factors, verbose=verbose)
    return Terms([Term(Numerator([list(numerator_dict.keys())], [list(map(int, numerator_dict.values()))], [1]),
                       Denominator(list(denominator_dict.keys()), list(map(int, denominator_dict.values()))))])
