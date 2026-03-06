import functools
import numpy
import sympy

from copy import copy
from collections import OrderedDict
from multiset import Multiset
from pyadic import ModP, PAdic
from pyadic.interpolation import Newton_polynomial_interpolation, Thiele_rational_interpolation, multivariate_Newton_polynomial_interpolation
from pycoretools import flatten

from ..terms.terms import Terms, Term, Numerator, Denominator


@functools.lru_cache(maxsize=1024)
def get_invariant_dict(possible_denominators, evaluator, keep_all_non_unique=False, field=None, verbose=True):
    if field is None:
        field = evaluator.field
    candidate_denoms_dict = {}
    for possible_denominator in possible_denominators:
        possible_denominator_of_t = evaluator(possible_denominator)
        if isinstance(possible_denominator_of_t, (ModP, PAdic, sympy.core.numbers.Number)):
            candidate_denoms_dict[possible_denominator] = possible_denominator_of_t
        elif isinstance(possible_denominator_of_t, numpy.ndarray):
            pass
        else:
            this_candidate = 4 * possible_denominator_of_t.expand()
            if this_candidate == 0 and field.name in ['finite field', 'Fp']:
                print(f"Warning: candidate denominator {possible_denominator} evaluates to zero on an Fp slice, skipping it.")
                continue
            this_candidate = sympy.poly(
                this_candidate, modulus=field.characteristic
            ).as_expr().factor(modulus=field.characteristic)
            if this_candidate == 0 and field.name in ['padic', 'Qp']:
                this_candidate = sympy.poly(
                    4 * possible_denominator_of_t.expand() / field.characteristic, modulus=field.characteristic
                ).as_expr().factor(modulus=field.characteristic)
            candidate_denoms_dict[possible_denominator] = this_candidate
    constant_denom_candidates = dict((key, val) for key, val in candidate_denoms_dict.items() if isinstance(val, (ModP, PAdic)) or val.free_symbols == set())
    if verbose and len(constant_denom_candidates) > 0:
        print(f"The following candidate demoninator factors are constant, they will be discarded:\n{constant_denom_candidates}")
    candidate_denoms_dict = {key: val for key, val in candidate_denoms_dict.items() if not isinstance(val, (ModP, PAdic))}
    if not set(flatten([list(val.as_powers_dict().values()) for key, val in candidate_denoms_dict.items()])) == {1}:
        # raise NotImplementedError(f"Expecting only factors to power 1 on a generic slice, dropping factors with higher powers.\
        # Got powers {set(flatten([list(val.as_powers_dict().values()) for key, val in candidate_denoms_dict.items()]))}")
        print("Warning: expected only factors to power 1 on a generic slice.", end=" ")
        print(f"Got powers {set(flatten([list(val.as_powers_dict().values()) for key, val in candidate_denoms_dict.items()]))}.", end=" ")
        print("Dropping factors with higher powers.")
        candidate_denoms_dict = {key: val for key, val in candidate_denoms_dict.items() if set(list(val.as_powers_dict().values())) == {1}}
    # check that no factors are shared
    candidate_denom_factors_dict = {key: [entry for entry in val.as_ordered_factors() if 't' in str(entry)] for key, val in candidate_denoms_dict.items()}
    candidate_unique_factors = [key for key, val in Multiset(flatten([val for key, val in candidate_denom_factors_dict.items()])).items() if val == 1]
    unique_candidate_denom_factors_dict = {key: [entry for entry in val if entry in candidate_unique_factors] for key, val in candidate_denom_factors_dict.items()}
    unique_candidate_denom_factors_dict = {key: val for key, val in unique_candidate_denom_factors_dict.items() if val != []}
    non_unique_candidate_denom_factors_dict = {key: val for key, val in candidate_denom_factors_dict.items() if key not in unique_candidate_denom_factors_dict.keys()}
    if not non_unique_candidate_denom_factors_dict == {}:
        print("Warning: non unique candidate denominator factors found. Are you sure they all generate distinct ideals at codimension 1?")
        print("All dicts:", candidate_denom_factors_dict)
        print("Unique:", unique_candidate_denom_factors_dict)
        print("Non unique:", non_unique_candidate_denom_factors_dict)
    if keep_all_non_unique:
        print("Keeping all non-unique candidate denominator factors.")
        return candidate_denom_factors_dict
    if not non_unique_candidate_denom_factors_dict == {}:
        print("Keeping one representative for each family. Warning: families need not be disjoint.")  # TO BE IMPROVED
        # Remove any entries with empty value lists
        non_unique_candidate_denom_factors_dict = {
            key: val for key, val in non_unique_candidate_denom_factors_dict.items() if val
        }

        # Sort by length of the key
        non_unique_candidate_denom_factors_dict = OrderedDict(
            sorted(non_unique_candidate_denom_factors_dict.items(), key=lambda kv: len(kv[0]))
        )

        seen = set()
        filtered_dict = OrderedDict()
        equivalence_classes = {}  # maps filtered-out keys to the key that claimed the exprs

        for k, v in non_unique_candidate_denom_factors_dict.items():
            if all(e not in seen for e in v):
                filtered_dict[k] = v
                seen.update(v)
            else:
                # Find the first key in filtered_dict whose values overlap
                for rep_key, rep_vals in filtered_dict.items():
                    if any(e in rep_vals for e in v):
                        equivalence_classes.setdefault(rep_key, []).append(k)
                        break

        # Print equivalence classes
        for rep, aliases in equivalence_classes.items():
            print(f"Equivalence class for {rep}:")
            for alias in aliases:
                print(f"  ↪ {alias}")

        # Store the final filtered result
        non_unique_candidate_denom_factors_dict = filtered_dict
    return unique_candidate_denom_factors_dict | non_unique_candidate_denom_factors_dict


def match_factors(polynomial, candidate_factors_dict, field, degree_bounds={}, assert_full_match=False, assert_factors=True, verbose=False):
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
    if verbose:
        print("[match factors] Polynomial:", polynomial_power_dict)
    key_val_pairs = list(candidate_factors_dict.items())
    key_val_pairs = sorted(candidate_factors_dict.items(),
                           key=lambda item: (0 if 's' in item[0] else 1, -len(item[1]), len(item[0])))
    for candidate, all_factors in key_val_pairs:
        if all([factor in unmatched_polynomial_power_dict.keys() for factor in all_factors]):  # check all factors appear in left over part
            if verbose:
                print(f"[match factors] {candidate} matches {all_factors}")
            # take the lowest power with which they appear
            all_factors_powers = [unmatched_polynomial_power_dict[factor] for factor in all_factors]
            min_power = min(all_factors_powers)
            reject = False
            # apply degree bounds
            for key, val in degree_bounds.items():
                if any(entry in candidate for entry in key):
                    if val <= 0:
                        if verbose:
                            print(f"[match factors] {key} was found in {candidate}, rejected due to left over degree bound {val}")
                        reject = True
                        break
                    else:
                        if verbose:
                            print(f"[match factors] {key} was found in {candidate}, applying degree bound {val}")
                        if val < min_power:
                            min_power = val
                        degree_bounds[key] = val - min_power  # for now assume it can only appear linearly within the candidate, to be improved if needed
            if reject:
                continue
            # Following unmatched poly factors check could be improved when one of the two Warnings above are raised
            unmatched_polynomial_power_dict = {
                key: (val if key not in all_factors else val - min_power)
                for key, val in unmatched_polynomial_power_dict.items()
                if (val if key not in all_factors else val - min_power) != 0
            }
            if verbose:
                print(f"[match factors] left over {unmatched_polynomial_power_dict}")
            matched_irreds[candidate] = min_power
            matched_degree += sum([factor.as_poly(modulus=field.characteristic).degree() for factor in all_factors]) * min_power
        # Turned off - to be revised
        if False and any([factor in polynomial_power_dict.keys() for factor in all_factors]):  # check if any matches
            all_factors_appear = all([factor in polynomial_power_dict.keys() for factor in all_factors])
            # check that all factors actually do appear
            if assert_factors:
                assert all_factors_appear, ("Assertion failed: not all factors appear:", candidate, all_factors)
            else:
                if not all_factors_appear:
                    print("Warning: would have failed assertion, some factors are missing:", candidate,
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
        print(f"[match factors] Matched {matched_degree} / {polynomial_degree}: {matched_irreds}")
    if assert_full_match:
        if polynomial_degree != matched_degree:
            print(f"[match factors] Unmatched factors: {unmatched_polynomial_power_dict}")
            raise AssertionError
    # sort back into the original order of candidates
    matched_irreds = {key: matched_irreds[key] for key in candidate_factors_dict if key in matched_irreds}
    return matched_irreds


def univariate_Newton_on_slice(oFunc, oSlice, verbose=False):
    field = oSlice.field

    def f(tval):
        oPoint = oSlice.copy()
        oPoint.subs({'t': tval})
        return ModP(oFunc(oPoint), oPoint.field.characteristic)

    rat_func_t = Newton_polynomial_interpolation(f, field.characteristic, verbose=verbose)
    return rat_func_t


def univariate_Thiele_on_slice(oFunc, oSlice, verbose=False):
    field = oSlice.field

    def f(tval):
        oPoint = oSlice.copy()
        oPoint.subs({'t': tval})
        return ModP(oFunc(oPoint), oPoint.field.characteristic)

    rat_func_t = Thiele_rational_interpolation(f, field.characteristic, verbose=verbose)
    return rat_func_t


def bivariate_Newton_on_slice(oFunc, oSlice, verbose=False):
    field = oSlice.field

    def f(tval1, tval2):
        oPoint = oSlice.copy()
        oPoint.subs({'t1': tval1, 't2': tval2})
        return ModP(oFunc(oPoint), oPoint.field.characteristic)

    rat_func_ts = multivariate_Newton_polynomial_interpolation(f, field.characteristic, verbose=verbose)
    return rat_func_ts


def univariate_Thiele_on_slice_given_LCD(oFunc, oTerms, oSlice, verbose=False):
    """univariate rational functions of t via Newton polynomial interpolation of BlackBox function, with common rational factor pulled out"""
    # this is guaranteed to be a polynomial, if full set of denominator factors is known
    tnum = univariate_Newton_on_slice(lambda oPs: oFunc(oPs) / oTerms(oPs), oSlice, verbose=verbose)
    if verbose:
        print()
    # this may be a rational function if common numerator factor is found
    tdenom = univariate_Thiele_on_slice(oTerms, oSlice, verbose=verbose)
    # univariate field of fraction of galois field
    FFGF = sympy.GF(oSlice.field.characteristic).frac_field(sympy.symbols('t'))
    return (FFGF(tnum) * FFGF(tdenom)).as_expr()
    # return FFGF(tnum) * FFGF(tdenom)


def bivariate_Thiele_on_slice_given_LCD(oFunc, oTerms, oSlice, verbose=False):
    """univariate rational functions of t via Newton polynomial interpolation of BlackBox function, with common rational factor pulled out"""
    # this is guaranteed to be a polynomial, if full set of denominator factors is known
    if verbose:
        print("[bivariate Thiele on slice given LCD] Obtaining unknown numerator:")
    tnum = bivariate_Newton_on_slice(lambda oPs: oFunc(oPs) / oTerms(oPs), oSlice, verbose=verbose)
    # this may be a rational function if common numerator factor is found - split it here
    if verbose:
        print("[bivariate Thiele on slice given LCD] Obtaining known denominator:")
    tdenom = bivariate_Newton_on_slice(oTerms[0].oDen.as_term(), oSlice, verbose=verbose)
    if verbose:
        print("[bivariate Thiele on slice given LCD] Obtaining known numerator:")
    tnumknown = bivariate_Newton_on_slice(oTerms[0].oNum.as_term(), oSlice, verbose=verbose)
    # univariate field of fraction of galois field
    if verbose:
        print("[bivariate Thiele on slice given LCD] Assembling result:")
    return tnum * tnumknown / tdenom
    FFGF = sympy.GF(oSlice.field.characteristic).frac_field(*sympy.symbols(['t1', 't2']))
    return (FFGF(tnum) * FFGF(tnumknown) / FFGF(tdenom)).as_expr()


def do_codimension_one_study(oFunc, oSlice, denominator_candidates, oSlice_for_invariants=None, oTermsDenom=None,
                             assert_factors=True, degree_bounds={}, keep_all_non_unique=False, verbose=False):
    """Returns the least common denominator (incl. common numerator factors)
    from a univariate slice 'oSlice' and a sufficient list of guesses 'denominator_candidates'.
    If 'oTermsDenom' is provided univariate_Thiele_on_slice_given_LCD is used, instead of univariate_Thiele_on_slice.
    'oTermsDenom' may be larger than the true LCD.
    """
    field = oSlice.field
    if oTermsDenom is None:
        rat_func_t = univariate_Thiele_on_slice(oFunc, oSlice, verbose=verbose)
    else:
        rat_func_t = univariate_Thiele_on_slice_given_LCD(oFunc, oTermsDenom, oSlice, verbose=verbose)
    if oSlice_for_invariants is None:
        oSlice_for_invariants = oSlice
    if verbose:
        print("\n", rat_func_t)
    numerator = rat_func_t.as_numer_denom()[0]
    denominator = rat_func_t.as_numer_denom()[1]
    invariant_dict = get_invariant_dict(tuple(denominator_candidates), oSlice_for_invariants, keep_all_non_unique=keep_all_non_unique)
    numerator_dict = match_factors(numerator, invariant_dict, field, assert_full_match=False, degree_bounds=degree_bounds, assert_factors=assert_factors, verbose=verbose)
    denominator_dict = match_factors(denominator, invariant_dict, field, assert_full_match=True, assert_factors=assert_factors, verbose=verbose)
    return Terms([Term(Numerator([list(numerator_dict.keys())], [list(map(int, numerator_dict.values()))], [1]),
                       Denominator(list(denominator_dict.keys()), list(map(int, denominator_dict.values()))))])
