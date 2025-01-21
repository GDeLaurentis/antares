#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Author: Giuseppe

import os
import sys
import functools
import operator
import pandas
import numpy
import mpmath

from fractions import Fraction as Q

from lips import Particles
from lips.invariants import Invariants
from lips.symmetries import phase_weights_compatible_symmetries

# from linac.tensor_function import tensor_function as _tensor_function

from pyadic.finite_field import ModP, rationalise
from pyadic.padic import padic_log

from antares.core.settings import settings
from antares.scalings.pair import pair_scalings

local_directory = os.path.dirname(os.path.abspath(__file__))
mpmath.mp.dps = 300

if sys.version_info.major > 2:
    unicode = str


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def composed(*decorators):
    """
    Apply an arbitrary number of decorators in a single line.

    The following are equivalent:

    Examples
    --------
    Using `composed`:

    >>> @composed(decorator1, decorator2)
    ... def func(...):
    ...     pass

    Using individual decorators:

    >>> @decorator1
    ... @decorator2
    ... def func(...):
    ...     pass

    Parameters
    ----------
    decorators : callable
        An arbitrary number of decorator functions to be applied.

    Returns
    -------
    callable
        A function decorated with all the provided decorators.
    """
    def composed_decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            f = func
            for decorator in reversed(decorators):
                f = decorator(f)
            return f(*args, **kwargs)
        return wrapper
    return composed_decorator


def as_scalar_if_scalar(func):
    """Turns numpy arrays with zero dimensions into 'real' scalars."""
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        res = func(*args, **kwargs)
        if res.shape == ():
            return res[()]  # pops the scalar out of array(scalar) or does nothing if array has non-trivial dimensions.
        elif functools.reduce(operator.mul, res.shape) == 1:
            return res.flatten()[0]
        else:
            return res
    return wrapper


def numpy_vectorized(*decorator_args, **decorator_kwargs):
    """Similar to numpy.vectorize when used as a decorator, but retains the __name__ of the decorated function and accepts args and kwargs."""
    if len(decorator_args) == 1 and len(decorator_kwargs) == 0 and callable(decorator_args[0]):
        func = decorator_args[0]

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            return numpy.vectorize(func)(*args, **kwargs)
        return wrapper
    else:
        def numpy_vectorized_decorator(func):
            @functools.wraps(func)
            def wrapper(*args, **kwargs):
                return numpy.vectorize(func, *decorator_args, **decorator_kwargs)(*args, **kwargs)
            return wrapper
        return numpy_vectorized_decorator


@numpy_vectorized(otypes='O')
def regulated_division(x, y):
    """Like division, but sends 0 / 0 to 1."""
    return 1 if x == y else x / y


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


class Numerical_Methods:

    # Mass Dimension Functions

    @property
    def internal_masses(self):
        if hasattr(self, '_internal_masses'):
            return self._internal_masses
        else:
            return set()

    @internal_masses.setter
    def internal_masses(self, temp_internal_masses):
        self._internal_masses = temp_internal_masses

    @staticmethod
    @composed(as_scalar_if_scalar, numpy_vectorized)
    def mpc_mass_dimension(ratio, z):
        return float(round(2 * mpmath.log(abs(ratio)) / mpmath.log(z ** 2)) / 2)  # round twice quantity to get half integers

    @staticmethod
    @composed(as_scalar_if_scalar, numpy_vectorized)
    def padic_mass_dimension(ratio, z):
        return float(rationalise(int(padic_log(ratio, z ** 2)), settings.field.characteristic ** settings.field.digits))

    @staticmethod
    @composed(as_scalar_if_scalar, numpy_vectorized(excluded=[1, 2]))
    def finite_field_mass_dimension(ratio, search_range, powers):
        if int(ratio) not in search_range:
            raise ValueError("Mass dimension not found in search range.")
        return powers[search_range.index(int(ratio))]

    @property
    def mass_dimension(self):
        """Returns the mass (a.k.a. energy) dimension. The value is expected to be a (half-)integer. Return type is float, vectorized depending on input."""
        oParticles = Particles(self.multiplicity, field=settings.field, seed=0, internal_masses=self.internal_masses)
        before = self(oParticles.copy())
        z = 3  # do not use 2! Different (small) powers have same values in finite fields.
        for oParticle in oParticles:
            oParticle.r_sp_u = oParticle.r_sp_u * z
            oParticle.l_sp_u = oParticle.l_sp_u * z
        for internal_mass in oParticles.internal_masses:
            mass_pow = int(internal_mass[-1]) if internal_mass[-1].isdigit() else 1
            oParticles.__setattr__(internal_mass, getattr(oParticles, internal_mass) * z ** 2 ** mass_pow)
        after = self(oParticles)
        if settings.field.name == "mpc":
            mass_dimension = self.mpc_mass_dimension(regulated_division(after, before), z)
        elif settings.field.name == "padic" and settings.field.digits > 1:
            mass_dimension = self.padic_mass_dimension(regulated_division(after, before), z)
        elif settings.field.name in ["padic", "finite field"]:
            search_start, search_stop = -500, 500
            powers = [entry / 2 for entry in range(search_start, search_stop, 1)]
            search_range = [int(ModP(Q(z) ** int(2 * power), settings.field.characteristic)) for power in powers]
            assert len(search_range) == len(set(search_range))  # there should be no duplicates (like e.g. for z = 2)
            mass_dimension = self.finite_field_mass_dimension(regulated_division(after, before), search_range, powers)
        return mass_dimension

    @staticmethod
    def mass_dimension_lParticles(multiplicity):
        oParticles = Particles(multiplicity, field=settings.field, seed=0)
        lParticles = [oParticles.copy()]
        z = 3  # do not use 2! Different (small) powers have same values in finite fields.
        for oParticle in oParticles:
            oParticle.r_sp_u = oParticle.r_sp_u * z
            oParticle.l_sp_u = oParticle.l_sp_u * z
        lParticles += [oParticles.copy()]
        return lParticles

    # Phase Weights Functions

    @staticmethod
    @composed(as_scalar_if_scalar, numpy_vectorized)
    def mpc_phase_weight(ratio, z):
        return int(round(2 * mpmath.log(abs(ratio)) / mpmath.log(z)) / 2)  # round twice quantity to get half integers

    @staticmethod
    @composed(as_scalar_if_scalar, numpy_vectorized)
    def padic_phase_weight(ratio, z):
        return int(rationalise(int(padic_log(ratio, z)), settings.field.characteristic ** settings.field.digits))

    @staticmethod
    @composed(as_scalar_if_scalar, numpy_vectorized(excluded=[1, 2]))
    def finite_field_phase_weight(ratio, search_range, powers):
        if int(ratio) not in search_range:
            raise ValueError("Phase weight not found in search range.")
        return powers[search_range.index(int(ratio))]

    @property
    def phase_weights(self):
        """Returns the phase weights (a.k.a. little group scalings). Return type is a list of int's, vectorized depending on input."""
        oParticles = Particles(self.multiplicity, field=settings.field, seed=0, internal_masses=self.internal_masses)
        phase_weights = []
        z = 3  # do not use 2! Different (small) powers have same values in finite fields.
        if settings.field.name in ["padic", "finite field"]:
            search_start, search_stop = -100, 100
            powers = range(search_start, search_stop, 1)
            search_range = [int(ModP(Q(z) ** power, settings.field.characteristic)) for power in powers]
        for oParticle in oParticles:
            before = self(oParticles.copy())
            oParticle.r_sp_u = oParticle.r_sp_u * z
            oParticle.l_sp_u = oParticle.l_sp_u / z
            after = self(oParticles)
            if settings.field.name == "mpc":
                temp_p_w = self.mpc_phase_weight(regulated_division(after, before), z)
            elif settings.field.name == "padic" and settings.field.digits > 1:
                temp_p_w = self.padic_phase_weight(regulated_division(after, before), z)
            elif settings.field.name in ["padic", "finite field"]:
                temp_p_w = self.finite_field_phase_weight(regulated_division(after, before), search_range, powers)
            phase_weights += [temp_p_w]
        phase_weights = numpy.moveaxis(numpy.array(phase_weights), 0, -1)
        if len(phase_weights.shape) == 1:
            phase_weights = phase_weights.tolist()
        return phase_weights

    @property
    def all_symmetries(self, possible_symmetries=None, field=settings.field):
        """Obtain the symmetries of the function, by numerical evaluation."""
        symmetries = []
        if possible_symmetries is None:
            possible_symmetries = phase_weights_compatible_symmetries(self.phase_weights)
        for sym in possible_symmetries:
            sym = sym[:2]  # make sure this sym does not have a sign yet
            oParticles = Particles(self.multiplicity, seed=0, field=field)
            # One may wish to obtain symmetries valid only on certain codimension X surfaces
            # may need to adjust tollerance for this, as the symmetry would be approximate
            # e.g.: oParticles._set("Δ_23|14|56", 10 ** -50,)
            oNewParticles = oParticles.image(sym)
            if abs(self(oParticles) - self(oNewParticles)) <= field.tollerance:
                symmetries += [sym + ("+", )]
            elif abs(self(oParticles) + self(oNewParticles)) <= field.tollerance:
                symmetries += [sym + ("-", )]
        return symmetries

    @staticmethod
    def phase_weights_lParticles(multiplicity):
        lParticles = []
        oParticles = Particles(multiplicity, field=settings.field, seed=0)
        z = 3  # do not use 2! Different (small) powers have same values in finite fields.
        lParticles += [oParticles.copy()]
        for oParticle in oParticles:
            oParticle.r_sp_u = oParticle.r_sp_u * z
            oParticle.l_sp_u = oParticle.l_sp_u / z
            lParticles += [oParticles.copy()]
        return lParticles

    @staticmethod
    def mass_dimension_and_phase_weights_lParticles(multiplicity):
        return list(set(Numerical_Methods.mass_dimension_lParticles(multiplicity) +
                        Numerical_Methods.phase_weights_lParticles(multiplicity)))

    # Rest

    @property
    def soft_weights(self):
        oParticles = Particles(self.multiplicity, seed=0)
        soft_weights = [[], []]
        for i, oParticle in enumerate(oParticles):
            a, b = self.multiplicity - 1, self.multiplicity
            if i in [a, b]:
                a, b = 1, 2
            # Angles
            oParticles.randomise_all()
            oParticle.r_sp_u = numpy.array([mpmath.mpc('1e-50', '2e-50'), mpmath.mpc('3e-50', '2e-50')])
            oParticles.fix_mom_cons(a, b)
            before = self(oParticles)
            oParticle.r_sp_u = oParticle.r_sp_u / 10
            oParticles.fix_mom_cons(a, b)
            after = self(oParticles)
            temp_p_w = mpmath.log((abs(after / before))) / mpmath.log(10)
            soft_weights[0] += [int(round(float(temp_p_w)))]
            # Squares
            oParticles.randomise_all()
            oParticle.l_sp_u = numpy.array([mpmath.mpc('1e-50', '2e-50'), mpmath.mpc('3e-50', '2e-50')])
            oParticles.fix_mom_cons(a, b)
            before = self(oParticles)
            oParticle.l_sp_u = oParticle.l_sp_u / 10
            oParticles.fix_mom_cons(a, b)
            after = self(oParticles)
            temp_p_w = mpmath.log((abs(after / before))) / mpmath.log(10)
            soft_weights[1] += [int(round(float(temp_p_w)))]
        return soft_weights

    @property
    def is_iterable(self):
        return hasattr(self, "__iter__")

    @property
    def is_zero(self):
        oParticles = Particles(self.multiplicity, field=settings.field, internal_masses=self.internal_masses)
        if abs(self(oParticles)) > settings.field.tollerance:
            return False
        else:
            return True

    def do_single_collinear_limits(self, invariants=None, verbose=True):
        from antares.terms.term import Term
        from antares.terms.terms import Terms
        self.oTermFromSingleScalings = Term.from_single_scalings(self, invariants=invariants, verbose=verbose)
        if hasattr(self, "save_original_unknown_call_cache"):
            self.save_original_unknown_call_cache()
        self.num_invs, _num_exps = self.oTermFromSingleScalings.oNum.llInvs[0], self.oTermFromSingleScalings.oNum.llExps[0]
        self.den_invs, _den_exps = self.oTermFromSingleScalings.oDen.lInvs, self.oTermFromSingleScalings.oDen.lExps
        self.num_exps = dict(zip(self.num_invs, _num_exps))
        self.den_exps = dict(zip(self.den_invs, _den_exps))
        self.oTermsFromSingleScalings = Terms([self.oTermFromSingleScalings])
        self.oTermsFromSingleScalings.oUnknown = self
        if verbose:
            print("")
            print("Mass dimension & phase weights:", end='')
            print("{},".format(self.mass_dimension), end='')
            print(self.phase_weights, end='')
            sys.stdout.flush()
            print(u" \u2192 ", end='')
            print("{},".format(self.oTermsFromSingleScalings.ansatze_mass_dimensions[0]), end='')
            print(self.oTermsFromSingleScalings.ansatze_phase_weights[0])
            print("")

    def do_double_collinear_limits(self, silent=False):
        oInvariants = Invariants(self.multiplicity, Restrict3Brackets=settings.Restrict3Brackets,
                                 Restrict4Brackets=settings.Restrict4Brackets, FurtherRestrict4Brackets=settings.FurtherRestrict4Brackets)
        _pair_invs, _pair_exps, _pair_friends = pair_scalings(self, self.den_invs, self.den_invs, oInvariants.full)
        if hasattr(self, "save_original_unknown_call_cache"):
            self.save_original_unknown_call_cache()
        partial_filter = functools.partial(filter, lambda x: x in self.den_invs or (x in self.spurious_poles if hasattr(self, 'spurious_poles') else False) or x == "F")
        _true_friends = list(map(list, map(partial_filter, _pair_friends)))
        if silent is False:
            print("\rFinished calculating pair scalings. They are:                         ")
        else:
            print("\r                                                                      ")
        if _pair_exps != []:
            for i, pair_exp in enumerate(_pair_exps):
                if pair_exp != 'F' and (pair_exp == int(pair_exp) or pair_exp * 2 == int(pair_exp * 2)):
                    _pair_exps[i] = abs(_pair_exps[i])
            # old printing
            if silent is False:
                col_width = max([len(pair_inv[0]) + len(pair_inv[1]) + 6 for pair_inv in _pair_invs])
                for i in range(len(_pair_invs)):
                    if hasattr(self, "recursively_extract_original_unknown"):
                        if tuple(_pair_invs[i]) in self.recursively_extract_original_unknown().pair_exps.keys():
                            original_scaling = self.recursively_extract_original_unknown().pair_exps[tuple(_pair_invs[i])]
                        else:
                            original_scaling = 0.0
                        original_scaling = " ({})".format(original_scaling)
                    else:
                        original_scaling = ""
                    print(("[" + _pair_invs[i][0] + ", " + _pair_invs[i][1] + "]:").ljust(col_width) + unicode(_pair_exps[i]) + original_scaling + ", " + (unicode(
                        len(_pair_friends[i]))).ljust(2) + u" \u2192 " + unicode(len(_true_friends[i])))
                print("")
        else:
            if silent is False:
                print("Non-existent")
        # Remove the failed ones.
        # __pair_invs = [entry for entry in _pair_invs]
        # __pair_exps = [entry for entry in _pair_exps]
        # for j, jpair in enumerate(__pair_invs):
        #     if __pair_exps[j] == "F":
        #         index = _pair_invs.index(jpair)
        #         _pair_invs.pop(index)
        #         _pair_exps.pop(index)
        #         _pair_friends.pop(index)
        self.pair_friends = dict(zip(list(map(tuple, _pair_invs + [pair[::-1] for pair in _pair_invs])), _pair_friends + _pair_friends))
        self.pair_exps = dict(zip(list(map(tuple, _pair_invs + [pair[::-1] for pair in _pair_invs])), _pair_exps + _pair_exps))
        self.true_friends = dict(zip(list(map(tuple, _pair_invs + [pair[::-1] for pair in _pair_invs])), _true_friends + _true_friends))

    @property
    def collinear_data(self):
        df = pandas.DataFrame(columns=self.den_invs, index=self.den_invs, data=[
            [str(self.pair_exps[inv1, inv2] if type(self.pair_exps[inv1, inv2]) in [str, unicode] else int(self.pair_exps[inv1, inv2])
                 if self.pair_exps[inv1, inv2].is_integer() else self.pair_exps[inv1, inv2]) +
             "/" + str(len(self.pair_friends[inv1, inv2])) + "/" + str(len(self.true_friends[inv1, inv2]))
             if (inv1, inv2) in self.pair_exps else self.den_exps[inv1] if inv1 == inv2 else None for inv2 in self.den_invs] for inv1 in self.den_invs])
        df.style.caption = "Collinear data.\nPower/Degeneracy of phase space/Degeneracy of restricted phase space."
        return df


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def memoized(*decorator_args, **decorator_kwargs):
    """Diskcaching decorator generator."""
    def memoized_decorator(func):

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            self = args[0]

            if not hasattr(self, 'diskcache'):
                return func(*args, **kwargs)

            @self.diskcache.memoize(*decorator_args, **decorator_kwargs)
            def diskcached_func(*args, **kwargs):
                return func(*args, **kwargs)
            return diskcached_func(*args, **kwargs)
        return wrapper
    return memoized_decorator


class num_func(Numerical_Methods, object):

    def __init__(self, evaluable_function, verbose=False):
        super().__init__()
        self.evaluable_function = evaluable_function
        self.multiplicity = evaluable_function.multiplicity
        if hasattr(evaluable_function, 'internal_masses'):
            self.internal_masses = evaluable_function.internal_masses
        self.__name__ = evaluable_function.__name__
        self.spurious_poles = []
        self.easy_boxes = []
        self.amppart = None
        if verbose:
            print("Multiplicity set to {}, name set to '{}'".format(self.multiplicity, self.__name__))

    @memoized(name='num_func.__call__', ignore={0})
    def __call__(self, oParticles):
        return self.evaluable_function(oParticles)


class _tensor_function(object):
    """Tensor function supporting indexing and iteration.
       For instance, the initializer 'callable_function' can be a function returning a numpy.array"""
    # Taken from linac, import from there once made public

    def __init__(self, callable_function):
        self.callable_function = callable_function
        if hasattr(callable_function, "__name__"):
            self.__name__ = callable_function.__name__

    def flatten(self):
        return tensor_function(lambda args: self(args).flatten())

    def __getitem__(self, index):
        return tensor_function(lambda args: self(args)[index])

    def __matmul__(self, other):
        assert isinstance(other, numpy.ndarray)
        return tensor_function(lambda args: self(args) @ other)

    @memoized(name='tensor_function.__call__', ignore={0})
    def __call__(self, *args, **kwargs):
        res = self.callable_function(*args, **kwargs)
        if hasattr(res, "shape") and not hasattr(self, "__shape__"):
            self.shape = res.shape
        elif isinstance(res, list):
            self.__size__ = len(res)
        return res

    def __len__(self):
        if hasattr(self, "shape"):
            return self.shape[0]
        elif hasattr(self, "__size__"):
            return self.__size__
        else:
            raise AttributeError("Length not known. Have you tried evaluating the function at least once?")

    @property
    def shape(self):
        if hasattr(self, "__shape__"):
            return self.__shape__
        else:
            raise AttributeError("Shape not known. Have you tried evaluating the function at least once?")

    @shape.setter
    def shape(self, value):
        self.__shape__ = value

    def __iter__(self):
        for entry in self.flatten():
            yield entry


class tensor_function(Numerical_Methods, _tensor_function):
    pass
