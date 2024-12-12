#!/usr/bin/env python
# -*- coding: utf-8 -*-

#   ___ _  _ ___      _      _
#  | _ ) || | _ \__ _| |_ __| |_
#  | _ \ __ |  _/ _` |  _/ _| ' \
#  |___/_||_|_| \__,_|\__\__|_||_|
#

# Author: Giuseppe

from __future__ import print_function

import sys
import numpy
import collections

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

old_print = print


def file_compatible_print(*args, **kwargs):
    if sys.stdout.isatty() is True or not hasattr(sys.stdout, "__str__") or "ipykernel.iostream.OutStream" in sys.stdout.__str__():
        old_print(*args, **kwargs)
        if "end" in kwargs and "\r" in kwargs["end"]:
            sys.stdout.flush()
    else:
        to_be_printed = args[0]
        if "\033[1A" in to_be_printed:
            to_be_printed = to_be_printed.replace("\033[1A", "")
            if len(file_compatible_print.refs) >= 2:
                sys.stdout.seek(file_compatible_print.refs[-2])
            elif len(file_compatible_print.refs) >= 1:
                sys.stdout.seek(file_compatible_print.refs[-1])
            file_compatible_print(to_be_printed, **kwargs)
        elif len(to_be_printed) > 0 and to_be_printed[0] == "\r":
            sys.stdout.seek(file_compatible_print.refs[-1])
            to_be_printed = to_be_printed.replace("\r", "")
            file_compatible_print(to_be_printed, **kwargs)
        elif "end" not in kwargs or ("end" in kwargs and "\n" in kwargs["end"]):
            file_compatible_print.refs.append(sys.stdout.tell())
            old_print(to_be_printed.encode('utf-8'), **kwargs)
            file_compatible_print.refs.append(sys.stdout.tell())
        elif "end" in kwargs and kwargs["end"] == "":
            file_compatible_print.ref.append(sys.stdout.tell())
            old_print(to_be_printed.encode('utf-8'), **kwargs)
        elif "end" in kwargs and "\r" in kwargs["end"]:
            old_print(to_be_printed.encode('utf-8'), **kwargs)
            sys.stdout.seek(file_compatible_print.refs[-1])
        else:
            raise Exception("Ending not supported by file_compatible_print.")


if sys.stdout.isatty() is False and hasattr(sys.stdout, "__str__") and "ipykernel.iostream.OutStream" not in sys.stdout.__str__():
    file_compatible_print.refs = collections.deque(maxlen=25)
# print = file_compatible_print

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

if sys.version_info[0] <= 2:
    try:
        import BH
        BH_found = True
    except ImportError:
        BH_found = False
    else:
        BH_found = True

    try:
        import gmpTools
    except ImportError:
        gmpTools_found = False
    else:
        gmpTools_found = True
else:
    BH_found = False
    gmpTools_found = False


if gmpTools_found and BH_found:
    def accuracy():
        return BH.RGMP_get_current_nbr_digits()
else:
    def accuracy():
        return 300


if gmpTools_found and BH_found:

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

    C = BH.CGMP
    R = BH.RGMP
    BH.to_int_prec = BH.RGMP(str('1e-6'))

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

    old_real = C.real

    @property
    def my_real(self):
        return old_real(self)

    C.real = my_real

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

    old_imag = C.imag

    @property
    def my_imag(self):
        return old_imag(self)

    C.imag = my_imag

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

    def log(x):
        if str(type(x)) in ["<type \'int\'>", "<type \'complex\'>", "<type \'numpy.complex128\'>", "<type \'numpy.float64\'>"]:
            return numpy.log(x)
        elif str(type(x)) == "<class \'BH.RGMP\'>":
            return BH.log(x)
        elif str(type(x)) in ["<type \'gmpTools.CGMP\'>", "<type \'gmpTools.RGMP\'>"]:
            return gmpTools.log(x)
        else:
            raise Exception("Log of type {} is not defined".format(type(x)))

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

    def to_CGMP(x):
        c1 = "<type \'complex\'>"
        c2 = "<type \'numpy.complex128\'>"
        if str(type(x)) == "<type \'int\'>":
            return C(R(x), R(0))
        elif str(type(x)) == c1 or str(type(x)) == c2:
            if x.real.is_integer() and x.imag.is_integer():
                return C(R(x.real), R(x.imag))
        elif str(type(x)) == "<class \'BH.RGMP\'>":
            return C(x, R(0))
        elif str(type(x)) == "<class \'BH.CGMP\'>":
            return x
        else:
            raise Exception("failed conversion to CGMP for type {}".format(str(type(x))))

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

    old_abs = abs

    def my_abs(num):
        if str(type(num)) == '<class \'BH.CGMP\'>':
            mod_num = num.real * num.real + num.imag * num.imag
            mod_num = BH.sqrt(mod_num)
            return mod_num
        elif str(type(num)) == '<class \'BH.RGMP\'>':
            mod_num = num * num
            mod_num = BH.sqrt(mod_num)
            return mod_num
        else:
            return old_abs(num)

    abs = my_abs

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

    def to_int(num):
        from core.my_warnings import oWarning as w
        from core.settings import settings
        rounded = to_int_inner(num)
        if str(type(num)) == '<class \'BH.RGMP\'>':
            diff = abs(R(rounded) - num)
            if diff > eval(settings.ToIntegerPrecision):
                _warn = "Tollerance exceeded in (RGMP) to_int. "
                _orig = "The original number was {}".format(repr(num))
                w.warn("{}{}".format(_warn, _orig), "Rounded precision")
            return rounded
        elif str(type(num)) == '<class \'gmpTools.RGMP\'>':
            diff = abs(gmpTools.RGMP(rounded) - num)
            if diff > gmpTools.RGMP(BH.to_int_prec):
                _warn = "Tollerance exceeded in (RGMP) to_int. "
                _orig = "The original number was {}".format(repr(num))
                w.warn("{}{}".format(_warn, _orig), "Rounded precision")
            return rounded
        else:
            diff = abs(rounded - num)
            if diff > float('1e-3'):
                _warn = "Tollerance exceeded in to_int. "
                _orig = "The original number was {}".format(num)
                w.warn("{}{}".format(_warn, _orig), "Rounded precision")
            return rounded

    def to_int_inner(num):
        if str(type(num)) == '<class \'BH.RGMP\'>':
            return BH.to_int(num)
        elif str(type(num)) == '<type \'gmpTools.RGMP\'>':
            return BH.to_int(num)
        else:
            return int(round(num))

    def __my_pow__(self, other):
        if other != int(other) and other * 2 != int(other * 2):
            raise Exception("Power implemented only for integers and half-integers.", "Not supported type in C.__pow__")
        elif other == int(other):
            result = to_CGMP(1)
            if other > 0:
                result = to_CGMP(1)
                for i in range(int(other)):
                    result = result * self
            elif other < 0:
                result = to_CGMP(1)
                for i in range(abs(int(other))):
                    result = result / self
            return result
        elif other * 2 == int(other * 2):
            result = to_CGMP(1)
            if other > 0:
                result = to_CGMP(1)
                for i in range(int(other * 2)):
                    result = result * self
            elif other < 0:
                result = to_CGMP(1)
                for i in range(abs(int(other * 2))):
                    result = result / self
            return C.sqrt(result)

    C.__pow__ = __my_pow__

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

    # Making SwigPyObjects Picklable

    # File "/home/gdl/Desktop/Programs/BHlib_Jul20//build/lib/blackhat/BH.py", line 69
    # getattr -> getattribute

    def __my_setstate__(self, state):
        # print "my_setstate", state[0], state[1]
        num = C(eval(state[0]), eval(state[1]))
        self.__init__(num)

    C.__setstate__ = __my_setstate__

    def __my_getstate__(self):
        real_repr = repr(self.real).replace("GMP", "").replace("(", "('").replace(")", "')")
        imag_repr = repr(self.imag).replace("GMP", "").replace("(", "('").replace(")", "')")
        # print "my_getstate", real_repr, imag_repr
        return real_repr, imag_repr

    C.__getstate__ = __my_getstate__

    old_init = C.__init__

    def __my_init__(self, *args):
        if str(type(args[0])) == "<class 'BH.CGMP'>":
            args = (args[0].real, args[0].imag)
            # print "my_init", args
        old_init(self, *args)

    C.__init__ = __my_init__
