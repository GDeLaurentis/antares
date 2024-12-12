#!/usr/bin/env python
# -*- coding: utf-8 -*-

#    ___             _____         _
#   / __|___ _ _ ___|_   _|__  ___| |___
#  | (__/ _ \ '_/ -_)_| |/ _ \/ _ \ (_-<
#   \___\___/_| \___(_)_|\___/\___/_/__/
#

# Author: Giuseppe

from __future__ import print_function
from __future__ import unicode_literals

import sys
import os
import subprocess
import re
import shelve
import math
import time
import random
import functools
import itertools
import multiprocessing
import multiprocessing.pool
import threading
import numpy
import types
import mpmath
import fractions
import copyreg

from antares.core.bh_patch import BH_found
if BH_found:
    from antares.core.bh_patch import BH

if sys.version_info[0] > 2:
    unicode = str

MainPythonDirectory = os.path.dirname(os.path.abspath(__file__))[:-5]


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def _pickle_method(m):
    if m.im_self is None:
        return getattr, (m.im_class, m.im_func.func_name)
    else:
        return getattr, (m.im_self, m.im_func.func_name)


copyreg.pickle(types.MethodType, _pickle_method)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def configuration_unpacker(sys_argv):
    help_msg = r"Usage: python numerical_to_analytical.py|analytical_to_analytical.py UniqueRunFileName HELCONF:[\+\-] " \
               r"AMPPART:[custom|tree|box|triangle|bubble|rational] AMPINDEX=[#] LOOPID:[G|nf]]"

    helconf_regex = re.compile(r'^(?:HELCONF[=|:])(?P<HELCONF>[\+\-qbpmy]+)$')
    amppart_regex = re.compile(r'^(?:AMPPART[=|:])(?P<AMPPART>custom$|external$|external_.*|tree$|rational$|box$|triangle$|bubble$)')
    ampindex_regex = re.compile(r'^(?:AMPINDEX[=|:])(?P<AMPINDEX>\d*)$')
    loopid_regex = re.compile(r'^(?:LOOPID[=|:])(?P<LOOPID>nf$|G$|RT$|LT$|LC$|SLC$|None$)$')

    if any([arg == "h" or arg == "help" for arg in sys_argv]):
        sys.exit(help_msg)

    helconf = amppart = ampindex = loopid = None
    for i in range(2, len(sys_argv)):
        arg = sys_argv[i]
        if helconf_regex.findall(arg) != []:
            helconf = helconf_regex.findall(arg)[0].replace("+", "p").replace("-", "m").replace(",", "").replace(" ", "")
        elif amppart_regex.findall(arg) != []:
            amppart = amppart_regex.findall(arg)[0]
        elif ampindex_regex.findall(arg) != []:
            ampindex = ampindex_regex.findall(arg)[0]
        elif loopid_regex.findall(arg) != []:
            loopid = loopid_regex.findall(arg)[0]
            if loopid == "None":
                loopid = None
        else:
            raise Exception(help_msg)

    if helconf is None:
        raise Exception("Helconf is a required parameter.")
    elif amppart is None:
        raise Exception("Amppart is a required parameter.")
    elif amppart in ["box", "triangle", "bubble"] and ampindex is None:
        raise Exception("Ampindex is a required parameter given the amppart: {}.".format(amppart))
    elif amppart in ["box", "triangle", "bubble", "rational"] and loopid is None:
        raise Exception("Loopid is a required paramrter given the amppart: {}.".format(amppart))

    return {'helconf': helconf, 'amppart': amppart, 'ampindex': ampindex, 'loopid': loopid}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


# Metrics and Pauli Matrices
LeviCivita = numpy.array([[0, 1], [-1, 0]])
MinkowskiMetric = numpy.diag([1, -1, -1, -1])
Pauli_zero = numpy.diag([1, 1])
Pauli_x = numpy.array([[0, 1], [1, 0]])
Pauli_y = numpy.array([[0, -1j], [1j, 0]])
Pauli_z = numpy.array([[1, 0], [0, -1]])
Pauli = numpy.array([Pauli_zero, Pauli_x, Pauli_y, Pauli_z])
Pauli_bar = numpy.array([Pauli_zero, -Pauli_x, -Pauli_y, -Pauli_z])

# Lorentz Invariants String Patterns
pSijk = re.compile(r'^(?:s_|S_)(\d+)$')
pd5 = re.compile(r'^δ5$')
ptr5 = re.compile(r'^(?:tr5_)(\d+)$')
pDijk = re.compile(r'(?:^Δ_(\d+)$)|(?:^Δ_(\d+)\|(\d+)\|(\d+)$)')
pOijk = re.compile(r'^(?:Ω_)(\d+)$')
pPijk = re.compile(r'^(?:Π_)(\d+)$')
pA2 = re.compile(r'^(?:⟨)([0-9])(?:\|)([0-9])(?:⟩)$')
pS2 = re.compile(r'^(?:\[)([0-9])(?:\|)([0-9])(?:\])$')
p3B = re.compile(r'^(?:⟨|\[)(\d)(?:\|\({0,1})([\d[\+|-]*]*)(?:\){0,1}\|)(\d)(?:⟩|\])$')
pNB = re.compile(r'^(?:⟨|\[)(?P<start>\d)(?:\|)(?P<middle>(?:(?:\([\d[\+|-]{1,}]{,1}\))|(?:[\d[\+|-]{1,}]{,1}))*)(?:\|)(?P<end>\d)(?:⟩|\])$')

# BlackHat variables
if BH_found:
    Lambda = BH.LambdaRGMP
    Lambdat = BH.LambdatRGMP
    C_Mom = BH.Cmomgmp
    Mom_Conf = BH.mcgmp


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def rand_frac():
    return mpmath.mpc(random.randrange(-100, 101)) / mpmath.mpc(random.randrange(1, 201))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def cgmp_to_mpc(cgmp_nbr):
    return mpmath.mpc(str(cgmp_nbr.real), str(cgmp_nbr.imag))


def mpc_to_cgmp(mpc_nbr):
    import gmpTools
    return gmpTools.CGMP(str(mpc_nbr.real), str(mpc_nbr.imag))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


class MyShelf(object):     # context manager shelf

    def __init__(self, filename, flag):
        self.filename = filename
        self.flag = flag

    def __enter__(self):
        if self.flag == "c":  # make sure path exists
            folder_path = "/".join(self.filename.split("/")[:-1])
            if not os.path.exists(folder_path):
                os.makedirs(folder_path)
        self.obj = shelve.open(self.filename, flag=self.flag)
        return self.obj

    def __exit__(self, exc_type, exc_value, traceback):
        self.obj.close()


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def retry(ExceptionsToCheck, max_tries=2, silent=False):
    def deco_retry(func):
        @functools.wraps(func)
        def f_retry(*args, **kwargs):
            current_try = 0
            while current_try < max_tries:
                try:
                    return_value = func(*args, **kwargs)
                except ExceptionsToCheck as e:
                    last_exception = e
                    if silent is False:
                        print("\r{}: {}...  ".format(type(e).__name__, str(e)[:75]), end="")
                        if type(args[-1]) in [list, tuple]:
                            print(", ".join(args[-1]), end="")
                        else:
                            print("Last arg: {}".format(args[-1]), end="")
                    current_try += 1
                else:
                    if current_try != 0 and silent is False:
                        print("\r{}: {}... ".format(type(last_exception).__name__, str(last_exception)[:75]), end="")
                        if type(args[-1]) in [list, tuple]:
                            print(", ".join(args[-1]), end="")
                        else:
                            print("Last arg: {}".format(args[-1]), end="")
                        print(" ~ Corrected on try {}!                                                               ".format(current_try + 1))
                    return return_value
            if silent is False:
                print("\r{}: {}... ".format(type(last_exception).__name__, str(last_exception)[:75]), end="")
                if type(args[-1]) in [list, tuple]:
                    print(", ".join(args[-1]), end="")
                else:
                    print("Last arg: {}".format(args[-1]), end="")
                print(" ~ Uncorrected after {} tries. :(                                                           ".format(max_tries))
            return None
        return f_retry
    return deco_retry


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


class fakeValue(object):

    def __init__(self, type_, init_value):
        if type_ == 'f':
            self.value = float(init_value)
        elif type_ == 'i':
            self.value = int(init_value)


class fakeManager(object):

    def __init__(self):
        pass

    @staticmethod
    def Value(type_, init_value):
        return fakeValue(type_, init_value)


class Progress(object):

    def __init__(self, maximum, UseParallelisation):
        if UseParallelisation is True:
            self.manager = multiprocessing.Manager()
        else:
            self.manager = fakeManager()
        self.time = self.manager.Value('f', 0)
        self.average_time = self.manager.Value('f', 0)
        self.current = self.manager.Value('i', -1)
        self.maximum = maximum
        self.hours = self.manager.Value('i', 0)
        self.minutes = self.manager.Value('i', 0)
        self.seconds = self.manager.Value('i', 0)

    def increment(self):
        self.current.value += 1

    def decrement(self):
        self.current.value -= 1

    def write(self):
        if self.current.value >= 1:
            previous_percentage = float(math.ceil(float(self.current.value - 1) / self.maximum * 10000)) / 100
        current_percentage = float(math.ceil(float(self.current.value) / self.maximum * 10000)) / 100
        previous_time = self.time.value
        self.time.value = time.time()
        if current_percentage != 0:
            msg = "{0:.2f}% completed.".format(current_percentage)
            left_percentage = 100 - current_percentage
            step_percentage = current_percentage - previous_percentage
            if step_percentage == 0:   # not pretty but avoids division by zero
                step_percentage = 0.01
            nbr_left_steps = left_percentage / step_percentage
            self.average_time.value = (
                (self.average_time.value * (self.current.value - 1) +
                 (self.time.value - previous_time)) / self.current.value)
            seconds = int((self.average_time.value) * nbr_left_steps)
            minutes, seconds = divmod(seconds, 60)
            hours, minutes = divmod(minutes, 60)
            delta_time = ((hours - self.hours.value) * 3600 +
                          (minutes - self.minutes.value) * 60 +
                          (seconds - self.seconds.value))
            if delta_time > 5 or delta_time < 0:
                ETA = "%d:%02d:%02d" % (hours, minutes, seconds)
                self.hours.value = hours
                self.minutes.value = minutes
                self.seconds.value = seconds
            else:
                ETA = "%d:%02d:%02d" % (self.hours.value, self.minutes.value,
                                        self.seconds.value)
            msg += " ETA: {}.".format(ETA)
        else:
            msg = "{0:.2f}% completed.".format(current_percentage)
        return msg


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


try:
    from pytest_cov.embed import cleanup_on_sigterm
except ImportError:
    pass
else:
    cleanup_on_sigterm()


class NoDaemonProcess(multiprocessing.Process):

    def _get_daemon(self):  # make 'daemon' attribute always return False
        return False

    def _set_daemon(self, value):
        pass

    daemon = property(_get_daemon, _set_daemon)


class NoDaemonProcessPool(multiprocessing.pool.Pool):
    # see https://stackoverflow.com/questions/52948447/error-group-argument-must-be-none-for-now-in-multiprocessing-pool
    def Process(self, *args, **kwargs):
        proc = super(NoDaemonProcessPool, self).Process(*args, **kwargs)
        proc.__class__ = NoDaemonProcess
        return proc


class MyProcessPool(object):      # context manager pool (spawns new processes, sidesteps the GIL)

    def __init__(self, processes=1, initializer=None, initargs=None):
        self.processes = processes
        self.initializer = initializer
        self.initargs = initargs

    def __enter__(self):
        self.obj = NoDaemonProcessPool(self.processes, self.initializer, self.initargs)
        # self.obj = multiprocessing.Pool(self.processes, self.initializer, self.initargs)
        return self.obj

    def __exit__(self, exc_type, exc_value, traceback):
        self.obj.close()
        self.obj.join()
        self.obj.terminate()


class MyThreadPool(object):      # context manager pool (no new processes, currently affected by GIL)

    def __init__(self, processes=1, initializer=None, initargs=None):
        self.processes = processes
        self.initializer = initializer
        self.initargs = initargs

    def __enter__(self):
        # self.obj = NoDaemonThreadPool(self.processes, self.initializer, self.initargs)
        self.obj = multiprocessing.pool.ThreadPool(self.processes, self.initializer, self.initargs)
        return self.obj

    def __exit__(self, exc_type, exc_value, traceback):
        self.obj.close()
        self.obj.join()
        self.obj.terminate()


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def progress_wrapper(func, *args, **kwargs):
    global lock, prog
    with lock if lock is not None else nullcontext():
        prog.increment()
        if type(args[-1]) in [list, tuple] and type(args[-1][0]) in [str, unicode] and len(", ".join(args[-1])) < 40:
            print("\r{} {} working on {}.                                   ".format(func.func.__name__ if hasattr(func.func, "__name__") else "function",
                                                                                     prog.write(), ", ".join(args[-1])), end="\r")
        elif len(unicode(args[-1]).replace("\n", "")) < 40:
            print("\r{} {} working on {}.                                   ".format(func.func.__name__ if hasattr(func.func, "__name__") else "function",
                                                                                     prog.write(), unicode(args[-1]).replace("\n", "")), end="\r")
        else:
            print("\r{} {}.                                                 ".format(func.func.__name__ if hasattr(func.func, "__name__") else "function",
                                                                                     prog.write()), end="\r")
        sys.stdout.flush()
    return func(*args, **kwargs)


_lambda_compatible_func = None


def worker(x):
    return _lambda_compatible_func(x)


def _init(l, p, func, ):
    global lock, prog, _lambda_compatible_func
    lock = l
    prog = p
    _lambda_compatible_func = func


def mapThreads(func, *args, **kwargs):
    # Default keyword arguments
    UseParallelisation = kwargs.pop("UseParallelisation", True)
    ParallelisationType = kwargs.pop("ParallelisationType", ('Thread', 'Process')[1])
    Cores = kwargs.pop("Cores", os.cpu_count() // 2)
    verbose = kwargs.pop("verbose", True)

    # map is applied to the last argument
    _func_partial = functools.partial(func, *args[:-1], **kwargs)

    if verbose:
        func_partial = functools.partial(progress_wrapper, _func_partial)
    else:
        func_partial = _func_partial

    if UseParallelisation is True:
        l = multiprocessing.Lock()
        p = Progress(len(args[-1]), UseParallelisation)
        if ParallelisationType == 'Process' or ParallelisationType == 'Processing':
            with MyProcessPool(Cores, initializer=_init, initargs=(l, p, func_partial,)) as pool:
                results = pool.map(worker, args[-1])
        elif ParallelisationType == 'Thread' or ParallelisationType == 'Threading':
            with MyThreadPool(Cores, initializer=_init, initargs=(l, p, func_partial,)) as pool:
                results = pool.map(worker, args[-1])
        else:
            raise ValueError("ParallelisationType must be either Process(ing) or Thread(ing).")
    else:
        global prog, lock
        lock = None
        prog = Progress(len(args[-1]), UseParallelisation)
        results = list(map(func_partial, args[-1]))

    if verbose:
        print("\r                                                                                                                   ", end="\r")
        sys.stdout.flush()
    return results


def filterThreads(lambda_func, iterable):
    lambda_func.__name__ = str("filterThreads")
    TrueOrFalseList = mapThreads(lambda_func, iterable)
    iterable = [entry for i, entry in enumerate(iterable) if TrueOrFalseList[i] is True]
    return iterable


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


class nullcontext(object):

    def __init__(self, enter_result=None):
        self.enter_result = enter_result

    def __enter__(self):
        return self.enter_result

    def __exit__(self, *excinfo):
        pass


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def output_needs_grabbing():
    return "ipykernel.iostream.OutStream" not in sys.stdout.__str__() if hasattr(sys.stdout, "__str__") else False


class OutputGrabber(object):
    """ Class used to grab standard output or another stream. """
    escape_char = "\b"

    def __init__(self, stream=None, threaded=False):
        if output_needs_grabbing():
            self.origstream = stream
            self.threaded = threaded
            if self.origstream is None:
                self.origstream = sys.stdout
            self.origstreamfd = self.origstream.fileno()
            self.capturedtext = ""
            # Create a pipe so the stream can be captured:
            self.pipe_out, self.pipe_in = os.pipe()
        else:
            pass

    def __enter__(self):
        if output_needs_grabbing():
            self.start()
            return self
        else:
            return None

    def __exit__(self, type, value, traceback):
        if output_needs_grabbing():
            self.stop()
        else:
            pass

    def start(self):
        """
        Start capturing the stream data.
        """
        self.capturedtext = ""
        # Save a copy of the stream:
        self.streamfd = os.dup(self.origstreamfd)
        # Replace the original stream with our write pipe:
        os.dup2(self.pipe_in, self.origstreamfd)
        if self.threaded:
            # Start thread that will read the stream:
            self.workerThread = threading.Thread(target=self.readOutput)
            self.workerThread.start()
            # Make sure that the thread is running and os.read() has executed:
            time.sleep(0.01)

    def stop(self):
        """
        Stop capturing the stream data and save the text in `capturedtext`.
        """
        # Print the escape character to make the readOutput method stop:
        self.origstream.write(self.escape_char)
        # Flush the stream to make sure all our data goes in before
        # the escape character:
        self.origstream.flush()
        if self.threaded:
            # wait until the thread finishes so we are sure that
            # we have until the last character:
            self.workerThread.join()
        else:
            self.readOutput()
        # Close the pipe:
        os.close(self.pipe_in)
        os.close(self.pipe_out)
        # Restore the original stream:
        os.dup2(self.streamfd, self.origstreamfd)
        # Close the duplicate stream:
        os.close(self.streamfd)

    def readOutput(self):
        """
        Read the stream data (one byte at a time)
        and save the text in `capturedtext`.
        """
        while True:
            char = os.read(self.pipe_out, 1)
            if not char or self.escape_char in char:
                break
            self.capturedtext += char


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


class log_linear_fit_Exception(Exception):
    pass


def log_linear_fit(xaxis, yaxis):
    """Fit slope of straight line in a log-log scatter plot. Slope has to be integer or half integer."""

    if not all([entry > 0 for entry in xaxis]):
        raise log_linear_fit_Exception("Negative entry in xaxis: {}".format(xaxis))
    if not all([entry > 0 for entry in yaxis]):
        raise log_linear_fit_Exception("Negative entry in yaxis: {}".format(yaxis))

    xaxis, yaxis = list(map(mpmath.log, xaxis)), list(map(mpmath.log, yaxis))
    slopes = [(y2 - y1) / (x2 - x1) for x1, x2, y1, y2 in zip(xaxis, xaxis[1:], yaxis, yaxis[1:])]
    slopes = list(map(float, slopes))
    rounded_slopes = [b / 2 for b in map(round, [2 * a for a in slopes])]

    if not abs(slopes[-1] - rounded_slopes[-1]) < 10 ** -6:
        raise log_linear_fit_Exception("Log linear fit harsh rounding: {} rounded to {}.".format(slopes[-1], rounded_slopes[-1]))

    return rounded_slopes[-1]


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def Write(lBrackets, lExponents=[], lCoefficients=[]):
    if lExponents == []:
        lExponents = [1 for brk in lBrackets]
    if lCoefficients == [] and lBrackets != []:
        lCoefficients = [1]
    # a few checks on the input
    if (not all(isinstance(item, list) for item in lBrackets)):
        lBrackets = [lBrackets]
    if (not all(isinstance(item, list) for item in lExponents)):
        lExponents = [lExponents]
    if (not all(len(brks) == len(exps) for (brks, exps) in zip(lBrackets, lExponents))):
        raise Exception("Brackets and exponents length mismatch.")
    if (not (len(lBrackets) == len(lExponents) and len(lBrackets) == len(lCoefficients))):
        raise Exception("lBrackets and lExp.s and lCoef.s length mismatch. {}, {}, {}".format(lBrackets, lExponents, lCoefficients))
    # proper write function
    final_string = ""
    for i, (brackets, exponents, coefficient) in enumerate(zip(lBrackets, lExponents, lCoefficients)):
        # write each term in a new line
        if i >= 1:
            final_string += "\n"
        # add the coefficient
        if coefficient == 1 and i == 0:
            pass
        elif coefficient == 1:
            final_string += "+"
        elif coefficient[0] == 0:
            final_string += "{}I".format(str(coefficient[1]))
        elif coefficient[1] == 0:
            final_string += "{}".format(str(coefficient[1]))
        else:
            final_string += "{}+{}I".format(str(coefficient[0]), str(coefficient[1]))
        # single fraction string construction
        _string = "("
        at_leat_one_in_numerator = False
        for i in range(len(brackets)):
            if exponents[i] > 0:
                at_leat_one_in_numerator = True
                if exponents[i] == 1:
                    _string += brackets[i] + " "
                else:
                    _string += brackets[i] + "^" + str(exponents[i]) + " "
        if at_leat_one_in_numerator is False:
            _string += "1 "
        _string = _string[:-1]
        _string += ")/("
        for i in range(len(brackets)):
            if exponents[i] < 0:
                if exponents[i] == -1:
                    _string += brackets[i] + " "
                else:
                    _string += brackets[i] + "^" + str(abs(exponents[i])) + " "
        _string = _string[:-1]
        _string += ")"
        final_string += _string
    return final_string


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def PartialCompute(lBrackets, lExponents=[], lCoefficients=[1]):
    return functools.partial(Compute, lBrackets, lExponents, lCoefficients)


def Compute(lBrackets, lExponents, lCoefficients, oParticles):
    if lExponents == []:
        lExponents = [1 for brk in lBrackets]
    if lCoefficients == []:
        lCoefficients = [1]
    # a few checks on the inputs
    if not all(isinstance(item, list) for item in lBrackets):
        lBrackets = [lBrackets]
    if not all(isinstance(item, list) for item in lExponents):
        lExponents = [lExponents]
    if not all(len(brks) == len(exps) for (brks, exps) in zip(lBrackets, lExponents)):
        raise Exception("Brackets and exponents length mismatch.")
    if not (len(lBrackets) == len(lExponents) and len(lBrackets) == len(lCoefficients)):
        raise Exception("lBrackets and lExp.s and lCoef.s length mismatch")
    # proper compute function
    final_result = 0
    for i, (brackets, exponents, coefficient) in enumerate(
            zip(lBrackets, lExponents, lCoefficients)):
        _result = 1
        for j, (bracket, exponent) in enumerate(zip(brackets, exponents)):
            if exponent != 0:
                _result = _result * oParticles.compute(bracket) ** exponent
        # if the coefficient is in fraction format needs to be unpacked
        if type(coefficient) is tuple:
            real_fraction = coefficient[0]
            imag_fraction = coefficient[1]
            real = real_fraction.numerator / real_fraction.denominator
            imag = imag_fraction.numerator / imag_fraction.denominator
            coefficient = real + 1j * imag
        final_result = final_result + coefficient * _result
    return final_result


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def Generate_LaTeX_and_PDF(message, res_path, partial=False, compile_tex_to_pdf=True):
    if partial is False:
        print("\nWriting LaTeX and PDF result file.")
    else:
        print("Writing LaTeX and PDF partial result file.\n")
    if partial is False:
        res_path += ".tex"
    else:
        res_path += "_partial.tex"
    os.makedirs(os.path.dirname(res_path), exist_ok=True)
    with open(res_path, "w") as oFile:
        oFile.write("\\documentclass[varwidth, border=5pt]{standalone}\n")
        oFile.write("\\usepackage[paperwidth=575cm, paperheight=575cm, margin=1in, landscape]{geometry}\n")
        oFile.write("\\DeclareMathSizes{1}{1}{1}{0.1}\n")
        oFile.write("\\usepackage{unicode-math}\n")
        oFile.write("\\usepackage{mathtools}\n")
        oFile.write("\\standaloneenv{my}\n\n")
        oFile.write("\\begin{document}\n")
        if sys.version_info[0] > 2:
            oFile.write(message)
        else:
            oFile.write(message.encode('utf8'))
        oFile.write("\\end{document}\n")
    if compile_tex_to_pdf:
        with open(os.devnull, 'wb') as devnull:
            subprocess.check_call([MainPythonDirectory + '/scripts/SpinorLatexCompiler.py', res_path], stdout=devnull, stderr=devnull)
        if partial is False:
            if os.path.isfile(res_path[:-4] + "_partial.tex"):
                os.system("rm " + res_path[:-4].replace("(", "(").replace(")", ")") + "_partial.tex")
            if os.path.isfile(res_path[:-4] + "_partial.pdf"):
                os.system("rm " + res_path[:-4].replace("(", "(").replace(")", ")") + "_partial.pdf")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def generate_BH_cpp_files(shelconf, PWRespath, convert_cut=True, convert_rational=True):
    from terms.terms import LoadResults
    from lips import Particles
    from core.unknown import Unknown
    import topologies.topology as topology
    multiplicity = len(shelconf)

    if convert_cut is True:
        # cut load templates
        with open(MainPythonDirectory + "/../C++/template_loop.cpp", "r") as templatefile:
            loop_BH_cpp_file = templatefile.read()
        with open(MainPythonDirectory + "/../C++/template_loop.hpp", "r") as templatefile:
            loop_BH_hpp_file = templatefile.read()
        with open(MainPythonDirectory + "/../C++/template_coefficient.cpp", "r") as templatefile:
            coefficient_BH_cpp_file = templatefile.read()
        with open(MainPythonDirectory + "/../C++/template_coefficient-impl.cpp", "r") as templatefile:
            coefficient_impl_BH_cpp_file = templatefile.read()
    if convert_rational is True:
        # rational load templates
        with open(MainPythonDirectory + "/../C++/template_rational.cpp", "r") as templatefile:
            rational_BH_cpp_file = templatefile.read()

    # BH handles
    p, m = BH.cvar.p, BH.cvar.m
    helconf = BH.process(*[p if character == "p" else m for character in shelconf])
    A1 = BH.One_Loop_Helicity_Amplitude(helconf, BH.glue)
    if convert_cut is True:
        with OutputGrabber():
            A1_cutpart = A1.cut_part().makeDarrenCutPart()
        nbr_boxes = int(A1_cutpart.nbr_boxes())
        nbr_triangles = int(A1_cutpart.nbr_triangles())
        nbr_bubbles = int(A1_cutpart.nbr_bubbles())
        tot_nbr_parts = nbr_boxes + nbr_triangles + nbr_bubbles

    # {0} for template_loop.cpp; {1} for template_rational.cpp  ---  Amplitude Name
    if convert_cut is True:
        NameAmpl = "{}g_{}_G".format(multiplicity, shelconf)
    if convert_rational is True:
        RNameAmpl = "{}g{}".format(multiplicity, int(shelconf.replace("p", "1").replace("m", "0")[::-1], 2))

    if convert_cut is True:
        # {1} for template_loop.cpp  ---  Corners
        allcorners = []
        for i in range(multiplicity - 2):
            allcorners += ["".join([str((j + k) % multiplicity + 1) for k in range(i + 1)]) for j in range(multiplicity)]
        printCorners = ""
        for corner in allcorners:
            printCorners += "         vector<int> c" + corner + "; "
            for i in corner:
                printCorners += "c" + corner + ".push_back(ind.at(" + str(int(i) - 1) + ")); "
            printCorners += "\n"

        # {2} for template_loop.cpp  ---  Integrals
        printIntegrals = ''
        for ii in range(1, nbr_boxes + 1):
            code = str(eval('A1_cutpart.{}({}).get_code()'.format('box', ii)))
            codelista = sorted([int(code[i]) for i in range(len(code))])
            codelist = codelista + [multiplicity + codelista[0]]
            codelistgr = ['c' + (''.join([str((codelist[i] + j - 1) % multiplicity + 1) for j in range(codelist[i + 1] - codelist[i])])) for i in range(len(code))]
            printIntegrals += ('CI_users.push_back(new Cached_Box_Integral_User(' + codelistgr[0] + ', ' + codelistgr[1] + ', ' + codelistgr[2] + ', ' + codelistgr[3] + '));\n')
        for ii in range(1, nbr_triangles + 1):
            code = str(eval(
                'A1_cutpart.{}({}).get_code()'.format('triangle', ii)))
            codelista = sorted([int(code[i]) for i in range(len(code))])
            codelist = codelista + [multiplicity + codelista[0]]
            codelistgr = ['c' + (''.join([str((codelist[i] + j - 1) % multiplicity + 1) for j in range(codelist[i + 1] - codelist[i])])) for i in range(len(code))]
            printIntegrals += ('CI_users.push_back(new Cached_Triangle_Integral_User(' + codelistgr[0] + ', ' + codelistgr[1] + ', ' + codelistgr[2] + '));\n')
        for ii in range(1, nbr_bubbles + 1):
            code = str(eval('A1_cutpart.{}({}).get_code()'.format('bubble', ii)))
            codelista = sorted([int(code[i]) for i in range(len(code))])
            codelist = codelista + [multiplicity + codelista[0]]
            codelistgr = ['c' + (''.join([str((codelist[i] + j - 1) % multiplicity + 1) for j in range(codelist[i + 1] - codelist[i])])) for i in range(len(code))]
            printIntegrals += ('CI_users.push_back(new Cached_Bubble_Integral_User(' + codelistgr[0] + ', ' + codelistgr[1] + '));\n')

    # {4} for template_loop.cpp; {3} for template_rational.cpp  ---  Coefficients
    all_invariants = {}
    all_coefficients = {}
    if convert_cut is True:
        coefficients_declarations = ""
        lZeroCoeffs = []
        lCoefficients = []
        lTopos = topology.independent_topologies(shelconf)
        for i, topos in enumerate(lTopos):
            for topo in topos:
                print("\rReading topology {} ".format(topo), end="")
                coefficients_declarations += '\n // Topology: {} \n'.format(topo)
                amppart, ampindices = topo.models
                print("({}({})).              ".format(amppart, ampindices[0]), end="")
                oCutModel = LoadResults(shelconf, amppart, ampindices[0], load_partial_results_only=False, check_loaded_results=True, silent=True)
                for ampindex in ampindices:
                    if amppart == "box":
                        coef_index = ampindex
                    elif amppart == "triangle":
                        coef_index = ampindex + nbr_boxes
                    else:
                        coef_index = ampindex + nbr_boxes + nbr_triangles
                    if oCutModel[1] is True:
                        raise Exception("Loaded result is partial while printing BH cpp files.")
                    elif oCutModel[0] == 0:
                        lZeroCoeffs += [coef_index]
                    else:
                        if ampindex == ampindices[0]:
                            oCut = oCutModel[0][0]
                            oCut = oCut.explicit_representation()
                        else:
                            conversion_rule = topology.conversion_rules(topology.get_label(shelconf, amppart, ampindices[0]), topology.get_label(shelconf, amppart, ampindex))[0]
                            oCut = oCutModel[0][0].Image(conversion_rule)
                            oCut = oCut.explicit_representation()
                        # check converted result
                        oParticles = Particles(len(shelconf))
                        oParticles.fix_mom_cons()
                        oNumericalCut = Unknown(shelconf, amppart, ampindex)
                        if str(oNumericalCut(oParticles)) != str(oCut(oParticles)):
                            raise Exception("Conversion failed for {}({})".format(amppart, ampindex))
                        # end of check
                        all_invariants = set(list(all_invariants) + list(oCut.invs_set))
                        all_coefficients = set(list(all_coefficients) + flatten(oCut.llCoefs))
                        coefficients_declarations += '\ntemplate <class T> Complex C' + NameAmpl + '_co' + str(coef_index) + ' (C' + NameAmpl + '_mom_conf_info<T>& mci);\n'
                        lCoefficients += ['{}template <class T> std::complex<T> C{NameAmpl}_co{CoefIndex} (C{NameAmpl}_mom_conf_info<T>& mci) [[return {};]]'.format(
                            *printBHcpp(oCut, NameAmpl, coef_index), NameAmpl=NameAmpl, CoefIndex=coef_index)]
            print("\r                                                                           \r", end="")

    # rational
    if convert_rational is True:
        oRational = LoadResults(shelconf, "rational", load_partial_results_only=False, check_loaded_results=True, silent=True)
        if oRational[1] is True:
            raise Exception("Loaded rational result is partial while printing BH cpp files.")
        elif oRational[0] == 0:
            raise Exception("Is the rational really zero?")
        else:
            oRational = oRational[0][0]
            oRational = oRational.explicit_representation()
            all_invariants = set(list(all_invariants) + list(oRational.invs_set))
            all_coefficients = set(list(all_coefficients) + flatten(oRational.llCoefs))
            printRational = '{}template <class T> std::complex<T> R{NameAmpl} (const eval_param<T>& ep,const mass_param_coll& mpc) '\
                '[[R{NameAmpl}_mom_conf_info<T> epi(ep); return {};]]'.format(
                    *printBHcpp(oRational.explicit_representation(), RNameAmpl), NameAmpl=RNameAmpl).replace("mci", "epi")

    # {3} for template_loop.cpp; {0} for template_rational.cpp  ---  Invariants
    printInvariants = ""
    # cache the coefficients in here as well
    all_numerators = set(flatten([[abs(complex_number[1].numerator) % complex_number[1].denominator,
                                   abs(complex_number[1].numerator) / complex_number[1].denominator] for complex_number in all_coefficients]))
    all_denominators = set([complex_number[1].denominator for complex_number in all_coefficients])
    for numerator in all_numerators:
        if numerator == 0:
            continue
        printInvariants += 'Complex n' + str(numerator) + 'i = Complex(0, T(' + str(numerator) + ')); '
    for denominator in all_denominators:
        printInvariants += 'Complex n_' + str(denominator) + ' = Complex(T(1)/T(' + str(denominator) + '), 0); '
    printInvariants = printInvariants[:-1] + '\n'
    # proper invariants
    for submult in range(1, multiplicity):   # SPA & SPB & Sij
        for i in range(submult + 1, multiplicity + 1):
            printInvariants += ('Complex spa' + str(submult) + str(i) + ' = SPA(' + str(submult) + ',' + str(i) + '); ')
            printInvariants += ('Complex spa' + str(i) + str(submult) + ' = -spa' + str(submult) + str(i) + '; ')
            printInvariants += ('Complex spb' + str(submult) + str(i) + ' = SPB(' + str(submult) + ',' + str(i) + '); ')
            printInvariants += ('Complex spb' + str(i) + str(submult) + ' = -spb' + str(submult) + str(i) + '; ')
            printInvariants += ('Complex s' + str(submult) + str(i) + ' = S(' + str(submult) + ',' + str(i) + '); ')
            printInvariants += ('Complex s' + str(i) + str(submult) + ' = s' + str(submult) + str(i) + ';\n')
            # powers for angle and square brackets
            for j in range(2, 10):
                printInvariants += ('Complex spa' + str(submult) + str(i) + '_' + str(j) + ' = pow(spa' + str(submult) + str(i) + ',' + str(j) + '); ')
                printInvariants += ('Complex spa' + str(i) + str(submult) + '_' + str(j) + ' = pow(spa' + str(i) + str(submult) + ',' + str(j) + '); ')
                printInvariants += ('Complex spb' + str(submult) + str(i) + '_' + str(j) + ' = pow(spb' + str(submult) + str(i) + ',' + str(j) + '); ')
                printInvariants += ('Complex spb' + str(i) + str(submult) + '_' + str(j) + ' = pow(spb' + str(i) + str(submult) + ',' + str(j) + ');\n')
    for submult in range(1, multiplicity):   # Sijk
        for i in range(submult + 1, multiplicity + 1):
            for j in range(i + 1, multiplicity + 1):
                printInvariants += ('Complex s' + str(submult) + str(i) + str(j) + ' = SS(' + str(submult) + ',' + str(i) + ',' + str(j) + '); ')
                printInvariants += ('Complex s' + str(j) + str(submult) + str(i) + ' = s' + str(submult) + str(i) + str(j) + '; ')
                printInvariants += ('Complex s' + str(i) + str(j) + str(submult) + ' = s' + str(submult) + str(i) + str(j) + ';\n')
    for invariant in all_invariants:
        if pA2.findall(invariant) == [] and pS2.findall(invariant) == [] and pSijk.findall(invariant) == []:
            invariant_name, invariant_expr = None, None
            spab = re.compile(r'⟨(\d)\|\((\d)\+(\d)\)\|(\d)\]')
            spba = re.compile(r'\[(\d)\|\((\d)\+(\d)\)\|(\d)⟩')
            if spab.findall(invariant) != []:     # SPAB
                invariant_name = spab.sub(r'spab\1\2\3\4', invariant)
                invariant_expr = spab.sub(r'spa\1\2*spb\2\4+spa\1\3*spb\3\4', invariant)
            elif spba.findall(invariant) != []:   # SPBA
                invariant_name = spba.sub(r'spba\1\2\3\4', invariant)
                invariant_expr = spba.sub(r'spb\1\2*spa\2\4+spb\1\3*spa\3\4', invariant)
            elif pDijk.findall(invariant) != [] or pOijk.findall(invariant) != [] or pPijk.findall(invariant) != []:  # Delta's, Omega's, Pi's
                oParticles = Particles(multiplicity)
                if pDijk.findall(invariant) != []:
                    ijk = map(int, pDijk.findall(invariant)[0])
                    nol = oParticles.ijk_to_3NonOverlappingLists(ijk)
                    K12xK12 = "T(1)/T(4)*pow((s{}{}+s{}{}+s{}{}+s{}{}),2)".format(nol[0][0], nol[1][0], nol[0][0], nol[1][1], nol[0][1], nol[1][0], nol[0][1], nol[1][1])
                    K11xK22 = "s{}{}*s{}{}".format(nol[0][0], nol[0][1], nol[1][0], nol[1][1])
                    invariant_name = pDijk.sub(r'Delta\1', invariant)
                    invariant_expr = K12xK12 + "-" + K11xK22
                elif pOijk.findall(invariant) != []:
                    ijk = map(int, pOijk.findall(invariant)[0])
                    nol = oParticles.ijk_to_3NonOverlappingLists(ijk)
                    invariant_name = pOijk.sub(r'Omega\1', invariant)
                    invariant_expr = "T(2)*s{}{}*s{}{}".format(nol[2][0], nol[2][1], nol[1][0], nol[1][1])
                    invariant_expr = invariant_expr + "-" + "(s{}{}+s{}{}-s{}{})*s{}{}{}".format(*nol[2] + nol[1] + nol[0] + nol[2] + [nol[0][0]])
                elif pPijk.findall(invariant) != []:
                    ijk = map(int, pPijk.findall(invariant)[0])
                    nol = oParticles.ijk_to_3NonOverlappingLists(ijk)
                    invariant_name = pPijk.sub(r'Pi\1', invariant)
                    invariant_expr = "s{}{}{}-s{}{}{}".format(*nol[2] + [nol[0][0]] + nol[2] + [nol[0][1]])
            if invariant_name is None or invariant_expr is None:
                raise Exception("Invariant {} not understood for c++ compilation.".format(invariant))
            else:
                printInvariants += 'Complex {} = {};\n'.format(invariant_name, invariant_expr)
    printInvariants = printInvariants[:-1]
    declarations = "".join([txt.split('=')[0] + ';' for txt in printInvariants.split(';') if txt])
    initialisations = "\n".join([txt.replace("Complex ", "") for txt in printInvariants.split('\n') if txt])

    # {5} for template_loop.cpp  ---  Result: sum of coefficients * integrals
    if convert_cut is True:
        if 1 not in lZeroCoeffs:
            printResult2 = ('\nSeriesC<T> result = C' + NameAmpl + '_co1(mci)*(*CI_users[0]->' + 'get_value(mc,ind,mu))')
        for i in range(1, tot_nbr_parts - 1):
            if i + 1 not in lZeroCoeffs:
                printResult2 += (' + C' + NameAmpl + '_co' + str(i + 1) + '(mci)*(*CI_users[' + str(i) + ']->get_value(mc,ind,mu))')
        if tot_nbr_parts not in lZeroCoeffs:
            printResult2 += (' + C' + NameAmpl + '_co' + str(tot_nbr_parts) + '(mci)*(*CI_users[' + str(tot_nbr_parts - 1) + ']->get_value(mc,ind,mu));\n')

    if convert_cut is True:
        # wCI.cpp master file
        outfile_txt = loop_BH_cpp_file.format(AmplitudeName=NameAmpl, CornerDefinitions=printCorners, IntegralDefinitions=printIntegrals,
                                              SumOfCoefficientsTimesIntegrals=printResult2).replace('[[', '{').replace(']]', '}')
        with open(PWRespath + "/cpp/C" + NameAmpl + "_wCI.cpp", 'w+') as oFile:
            oFile.write(outfile_txt)
        # wCI.hpp header file
        outfile_txt = loop_BH_hpp_file.format(AmplitudeName=NameAmpl, InvariantsDeclarations=declarations, InvariantsInitialisations=initialisations,
                                              CoefficientsDeclarations=coefficients_declarations).replace('[[', '{').replace(']]', '}')
        with open(PWRespath + "/cpp/C" + NameAmpl + "_wCI.hpp", 'w+') as oFile:
            oFile.write(outfile_txt)
        # co#.hpp puppet files ~~~~ this also prints the files to be added in the split.inc file
        print("Split files:")
        print("C" + NameAmpl + "_wCI.cpp ", end="")
        for coef in lCoefficients:
            coef_index = int(re.findall(r"_co(\d*)", coef)[0])
            file_name = "C" + NameAmpl + "_co" + str(coef_index)
            if len(coef) < 10 ** 4:  # don't split it
                outfile_txt = coefficient_BH_cpp_file.format(CoefficientDefinition=coef, AmplitudeName=NameAmpl, CoefficientIndex=coef_index).replace('[[', '{').replace(']]', '}')
                with open(PWRespath + "/cpp/" + file_name + ".cpp", 'w+') as oFile:
                    oFile.write(outfile_txt)
                print(file_name + ".cpp ", end="")
            else:  # split it
                with open(PWRespath + "/cpp/" + file_name + ".cpp", 'w+') as oFile:
                    oFile.write(coef.replace('[[', '{').replace(']]', '}'))
                for _type in ["R", "RHP", "RVHP", "RGMP"]:
                    outfile_txt = coefficient_impl_BH_cpp_file.format(_type=_type, AmplitudeName=NameAmpl, CoefficientIndex=coef_index).replace('[[', '{').replace(']]', '}')
                    with open(PWRespath + "/cpp/" + file_name + "-impl-{}.cpp".format(_type), 'w+') as oFile:
                        oFile.write(outfile_txt)
                        print(file_name + "-impl-{}.cpp ".format(_type), end="")
        print("")

    if convert_rational is True:
        Rinits = re.sub(r"\((\d)\,(\d)\)", r"(\1-1,\2-1)", initialisations)
        Rinits = re.sub(r"\((\d)\,(\d)\,(\d)\)", r"(\1-1,\2-1,\3-1)", Rinits)
        outfile_txt = rational_BH_cpp_file.format(RNameAmpl.split("g")[0] + 'g', printRational=printRational, AmplitudeName=RNameAmpl, InvariantsDeclarations=declarations,
                                                  InvariantsInitialisations=Rinits).replace('[[', '{').replace(']]', '}')
        with open(PWRespath + '/cpp/' + 'R_' + RNameAmpl.split("g")[0] + 'g_eval_' + RNameAmpl.split("g")[1] + '-template.cpp', 'w+') as oFile:
            oFile.write(outfile_txt)

    return

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def printBHcpp(oTerms, NameAmpl, coef_index=None, size_limit=20):
    from terms.term import Numerator, Denominator
    if type(size_limit) is not int:
        raise Exception("size_limit is neither int nor None in printBHcpp")
    if coef_index is None:
        UniqueName = 'R' + NameAmpl
        NameAmpl = UniqueName
    else:
        UniqueName = 'C' + NameAmpl + '_co' + str(coef_index) + '_'
        NameAmpl = 'C' + NameAmpl
    splitted_terms, sum_of_splitted_terms = "", ""
    for i, iTerm in enumerate(oTerms):
        if iTerm.oNum.lCommonInvs != []:
            numerator_common = "template <class T> std::complex<T> " + UniqueName + "numerator_common_{0} ({NameAmpl}_mom_conf_info<T>& mci) [[return {1};]]\n".format(
                i, printBHcpp_inner(unicode(Numerator([iTerm.oNum.lCommonInvs], [iTerm.oNum.lCommonExps]))), NameAmpl=NameAmpl)
        denominator = "template <class T> std::complex<T> " + UniqueName + "denominator_{0} ({NameAmpl}_mom_conf_info<T>& mci) [[return {1};]]\n".format(
            i, printBHcpp_inner(unicode(Denominator(iTerm.oDen.lInvs, iTerm.oDen.lExps))), NameAmpl=NameAmpl)
        numerators = []
        for chunk in chunks(zip(iTerm.oNum.llInvs, iTerm.oNum.llExps, iTerm.oNum.lCoefs), size_limit):
            chunk_of_llInvs, chunk_of_llExps, chunk_of_lCoefs = zip(*chunk)
            numerators += [printBHcpp_inner(unicode(Numerator(chunk_of_llInvs, chunk_of_llExps, chunk_of_lCoefs)))]
        numerator = "".join(["template <class T> std::complex<T> " + UniqueName + "numerator_{0}_{1} ({NameAmpl}_mom_conf_info<T>& mci) [[return {2};]]\n".format(
            i, j, numerator, NameAmpl=NameAmpl) for j, numerator in enumerate(numerators)])
        if iTerm.oNum.lCommonInvs != []:
            term = numerator_common + denominator + numerator
            term += ("template <class T> std::complex<T> " + UniqueName +
                     "term_{0} ({NameAmpl}_mom_conf_info<T>& mci) [[return {UniqueName}numerator_common_{0}(mci)*({1})/{UniqueName}denominator_{0}(mci);]]\n\n".format(
                         i, "+".join([UniqueName + "numerator_{}_{}(mci)".format(i, j) for j in range(len(numerators))]), UniqueName=UniqueName, NameAmpl=NameAmpl))
        else:
            term = denominator + numerator
            if len(numerators) == 1:
                term += "template <class T> std::complex<T> {UniqueName}term_{0} ({NameAmpl}_mom_conf_info<T>& mci) [[return {1}/{UniqueName}denominator_{0}(mci);]]\n\n".format(
                    i, "+".join([UniqueName + "numerator_{}_{}(mci)".format(i, j) for j in range(len(numerators))]), UniqueName=UniqueName, NameAmpl=NameAmpl)
            else:
                term += ("template <class T> std::complex<T> {UniqueName}term_{0} ({NameAmpl}_mom_conf_info<T>& mci) [[return ({1})/{UniqueName}denominator_{0}(mci);]]\n\n".format(
                    i, "+".join([UniqueName + "numerator_{}_{}(mci)".format(i, j) for j in range(len(numerators))]), UniqueName=UniqueName, NameAmpl=NameAmpl))
        splitted_terms += term
    sum_of_splitted_terms = "+".join(UniqueName + "term_{}(mci)".format(i) for i in range(len(oTerms)))
    return splitted_terms, sum_of_splitted_terms


def printBHcpp_inner(BHcpp):
    BHcpp = re.sub(r'\n$', r'', BHcpp)
    BHcpp = re.sub(r'⟨(\d)\|(\d)⟩', r'mci.spa\1\2', BHcpp)
    BHcpp = re.sub(r'⟨(\d)\|\((\d)\+(\d)\)\|(\d)\]', r'mci.spab\1\2\3\4', BHcpp)
    BHcpp = re.sub(r'\[(\d)\|\((\d)\+(\d)\)\|(\d)⟩', r'mci.spba\1\2\3\4', BHcpp)
    BHcpp = re.sub(r'\[(\d)\|(\d)\]', r'mci.spb\1\2', BHcpp)
    BHcpp = re.sub(r'Δ_(\d)(\d)(\d)', r'mci.Delta\1\2\3', BHcpp)
    BHcpp = re.sub(r'Ω_(\d)(\d)(\d)', r'mci.Omega\1\2\3', BHcpp)
    BHcpp = re.sub(r'Π_(\d)(\d)(\d)', r'mci.Pi\1\2\3', BHcpp)
    BHcpp = re.sub(r's_(\d)(\d)(\d)', r'mci.s\1\2\3', BHcpp)
    BHcpp = re.sub(r's_(\d)(\d)', r'mci.s\1\2', BHcpp)
    BHcpp = re.sub(r'²', r'^2', BHcpp)
    BHcpp = re.sub(r'³', r'^3', BHcpp)
    BHcpp = re.sub(r'⁴', r'^4', BHcpp)
    BHcpp = re.sub(r'(?P<nbr>\d)(?P<letter>[A-HJ-Za-z])', r'\g<nbr>*\g<letter>', BHcpp)
    BHcpp = re.sub(r'(?P<parentesis>\))(?P<letter>[A-HJ-Za-z])', r'\g<parentesis>*\g<letter>', BHcpp)
    BHcpp = re.sub(r'\)\(', r')*(', BHcpp)
    BHcpp = re.sub(r'-(\d*)/(\d*)I', r'+Complex(0,-T(\1)/T(\2))*', BHcpp)
    BHcpp = re.sub(r'\+(\d*)/(\d*)I', r'+Complex(0,T(\1)/T(\2))*', BHcpp)
    BHcpp = re.sub(r'(\d*)/(\d*)I', r'+Complex(0,T(\1)/T(\2))*', BHcpp)
    BHcpp = re.sub(r'-(\d*)I', r'+Complex(0,-T(\1))*', BHcpp)
    BHcpp = re.sub(r'\+(\d*)I', r'+Complex(0,T(\1))*', BHcpp)
    BHcpp = re.sub(r'(\d*)I', r'+Complex(0,T(\1))*', BHcpp)
    BHcpp = re.sub(r'(?P<nbr>\d)\(', r'\g<nbr>*(', BHcpp)
    BHcpp = re.sub(r'\(\+', r'(', BHcpp)
    BHcpp = re.sub(r'/(?!T)', r')/(', BHcpp)
    BHcpp = re.sub(r'^', r'(', BHcpp)
    BHcpp = re.sub(r'\n', r')+\n(', BHcpp)
    BHcpp = re.sub(r'$', r')', BHcpp)
    BHcpp = re.sub(r'\(\+', r'(', BHcpp)
    BHcpp = re.sub(r'(mci.spa\d\d)\^(\d)', r'\1_\2', BHcpp)
    BHcpp = re.sub(r'(mci.spb\d\d)\^(\d)', r'\1_\2', BHcpp)
    BHcpp = re.sub(r'([0-9A-Za-z\.]*)\^(\d)', r'pow(\1,\2)', BHcpp)
    BHcpp = re.sub(r"T\((?P<num>\d*)\)/T\((?P<den>\d*)\)\)", fraction_to_proper_fraction, BHcpp)
    BHcpp = re.sub(r"T\((?P<num>\d*)\)\)", r"mci.n\g<num>i", BHcpp)
    BHcpp = BHcpp.replace("+Complex(0,-", "-")
    BHcpp = BHcpp.replace("+Complex(0,", "+")
    BHcpp = BHcpp.replace("Complex(0,", "")
    return BHcpp

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def fraction_to_proper_fraction(match):
    num = int(match.group("num"))
    den = int(match.group("den"))
    if num > den:
        integer_part = num / den
        proper_num = num % den
        result = "(mci.n{}i+mci.n{}i*mci.n_{})".format(integer_part, proper_num, den)
    else:
        result = "mci.n{}i*mci.n_{}".format(num, den)
    return result


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def LaTeXToPython(res_path, partial_only=False):
    from ..terms.terms import Terms

    res_path = str(res_path)

    partial = False
    if os.path.exists(res_path + ".tex") and partial_only is False:
        with open(res_path + ".tex", 'r') as content_file:
            latex_result = content_file.read()
    elif os.path.exists(res_path + "_partial.tex"):
        partial = True
        with open(res_path + "_partial.tex", 'r') as content_file:
            latex_result = content_file.read()
    else:
        return None, None

    latex_result = latex_result.replace(" ", "")
    latex_terms_pattern_selector = re.compile(r'\\begin{my}\n\$\\begin{gathered}(.*?)\\end{gathered}\$\n\\end{my}\n', re.DOTALL)
    python_results = latex_terms_pattern_selector.findall(latex_result)
    for i, python_result in enumerate(python_results):
        python_result = re.sub(r"\n$", r"", python_result)
        python_result = re.sub(r"^\n", r"", python_result)
        python_result = python_result.replace("\\scriptscriptstyle", "").replace("\\phantom{+}", "").replace("+\\\\", "+")
        if python_result == r"\text{Therequestedquantityisidenticallyzero.}":
            return 0, partial
        python_result = python_result.replace("\\langle", "⟨").replace("\\rangle", "⟩").replace("\\frac{", "(")
        python_result = re.sub(r"{([\d\|]+)}", r"\1", python_result)
        python_result = re.sub(r"{(\d).(\d)}", r"\1.\2", python_result)
        python_result = python_result.replace("}{", ")\\(").replace("}+", ")+")
        if python_result[-1] == "}":
            python_result = python_result[:-1] + ")"
        python_result = python_result.replace("}", "").replace("{", "")
        python_result = python_result.replace(")+\n", ")\n+")
        python_result = re.sub(r"\|(?P<nbr>\d)\+", r"|(\g<nbr>+", python_result)
        python_result = re.sub(r"\|(?P<nbr>\d)\-", r"|(\g<nbr>-", python_result)
        python_result = re.sub(r"\+(?P<nbr>\d)\|", r"+\g<nbr>)|", python_result)
        python_result = re.sub(r"\-(?P<nbr>\d)\|", r"-\g<nbr>)|", python_result)
        python_result = python_result.replace("\\", "/")
        python_result = re.sub(r"(\d+),/;/text", r"'\1',", python_result)
        python_result = python_result.replace(",/;-", ",'-'").replace(",/;+", ",'+'")
        python_result = python_result.replace("False)", "False,'+')").replace("True)", "True,'+')")
        python_result = re.sub(r"⟨(\d)(\d)⟩", r"⟨\1|\2⟩", python_result)
        python_result = re.sub(r"\[(\d)(\d)\]", r"[\1|\2]", python_result)
        python_result = python_result.replace("/Delta", "Δ").replace("/Omega", "Ω").replace("/Pi", "Π")
        # print(python_result, end="\n\n")
        python_results[i] = Terms(python_result)
    return python_results, partial


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def Sum_Abs(lIntegers):
    return sum([abs(Int) for Int in lIntegers])


def Positive_Entries(lIntegers):
    return sum([1 for entry in lIntegers if entry > 0])


def Negative_Entries(lIntegers):
    return sum([1 for entry in lIntegers if entry < 0])


def get_max_denominator(iterable):
    return max([entry.denominator for entry in iterable])


def get_max_abs_numerator(iterable):
    return max([abs(entry).numerator for entry in iterable])


def chained_gcd(array, outlier_fraction=0):
    """GCD of array, exluding outliers from both ends of the sorted array."""
    number_of_outliers = int(len(array) * outlier_fraction)
    if number_of_outliers > 0:
        array_without_outliers = sorted(array)[number_of_outliers:-number_of_outliers]
    else:
        array_without_outliers = array
    return functools.reduce(math.gcd, array_without_outliers)


def get_common_Q_factor(array, outlier_fraction=0):
    numer_gcd = chained_gcd([entry.numerator for entry in array.flatten() if entry != 0], outlier_fraction)
    denom_gcd = chained_gcd([entry.denominator for entry in array.flatten() if entry != 0], outlier_fraction)
    return fractions.Fraction(numer_gcd, denom_gcd)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


# order first list wrt number of forbidden pairs with second one
def forbidden_ordering(list_one, list_two, oUnknown):
    list_one.sort(key=lambda iInv: sum(
        [1 if oUnknown.pair_exps[(iInv, jInv)] == min([oUnknown.den_exps[iInv], oUnknown.den_exps[jInv]]) and oUnknown.pair_friends[(iInv, jInv)] == 2 else 0 for jInv in list_two]))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def flatten(temp_list, recursion_level=0, treat_list_subclasses_as_list=True, treat_tuples_as_lists=False, max_recursion=None):
    from sympy.matrices.dense import MutableDenseMatrix
    from numpy import ndarray
    flat_list = []
    for entry in temp_list:
        if type(entry) is list and (max_recursion is None or recursion_level < max_recursion):
            flat_list += flatten(entry, recursion_level=recursion_level + 1, treat_list_subclasses_as_list=treat_list_subclasses_as_list,
                                 treat_tuples_as_lists=treat_tuples_as_lists, max_recursion=max_recursion)
        elif ((issubclass(type(entry), list) or type(entry) in [MutableDenseMatrix, ndarray]) and
              treat_list_subclasses_as_list is True and (max_recursion is None or recursion_level < max_recursion)):
            flat_list += flatten(entry, recursion_level=recursion_level + 1, treat_list_subclasses_as_list=treat_list_subclasses_as_list,
                                 treat_tuples_as_lists=treat_tuples_as_lists, max_recursion=max_recursion)
        elif (type(entry) is tuple and treat_tuples_as_lists is True and (max_recursion is None or recursion_level < max_recursion)):
            flat_list += flatten(entry, recursion_level=recursion_level + 1, treat_list_subclasses_as_list=treat_list_subclasses_as_list,
                                 treat_tuples_as_lists=treat_tuples_as_lists, max_recursion=max_recursion)
        else:
            flat_list += [entry]
    return flat_list


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

def all_non_empty_subsets(iterable):
    return itertools.chain(*map(lambda x: itertools.combinations(iterable, x), range(1, len(iterable) + 1)))
