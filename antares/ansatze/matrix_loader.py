#   __  __      _       _     _                 _
#  |  \/  |__ _| |_ _ _(_)_ _| |   ___  __ _ __| |___ _ _
#  | |\/| / _` |  _| '_| \ \ / |__/ _ \/ _` / _` / -_) '_|
#  |_|  |_\__,_|\__|_| |_/_\_\____\___/\__,_\__,_\___|_|

# Author: Giuseppe
# Created: 10/07/2018

import os
import numpy
import operator
import copy
import functools

from lips import Particles
from linac import timeit
from pycoretools import mapThreads, flatten

from ..core.settings import settings
from ..core.tools import MyShelf
from .interface import SymmetriseAnsatze

local_directory = os.path.dirname(os.path.abspath(__file__))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


@timeit
def load_matrix(oTerms, oAnsatz, nRows=None, bAugmentedMatrix=True):
    """Load the matrix given a oTerms object (~denominator ansatz) and oAnsatz object (~numerator ansatz)."""

    if nRows is None:
        nRows = oAnsatz.nInput
    nColumns = oAnsatz.nInput + bAugmentedMatrix
    print("\rThe matrix size is {}x{}.                                                       ".format(nRows, nColumns), end="")

    if any([oTerm.am_I_a_symmetry for oTerm in oTerms]):
        oAnsatz = SymmetriseAnsatze(oAnsatz, oTerms)

    Matrix = []

    lphase_space_points = phase_space_points(multiplicity=oTerms.multiplicity, nbr_points=nRows,
                                             small_invs=oTerms.oFittingSettings.lSmallInvs,
                                             small_invs_exps=oTerms.oFittingSettings.lSmallInvsExps,
                                             UseParallelisation=settings.UseParallelisation, Cores=settings.Cores)
    oTermsExplicit = oTerms.explicit_representation(raw=True)

    basis = oAnsatz.lInvariants + list(oTermsExplicit)
    lindices = numpy.zeros((oAnsatz.nInput, max([len(term) + 1 for term in flatten(oAnsatz, max_recursion=1) if term != 1])), dtype=numpy.int32)
    for i, ansatz_term in enumerate(oAnsatz):
        for j, jInvs in enumerate(ansatz_term):
            row_index = sum([len(entry) for entry in oAnsatz[:i]]) + j
            lindices[row_index, 0] = basis.index(oTermsExplicit[i])
            for k, kInv in enumerate(jInvs):
                lindices[row_index, k + 1] = basis.index(kInv)
    bases = [[entry if type(entry) is int else oParticles.compute(entry) if isinstance(entry, str) else entry(oParticles) for entry in basis]
             for oParticles in lphase_space_points]

    if settings.UseGpu is True and settings.field.name == "mpc":
        bases = numpy.array([list(map(complex, _basis)) for _basis in bases])
    elif settings.UseGpu is True and settings.field.name == "finite field":
        bases = numpy.array([list(map(int, _basis)) for _basis in bases])

    # Load "Square" Matrix --- Might be rectangular if there are symmetries, but it will be folded back to square
    if settings.UseGpu is True and settings.field.name in ["mpc", "finite field"]:
        print("\rLoading square matrix on gpu.                             ", end="\n")
        rectangular_matrix = cuda_load_square_matrix(bases, lindices, (nRows, oAnsatz.nInput))
        rectangular_matrix = rectangular_matrix.tolist()
    else:
        print("\rLoading square matrix on cpu.                             ", end="\n")
        rectangular_matrix = mapThreads(evaluate_row, lindices, bases, UseParallelisation=settings.UseParallelisation, Cores=settings.Cores)
    print("\rFinished loading the square matrix.                  ", end="")

    # Fold it if symmetries were used
    if oAnsatz.nInput > nColumns - bAugmentedMatrix:
        lSymmetries = [oTerm for oTerm in oTerms if oTerm.am_I_a_symmetry is True]
        for i, irow in enumerate(rectangular_matrix):
            for j in range(len(irow)):
                if j % (nColumns - bAugmentedMatrix) == j:
                    continue
                elif lSymmetries[j // (nColumns - bAugmentedMatrix) - 1].tSym[2] == "+":
                    irow[j % (nColumns - bAugmentedMatrix)] += irow[j]
                elif lSymmetries[j // (nColumns - bAugmentedMatrix) - 1].tSym[2] == "-":
                    irow[j % (nColumns - bAugmentedMatrix)] -= irow[j]
                else:
                    raise Exception("Sign of symmetry unspecfied.")
            rectangular_matrix[i] = irow[0:(nColumns - bAugmentedMatrix)]
    square_matrix = rectangular_matrix

    # Load Augmented Column
    if bAugmentedMatrix is True:
        print("\rLoading augmented column.                           ", end="")
        augmented_column = mapThreads(oTerms.oUnknown, lphase_space_points, UseParallelisation=settings.UseParallelisation, Cores=settings.Cores)
        oTerms.oUnknown.save_original_unknown_call_cache()
        if settings.UseGpu is True:
            if settings.field.name == "mpc":
                augmented_column = list(map(complex, augmented_column))
            elif settings.field.name == "finite field":
                augmented_column = list(map(int, augmented_column))

        # Stitch Them Together
        Matrix = [row + [augmented_column[i]] for i, row in enumerate(square_matrix)]
    else:
        Matrix = square_matrix
    print("\rFinished loading augmented column.                                                                            ", end="")

    # old attempt at increasing matrix row reduction stability by recomputing rows with large scale differences (not needed?)
    if False and settings.InversionUsesGMPs is False:
        max_rows = [max(abs(entry) for entry in row) for row in Matrix]
        min_rows = [min(abs(entry) for entry in row) for row in Matrix]
        ratios = [max_rows[i] / min_rows[i] if min_rows[i] != 0 else 1 for i in range(len(max_rows))]   # note: cheeky fix for min_row == 0 case, possibly due to smallinv in ansatz
        ToBeRejected = [i for (i, entry) in enumerate(ratios) if entry > 10 ** 6]
        if len(ToBeRejected) > 0 and len(ToBeRejected) < len(ratios):
            print("\033[1A\r                                                                                                         ", end="")
            print("\rRejecting {}/{}.".format(len(ToBeRejected), len(ratios)))
            NewMatrix = [entry for (i, entry) in enumerate(Matrix) if i not in ToBeRejected]
            NewMatrix += load_matrix(oTerms, oAnsatz, (len(ToBeRejected), nColumns))[1]
        else:
            NewMatrix = [entry for (i, entry) in enumerate(Matrix)]
        return NewMatrix

    return numpy.array(Matrix)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def cuda_load_square_matrix(bases, lindices, nRows_nColumns):
    nRows, nColumns = nRows_nColumns

    import pycuda.driver as cuda
    import pycuda.autoinit                    # noqa
    from pycuda.compiler import SourceModule  # noqa
    from linac.pycuda_tools import cuda_set_vars_and_get_funcs, folded_number_of_columns

    # Compile Cuda Code
    exec(str(cuda_set_vars_and_get_funcs(path_to_cuda_script=local_directory + "/matrix_loader.cu", BASIS_LENGTH=len(bases[0]), NBR_ROWS=nRows,
                                         NBR_COLUMNS=nColumns, MASS_DIMENSION=len(lindices[0]))), locals(), globals())

    # Push Bases And Indices To Device
    bases_gpu = cuda.mem_alloc(bases.size * bases.dtype.itemsize)
    cuda.memcpy_htod(bases_gpu, bases)
    lindices_gpu = cuda.mem_alloc(lindices.size * lindices.dtype.itemsize)
    cuda.memcpy_htod(lindices_gpu, lindices)

    # Push Empty Matrix To Device
    matrix = numpy.zeros((nRows, nColumns), dtype=numpy.complex128)
    matrix_gpu = cuda.mem_alloc(matrix.size * matrix.dtype.itemsize)
    cuda.memcpy_htod(matrix_gpu, matrix)

    # Load Matrix On GPGPU
    CudaLoadMatrix(matrix_gpu, bases_gpu, lindices_gpu, block=(folded_number_of_columns(nColumns), 1, 1), grid=(nRows, 1))  # noqa

    # Pull Loaded Matrix Back To Host
    cuda.memcpy_dtoh(matrix, matrix_gpu)

    return matrix


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def evaluate_row(lindices, basis):
    return [functools.reduce(operator.mul, [basis[index] for index in indices]) for indices in lindices]


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


@timeit
def phase_space_points(multiplicity=None, nbr_points=None, small_invs=None, small_invs_exps=None, UseParallelisation=True, Cores=6):
    """Returns phase space points (Particles objects) of given multiplicity in collinear limit described by small_invs & small_invs_exps."""
    lParticles = mapThreads(phase_space_point, multiplicity, small_invs, small_invs_exps, range(nbr_points), UseParallelisation=UseParallelisation, Cores=Cores)
    return lParticles


def phase_space_point(multiplicity, small_invs, small_invs_exps, nbr_point):
    while True:  # avoid field extension with padics in phase space construction
        if settings.field.name == "mpc":
            seed = nbr_point
        else:
            seed = None  # can't be seeded if need to retry
        oParticles = Particles(multiplicity, field=settings.field, seed=seed)
        if small_invs is None or len(small_invs) == 0:
            pass
        elif len(small_invs) == len(small_invs_exps) == 1:
            if settings.field.name == "mpc":
                oParticles.variety(small_invs, (10 ** -int(0.80 * 300 / 5), ))
            elif settings.field.name == "padic":
                oParticles.variety(small_invs, (1, ))
        elif len(small_invs) == len(small_invs_exps) == 2:
            if settings.field.name == "mpc":
                oParticles.variety(small_invs, (10 ** -int(0.80 * 300 / 10), 10 ** - int(2 * 0.80 * 300 / 10), ))
            elif settings.field.name == "padic":
                oParticles.variety(small_invs, (1, 1, ))
        else:
            raise Exception("Bad format for small_invs and small_invs_exps in phase_space_points.")
        if oParticles.spinors_are_in_field_extension:
            continue
        else:
            break
    return oParticles


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# Deprecated

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def LoadMatrixInfoFromCache(terms, nInput):
    print("\rLoading Matrix Info from cache.")
    with MyShelf(settings.PWMatrixInfo, 'r') as sMatrixInfo:
        dropped_redundant = sMatrixInfo[str("dropped_redundant")]
        dropped_zero = sMatrixInfo[str("dropped_zero")]
    _terms = copy.deepcopy(terms)
    for i, term in enumerate(terms):
        for j, invariants in enumerate(term):
            index = j + sum(len(entry) for entry in terms[0:i])
            if index in dropped_redundant + dropped_zero:
                _terms[i].remove(invariants)
    terms = _terms
    nInput = sum(len(entry) for entry in terms)
    print("\rThe new matrix size is {}x{}.                                                       ".format(nInput, nInput + 1), end="")
    return terms, nInput


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def AnalyseNumerics(Matrix, ID):
    import matplotlib.pyplot as plt
    import math
    AbsMatrixEntries = [abs(entry) for row in Matrix for entry in row]
    plt.hist([entry / math.log(10) for entry in map(math.log, filter(lambda x: x > 0, AbsMatrixEntries))])
    plt.title("AbsMatrixEntries")
    plt.xlabel("Abs(Entry)")
    plt.ylabel("Frequency")
    plt.savefig("AbsMatrixEntries{}.png".format(ID))
    plt.clf()
    max_rows = [max(abs(entry) for entry in row) for row in Matrix]
    min_rows = [min(abs(entry) for entry in row) for row in Matrix]
    ratios = [max_rows[i] / min_rows[i] if min_rows[i] != 0 else 1 for i in range(len(max_rows))]
    plt.hist([entry / math.log(10) for entry in map(math.log, filter(lambda x: x > 0, ratios))])
    plt.title("RatiosByRows")
    plt.xlabel("max_rows/min_rows")
    plt.ylabel("Frequency")
    plt.savefig("RatiosByRows{}.png".format(ID))
    # print ratios
    # ToBeRejected = [i for (i, entry) in enumerate(ratios) if entry > 10 ** 6]
    # print "Hypothetically reject {}/{}.".format(len(ToBeRejected), len(ratios))


# if False and settings.InversionUsesGMPs is False and "SSH_CONNECTION" not in os.environ:
#     AnalyseNumerics(Matrix, 1)
#     AnalyseNumerics(NewMatrix, 2)
#     print("")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


# def LoadUniqueMatrix(nInput, dropped, GMPMatrix):
#     # reconstruct the unique matrix
#     nUnique = nInput - len(dropped)
#     UniqueMatrix = gmpTools.GMPCmatrix(nUnique, nUnique + 1)
#     for i in range(nUnique):
#         # square part
#         vs = [GMPMatrix[i, j] for j in range(nInput) if j not in dropped]
#         for j, v in enumerate(vs):
#             UniqueMatrix[i, j] = v
#         # add back the last column
#         UniqueMatrix[i, nUnique] = GMPMatrix[i, nInput]
#     return UniqueMatrix
