import pytest

from multiprocessing import Pool
from antares.core.tools import mapThreads


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def test_mapThreads_simple():
    def shift(x):
        return x + 1
    assert mapThreads(lambda x: x + 1, range(100), Cores=4) == mapThreads(shift, range(100), Cores=4) == list(range(1, 101))


def test_mapThreads_multiple_arguments():
    def shift(x, y):
        return x + y
    assert mapThreads(shift, 5, range(100), Cores=4) == list(range(5, 105))


def test_mapThreads_nodaemon():
    def shift(xs, y):
        return sum(mapThreads(lambda x: x + 1, xs, Cores=4)) + y
    assert mapThreads(shift, range(100), range(100), Cores=4) == list(range(5050, 5150))


def f(x):
    return sum(mapThreads(lambda x: x + 1, range(x), UseParallelisation=False))


def g(x):
    return sum(mapThreads(lambda x: x + 1, range(x), UseParallelisation=True))


def test_daemonic_parallelisation_off():
    with Pool() as pool:
        res = pool.map(f, range(10))
    assert sum(res) == 165


def test_daemonic_parallelisation_on():
    with pytest.raises(AssertionError):
        with Pool() as pool:
            pool.map(g, range(10))
