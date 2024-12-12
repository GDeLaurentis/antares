import pathlib

from lips import Particles
from syngular import Field

from antares.terms.lterms import TermsList
from antares.terms.terms import Terms


test_folder = pathlib.Path(__file__).parent


def test_basis_eval():
    basis = TermsList(test_folder / 'test_basis', multiplicity=6, verbose=True)
    oPs = Particles(6, field=Field("finite field", 2 ** 31 - 1, 1))
    numerical_basis = basis(oPs)
    print(numerical_basis)


def test_str_and_rstr():
    basis = TermsList(test_folder / 'test_basis', multiplicity=6, verbose=True)
    oPs = Particles(6, field=Field("finite field", 2 ** 31 - 1, 1))
    for i, oTerms in enumerate(basis):
        if isinstance(oTerms, Terms):
            print(f"checking {i}:\n {oTerms} \n\n {Terms(str(oTerms))} \n\n")
            # print(hash(oTerms), hash(Terms(str(oTerms))))  # perhaps eventually check hash here
            assert oTerms(oPs) == Terms(str(oTerms))(oPs)


def test_basis_explicit_representation():
    basis = TermsList(test_folder / 'test_basis', multiplicity=6, verbose=True)
    oPs = Particles(6, field=Field("finite field", 2 ** 31 - 1, 1))
    assert basis(oPs) == basis.explicit_representation()(oPs)
