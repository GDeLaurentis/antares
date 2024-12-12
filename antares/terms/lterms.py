import functools
import numpy
import pathlib

from lips.symmetries import inverse

from antares.core.tools import Generate_LaTeX_and_PDF
from antares.terms.terms import LoadResults, Terms
from ..core.numerical_methods import Numerical_Methods


class TermsList(Numerical_Methods, list):

    def __init__(self, list_of_terms, multiplicity, verbose=False):
        list.__init__(self)
        if isinstance(list_of_terms, list):
            self.extend(list_of_terms)
            self.multiplicity = multiplicity
        elif isinstance(list_of_terms, (str, pathlib.Path)):
            path = list_of_terms
            path = pathlib.Path(path)
            path.resolve(strict=False)
            try:
                with open(path / "basis.txt", "r") as file:
                    content = file.read()
            except FileNotFoundError:
                if verbose:
                    print("\rNo basis found, returning empty basis.                  ")
                return
            basis = eval(content)
            for index, basis_entry_file_or_symmetry in enumerate(basis):
                if verbose:
                    print(f"\r @ {index}", end="")
                if isinstance(basis_entry_file_or_symmetry, str):
                    basis[index] = LoadResults(path / basis_entry_file_or_symmetry)[0][0]
                    basis[index].multiplicity = multiplicity
            if verbose:
                print(f"\rLoaded basis of size {len(basis)}                           ")
            self.__init__(basis, multiplicity, verbose)
        else:
            raise TypeError("Expected a list or a path as input.")

    def __hash__(self):
        return hash(tuple(map(hash, self)))

    def __getitem__(self, item):
        if isinstance(item, slice):
            return TermsList(super().__getitem__(item), self.multiplicity)
        else:
            return super().__getitem__(item)

    @functools.lru_cache(maxsize=256)
    def __call__(self, oPs):
        numerical_basis, last_coeff = [], None
        for basis_element in self:
            if isinstance(basis_element, tuple):
                numerical_basis += [last_coeff(oPs.image(basis_element)) if isinstance(last_coeff, Terms) else oPs.image(basis_element)(last_coeff)]
            else:
                last_coeff = basis_element
                numerical_basis += [last_coeff(oPs) if isinstance(last_coeff, Terms) else oPs(last_coeff)]
        if isinstance(self, numpy.ndarray):
            return numpy.array(numerical_basis)
        else:
            return numerical_basis

    def save(self, result_path, naming_convention=["dense", "sparse"][0], overwrite_basis=True):
        assert naming_convention in ["dense", "sparse"]
        with open(result_path + "basis.txt", "w") as f:
            f.write("[" + ",\n ".join(map(str, [entry if isinstance(entry, tuple) else
                                                f"\'coeff_{i if naming_convention == 'sparse' else sum([1 for _entry in self[:i] if isinstance(_entry, Terms)])}\'"
                                                for i, entry in enumerate(self)])) + "]")
        for i, entry in enumerate(self):
            index = (i if naming_convention == 'sparse' else sum([1 for _entry in self[:i] if isinstance(_entry, Terms)]))
            if isinstance(entry, Terms) and (overwrite_basis or not pathlib.Path(result_path + f"coeff_{index}.pdf").is_file()):
                Generate_LaTeX_and_PDF(entry.Write_LaTex(), result_path + f"coeff_{index}")

    def explicit_representation(self):
        basis_explicit, last_coeff = [], None
        for basis_element in self:
            if isinstance(basis_element, tuple):
                basis_explicit += [last_coeff.Image(inverse(basis_element))]
            else:
                last_coeff = basis_element.explicit_representation()
                basis_explicit += [last_coeff]
        return TermsList(basis_explicit, self.multiplicity)
