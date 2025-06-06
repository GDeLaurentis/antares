import functools
import numpy
import pathlib

from lips.symmetries import inverse

from ..core.tools import Generate_LaTeX_and_PDF
from ..core.numerical_methods import Numerical_Methods
from ..core.settings import settings
from .terms import LoadResults, Terms


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

    @property
    def poles_and_orders(self):
        poles_and_orders = {}
        for i, oTerms in enumerate(self):
            these_poles_and_orders = set(zip(oTerms[0].oDen.lInvs, oTerms[0].oDen.lExps))
            for pole, order in these_poles_and_orders:
                if pole in poles_and_orders.keys():
                    poles_and_orders[pole] = max(poles_and_orders[pole], order)
                else:
                    poles_and_orders[pole] = order
        if settings.invariants is not None:
            poles_and_orders = dict(sorted(poles_and_orders.items(), key=lambda x: settings.invariants.index(x[0])
                                           if x[0] in settings.invariants else 99))
        return poles_and_orders

    @property
    def max_sizes_poles_vector_spaces(self):
        max_size_vector_spaces = {}
        for key, val in self.poles_and_orders.items():
            counter_dict = {}
            for oTerms in self:
                if key in oTerms[0].oDen.lInvs:
                    exp = oTerms[0].oDen.lExps[oTerms[0].oDen.lInvs.index(key)]
                    if exp in counter_dict.keys():
                        counter_dict[exp] += 1
                    else:
                        counter_dict[exp] = 1
            max_size_vector_spaces[key] = counter_dict
        return max_size_vector_spaces
    
    @staticmethod
    def cumulative_pole_orders(pole_dict):
        orders = sorted(pole_dict.keys(), reverse=True)
        
        cumulative = {}
        running_total = 0
        for order in orders:
            running_total += pole_dict[order]
            cumulative[order] = running_total

        return cumulative

    @property
    def comulative_max_sizes_poles_vector_spaces(self):
        return {key: self.cumulative_pole_orders(val) for key, val in self.max_sizes_poles_vector_spaces.items()}

    @property
    def max_size_of_all_poles_vector_spaces(self):
        return max([max([val2 for _, val2 in val.items()]) for _, val in self.comulative_max_sizes_poles_vector_spaces.items()]) + 1
