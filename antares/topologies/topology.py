#!/usr/bin/env python
# -*- coding: utf-8 -*-

#   _____              _           _
#  |_   _|__ _ __  ___| |___  __ _(_)___ ___
#    | |/ _ \ '_ \/ _ \ / _ \/ _` | / -_|_-<
#    |_|\___/ .__/\___/_\___/\__, |_\___/__/
#           |_|              |___/

# Author: Giuseppe

from __future__ import unicode_literals
from __future__ import print_function

import sys
import re
import itertools

from lips import Particles
from lips.tools import pNB as pNB_internal
from lips.particles_eval import pNB as pNB_overall

from antares.core.tools import flatten

if sys.version_info[0] > 2:
    unicode = str


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


label_to_list_pattern = re.compile(r"([E|I|\d]q[2-4]*p|[E|I|\d]q[2-4]*m|[E|I|\d]qb[2-4]*p|[E|I|\d]qb[2-4]*m|[E|I|\d]gp|[E|I|\d]gm|[E|I|\d]p|[E|I|\d]m|[E|I|\d]yp|[E|I|\d]ym)")


def get_label(oBHUnknown):
    if oBHUnknown.amppart in ['tree', 'rational']:
        return "".join(map("".join, zip(map(str, range(1, len(oBHUnknown.helconf) + 1)), oBHUnknown.helconf))).replace("+", "p").replace("-", "m")
    if oBHUnknown.amppart == 'box':
        corners = [1, 2, 3, 4]
    elif oBHUnknown.amppart == 'triangle':
        corners = [1, 2, 3]
    elif oBHUnknown.amppart == 'bubble':
        corners = [1, 2]
    full_corners = [map(int, [str(eval('oBHUnknown.A1_cutpart.{}({}).corner_ind({}, {})'.format(oBHUnknown.amppart, oBHUnknown.ampindex, corner, corner_ind)))
                              for corner_ind in range(1, eval('oBHUnknown.A1_cutpart.{}({}).corner_size({})'.format(oBHUnknown.amppart, oBHUnknown.ampindex, corner)) + 1)])
                    for corner in corners]
    links = []
    for corner in corners:
        link = str(eval('oBHUnknown.A1_cutpart.{}({}).get_process({})'.format(
            oBHUnknown.amppart, oBHUnknown.ampindex, corner))).replace("(", "").replace(")", "").replace("-", "m").replace("+", "p").replace("101", "").replace("11", "").replace(
                "21", "").replace("gp", "p").replace("gm", "m").split(",")
        links += [[link[0], link[-1]]]
    for i, link in enumerate(links):
        links[i] = ["I" + link[0], "I" + link[-1]]
    label = ""
    for i, full_corner in enumerate(full_corners):
        label += "{}".format(links[i][0]).replace("-", "m").replace("+", "p")
        for particle in full_corner:
            label += "{}{}".format(particle, oBHUnknown.helconf[particle - 1])
        label += "{}".format(links[i][1]).replace("-", "m").replace("+", "p")
    # compactify quark numbering
    for i in range(1, 4):
        index1 = "" if i == 1 else str(i)
        index2 = str(i + 1)
        if (any(["qb{}p".format(index2) in entry or "qb{}m".format(index2) in entry or "q{}p".format(index2) in entry or "q{}m".format(index2) in entry
                 for entry in label_to_list_pattern.findall(label)]) and not
            any(["qb{}p".format(index1) in entry or "qb{}m".format(index1) in entry or "q{}p".format(index1) in entry or "q{}m".format(index1) in entry
                 for entry in label_to_list_pattern.findall(label)])):
            label = label.replace("qb{}p".format(index2), "qb{}p".format(index1))
            label = label.replace("qb{}m".format(index2), "qb{}m".format(index1))
            label = label.replace("q{}p".format(index2), "q{}p".format(index1))
            label = label.replace("q{}m".format(index2), "q{}m".format(index1))
    assert get_nbr_quark_lines(label) == "".join(oBHUnknown.helconf).count("q") / 2 + (oBHUnknown.loopid == "nf")
    assert len(get_external_quarks(label)) == "".join(oBHUnknown.helconf).count("q") / 2
    return label


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def corners_from_label(label):
    label_list = label_to_list_pattern.findall(label)
    internal_indices = [i for i, entry in enumerate(label_list) if entry[0] == "I"]
    corners = [label_list[internal_indices[0 + 2 * i]:internal_indices[1 + 2 * i] + 1] for i in range(len(internal_indices) // 2)]
    return corners


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def helconf_from_label(label):
    helconf = label
    helconf = helconf.replace("Im", "").replace("Ip", "")
    helconf = "".join([entry for entry in helconf[helconf.index("1"):] + helconf[:helconf.index("1")] if entry == "p" or entry == "m"])
    return helconf


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


class Topology(str):

    def __new__(cls, helconf_or_BHUnknown):
        if type(helconf_or_BHUnknown) in [str, unicode]:
            label = helconf_or_BHUnknown
        else:
            label = get_label(helconf_or_BHUnknown)
        topology = label_to_list_pattern.findall(label)
        for i, entry in enumerate(topology):
            if topology[i][0].isdigit():
                topology[i] = "E" + topology[i][1:]
        topology = "".join(topology)
        return super(Topology, cls).__new__(cls, topology)

    def __init__(self, helconf_or_BHUnknown):
        if type(helconf_or_BHUnknown) in [str, unicode]:
            self.label = helconf_or_BHUnknown
        else:
            self.label = get_label(helconf_or_BHUnknown)
            self.oBHUnknown = helconf_or_BHUnknown

    @property
    def topology(self):
        nbr_internal_lines = self.label.count("I") / 2
        if nbr_internal_lines == 4:
            return 'box'
        elif nbr_internal_lines == 3:
            return 'triangle'
        elif nbr_internal_lines == 2:
            return 'bubble'
        else:
            raise Exception("Could not identify topology of {}.".format(self.label))

    @property
    def nbr_external_quarks(self):
        return len(re.findall(pattern="Eqp|Eqm", string=self))

    @property
    def nbr_external_antiquarks(self):
        return len(re.findall(pattern="Eqbp|Eqbm", string=self))

    @property
    def nbr_external_gluons(self):
        return len(re.findall(pattern="Ep|Em", string=self))

    @property
    def nbr_external_photons(self):
        return len(re.findall(pattern="yp|ym", string=self))

    @property
    def equivalence_class(self):
        if not hasattr(self, "_equivalence_class"):
            self._equivalence_class = equivalent_strings(self)
            return self._equivalence_class
        else:
            return self._equivalence_class

    @property
    def models(self):
        if self.topology == 'box':
            return [j + 1 for j, topo in enumerate(all_topologies(self.oBHUnknown)) if topo == self]
        if self.topology == 'triangle':
            nbr_boxes = int(self.oBHUnknown.A1_cutpart.nbr_boxes())
            return [j + 1 - nbr_boxes for j, topo in enumerate(all_topologies(self.oBHUnknown)) if topo == self]
        if self.topology == 'bubble':
            nbr_boxes = int(self.oBHUnknown.A1_cutpart.nbr_boxes())
            nbr_triangles = int(self.oBHUnknown.A1_cutpart.nbr_triangles())
            return [j + 1 - nbr_boxes - nbr_triangles for j, topo in enumerate(all_topologies(self.oBHUnknown)) if topo == self]

    def __eq__(self, other):
        if (hasattr(self, "oBHUnknown") and hasattr(other, "oBHUnknown") and hasattr(self, "_hash") and hasattr(other, "_hash") and
           self.oBHUnknown.process_name == other.oBHUnknown.process_name):
            return self._hash == other._hash
        elif other in self.equivalence_class:
            return True
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):    # returns a value, it is unique for each topology in a give HELCONF
        if hasattr(self, "_hash"):
            return self._hash
        else:
            topologies = all_topologies(self.oBHUnknown)
            for i, topology in enumerate(topologies):
                if self == topology:
                    self._hash = i + 1
                    return i + 1
            raise Exception("Failed hashing of {}.".format(self))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def flip_helicities(string):
    return string.replace("p", "X").replace("m", "p").replace("X", "m")


def flip_quark_line(string, quark_line_nbr):
    if quark_line_nbr == 0:
        return string
    elif quark_line_nbr == 1:
        string = re.sub(r"(qb)(?!\d)", "X", string)
        string = re.sub(r"(q)(?!\d|b)", "qb", string)
        return string.replace("X", "q")
    else:
        quark = "q{}".format(quark_line_nbr)
        anti_quark = "qb{}".format(quark_line_nbr)
        return string.replace(anti_quark, "X").replace(quark, anti_quark).replace("X", quark)


def get_nbr_quark_lines(string):
    if list(set(map(int, re.findall(r"q(\d)", string)))) != []:
        nbr_quark_lines = max(list(set(map(int, re.findall(r"q(\d)", string)))))
    elif "q" in string:
        nbr_quark_lines = 1
    else:
        nbr_quark_lines = 0
    return nbr_quark_lines


def get_external_quarks(label):
    external_quark_pattern = re.compile(r"(?<!I)q(\d*)")
    external_quarks = list(set(external_quark_pattern.findall(label)))
    if "" in external_quarks:
        external_quarks[external_quarks.index("")] = 1
    return external_quarks


def equivalent_strings(string):
    nbr_quark_lines = get_nbr_quark_lines(string)
    self_as_list = label_to_list_pattern.findall(string)
    equivalence_class = []
    for i in range(len(self_as_list)):
        # read forewards
        permutation_as_list = self_as_list[i:] + self_as_list[:i]
        permutation = "".join(permutation_as_list)
        if not ("I" not in permutation or (permutation_as_list[0][0] == "I" and permutation_as_list[-1][0] == "I")):
            continue
        # quark lines flips
        for combination in flatten([[entry for entry in itertools.combinations(range(1, nbr_quark_lines + 1), i)] for i in range(nbr_quark_lines + 2)]):
            new_permutation = permutation
            for quark_line_index in combination:
                new_permutation = flip_quark_line(new_permutation, quark_line_index)
            equivalence_class += [new_permutation]
        permutation = flip_helicities(permutation)  # chirality flip
        # quark lines flips
        for combination in flatten([[entry for entry in itertools.combinations(range(1, nbr_quark_lines + 1), i)] for i in range(nbr_quark_lines + 2)]):
            new_permutation = permutation
            for quark_line_index in combination:
                new_permutation = flip_quark_line(new_permutation, quark_line_index)
            equivalence_class += [new_permutation]
    for i in range(len(self_as_list)):
        # read backwards
        permutation_as_list = self_as_list[::-1][i:] + self_as_list[::-1][:i]
        permutation = "".join(permutation_as_list)
        if not ("I" not in permutation or (permutation_as_list[0][0] == "I" and permutation_as_list[-1][0] == "I")):
            continue
        # quark lines flips
        for combination in flatten([[entry for entry in itertools.combinations(range(1, nbr_quark_lines + 1), i)] for i in range(nbr_quark_lines + 2)]):
            new_permutation = permutation
            for quark_line_index in combination:
                new_permutation = flip_quark_line(new_permutation, quark_line_index)
            equivalence_class += [new_permutation]
        permutation = flip_helicities(permutation)  # chirality flip
        # quark lines flips
        for combination in flatten([[entry for entry in itertools.combinations(range(1, nbr_quark_lines + 1), i)] for i in range(nbr_quark_lines + 2)]):
            new_permutation = permutation
            for quark_line_index in combination:
                new_permutation = flip_quark_line(new_permutation, quark_line_index)
            equivalence_class += [new_permutation]
    return equivalence_class


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def conversion_rules(source_label_or_bhunknown, target_label_or_bhunknown):
    source_topo, target_topo = Topology(source_label_or_bhunknown), Topology(target_label_or_bhunknown)
    if source_topo != target_topo:
        raise Exception("Source and target are of different topologies. No conversion rule is known.")
    parity_of_process = (-1) ** source_topo.nbr_external_gluons * (-1) ** source_topo.nbr_external_photons * (-1) ** source_topo.nbr_external_quarks
    nbr_quark_lines = get_nbr_quark_lines(source_topo)
    # external_quarks = get_external_quarks(source_topo)
    quark_lines_flips = flatten([[entry for entry in itertools.combinations(range(1, nbr_quark_lines + 1), i)] for i in range(nbr_quark_lines + 2)])
    len_without_chirality_flip = len(quark_lines_flips)
    all_permutations = [[entry[0] for entry in label_to_list_pattern.findall(equivalent_string) if entry[0].isdigit()] for equivalent_string in equivalent_strings(target_topo.label)]
    base_permutation = [entry[0] for entry in label_to_list_pattern.findall(source_topo.label) if entry[0].isdigit()]
    # len_cycle = source_topo.count("E") if "I" not in source_topo else source_topo.count("I") / 2  # multiplicity: tree/rational, 4: box, 3: triangle, 2:bubble.
    permutation_indices, complex_conjugates, signs = [], [], []
    for i, equivalent_topology in enumerate(target_topo.equivalence_class):
        if equivalent_topology == str(source_topo):
            # sign from flipping an external quark line
            # if len(external_quarks) > 0:
            #     flip_index = (i % (len(all_permutations) / 2)
            #                   % ((len(all_permutations) / 2) / len_cycle)
            #                   % (((len(all_permutations) / 2) / len_cycle) / len_without_chirality_flip))
            #     flip_tuple = quark_lines_flips[flip_index]
            #     if len(external_quarks) == 1:
            #         sign = -1 if external_quarks[0] in flip_tuple else 1
            #     else:
            #         raise Exception("Multiple external lines not supported yet!")
            # else:
            #     sign = +1
            # sign from backwards permutation, i.e. parity
            backwards = 1 if i < len(all_permutations) / 2 else parity_of_process
            sign = backwards
            permutation_indices += [i]
            if i % (2 * len_without_chirality_flip) >= len_without_chirality_flip:
                complex_conjugates += [True]
            else:
                complex_conjugates += [False]
            signs += ["+" if sign == +1 else "-"]
    symmetries = list(set(zip(["".join(map(str, entry)) for i, entry in enumerate(all_permutations) if i in permutation_indices], complex_conjugates, signs)))
    symmetries.sort()
    if base_permutation != map(str, range(1, len(all_permutations[0]) + 1)):
        for i, symmetry in enumerate(symmetries):
            rule_on_base_permutation = symmetry[0]
            rule_on_ordered = "".join([rule_on_base_permutation[base_permutation.index(entry)] for entry in "".join(map(str, range(1, len(all_permutations[0]) + 1)))])
            symmetries[i] = (rule_on_ordered, symmetry[1], symmetry[2])
    return symmetries


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def internal_symmetry(oBHUnknown):
    return [conversion_rule for conversion_rule in conversion_rules(oBHUnknown, oBHUnknown) if conversion_rule != ("".join(map(str, range(1, oBHUnknown.multiplicity + 1))), False, "+")]


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def box_topologies(oBHUnknown):
    original_amppart, original_ampindex = oBHUnknown.amppart, oBHUnknown.ampindex
    topologies = []
    for i in range(1, oBHUnknown.A1_cutpart.nbr_boxes() + 1):
        oBHUnknown.amppart, oBHUnknown.ampindex = "box", i
        topologies += [Topology(oBHUnknown)]
    oBHUnknown.amppart, oBHUnknown.ampindex = original_amppart, original_ampindex
    return topologies


def triangle_topologies(oBHUnknown):
    original_amppart, original_ampindex = oBHUnknown.amppart, oBHUnknown.ampindex
    topologies = []
    for i in range(1, oBHUnknown.A1_cutpart.nbr_triangles() + 1):
        oBHUnknown.amppart, oBHUnknown.ampindex = "triangle", i
        topologies += [Topology(oBHUnknown)]
    oBHUnknown.amppart, oBHUnknown.ampindex = original_amppart, original_ampindex
    return topologies


def bubble_topologies(oBHUnknown):
    original_amppart, original_ampindex = oBHUnknown.amppart, oBHUnknown.ampindex
    topologies = []
    for i in range(1, oBHUnknown.A1_cutpart.nbr_bubbles() + 1):
        oBHUnknown.amppart, oBHUnknown.ampindex = "bubble", i
        topologies += [Topology(oBHUnknown)]
    oBHUnknown.amppart, oBHUnknown.ampindex = original_amppart, original_ampindex
    return topologies


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


all_topologies_cache = {}


def all_topologies(oBHUnknown):
    if oBHUnknown.process_name in all_topologies_cache.keys():  # read from cache
        return all_topologies_cache[oBHUnknown.process_name]
    else:
        topologies = box_topologies(oBHUnknown) + triangle_topologies(oBHUnknown) + bubble_topologies(oBHUnknown)
        all_topologies_cache[oBHUnknown.process_name] = topologies
        for i, topology in enumerate(topologies):
            hash(topology)
            all_topologies_cache[oBHUnknown.process_name] = topologies
    return topologies


def independent_topologies(oBHUnknown):
    inde_topos = list(set(all_topologies(oBHUnknown)))
    independent_boxes = [topo for topo in inde_topos if topo.topology == 'box']
    independent_boxes.sort(key=lambda x: hash(x))
    independent_triangles = [topo for topo in inde_topos if topo.topology == 'triangle']
    independent_triangles.sort(key=lambda x: hash(x))
    independent_bubbles = [topo for topo in inde_topos if topo.topology == 'bubble']
    independent_bubbles.sort(key=lambda x: hash(x))
    return independent_boxes, independent_triangles, independent_bubbles


def topology_info(oBHUnknown):
    independent_boxes, independent_triangles, independent_bubbles = independent_topologies(oBHUnknown)
    topology_info = ""
    for i, independent_box in enumerate(independent_boxes):
        topology_info += "box_topology({}): {}".format(i + 1, independent_box.models) + "\n"
    for i, independent_triangle in enumerate(independent_triangles):
        topology_info += "triangle_topology({}): {}".format(i + 1, independent_triangle.models) + "\n"
    for i, independent_bubble in enumerate(independent_bubbles):
        topology_info += "bubble_topology({}): {}".format(i + 1, independent_bubble.models) + "\n"
    return topology_info[:-1]


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def convert_invariant(Invariant, Rule):
    if type(Invariant) not in [unicode, str]:
        return Invariant
    if Invariant == "1":
        return "1"
    # Rule Parsing
    OneLinePermutation = Rule[0]
    ComplexConjugation = Rule[1]
    n = len(OneLinePermutation)
    RuleDict = dict(zip([unicode(entry) for entry in range(1, n + 1)], [entry for entry in OneLinePermutation]))
    # Special Variables
    if "Δ" in Invariant or "Ω" in Invariant or "Π" in Invariant:
        from antares.core.tools import pDijk, pOijk, pPijk
        oParticles = Particles(n)
        if "Δ" in Invariant:
            match_list = pDijk.findall(Invariant)[0]
            if match_list[0] == '':
                # adjacent = False  # not used?
                NonOLists = [list(map(int, corner)) for corner in match_list[1:]]
            else:
                # adjacent = True  # not used?
                NonOLists = oParticles.ijk_to_3NonOverlappingLists(list(map(int, match_list[0])))
        elif "Ω" in Invariant:
            ijk = list(map(int, pOijk.findall(Invariant)[0]))
            NonOLists = oParticles.ijk_to_3NonOverlappingLists(ijk)
        elif "Π" in Invariant:
            ijk = list(map(int, pPijk.findall(Invariant)[0]))
            NonOLists = oParticles.ijk_to_3NonOverlappingLists(ijk)
        NonOLists = [sorted([int(RuleDict[str(entry)]) for entry in NonOList]) for NonOList in NonOLists]
        NonOLists = [entry if 1 not in entry and n not in entry else
                     entry[[i for i, _entry in enumerate(entry) if _entry != entry[i - 1] + 1][-1]:] +
                     entry[:[i for i, _entry in enumerate(entry) if _entry != entry[i - 1] + 1][-1]] for entry in NonOLists]
    # Convert Invariant
    Invariant = "".join([character if character not in RuleDict or (character == "5" and Invariant[i - 2:i] == "tr")
                         else RuleDict[character] for i, character in enumerate(Invariant)])
    if ("Δ" in Invariant and Invariant.count("|") == 1) or "Ω" in Invariant or "Π" in Invariant:
        if "Δ" in Invariant:
            ijk = list(map(int, pDijk.findall(Invariant)[0][0]))
        elif "Ω" in Invariant:
            ijk = list(map(int, pOijk.findall(Invariant)[0]))
        elif "Π" in Invariant:
            ijk = list(map(int, pPijk.findall(Invariant)[0]))
        NonOListIndices = [i for entry in ijk for i, NonOList in enumerate(NonOLists) if entry in NonOList]
        if "Δ" in Invariant:
            ijk = [NonOLists[NonOListIndices[0]][0], NonOLists[NonOListIndices[1]][0], NonOLists[NonOListIndices[2]][0]]
            ijk.sort()
            if ijk[0] != 1:
                ijk = ijk[1:] + ijk[:1]
            Invariant = "Δ_" + "".join(list(map(str, ijk)))
        elif "Ω" in Invariant or "Π" in Invariant:
            a = NonOLists[NonOListIndices[0]][0]
            b = NonOLists[NonOListIndices[1]][0]
            c = NonOLists[NonOListIndices[2]][0]
            if b > c:         # make sure b < c
                b, c = c, b
            if a < b:
                pass
            elif a > c:
                pass
            else:
                b, c = c, b
            if "Ω" in Invariant:
                Invariant = "Ω_" + str(a) + str(b) + str(c)
            elif "Π" in Invariant:
                Invariant = "Π_" + str(a) + str(b) + str(c)
    if ComplexConjugation is True:
        subs_dict = {"⟨": "]", "⟩": "[", "]": "⟨", "[": "⟩"}

        def conjugate_pNB(match):
            abcd = pNB_internal.search(match.group())
            a = int(abcd.group('start'))
            bcs = abcd.group('middle').replace("(", "").replace(")", "").split("|")
            d = int(abcd.group('end'))
            for i, bc in enumerate(bcs):
                bcs[i] = re.findall(r"([+|-]{0,1}\d+)", bc)
            begin = re.sub(r"(⟨|⟩|\]|\[)", lambda match_inner: subs_dict[match_inner.group()], match.group()[-1])
            end = re.sub(r"(⟨|⟩|\]|\[)", lambda match_inner: subs_dict[match_inner.group()], match.group()[0])
            if "(" in match.group() and ")" in match.group():
                return f"{begin}{d}|({')|('.join(['+'.join(bc) for bc in bcs[::-1]])})|{a}{end}".replace("++", "+").replace("+-", "-")
            else:
                return f"{begin}{d}|{'|'.join(['+'.join(bc) for bc in bcs[::-1]])}|{a}{end}".replace("++", "+").replace("+-", "-")

        if re.findall(r"\[(\d)\|(\d)\]", Invariant) != []:
            Invariant = re.sub(r"\[(\d)\|(\d)\]", r"⟨\2|\1⟩", Invariant)
        elif re.findall(r"⟨(\d)\|(\d)⟩", Invariant) != []:
            Invariant = re.sub(r"⟨(\d)\|(\d)⟩", r"[\2|\1]", Invariant)
        elif pNB_overall.findall(Invariant) != []:
            Invariant = pNB_overall.sub(lambda match: conjugate_pNB(match), Invariant)
        elif "⟨" in Invariant or "⟩" in Invariant or "[" in Invariant or "]" in Invariant:
            raise Exception(f"Could not conjugate {Invariant}.")
    return Invariant


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def sign_of_permutation(permutation):
    permutation = map(int, list(permutation))
    permutation_counter = 0
    while True:
        switch_happend = False
        for i, entry in enumerate(permutation):
            if i == len(permutation) - 1:
                continue
            if permutation[i] > permutation[i + 1]:
                permutation[i], permutation[i + 1] = permutation[i + 1], permutation[i]
                switch_happend = True
                permutation_counter += 1
        if switch_happend is False:
            break
    return permutation_counter % 2
