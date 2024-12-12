#!/usr/bin/env python
# -*- coding: utf-8 -*-

#   ___ _  _   _   _      _
#  | _ ) || | | | | |_ _ | |___ _  _____ __ ___ _
#  | _ \ __ | | |_| | ' \| / / ' \/ _ \ V  V / ' \
#  |___/_||_|  \___/|_||_|_\_\_||_\___/\_/\_/|_||_|


# Author: Giuseppe

import sys
import os
import re
import subprocess
import functools
import mpmath
import multiprocessing
import antares

from antares.core.settings import settings
from antares.core.bh_patch import BH_found, gmpTools_found
from antares.core.tools import OutputGrabber, MyShelf, mpc_to_cgmp, cgmp_to_mpc, flatten
from antares.core.numerical_methods import Numerical_Methods
from lips.invariants import Invariants

if BH_found:
    from antares.core.bh_patch import BH
if gmpTools_found:
    import gmpTools

mpmath.mp.dps = 300

if sys.version_info.major > 2:
    unicode = str

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


class BHUnknown(Numerical_Methods, object):

    def __init__(self, helconf, loopid=None, amppart=None, ampindex=None):
        self.helconf = helconf
        self.loopid = loopid
        self.amppart = amppart
        self.ampindex = ampindex
        if loopid is not None:
            self.print_graphs()
            self.print_topology_info()
            self.get_topology_info()
        self.get_call_cache()

    @property
    def helconf(self):
        """Helicity configuration. Accepts string or list with following words: (qp|qm|qbp|qbm|gp|gm|p|m). Returns a list."""
        return self._helconf

    @property
    def shelconf(self):
        return "".join(self.helconf)

    @helconf.setter
    def helconf(self, helconf):
        pParticlesIDs = re.compile("([q|G][2-4]*p|[q|G][2-4]*m|[q|G]b[2-4]*p|[q|G]b[2-4]*m|gp|gm|p|m|yp|ym)")
        if type(helconf) in [str, unicode]:
            self._helconf = pParticlesIDs.findall(helconf)
            assert "".join(self._helconf) == helconf
        elif type(helconf) is list:
            self._helconf = helconf
        else:
            raise Exception("Unsupported type for helconf: {}.".format(helconf))
        self._helconf = [entry.replace("gp", "p").replace("gm", "m") for entry in self._helconf]
        assert all([pParticlesIDs.findall(entry) is not [] for entry in self._helconf])
        self.A0 = BH.TreeHelAmpl(self.BH_process)

    @property
    def BH_process(self):
        return BH.process(*map(eval, ["BH.cvar." + entry for entry in self.helconf]))

    @property
    def BH_vectori(self):
        return BH.vectori([i for i in range(1, len(self.helconf) + 1)])

    @property
    def loopid(self):
        """Id of the particle in the loop. An helicity configuration needs to be set beforehand."""
        return self._loopid

    @loopid.setter
    def loopid(self, loopid):
        assert loopid in ["G", "nf", "RT", "RRT", "LRT", "LT", "LLT", "LC", "SLC", None]
        self._loopid = loopid
        if self.loopid is not None:
            BH.G = BH.glue
            BH.LC = BH.leading_color
            BH.SLC = BH.sub_leading_color
            self.A1 = BH.One_Loop_Helicity_Amplitude(self.BH_process, eval("BH.{}".format(self.loopid)))
            self.A1_rational = self.A1.rational_part()
            with OutputGrabber():
                self.A1_cutpart = self.A1.cut_part().makeDarrenCutPart()

    @property
    def helconf_and_loopid(self):
        if self.loopid is None:
            return self.shelconf
        else:
            return self.shelconf + "_" + self.loopid

    @property
    def multiplicity(self):
        return len(self.helconf)

    @property
    def amppart(self):
        return self._amppart

    @amppart.setter
    def amppart(self, amppart):
        assert amppart in ["tree", "box", "triangle", "bubble", "rational", None]
        self._amppart = amppart
        if self.amppart == "box":
            self.max_ampindex = int(self.A1_cutpart.nbr_boxes())
        elif self.amppart == "triangle":
            self.max_ampindex = int(self.A1_cutpart.nbr_triangles())
        elif self.amppart == "bubble":
            self.max_ampindex = int(self.A1_cutpart.nbr_bubbles())
        if amppart in ["box", "triangle", "bubble", None]:
            self.max_box_index = int(self.A1_cutpart.nbr_boxes())
            self.max_triangle_index = int(self.A1_cutpart.nbr_triangles())
            self.max_bubble_index = int(self.A1_cutpart.nbr_bubbles())

    @property
    def ampindex(self):
        return self._ampindex

    @ampindex.setter
    def ampindex(self, ampindex):
        if type(ampindex) is int:
            self._ampindex = ampindex
        elif type(ampindex) in [str, unicode]:
            assert ampindex.isdigit()
            self._ampindex = int(ampindex)
        elif ampindex is None:
            self._ampindex = None
        else:
            raise Exception("Unsupported ampindex: {}.".format(ampindex))
        if hasattr(self, "max_ampindex") and (self.ampindex > self.max_ampindex or self.ampindex <= 0) and self.ampindex is not None:
            raise Exception("Amplitude index out of bounds ({}/{}).".format(self.ampindex, self.max_ampindex))

    @property
    def call_cache_path(self):
        return settings.base_cache_path + "call_caches"

    @property
    def amppart_and_ampindex(self):
        return "{}({})".format(self.amppart, self.ampindex) if self.ampindex is not None and self.amppart is not None else self.amppart if self.amppart is not None else 'None'

    @property
    def count_gluons(self):
        return sum([1 for entry in self.helconf if entry in ["gp", "gm", "p", "m"]])

    @property
    def count_quarks(self):
        return len(re.findall(pattern="[q|G][2-4]*p|[q|G][2-4]*m|[q|G]b[2-4]*p|[q|G]b[2-4]*m", string="".join(self.helconf)))

    @property
    def count_photons(self):
        return sum([1 for entry in self.helconf if entry in ["yp", "ym"]])

    @property
    def short_process_name(self):
        short_process_name = ""
        if self.count_photons > 0:
            short_process_name += "{}y".format(self.count_photons)
        if self.count_quarks > 0:
            short_process_name += "{}q".format(self.count_quarks)
        if self.count_gluons > 0:
            short_process_name += "{}g".format(self.count_gluons)
        return short_process_name

    @property
    def process_name(self):
        return "_".join([self.short_process_name, self.helconf_and_loopid])

    @property
    def __name__(self):
        return "/".join([self.short_process_name, self.helconf_and_loopid, self.amppart_and_ampindex])

    @property
    def upper_res_path(self):
        upper_path = settings.base_res_path + "/".join([self.short_process_name, self.helconf_and_loopid])
        if not os.path.exists(upper_path + "/cpp"):
            os.makedirs(upper_path + "/cpp")
        if not os.path.exists(upper_path + "/math"):
            os.makedirs(upper_path + "/math")
        return upper_path

    @property
    def res_path(self):
        return settings.base_res_path + self.__name__

    def _evaluate(self, amppart, ampindex, oParticles):
        # look up in call cache
        if amppart in ["tree", "rational"] and amppart + "_" + str(hash(oParticles)) in self.dCallCache:
            return self.dCallCache[str(amppart + "_" + str(hash(oParticles)))]
        elif amppart in ["box", "triangle", "bubble"] and amppart + "(" + str(ampindex) + ")" + "_" + str(hash(oParticles)) in self.dCallCache:
            return self.dCallCache[str(amppart + "(" + str(ampindex) + ")" + "_" + str(hash(oParticles)))]
        # else compute and save to cache
        else:
            mom_conf = Upload_Momentum_Configuration(oParticles)
        if amppart == "tree":
            res = cgmp_to_mpc(gmpTools.CGMP(self.A0.eval(mom_conf, self.BH_vectori)))
            self.dCallCache[str(amppart + "_" + str(hash(oParticles)))] = res
            return res
        elif amppart in ["box", "triangle", "bubble"]:
            res = cgmp_to_mpc(gmpTools.CGMP(getattr(self.A1_cutpart, amppart)(int(ampindex)).eval(mom_conf, self.BH_vectori)))
            self.dCallCache[str(amppart + "(" + str(ampindex) + ")" + "_" + str(hash(oParticles)))] = res
            return res
        elif amppart == "rational":
            res = cgmp_to_mpc(gmpTools.CGMP(self.A1_rational.eval(mom_conf, self.BH_vectori)))
            self.dCallCache[str(amppart + "_" + str(hash(oParticles)))] = res
            return res
        else:
            raise Exception("Invalid amppart in evaluate: {}".format(amppart))

    def __call__(self, oParticles):
        # lookup in call cache
        if self.amppart_and_ampindex + "_" + str(hash(oParticles)) in self.dCallCache:
            return self.dCallCache[str(self.amppart_and_ampindex + "_" + str(hash(oParticles)))]
        # else compute all, save to cache, return value
        else:
            if self.amppart in ["box", "triangle", "bubble", None]:
                mom_conf = Upload_Momentum_Configuration(oParticles)
                self.A1_cutpart.eval(mom_conf, self.BH_vectori)
                # save all to cache!
                if self.amppart is None:
                    self.dCallCache[str(self.amppart_and_ampindex + "_" + str(hash(oParticles)))] = "Cut evaluated"
                # cache boxes
                for i in self.independent_boxes:
                    self._evaluate("box", i, oParticles)
                # cache triangles
                for i in self.independent_triangles:
                    self._evaluate("triangle", i, oParticles)
                # cache bubbles
                for i in self.independent_bubbles:
                    self._evaluate("bubble", i, oParticles)
                if self.amppart is not None and str(self.amppart_and_ampindex + "_" + str(hash(oParticles))) in self.dCallCache:
                    return self.dCallCache[str(self.amppart_and_ampindex + "_" + str(hash(oParticles)))]
                elif self.amppart is not None:
                    # not in cache: not one of the independent topology models  ---  not supposed to run these under production
                    print("Warning: Calling uncached cut part!")
                    return cgmp_to_mpc(gmpTools.CGMP(getattr(self.A1_cutpart, self.amppart)(int(self.ampindex)).eval(mom_conf, self.BH_vectori)))
            if self.amppart in ["tree", "rational"]:
                return self._evaluate(self.amppart, self.ampindex, oParticles)

    def reload_call_cache(self):
        if not os.path.exists(self.call_cache_path):
            os.makedirs(self.call_cache_path)
        with MyShelf(self.call_cache_path + "/" + self.process_name, 'c') as persistentCallCache:
            dCallCache = multiprocessing.Manager().dict(dict(persistentCallCache))
        setattr(antares.core.bh_unknown, "dCallCache_" + self.process_name, dCallCache)
        self.dCallCache = dCallCache

    def save_call_cache(self):
        if len(self.dCallCache.keys()) > self.len_keys_call_cache_at_last_save:
            self.len_keys_call_cache_at_last_save = len(self.dCallCache.keys())
            with MyShelf(self.call_cache_path + "/" + self.process_name, 'c') as persistentCallCache:
                persistentCallCache.update(self.dCallCache)

    def __setstate__(self, helconf_loopid_amppart_ampindex):
        self.__init__(*helconf_loopid_amppart_ampindex)

    def __getstate__(self):
        return self.helconf, self.loopid, self.amppart, self.ampindex

    def __getattr__(self, attr):
        splitted = attr.split("_")
        if len(splitted) == 2 and splitted[0] in ["tree", "box", "triangle", "bubble", "rational"] and splitted[1].isdigit():
            return functools.partial(self._evaluate, splitted[0], splitted[1])
        else:
            raise AttributeError("BHUnknown object has no attribute {}".format(attr))

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

    def print_graphs(self):
        if not os.path.exists(self.upper_res_path + "/graphs"):
            os.makedirs(self.upper_res_path + "/graphs")
            BH.print_cut_part_graph(self.A1_cutpart, str(self.upper_res_path + "/graphs"))
            try:
                with open(os.devnull, 'wb') as devnull:
                    output = subprocess.check_call(map(str, ['make', '--directory', self.upper_res_path + "/graphs", 'all']), stdout=devnull, stderr=subprocess.STDOUT)
            except subprocess.CalledProcessError as e:
                output = e.output
                print(output)

    def print_topology_info(self):
        if not os.path.exists(self.upper_res_path + "/all_topology_info"):
            from antares.topologies.topology import topology_info
            topology_info = topology_info(self)
            if topology_info != "":
                with open(self.upper_res_path + "/all_topology_info", "w+") as oFile:
                    oFile.write(topology_info)

    def get_topology_info(self):
        self.independent_boxes = []
        self.independent_triangles = []
        self.independent_bubbles = []
        if os.path.isfile(self.upper_res_path + "/all_topology_info"):
            with open(self.upper_res_path + "/all_topology_info") as all_topology_info:
                for line in all_topology_info.readlines():
                    if line.split("_")[0] == "box":
                        self.independent_boxes += [line.split("[")[1].split("]")[0].split(",")[0]]
                    elif line.split("_")[0] == "triangle":
                        self.independent_triangles += [line.split("[")[1].split("]")[0].split(",")[0]]
                    elif line.split("_")[0] == "bubble":
                        self.independent_bubbles += [line.split("[")[1].split("]")[0].split(",")[0]]

    def get_call_cache(self):
        if "dCallCache_" + self.process_name not in dir(antares.core.bh_unknown):
            self.reload_call_cache()
        else:
            self.dCallCache = getattr(antares.core.bh_unknown, "dCallCache_" + self.process_name)
        self.len_keys_call_cache_at_last_save = len(self.dCallCache.keys())

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

    def get_analytical_boxes(self):
        if not hasattr(self, "_all_analytical_boxes"):
            from antares.topologies.topology import independent_topologies, conversion_rules
            from antares.terms.terms import LoadResults
            oBHUnknown = BHUnknown(helconf=self.helconf, loopid=self.loopid if self.loopid is not None else "G")
            all_analytical_boxes = []
            for box_topo in independent_topologies(oBHUnknown)[0]:
                for i, ampindex in enumerate(box_topo.models):
                    if i == 0:
                        oBH1 = BHUnknown(helconf=self.helconf, loopid=self.loopid if self.loopid is not None else "G", amppart="box", ampindex=ampindex)
                        loaded_res = LoadResults(oBH1.res_path, load_partial_results_only=False, silent=True)[0]
                        if type(loaded_res) is list and len(loaded_res) > 0:
                            oAnalyticalBox1 = loaded_res[0]
                            oAnalyticalBox1.multiplicity = oBH1.multiplicity
                            all_analytical_boxes += [oAnalyticalBox1]
                        else:
                            print("Warning: missing box topology model.")
                            break
                    else:
                        oBH2 = BHUnknown(helconf=self.helconf, loopid=self.loopid if self.loopid is not None else "G", amppart="box", ampindex=ampindex)
                        symmetries = conversion_rules(oBH1, oBH2)
                        oAnalyticalBox2 = oAnalyticalBox1.Image(symmetries[0])
                        all_analytical_boxes += [oAnalyticalBox2]
            self._all_analytical_boxes = all_analytical_boxes
        return self._all_analytical_boxes

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

    @property
    def easy_boxes(self):
        return [box for box in self.get_analytical_boxes() if len(box) == 1]

    @property
    def basis_functions(self):
        if not hasattr(self, "_basis_functions"):
            basis_functions = []
            # read boxes
            oBoxes = self.get_analytical_boxes()
            for oBox in oBoxes:
                stripped_den = set([jInv for j, jInv in enumerate(oBox[0].oDen.lInvs) if oBox[0].oDen.lExps[j] == 1])
                spurious_poles = set([inv for inv in stripped_den if "(" in inv])
                for j, basis_function in enumerate(basis_functions):
                    spurious_poles_of_basis_function = set([inv for inv in basis_function if "(" in inv])
                    if spurious_poles == spurious_poles_of_basis_function and any([inv in basis_function and inv not in spurious_poles for inv in stripped_den]):
                        basis_functions[j] = set(list(basis_functions[j]) + list(stripped_den))
                        break
                else:
                    basis_functions += [stripped_den]
            # clean and order
            oInvariants = Invariants(self.multiplicity, Restrict3Brackets=settings.Restrict3Brackets,
                                     Restrict4Brackets=settings.Restrict4Brackets, FurtherRestrict4Brackets=settings.FurtherRestrict4Brackets)
            basis_functions = map(list, basis_functions)
            basis_functions = [sorted(basis_function, key=lambda t: oInvariants.full.index(t)) for basis_function in basis_functions]
            basis_functions = [ibasis_function for i, ibasis_function in enumerate(basis_functions) if not any(
                [all([inv in jbasis_function for inv in ibasis_function]) for j, jbasis_function in enumerate(basis_functions) if i != j])]
            self._basis_functions = basis_functions
        return self._basis_functions

    @property
    def basis_functions_invs(self):
        return list(set(flatten(self.basis_functions)))

    @property
    def spurious_poles(self):
        if not hasattr(self, "_spurious_poles"):
            _spurious_poles = []
            for basis_function in self.basis_functions:
                _spurious_poles_of_basis_function = list(set([inv for inv in basis_function if "(" in inv]))
                _spurious_poles = list(set(_spurious_poles + _spurious_poles_of_basis_function))
            if _spurious_poles == [] and self.basis_functions == []:
                oInvariants = Invariants(self.multiplicity, Restrict3Brackets=settings.Restrict3Brackets,
                                         Restrict4Brackets=settings.Restrict4Brackets, FurtherRestrict4Brackets=settings.FurtherRestrict4Brackets)
                _spurious_poles = oInvariants.invs_3
            self._spurious_poles = _spurious_poles
        return self._spurious_poles


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def Upload_Momentum_Configuration(oParticles):
    Lambda, Lambdat = BH.LambdaRGMP, BH.LambdatRGMP
    C_Mom, Mom_Conf = BH.Cmomgmp, BH.mcgmp
    n = len(oParticles)
    cms = []
    for j in range(1, n + 1):
        a = Lambda(gmpTools.to_BH_CGMP(mpc_to_cgmp(oParticles[j].r_sp_d[0, 0])), gmpTools.to_BH_CGMP(mpc_to_cgmp(oParticles[j].r_sp_d[1, 0])))
        b = Lambdat(gmpTools.to_BH_CGMP(mpc_to_cgmp(oParticles[j].l_sp_d[0, 0])), gmpTools.to_BH_CGMP(mpc_to_cgmp(oParticles[j].l_sp_d[0, 1])))
        cms += [C_Mom(a, b)]
        mom_conf = Mom_Conf(*cms)
    return mom_conf


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
