#!/usr/bin/env python
# -*- coding: utf-8 -*-

#   ___ ___   _   _      _
#  / __| __| | | | |_ _ | |___ _  _____ __ ___ _
#  \__ \ _|  | |_| | ' \| / / ' \/ _ \ V  V / ' \
#  |___/___|  \___/|_||_|_\_\_||_\___/\_/\_/|_||_|


# Author: Giuseppe

import subprocess
import re
import mpmath

mpmath.mp.dps = 300


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def upload_mom_conf_joe(oParticles):
    m = "[["
    for oParticle in oParticles:
        m += "[[[[{}, {}]],[[{}, {}]]]],".format(complex(oParticle.r_sp_d[0, 0]), complex(oParticle.r_sp_d[1, 0]),
                                                 complex(oParticle.l_sp_d[0, 0]), complex(oParticle.l_sp_d[0, 1]))
    m = m[:-1] + "]]"
    m = m.replace("[[", "{").replace("]]", "}").replace("(", "").replace(")", "").replace("j", "*I").replace("e", "*^")
    return m


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


class SEUnknown(object):

    def __init__(self, se_arguement):
        self.argument = se_arguement

    def __setstate__(self):
        self.__init__()

    def __getstate__(self):
        return

    def __call__(self, oParticles, verbose=False):
        output = subprocess.check_output(["/home/gdl/Programs/JosephsSEsolver/Package/SEscript.m",
                                          "m={}".format(upload_mom_conf_joe(oParticles)),
                                          "argument=Amplitude[{}]".format(self.argument)])
        if verbose is True:
            print(output)
        result = re.findall(r"([ \+\-\d\.\*\^I]+)", output.split("Numerical result:\n")[1])[-1:][0].replace(" ", "")
        if result == "":
            return 0
        else:
            result = eval(result.replace("*^", "e").replace("*I", "*1j"))
            return mpmath.mpc(result.real, result.imag)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
