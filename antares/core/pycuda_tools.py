#   ___       ___        _     _____         _
#  | _ \_  _ / __|  _ __| |__ |_   _|__  ___| |___
#  |  _/ || | (_| || / _` / _` || |/ _ \/ _ \ (_-<
#  |_|  \_, |\___\_,_\__,_\__,_||_|\___/\___/_/__/
#       |__/

# Author: Giuseppe
# Created: 24/01/2019

import math
import re


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def cuda_set_vars_and_get_funcs(path_to_cuda_script=None, **kwargs):
    # Set Variables
    oFile = open(path_to_cuda_script, "r")
    sFile = oFile.read()
    for var_name in kwargs:
        sFile = sFile.replace("{" + var_name + "}", str(kwargs[var_name]))
    oFile.close()

    # Get Exectuable
    sCommand = "mod = SourceModule(str(\"\"\"{}\"\"\"))\n\n".format(sFile)
    GlobalFunctions = sCommand.split("GLOBAL FUNCTIONS")[1].split(r"/*!!!  IMPLEMENTATION ")[0]
    last_preprocessor_condition = True
    for line in GlobalFunctions.split("\n"):
        if line == "":
            continue
        elif line[0:3] == "#if":
            condition = line.replace("#if", "")
            last_preprocessor_condition = eval(re.sub("([A-Z_]+)", lambda match: str(kwargs.get(match.group())), condition))
        elif line == "#else":
            last_preprocessor_condition = not last_preprocessor_condition
        elif line == "#endif":
            last_preprocessor_condition = True
        elif last_preprocessor_condition is True:
            FunctionName = re.findall(r"(?:int)|(?:void) (.*) (?:\(.*)", line)[0]
            sCommand += "Cuda{} = mod.get_function(str(\"{}\"))\n".format(FunctionName, FunctionName)

    return sCommand


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def number_of_foldings(x, FoldingMaxLength=1024):
    return int(math.ceil(x / float(FoldingMaxLength)))


def folded_number_of_columns(x, FoldingMaxLength=1024):
    return int(math.ceil((x / float(number_of_foldings(x, FoldingMaxLength)))))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
