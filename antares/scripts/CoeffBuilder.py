#!/usr/bin/env python
# -*- coding: utf-8 -*-

from sympy import pprint, factorint  # noqa
from collections import defaultdict
from fractions import Fraction
from math import log
import matplotlib.pyplot as plt
import codecs
import re
import sys


possible_denominators = [1, 2, 3, 4, 6, 8, 9, 12, 16, 18, 24, 27, 32, 36, 48, 54, 64, 72, 81, 96, 108, 128, 144, 162, 192, 216, 243, 256, 288, 324, 384, 432, 486, 512, 576,
                         648, 729, 768, 864, 972, 1024, 1152, 1296, 1458, 1536, 1728, 1944, 2048, 2187, 2304, 2592, 2916, 3072, 3456, 3888, 4096, 4374, 4608, 5184, 5832, 6144,
                         6561, 6912, 7776, 8192, 8748, 9216, 10368, 11664, 12288, 13122, 13824, 15552, 16384, 17496, 18432, 19683, 20736, 23328, 24576, 26244, 27648, 31104,
                         32768, 34992, 36864, 39366, 41472, 46656, 49152, 52488, 55296, 59049, 62208, 65536, 69984, 73728, 78732, 82944, 93312, 98304, 104976, 110592,
                         118098, 124416, 131072, 139968, 147456, 157464, 165888, 177147, 186624, 196608, 209952, 221184, 236196, 248832, 262144, 279936, 294912, 314928,
                         331776, 354294, 373248, 393216, 419904, 442368, 472392, 497664, 524288, 531441, 559872, 589824, 629856, 663552, 708588, 746496, 786432, 839808,
                         884736, 944784, 995328, 1048576, 1062882, 1119744, 1179648, 1259712, 1327104, 1417176, 1492992, 1572864, 1594323, 1679616, 1769472, 1889568, 1990656,
                         2097152, 2125764, 2239488, 2359296, 2519424, 2654208, 2834352, 2985984, 3145728, 3188646, 3359232, 3538944, 3779136, 3981312, 4194304, 4251528, 4478976,
                         4718592, 4782969, 5038848, 5308416, 5668704, 5971968, 6291456, 6377292, 6718464, 7077888, 7558272, 7962624, 8388608, 8503056, 8957952, 9437184, 9565938]


def view_by_frequency(dictionary):
    dictionary_view = [(v, k) for k, v in dictionary.iteritems()]
    dictionary_view.sort(reverse=True)
    for v, k in dictionary_view:
        print("%s: %d" % (k, v))


if __name__ == "__main__":

    if False:
        with open("/home/gdl/Desktop/tmux_clipboard_1", 'r') as content_file:
            replacement_content = str(content_file.read().decode('utf8'))
        replacement_content = replacement_content.split("\n")
        for i, entry in enumerate(replacement_content):
            if entry == "":
                replacement_content[i] = None
            elif entry[0:2] == "In":
                replacement_content[i] = None
            elif entry[0:2] == "Al":
                replacement_content[i] = "0"
            else:
                entry = entry.replace("Coeff. of this whole term is: ", "")
                replacement_content[i] = entry
        replacement_content = filter(None, replacement_content)
        # print replacement_content, len(replacement_content)
    else:
        replacement_content = None

    with open("/home/gdl/Desktop/tmux_clipboard", 'r') as content_file:
        content = str(content_file.read().decode('utf8'))

    # SPLIT SPINOR STRING FROM COEFFICIENTS

    content = content.replace("Coeff. of ", "")
    content = content.replace(" ", "")
    content = content.split("\n")

    for i, entry in enumerate(content):
        content[i] = entry.split(":")

    # REFINE COEFFICIENTS
    failed = -1

    numerators, denominators, to_be_refitted_terms, to_be_refitted_exps = [], [], [], []
    for i, entry in enumerate(content):
        coefficient = entry[-1]
        coefficient = re.split("(?<=...)([-,+])", coefficient)
        if len(coefficient) == 3:
            coefficient = [coefficient[0], coefficient[1] + coefficient[2]]
        if "*I" in coefficient[0]:
            coefficient = coefficient[0]
        else:
            real_part = coefficient[0].split("/")
            real_part = float(real_part[0]) / float(real_part[1])
            if real_part > 10 ** -6:
                # raise Exception("Real part isn't zero.")
                pass
            coefficient = coefficient[1]
            content[i][-1] = coefficient
        coefficient = coefficient.replace("*I", "")
        numerator, denominator = coefficient.split("/") if len(coefficient.split("/")) == 2 else (coefficient, 1)
        numerator, denominator = int(numerator), int(denominator)

        criterion = any([key > 3 for key in factorint(denominator).keys()])
        if criterion is True:
            new_fraction = Fraction(int(round(numerator / float(denominator) * 12288)), 12288)
            # distance between old and new fraction
            error1 = abs(new_fraction - (numerator / float(denominator)))
            # how far the new numerator is from being an integer
            error2 = abs(abs(numerator / float(denominator) * new_fraction.denominator) - abs(round(numerator / float(denominator) * new_fraction.denominator)))
            print("Re fractioning:", numerator, "/", denominator, " to:", new_fraction, " with error:", max([error1, error2]), end="")
            if max([error1, error2]) > 10 ** -5:
                new_fraction = Fraction(int(round(numerator / float(denominator) * 147456)), 147456)
                # distance between old and new fraction
                error1 = abs(new_fraction - (numerator / float(denominator)))
                # how far the new numerator is from being an integer
                error2 = abs(abs(numerator / float(denominator) * new_fraction.denominator) - abs(round(numerator / float(denominator) * new_fraction.denominator)))
                print("\rRe fractioning:", numerator, "/", denominator, " to:", new_fraction, " with error:", max([error1, error2]), end="")
                if max([error1, error2]) > 10 ** -5:
                    if replacement_content is not None:
                        failed += 1
                        coefficient = replacement_content[failed].replace("*I", "")
                        numerator, denominator = coefficient.split("/") if len(coefficient.split("/")) == 2 else (coefficient, 1)
                        numerator, denominator = int(numerator), int(denominator)
                        new_fraction = Fraction(int(numerator), int(denominator))
                        criterion = any([key > 3 for key in factorint(denominator).keys()])
                        print("\rTrying from replacement_content: {}.                                                     ".format(new_fraction), end="")
                    else:
                        criterion = True
                    if criterion is True:
                        print("Needs to be refitted.             ")
                        term = content[i][0]
                        ListOfSpinors = [term[k:k + 5] for k in range(0, len(term), 5)]
                        compact_terms = list(set(ListOfSpinors))
                        compact_exps = [ListOfSpinors.count(entry) for entry in compact_terms]
                        to_be_refitted_terms += [["\"" + entry + "\"" for entry in compact_terms]]
                        to_be_refitted_exps += [compact_exps]
                        content[i] = None
                    else:
                        print("                                   ")
                        content[i][-1] = str(new_fraction) + "*I"
                else:
                    print("                                   ")
                    content[i][-1] = str(new_fraction) + "*I"
            else:
                print("                                   ")
                content[i][-1] = str(new_fraction) + "*I"

        numerators += [numerator]
        denominators += [denominator]

    for term in to_be_refitted_terms:
        print("[[", end="")
        for entry in term:
            sys.stdout.write(entry + ", ",)
        print("\033[2D" + "]],")
    for entry in to_be_refitted_exps:
        print(str([entry]) + ",")
    content = filter(None, content)

    denominator_frequencies = defaultdict(int)
    primes_frequencies = defaultdict(int)
    for denominator in denominators:
        denominator_frequencies[denominator] += 1
        prime_factors = factorint(denominator)
        for prime_factor in prime_factors.keys():
            primes_frequencies[prime_factor] += prime_factors[prime_factor]

    ys = [denominator_frequencies[denominator] for denominator in denominators]
    xs = [min([primes_frequencies[prime] for prime in factorint(denominator)]) if len(factorint(denominator)) > 0 else 1 for denominator in denominators]
    xs = map(log, xs)
    plt.scatter(xs, ys)
    plt.xlabel('min(primes frequencies)')
    plt.ylabel('denominator frequencies')
    for i, denominator in enumerate(denominators):
        plt.annotate(denominator, (xs[i], ys[i]))
    plt.show()

    # BUILD THE RESULT

    for i, entry in enumerate(content):
        if "I" in entry[-1]:
            content[i] = [entry[-1]] + entry[:-1]

    for i, entry in enumerate(content):
        content[i] = "".join(entry)

    content = "+".join(content)
    content = content.replace("*I", "i")
    content = content.replace("+-", "-")
    content = content.replace("++", "+")
    content = content.replace("|", "")

    i = 0
    while i < len(content):
        char = content[i]
        if char == "⟨" or char == "[":
            # pprint(content)
            current_variable = content[i:i + 4]
            power_counter = 1
            while True:
                if content[i + 4 * power_counter:i + 4 * (power_counter + 1)] == current_variable:
                    power_counter += 1
                else:
                    break
            if power_counter > 1:
                content = content[0:i] + current_variable + "^" + str(power_counter) + content[i + 4 * power_counter:]
        i += 1

    file = codecs.open("/home/gdl/Desktop/tmux_clipboard_new", "w", "utf-8")
    file.write(content)
    file.close()
