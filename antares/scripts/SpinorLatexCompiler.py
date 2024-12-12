#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
import subprocess
import os
import tempfile
import shutil
import argparse


def generate_pdf(tex, texfile, output_dir=None, interaction=None, verbose=False):
    """Genertates the pdf from string"""

    current = os.getcwd()
    if verbose:
        print("Running from: ", current)

    # Create temporary directory for LaTeX compilation
    temp = tempfile.mkdtemp()
    os.chdir(temp)

    if shutil.which('lualatex') is None:
        raise Exception("lualatex not installed")

    f = open('temp.tex', 'w')
    f.write("\\RequirePackage{luatex85}\n" + tex)
    f.close()

    # Prepare the lualatex command
    cmd = ['lualatex']
    if interaction:
        cmd.extend(['-interaction', interaction])
    cmd.append('temp.tex')

    if verbose:
        print("Running command: ", ' '.join(cmd))

    # Execute lualatex
    if verbose:
        proc = subprocess.Popen(cmd)
    else:
        proc = subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    proc.communicate()
    code = proc.wait()
    if verbose:
        print("Exit code: ", code)

    # Ensure output directory exists or default to current dir
    if output_dir is None:
        output_dir = current
    else:
        output_dir = os.path.join(current, output_dir)  # Make relative to where the script was run
        os.makedirs(output_dir, exist_ok=True)

    # Rename and move the generated PDF
    pdfname = texfile[:-4] + ".pdf"
    if "/" in pdfname:
        pdfname = re.split(r'/', pdfname)[-1]

    shutil.move('temp.pdf', os.path.join(output_dir, pdfname))
    if verbose:
        print(f"PDF saved to: {os.path.join(output_dir, pdfname)}")

    # Clean up the temporary directory
    shutil.rmtree(temp)


def main():

    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Compile LaTeX files to PDF.')
    parser.add_argument('texfile', type=str, help='The LaTeX file to compile.')
    parser.add_argument('-o', '--output-directory', type=str, help='Output directory for the PDF.')
    parser.add_argument('-i', '--interaction', type=str, help='Set interaction mode for lualatex.', default='errorstopmode')
    parser.add_argument('-v', '--verbose', action='store_true', help='Print extra compilation info.')

    args = parser.parse_args()

    if not os.path.isfile(args.texfile):
        raise Exception(f"File {args.texfile} does not exist.")

    result = ""

    with open(args.texfile, "r") as file:
        for line in file:
            line = line.replace("*", "")  # suppress * for products
            if r"\scriptscriptstyle" not in line:
                result += line
            else:
                splitted = re.split(r'(?<!\(|{|}|\))(\+|\-)(?=\d+[/s⟨\[]|\d+i)', line)
                # pprint(splitted)
                if len(splitted) < 9:
                    result += line
                    continue
                splitted = splitted[:-2] + [splitted[-2] + splitted[-1]]
                nbr_terms = len(splitted[1:-1]) / 2
                result += splitted[0] + "\\,...⟪" + str(int(nbr_terms)) + r"\,\text{terms}" + "⟫...\\, " + splitted[-1]

    # symmetries display
    pSym_fp = re.compile(r"\((\d*),\\;\\text{False}\)")
    pSym_fm = re.compile(r"\((\d*),\\;\\text{False},\\;-\)")
    pSym_tp = re.compile(r"\((\d*),\\;\\text{True}\)")
    pSym_tm = re.compile(r"\((\d*),\\;\\text{True},\\;-\)")

    if pSym_fp.findall(result) != [] or pSym_fm.findall(result) != [] or pSym_tp.findall(result) != [] or pSym_tm.findall(result) != []:
        if len(pSym_fp.findall(result)) > 0:
            multiplicity = len(pSym_fp.findall(result)[0])
        elif len(pSym_fm.findall(result)) > 0:
            multiplicity = len(pSym_fm.findall(result)[0])
        elif len(pSym_tp.findall(result)) > 0:
            multiplicity = len(pSym_tp.findall(result)[0])
        elif len(pSym_tm.findall(result)) > 0:
            multiplicity = len(pSym_tm.findall(result)[0])
        base_permutation = "".join(map(str, range(1, multiplicity + 1)))

        result = pSym_fp.sub(r"(" + str(base_permutation) + r"\\; \\rightarrow \\; \1)", result)
        result = pSym_fm.sub(r"(" + str(base_permutation) + r"\\; \\rightarrow \\; -\1)", result)
        result = pSym_tp.sub(r"(" + str(base_permutation) + r"\\; \\rightarrow \\; \\overline{\1})", result)
        result = pSym_tm.sub(r"(" + str(base_permutation) + r"\\; \\rightarrow \\; -\\overline{\1})", result)

    if args.verbose:
        print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        print(f"compiling: {args.texfile}")
        print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        print("")
        print(result)
        print("")
        print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    else:
        print(f"compiling: {args.texfile}")

    file.close()

    generate_pdf(result, args.texfile, output_dir=args.output_directory, interaction=args.interaction, verbose=args.verbose)


if __name__ == "__main__":

    main()
