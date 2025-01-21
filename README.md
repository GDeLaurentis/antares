# <div align="center">Antares</div>
### <div align="center"><em>Automated Numerical To Analytical REconstruction Software</em></div>

<div align="center">

[![CI Lint](https://github.com/GDeLaurentis/antares-dev/actions/workflows/ci_lint.yml/badge.svg)](https://github.com/GDeLaurentis/antares-dev/actions/workflows/ci_lint.yml)
[![CI Test](https://github.com/GDeLaurentis/antares-dev/actions/workflows/ci_test.yml/badge.svg)](https://github.com/GDeLaurentis/antares-dev/actions/workflows/ci_test.yml)
[![Coverage](https://img.shields.io/badge/Coverage-37%25-red?labelColor=2a2f35)](https://github.com/GDeLaurentis/antares/actions)
[![Docs](https://github.com/GDeLaurentis/antares-dev/actions/workflows/cd_docs.yml/badge.svg?label=Docs)](https://gdelaurentis.github.io/antares-dev/)
[![DOI](https://zenodo.org/badge/902351393.svg)](https://doi.org/10.5281/zenodo.14501989)
[![Python](https://img.shields.io/pypi/pyversions/antares?label=Python)](https://pypi.org/project/antares/)

</div>

<div align="center">
The current version is meant to load and evaluate results only.<br>
Reconstruction routines coming soon!
</div>

<br>

Antares will provide tools to reconstruct (or simplify) analytical expressions for functions in the field of fractions over polynomial quotient rings, such as loop integral coefficients, from numerical evaluations. It is based on the study and fitting of pole residues in complexified momentum space. The divergence of expressions near the poles is handled via arbitrary-precision floating-point arithmetics or $p$-adic numbers. Alternatively finite fields can also be used (may involve semi-numerical slices).

## Installation
The package is available on the [Python Package Index](https://pypi.org/project/antares/)
```console
pip install antares-hep
```
Alternativelty, it can be installed by cloning the repo
```console
git clone https://github.com/GDeLaurentis/antares.git path/to/repo
pip install -e path/to/repo
```
