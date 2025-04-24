# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

### Changed

### Fixed

### Deprecated


## [0.7.0] - 2025-04-24

### Added

- Support for evaluating `Terms` that return a numpy array, e.g. for amplitude coefficients with massive fermions (open spin indices).
- Slicing as alternative to scalings, works for arbirary q-ring (not yet available in public version).

### Changed

- Refactored `Terms`, `Term`, `Numerator`, `Denominator` to be based on `Monomial` and `Polynomial` from syngular.
- Scalings functions accept invariants list as optional kwarg.
- Mass dimension and Phase weights can now be manually set (changed from property to setter/getter).
- `LoadResults` can now be passed `multiplicity` as kwarg.
- Simplified several functions related to `Terms` by using `__str__` and `__rstr__` rather than reimplementing them.

### Fixed

- Invariant patterns are provided by lips.


## [0.6.3] - 2025-01-30

### Changed

- `single_scalings` and `pair_scalings` raise an exception if run with finite fields.

### Fixed

- monomial parser for `Term` numerator no longer fails with linear masses.


## [0.6.2] - 2025-01-21

### Changed

- `Numerical_Methods.do_single_collinar_limits` and related functions now accept a list of invariants as optional `kwarg`.

### Fixed

- Fixed issue where `Terms.__rstr__` would fail to properly load a numerator with mass dependence.


## [0.6.1] - 2024-12-16

### Added

- `TermsList` class now replaces `basis` related functions.

### Fixed

- Fixed issue where `Particles` objects used in `Numerical_Methods.mass_dimensions` and `.phase_weights` would not all trigger an `lru_cache` hit, likely due to mutable nature of the object (despite hasing and equality being correctly implemented).


## [0.6.0] - 2024-06-18


[unreleased]: https://github.com/GDeLaurentis/antares/compare/v0.7.0...HEAD
[0.7.0]: https://github.com/GDeLaurentis/antares/releases/compare/v0.6.3...v0.7.0
[0.6.3]: https://github.com/GDeLaurentis/antares/releases/compare/v0.6.2...v0.6.3
[0.6.2]: https://github.com/GDeLaurentis/antares/releases/compare/v0.6.1...v0.6.2
[0.6.1]: https://github.com/GDeLaurentis/antares/releases/compare/v0.6.0...v0.6.1
[0.6.0]: https://github.com/GDeLaurentis/antares/releases/tag/v0.6.0
