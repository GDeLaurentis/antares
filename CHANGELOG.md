# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- `TermsList` class now replaces `basis` related functions.

### Changed

### Fixed

- Fixed issue where `Particles` objects used in `Numerical_Methods.mass_dimensions` and `.phase_weights` would not all trigger an `lru_cache` hit, likely due to mutable nature of the object (despite hasing and equality being correctly implemented).

### Deprecated


## [0.6.0] - 2024-06-18


[unreleased]: https://github.com/GDeLaurentis/antares/compare/v0.6.0...HEAD
[0.6.0]: https://github.com/GDeLaurentis/antares/releases/tag/v0.6.0
