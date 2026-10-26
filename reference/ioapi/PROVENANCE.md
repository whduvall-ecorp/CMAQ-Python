# Vendored reference: I/O API `POLY`

`poly.f` in this folder is the unchanged EDSS/Models-3 I/O API routine
`REAL FUNCTION POLY (XPT, XPTS, YPTS, NDEG)`, which CMAQ's
`CCTM/src/plrise/smoke/preplm.f` calls to interpolate ambient temperature and
wind to stack top. It is not part of the CMAQ source tree, so it is vendored
here solely as the differential-comparison reference for the Python port in
`cmaq_port/interpolation.py`.

## Provenance

- Upstream file: `ioapi/poly.f` in the I/O API source.
- Source URL: https://raw.githubusercontent.com/cjcoats/ioapi-3.2/master/ioapi/poly.f
- Embedded version string: `$Id: poly.f 219 2015-08-17 18:05:54Z coats $`.
- Copyright (C) 1992-2002 MCNC and Carlie J. Coats, Jr.;
  (C) 2003-2010 Baron Advanced Meteorological Systems.
- License: GNU Lesser General Public License version 2.1 (per the file header).

## Fidelity note

The file content is reproduced verbatim from the upstream source listed above.
Because it was transcribed through a text fetch rather than copied from a local
Git checkout, line endings may be normalized to this repository's convention;
the algorithm, statements, and ordering are unchanged. `cmaq_port.poly_validation`
records this file's SHA-256 in each run report. If/when the I/O API is built in
the Linux reference environment (workflow-plan Stage 1), re-point the validation
at that checkout's `poly.f` and compare hashes.

The routine is self-contained (`IMPLICIT NONE`, only REAL arguments, no
includes or modules), so it compiles standalone with the project-local GNU
Fortran toolchain used for the other CMAQ kernel comparisons.
