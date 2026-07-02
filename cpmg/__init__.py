###
# Copyright (c) 2025-2026 Prof. Dr. Donghan Lee, Korea Basic Science Institute (KBSI)
# Original code: Adam Mazur et al., MPI Goettingen, 2012
###
"""
cpmg — CPMG dispersion curve analysis package.

This package provides data structures and fitting models for analysing
NMR CPMG (Carr-Purcell-Meiboom-Gill) R2 relaxation-dispersion experiments.
It supports three exchange models:

  * Matrix  — exact Liouvillian matrix-exponential solution (any exchange regime)
  * London  — analytical slow-exchange approximation (London et al.)
  * Meiboom — analytical fast-exchange approximation (Meiboom et al.)

Modules
-------
cpmg          : Data containers (DispersionCurve, Residue, CPMGDataSet)
model_2state  : Fitting engine and R2 calculation models (CPMG_model)

Authors
-------
Prof. Dr. Donghan Lee
Korea Basic Science Institute (KBSI)

Original code:
Adam Mazur, Bjoern Hammesfahr, Christian Griesinger,
Donghan Lee, Martin Kollmar
Max-Planck-Institute for Biophysical Chemistry, Goettingen, 2012

Copyright (c) 2025-2026 Prof. Dr. Donghan Lee, Korea Basic Science Institute (KBSI)
"""

__version__ = '1.2.0'
