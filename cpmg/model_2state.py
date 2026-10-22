###
# Copyright (c) 2025-2026 Prof. Dr. Donghan Lee, Korea Basic Science Institute (KBSI)
# Original code: Adam Mazur et al., MPI Goettingen, 2012
###
"""
model_2state.py — Two-state exchange fitting models for NMR CPMG data.

This module implements the CPMG_model class, which wraps:
  * Three R2 calculation models (Matrix, London, Meiboom)
  * A nonlinear least-squares fitting engine (scipy.optimize.leastsq)
  * Log-file and PDF report generation
  * Conversion of Matrix/London results to Meiboom representation

Exchange Models
---------------
Matrix  : Exact Liouvillian matrix-exponential solution valid for any
          exchange regime.  Uses scipy.linalg.expm.
London  : Analytical slow-exchange approximation (London et al.).
Meiboom : Analytical fast-exchange approximation (Meiboom et al.);
          fits phi = pB*(1-pB)*delta_omega^2 instead of delta_omega.

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

import sys
import numpy as np
from numpy import array, sqrt, pi, dot, real, log, sin, sinh, cos, cosh, exp, zeros, sort, tanh, isfinite
from scipy.optimize import leastsq

from . import __version__
from . import report
from .cpmg import CPMGDataSet
from . import models
from . import fitting
from .fitting import _warn_if_not_converged
from .analysis import derived_std as _derived_std
from .models import _echo_count
from . import analysis
from scipy.linalg import expm          # matrix exponential
from numpy.linalg import matrix_power  # integer matrix power
from numpy import diag




class CPMG_model:
    """Two-state chemical-exchange model for CPMG R2 dispersion fitting.

    This class manages the full analysis pipeline:
      1. Hold a CPMGDataSet with all experimental dispersion curves.
      2. Provide R2 calculation functions for Matrix, London, and Meiboom models.
      3. Run nonlinear least-squares optimisation (scipy leastsq).
      4. Report results as a formatted log string, a multi-page PDF, and JSON.

    Global exchange-rate parameters (kAB, kBA, kex) are shared across all
    residues.  Per-residue parameters (dd / phi, R2_0) are stored on each
    Residue object.

    Attributes
    ----------
    kAB, kBA   : float  — forward/reverse rate constants (s^-1), Matrix/London
    kAB_std, kBA_std : float — their standard deviations
    kex, kex_std     : float — total exchange rate and std for Meiboom model
    dataset    : CPMGDataSet — all loaded experimental data
    model      : str  — 'Matrix', 'London', or 'Meiboom'
    verbose    : bool — print chi² progress during fitting when True
    """

    def __init__(self):
        # Program identification banner shown in the log file
        self.programName = ('ShereKhan\n'
                            'CPMG program\n'
                            'ver. ' + __version__ + '\n'
                            '\nProf. Dr. Donghan Lee\n'
                            'Korea Basic Science Institute (KBSI)\n'
                            '\nOriginal code: Adam Mazur, Bjoern Hammesfahr, '
                            'Christian Griesinger, Donghan Lee, Martin Kollmar\n'
                            'Max-Planck-Institute for Biophysical Chemistry, '
                            'Goettingen, 2012\n'
                            '\n\nPlease cite: Mazur A, Hammesfahr B, Griesinger C, Lee D, Kollmar M. '
                            'ShereKhan - calculating exchange parameters in relaxation dispersion '
                            'data from CPMG experiments. Bioinformatics 29(14):1819-1820 (2013). '
                            'doi:10.1093/bioinformatics/btt286\n')

        # Global exchange-rate parameters shared across all residues
        self.kAB = 0.0      # A→B rate constant (s^-1)
        self.kAB_std = 0.0  # standard deviation of kAB
        self.kBA = 0.0      # B→A rate constant (s^-1)
        self.kBA_std = 0.0  # standard deviation of kBA
        self.kex = 0.0      # total exchange rate kAB + kBA (Meiboom model)
        self.kex_std = 0.0  # standard deviation of kex
        self.covar = None   # covariance matrix of the last global fit

        self.dataset = CPMGDataSet()  # container for all experimental data

        self.model = ''      # 'Matrix', 'London', or 'Meiboom'
        self.verbose = False # set True to print chi² at each function evaluation

    # -------------------------------------------------------------------------
    # Parameter initialisation helpers
    # -------------------------------------------------------------------------

    def init_dd(self, dd_AB):
        return fitting.init_dd(self, dd_AB)

    def init_R2_0(self):
        return fitting.init_R2_0(self)

    # -------------------------------------------------------------------------
    # Parameter vector packing / unpacking
    # -------------------------------------------------------------------------

    def getLocalParams(self, p):
        return fitting.get_local_params(self, p)

    def setLocalParams(self, p, covar=None):
        return fitting.set_local_params(self, p, covar)

    # -------------------------------------------------------------------------
    # R2 calculation models
    # -------------------------------------------------------------------------

    def R2_clc_matrix(self, kAB, kBA, domega, R2_0, vCPMG, tcp):
        return models.r2_matrix(kAB, kBA, domega, R2_0, vCPMG, tcp)

    def R2_clc_matrix2(self, kAB, kBA, domega, R2_0, vCPMG, tcp):
        return models.r2_matrix_eig(kAB, kBA, domega, R2_0, vCPMG, tcp)

    def R2_clc_london(self, kAB, kBA, domega, R2_0, vCPMG):
        return models.r2_london(kAB, kBA, domega, R2_0, vCPMG)

    def R2_clc_meiboom(self, kex, phi, R2_0, vCPMG, field):
        return models.r2_meiboom(kex, phi, R2_0, vCPMG, field)

    def _cache_arrays(self):
        return fitting.cache_arrays(self)

    def _build_batch_plan(self):
        return fitting.build_batch_plan(self)

    # -------------------------------------------------------------------------
    # Vectorised R2 calculation models (accept v_array, return array)
    # -------------------------------------------------------------------------

    def R2_clc_matrix2_batch(self, kAB, kBA, domega_arr, R2_0_arr,
                             v_array, delta_arr, n_arr, tcp):
        return models.r2_matrix_eig_batch(kAB, kBA, domega_arr, R2_0_arr, v_array, delta_arr, n_arr, tcp)

    def R2_clc_london_vec(self, kAB, kBA, domega, R2_0, v_array):
        return models.r2_london_vec(kAB, kBA, domega, R2_0, v_array)

    def R2_clc_meiboom_vec(self, kex, phi, R2_0, v_array, field):
        return models.r2_meiboom_vec(kex, phi, R2_0, v_array, field)

    def _calc_r2(self, r, dspCurve, R2_0, v_array):
        return report.calc_r2(self, r, dspCurve, R2_0, v_array)

    # -------------------------------------------------------------------------
    # Error function (called by leastsq)
    # -------------------------------------------------------------------------

    def errFunc(self, p):
        return fitting.err_func(self, p)

    # -------------------------------------------------------------------------
    # Fitting routines
    # -------------------------------------------------------------------------

    def fit(self):
        return fitting.fit(self)

    # -------------------------------------------------------------------------
    # Individual fitting and AIC-based model comparison
    # -------------------------------------------------------------------------

    def _build_p0_current(self):
        return fitting.build_p0_current(self)

    def fit_individual(self):
        return analysis.fit_individual(self)

    _aic = staticmethod(analysis.aic)

    def compareModelsAIC(self):
        return analysis.compare_models_aic(self)

    def _refresh_global_stats(self):
        return analysis.refresh_global_stats(self)

    def _leave_one_out_kex(self):
        return analysis.leave_one_out_kex(self)

    def _per_residue_global_chi2(self):
        return analysis.per_residue_global_chi2(self)

    def jackknifeGlobal(self):
        return analysis.jackknife_global(self)

    def jackknifeReport(self, jk):
        return report.jackknife_report(jk)

    def aicReport(self, cmp):
        return report.aic_report(self, cmp)

    def _build_grid_p0(self, kex_val, pB_val, csd_val):
        return fitting.build_grid_p0(self, kex_val, pB_val, csd_val)

    def initGuessAll(self):
        return fitting.init_guess_all(self)

    # -------------------------------------------------------------------------
    # Reporting
    # -------------------------------------------------------------------------

    def _kex_pB(self):
        return analysis.kex_pB(self)

    def getLogBuffer(self, config_path=None):
        return report.log_buffer(self, config_path)

    def pdf(self, pdfFileName):
        return report.write_pdf(self, pdfFileName)

    # -------------------------------------------------------------------------
    # Model conversion
    # -------------------------------------------------------------------------

    def toMeiboom(self):
        return analysis.to_meiboom(self)

    # -------------------------------------------------------------------------
    # JSON result export
    # -------------------------------------------------------------------------

    def reportAllValues(self):
        return report.report_all_values(self)
