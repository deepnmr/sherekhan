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
from .models import _echo_count
from scipy.linalg import expm          # matrix exponential
from numpy.linalg import matrix_power  # integer matrix power
from numpy import diag


def _derived_std(kAB, kBA, dd, C3):
    """Propagate the fitted covariance to kex, pB and phi = pB*(1-pB)*dd^2.

    Args:
        kAB, kBA (float): Fitted exchange rates.
        dd       (float): Fitted chemical shift difference of the residue.
        C3       (3x3 array): Covariance of (kAB, kBA, dd).

    Returns:
        tuple: (kex_std, pB_std, phi_std), each clamped to >= 0 before the root.
    """
    C3 = np.asarray(C3, dtype=float)
    C2 = C3[:2, :2]
    k = kAB + kBA
    p = kAB / k
    g_kex = np.array([1.0, 1.0])
    g_pB = np.array([kBA / k ** 2, -kAB / k ** 2])
    g_phi = np.array([dd ** 2 * (1.0 - 2.0 * p) * g_pB[0],
                      dd ** 2 * (1.0 - 2.0 * p) * g_pB[1],
                      2.0 * p * (1.0 - p) * dd])
    kex_var = np.dot(g_kex, np.dot(C2, g_kex))
    pB_var = np.dot(g_pB, np.dot(C2, g_pB))
    phi_var = np.dot(g_phi, np.dot(C3, g_phi))
    return sqrt(max(kex_var, 0.0)), sqrt(max(pB_var, 0.0)), sqrt(max(phi_var, 0.0))


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
        """Fit every active residue independently (its own exchange parameters).

        In the global fit all residues share one exchange process (single kex /
        kAB,kBA).  Here each residue is fitted alone, so it gets its own
        exchange rate in addition to its own dd and R2_0.  Used as the
        alternative model in the AIC comparison (compareModelsAIC).

        The current global solution is used to warm-start each individual fit
        and is fully restored before returning, so global reporting downstream
        is unaffected.

        Returns:
            list of dict: one entry per active residue with keys
                'label', 'chi2', 'nvar', 'npar', and the fitted exchange
                parameters ('kex', 'pB', 'dd' for Matrix/London; 'kex', 'phi'
                for Meiboom) plus 'R2_0' (list).
        """
        rsds = self.dataset.rsds

        # --- Snapshot global state so we can restore it afterwards ---
        snap_active = [r.active for r in rsds]
        snap_glob   = (self.kAB, self.kBA, self.kex)
        snap_local  = [(r.dd, getattr(r, 'dd_std', 0.0),
                        list(r.R2_0), list(r.R2_0_std)) for r in rsds]
        snap_verbose = self.verbose
        self.verbose = False   # individual fits would otherwise spam chi² lines

        active_idx = [i for i, r in enumerate(rsds) if snap_active[i]]
        results = []

        for i in active_idx:
            # Activate only residue i
            for k, r in enumerate(rsds):
                r.active = (k == i)
            self._cache_arrays()   # rebuild batch plan for the single residue

            p0    = self._build_p0_current()
            out   = leastsq(self.errFunc, x0=p0, full_output=1)
            p1    = out[0]
            covar = out[1]                     # None when the fit is singular
            _warn_if_not_converged(out)
            fvec  = out[2]['fvec']             # residuals at the solution
            chi2  = float(np.dot(fvec, fvec))
            npar  = len(p1)
            nvar  = len(fvec)

            # Decode fitted parameters for this residue
            r = rsds[i]
            rec = {'label': r.label, 'chi2': chi2, 'nvar': nvar, 'npar': npar}
            if self.model == 'Meiboom':
                kex = p1[0]
                rec['kex'] = kex
                rec['kex_std'] = sqrt(covar[0][0]) if covar is not None else 0.0
                rec['phi'] = p1[1]
                rec['R2_0'] = list(p1[2:])
            else:
                kAB, kBA = p1[0], p1[1]
                kex = kAB + kBA
                rec['kex'] = kex
                # kex = kAB + kBA → var(kex) = var(kAB) + var(kBA) + 2 cov(kAB,kBA)
                if covar is not None:
                    rec['kex_std'] = sqrt(max(
                        covar[0][0] + covar[1][1] + 2.0 * covar[0][1], 0.0))
                else:
                    rec['kex_std'] = 0.0
                rec['pB']  = kAB / kex if kex != 0 else float('nan')
                rec['dd']  = p1[2]
                rec['R2_0'] = list(p1[3:])
            results.append(rec)

        # --- Restore global state ---
        for k, r in enumerate(rsds):
            r.active = snap_active[k]
            dd, dd_std, r2_0, r2_0_std = snap_local[k]
            r.dd = dd
            r.dd_std = dd_std
            r.R2_0 = r2_0
            r.R2_0_std = r2_0_std
        self.kAB, self.kBA, self.kex = snap_glob
        self.verbose = snap_verbose
        self._cache_arrays()   # rebuild plan for the full active set

        return results

    @staticmethod
    def _aic(chi2, k, n):
        """Akaike Information Criterion for a weighted least-squares fit.

        With residuals already divided by their standard deviations, chi² is
        -2 ln L up to an additive constant, so AIC = chi² + 2k.  The corrected
        AICc adds the small-sample penalty 2k(k+1)/(n-k-1); it falls back to
        +inf when n - k - 1 <= 0.

        Args:
            chi2 (float): sum of squared weighted residuals.
            k    (int)  : number of fitted parameters.
            n    (int)  : number of data points.

        Returns:
            tuple: (AIC, AICc).
        """
        aic = chi2 + 2.0 * k
        denom = n - k - 1
        aicc = aic + (2.0 * k * (k + 1)) / denom if denom > 0 else float('inf')
        return aic, aicc

    def compareModelsAIC(self):
        """Compare the global and individual fits via the AIC.

        Assumes a global fit (fit()) has already been run so the current state
        holds the global solution.  Runs fit_individual() for the alternative
        model, then scores both with AIC / AICc and computes Akaike weights.

        The global solution is left intact on return.

        Returns:
            dict: comparison summary with 'global', 'individual' sub-dicts
                  (each: chi2, k, n, aic, aicc), 'preferred', 'delta_aic',
                  'delta_aicc', Akaike weights, and per-residue individual
                  parameters under 'per_residue'.
        """
        # --- Global model score (current state = global solution) ---
        self._cache_arrays()
        p_glob = self._build_p0_current()
        res_g  = self.errFunc(p_glob)
        chi2_g = float(np.dot(res_g, res_g))
        n      = len(res_g)
        k_g    = len(p_glob)
        aic_g, aicc_g = self._aic(chi2_g, k_g, n)

        # --- Individual model score ---
        indiv = self.fit_individual()
        chi2_i = sum(d['chi2'] for d in indiv)
        k_i    = sum(d['npar'] for d in indiv)
        aic_i, aicc_i = self._aic(chi2_i, k_i, n)

        # --- Per-residue global-vs-individual preference ---
        # Attribute to each residue its own chi2 contribution under the shared
        # global solution, then score global vs individual for that residue.
        # Parameter-count convention: the individual model pays for its extra
        # exchange parameters (nExch), the global model borrows the shared
        # exchange rate for free (local params only).  So the AIC per residue
        # answers "is this residue's own exchange rate worth its cost?".
        nExch      = 1 if self.model == 'Meiboom' else 2
        kex_global = self.kex if self.model == 'Meiboom' else self.kAB + self.kBA
        pr_glob    = self._per_residue_global_chi2()  # {label: (chi2, n)}
        # Leave-one-out kex per residue: the unbiased reference for z(kex).
        # Scoring a residue against the full global kex is circular (that kex
        # was fitted using the residue itself), so use the kex fitted from the
        # other residues only.  Falls back to the full global kex when LOO is
        # undefined (a single active residue).
        loo_kex = self._leave_one_out_kex()  # {label: kex_without_residue} or {}

        for d in indiv:
            chi2_g_r, n_r = pr_glob[d['label']]
            k_i_r = d['npar']
            k_g_r = k_i_r - nExch          # global shares the exchange param(s)
            aic_g_r, aicc_g_r = self._aic(chi2_g_r, k_g_r, n_r)
            aic_i_r, aicc_i_r = self._aic(d['chi2'], k_i_r, n_r)

            # Prefer the corrected AICc; fall back to AIC if AICc is undefined
            # (n - k - 1 <= 0 for a data-starved residue).
            if aicc_g_r == float('inf') or aicc_i_r == float('inf'):
                sel_g, sel_i = aic_g_r, aic_i_r
            else:
                sel_g, sel_i = aicc_g_r, aicc_i_r
            preferred_r = 'global' if sel_g <= sel_i else 'individual'

            # Per-residue Akaike weights (from the selection metric)
            m_min = min(sel_g, sel_i)
            wg_r  = exp(-0.5 * (sel_g - m_min))
            wi_r  = exp(-0.5 * (sel_i - m_min))
            ws    = wg_r + wi_r
            wg_r, wi_r = wg_r / ws, wi_r / ws

            kstd = d.get('kex_std', 0.0)
            d['n']            = n_r
            d['chi2_global']  = chi2_g_r
            d['k_global']     = k_g_r
            d['k_indiv']      = k_i_r
            d['aic_global']   = aic_g_r
            d['aic_indiv']    = aic_i_r
            d['aicc_global']  = aicc_g_r
            d['aicc_indiv']   = aicc_i_r
            d['delta_aicc']   = sel_i - sel_g     # > 0 → global preferred
            d['weight_global'] = wg_r
            d['weight_indiv']  = wi_r
            d['preferred']    = preferred_r
            d['kex_global']   = kex_global
            # z-score of the individual kex against the leave-one-out reference
            # (kex fitted from the other residues); unbiased for outlier testing.
            kex_ref = loo_kex.get(d['label'], kex_global)
            d['kex_ref'] = kex_ref
            d['z_kex'] = (d['kex'] - kex_ref) / kstd if kstd > 0 else float('nan')

        # Akaike weights from the (small-sample) AICc
        aicc_min = min(aicc_g, aicc_i)
        wg = exp(-0.5 * (aicc_g - aicc_min))
        wi = exp(-0.5 * (aicc_i - aicc_min))
        wsum = wg + wi
        wg, wi = wg / wsum, wi / wsum

        preferred = 'global' if aicc_g <= aicc_i else 'individual'

        n_pref_global = sum(1 for d in indiv if d['preferred'] == 'global')
        n_pref_indiv  = len(indiv) - n_pref_global

        # Restore the reported global-fit statistics clobbered by the sub-fits.
        self._refresh_global_stats()

        return {
            'global':     {'chi2': chi2_g, 'k': k_g, 'n': n,
                           'aic': aic_g, 'aicc': aicc_g},
            'individual': {'chi2': chi2_i, 'k': k_i, 'n': n,
                           'aic': aic_i, 'aicc': aicc_i},
            'delta_aic':   aic_i - aic_g,
            'delta_aicc':  aicc_i - aicc_g,
            'weight_global':     wg,
            'weight_individual': wi,
            'preferred':   preferred,
            'n_pref_global':     n_pref_global,
            'n_pref_individual': n_pref_indiv,
            'per_residue': indiv,
        }

    def _refresh_global_stats(self):
        """Recompute self.chi2 / dof / npar / nvar for the full global solution.

        Individual fitting, the per-residue comparison, and the jackknife all
        make intermediate leastsq / errFunc calls that overwrite these scalars
        with sub-fit values.  Call this after such analyses so getLogBuffer
        reports the true full-active global-fit statistics.  Per-residue and
        global parameters are already restored by those methods; this only
        re-evaluates errFunc on the full active set to refresh the scalars.
        """
        self._cache_arrays()
        self.errFunc(self._build_p0_current())

    def _leave_one_out_kex(self):
        """Per-residue leave-one-out global exchange rate.

        Refits the global model with each active residue removed in turn and
        records the shared kex fitted from the *other* residues only.  This is
        the unbiased reference for judging whether a residue disagrees with the
        rest (the full global kex is contaminated by the residue itself), and it
        is exactly the per-residue quantity the jackknife needs.  The result is
        cached per active-residue set so the comparison and the jackknife share
        the N refits.  The full global solution is restored before returning.

        Returns:
            dict: {residue_label: kex_without_that_residue} for each active
            residue.  Empty dict when fewer than 2 residues are active
            (leave-one-out is ill-posed — nothing left to fit).
        """
        rsds = self.dataset.rsds
        active_labels = tuple(r.label for r in rsds if r.active)
        if len(active_labels) < 2:
            return {}

        cache = getattr(self, '_loo_cache', None)
        if cache is not None and cache[0] == active_labels:
            return dict(cache[1])

        # Snapshot the full global solution for warm-start and restore.
        snap_active = [r.active for r in rsds]
        snap_local  = [(r.dd, getattr(r, 'dd_std', 0.0),
                        list(r.R2_0), list(r.R2_0_std)) for r in rsds]
        snap_glob   = (self.kAB, self.kBA, self.kex)
        snap_std    = (self.kAB_std, self.kBA_std, self.kex_std)
        snap_verbose = self.verbose

        def restore_full():
            for k, r in enumerate(rsds):
                r.active = snap_active[k]
                dd, dd_std, r2_0, r2_0_std = snap_local[k]
                r.dd = dd
                r.dd_std = dd_std
                r.R2_0 = list(r2_0)
                r.R2_0_std = list(r2_0_std)
            self.kAB, self.kBA, self.kex = snap_glob
            self.kAB_std, self.kBA_std, self.kex_std = snap_std

        self.verbose = False   # keep the N refits quiet
        out = {}
        for i, r in enumerate(rsds):
            if not snap_active[i]:
                continue
            restore_full()          # each drop warm-starts from the full solution
            r.active = False
            self.fit()              # refit global on the remaining residues
            out[r.label] = self.kex if self.model == 'Meiboom' else self.kAB + self.kBA

        restore_full()
        self.verbose = snap_verbose
        self._refresh_global_stats()   # restore full-active global-fit statistics
        self._loo_cache = (active_labels, dict(out))
        return dict(out)

    def _per_residue_global_chi2(self):
        """Per-residue chi2 contribution under the current global solution.

        Evaluates each active residue in isolation at the shared global
        exchange parameters (and that residue's global-fit dd / R2_0 values),
        giving the residue's own contribution to the global chi2.  No refitting
        is done; the global state is restored before returning.

        Returns:
            dict: {residue_label: (chi2, n_points)} for every active residue.
        """
        rsds = self.dataset.rsds
        snap_active  = [r.active for r in rsds]
        snap_verbose = self.verbose
        self.verbose = False

        out = {}
        for i, r in enumerate(rsds):
            if not snap_active[i]:
                continue
            for k, rr in enumerate(rsds):
                rr.active = (k == i)
            self._cache_arrays()
            # _build_p0_current uses the (unchanged) global exchange rate and
            # this residue's global-fit dd / R2_0, so errFunc returns exactly
            # the residue's contribution to the global chi2.
            p   = self._build_p0_current()
            res = self.errFunc(p)
            out[r.label] = (float(np.dot(res, res)), len(res))

        for k, rr in enumerate(rsds):
            rr.active = snap_active[k]
        self.verbose = snap_verbose
        self._cache_arrays()
        return out

    def jackknifeGlobal(self):
        """Leave-one-residue-out jackknife of the shared global exchange rate.

        Refits the global model with each active residue removed in turn and
        tracks the shared exchange rate kex.  A stable kex across every
        leave-one-out refit means no single residue drives the global fit
        (robust global fit); a large swing when one residue is dropped flags
        that residue as influential / an outlier.

        Assumes fit() has already produced the full global solution.  The full
        global solution is restored before returning.

        Returns:
            dict or None: None if fewer than 2 active residues (jackknife
            undefined).  Otherwise a summary with keys:
              * 'kex_full'      — kex from the full global fit
              * 'kex_jack_mean' — mean of the leave-one-out kex values
              * 'kex_jack_est'  — bias-corrected jackknife estimate
              * 'bias'          — jackknife bias estimate
              * 'se'            — jackknife standard error of kex
              * 'n'             — number of active residues resampled
              * 'per_residue'   — list of {label, kex_drop, delta_kex,
                                  pseudo_value}
        """
        rsds = self.dataset.rsds
        # Ordered active labels; leave-one-out kex comes from the shared helper
        # (also used by compareModelsAIC, so the N refits happen only once).
        active_labels = [r.label for r in rsds if r.active]
        N = len(active_labels)
        if N < 2:
            return None

        kex_full = self.kex if self.model == 'Meiboom' else self.kAB + self.kBA
        loo   = self._leave_one_out_kex()               # {label: kex_without_residue}
        drops = [(label, loo[label]) for label in active_labels]

        # --- Jackknife statistics ---
        kex_vals = [k for _, k in drops]
        mean     = sum(kex_vals) / N
        bias     = (N - 1) * (mean - kex_full)
        # Pseudo-values: N*theta_full - (N-1)*theta_(-i)
        pseudo   = [N * kex_full - (N - 1) * k for k in kex_vals]
        jack_est = sum(pseudo) / N
        se       = sqrt((N - 1) / N * sum((k - mean) ** 2 for k in kex_vals))

        per_residue = []
        for (label, kex_i), ps in zip(drops, pseudo):
            per_residue.append({
                'label':        label,
                'kex_drop':     kex_i,
                'delta_kex':    kex_full - kex_i,
                'pseudo_value': ps,
                'influential':  bool(se > 0 and abs(kex_full - kex_i) > 2.0 * se),
            })

        return {
            'kex_full':      kex_full,
            'kex_jack_mean': mean,
            'kex_jack_est':  jack_est,
            'bias':          bias,
            'se':            se,
            'n':             N,
            'per_residue':   per_residue,
        }

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
        """Derive (kex, kex_std, pB, pB_std) from kAB/kBA and their stds."""
        kex     = self.kAB + self.kBA
        if self.covar is not None:
            # dd entries do not enter kex or pB; dd=0 only zeroes the phi part
            kex_std, pB_std, _ = _derived_std(self.kAB, self.kBA, 0.0, self.covar[:3, :3])
            return kex, kex_std, self.kAB / kex, pB_std
        kex_std = sqrt(self.kAB_std ** 2 + self.kBA_std ** 2)
        pB      = self.kAB / kex
        pB_std  = sqrt(self.kAB_std ** 2 * (1.0 / kex - self.kAB / kex ** 2) ** 2
                       + self.kBA_std ** 2 * (self.kAB / kex ** 2) ** 2)
        return kex, kex_std, pB, pB_std

    def getLogBuffer(self, config_path=None):
        return report.log_buffer(self, config_path)

    def pdf(self, pdfFileName):
        return report.write_pdf(self, pdfFileName)

    # -------------------------------------------------------------------------
    # Model conversion
    # -------------------------------------------------------------------------

    def toMeiboom(self):
        """Convert Matrix/London fit results to Meiboom (fast-exchange) representation.

        After fitting with the Matrix or London model in a fast-exchange regime,
        call this method to replace the per-residue chemical shift difference dd
        (ppm) with the Meiboom effective parameter::

            phi = pB * (1 - pB) * dd^2

        The original dd value is preserved in r.dd_copy so that _calc_r2() and
        pdf() can reconstruct the correct R2 curves even after the conversion.

        Uncertainty propagation::

            phi_std^2 = (dd_std * 2*pB*(1-pB)*dd)^2
                      + (pB_std * (1-2*pB)*dd^2)^2
        """
        # Total exchange rate, its uncertainty, and minor-state population pB
        # (with error propagation) — also stores kex / kex_std on self.
        self.kex, self.kex_std, pB, pB_std = self._kex_pB()

        for r in self.dataset.rsds:
            if r.active:
                # Save original csd before overwriting r.dd with phi
                dd_orig     = r.dd
                dd_std_orig = r.dd_std
                r.dd_copy   = dd_orig        # kept for plot/report reconstruction

                # Replace dd with phi = pB*(1-pB)*dd^2
                r.dd = pB * (1.0 - pB) * dd_orig ** 2
                # Propagate uncertainty through phi = pB*(1-pB)*dd^2
                if self.covar is not None:
                    idx = [0, 1, r._dd_idx]
                    _, _, r.dd_std = _derived_std(self.kAB, self.kBA, dd_orig,
                                                  self.covar[np.ix_(idx, idx)])
                else:
                    r.dd_std = sqrt(
                        (dd_std_orig * 2.0 * pB * (1.0 - pB) * dd_orig) ** 2
                        + (pB_std * (1.0 - 2.0 * pB) * dd_orig ** 2) ** 2
                    )

    # -------------------------------------------------------------------------
    # JSON result export
    # -------------------------------------------------------------------------

    def reportAllValues(self):
        return report.report_all_values(self)
