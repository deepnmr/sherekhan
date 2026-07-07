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

from matplotlib import use
use('Agg')  # non-interactive backend — must be set before importing pyplot

from sys import argv
import numpy as np
from numpy import array, sqrt, pi, dot, real, log, sin, sinh, cos, cosh, matrix, exp, zeros, sort, tanh
from scipy.optimize import leastsq
from matplotlib.pyplot import plot, errorbar, figure, title, legend, xlabel, ylabel, grid, xlim, ylim, close, \
    subplot
from matplotlib.backends.backend_pdf import PdfPages

from .cpmg import CPMGDataSet
from scipy.linalg import expm          # matrix exponential
from numpy.linalg import matrix_power  # integer matrix power
from numpy import diag
from os import getlogin, uname
from time import ctime


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
                            'ver. 1.2.0\n'
                            '\nProf. Dr. Donghan Lee\n'
                            'Korea Basic Science Institute (KBSI)\n'
                            '\nOriginal code: Adam Mazur, Bjoern Hammesfahr, '
                            'Christian Griesinger, Donghan Lee, Martin Kollmar\n'
                            'Max-Planck-Institute for Biophysical Chemistry, '
                            'Goettingen, 2012\n'
                            '\n\nPlease cite: xxx\n')

        # Global exchange-rate parameters shared across all residues
        self.kAB = 0.0      # A→B rate constant (s^-1)
        self.kAB_std = 0.0  # standard deviation of kAB
        self.kBA = 0.0      # B→A rate constant (s^-1)
        self.kBA_std = 0.0  # standard deviation of kBA
        self.kex = 0.0      # total exchange rate kAB + kBA (Meiboom model)
        self.kex_std = 0.0  # standard deviation of kex

        self.dataset = CPMGDataSet()  # container for all experimental data

        self.model = ''      # 'Matrix', 'London', or 'Meiboom'
        self.verbose = False # set True to print chi² at each function evaluation

    # -------------------------------------------------------------------------
    # Parameter initialisation helpers
    # -------------------------------------------------------------------------

    def init_dd(self, dd_AB):
        """Set the chemical shift difference (or phi) for all residues.

        Args:
            dd_AB (float): Chemical shift difference in ppm (Matrix/London)
                           or phi = pB*(1-pB)*dd^2 (Meiboom).
        """
        for r in self.dataset.rsds:
            r.dd = dd_AB
            r.dd_std = 0.0

    def init_R2_0(self):
        """Initialise R2_0 for every residue from the last (highest nu) data point.

        The last CPMG frequency typically gives the smallest Rex contribution,
        so R2exp at that point is a good proxy for the intrinsic R2_0.
        """
        for r in self.dataset.rsds:
            # highest nu_CPMG (last point) → smallest Rex → best R2_0 proxy
            r.R2_0 = [dsp.R2exp[-1] for dsp in r.dspCurves]
            r.R2_0_std = [0.0] * len(r.dspCurves)

    # -------------------------------------------------------------------------
    # Parameter vector packing / unpacking
    # -------------------------------------------------------------------------

    def getLocalParams(self, p):
        """Append per-residue parameters (dd, R2_0 values) to parameter list p.

        The parameter vector layout for active residues is::

            [global params, dd_r0, R2_0_r0_f0, R2_0_r0_f1, ...,
                            dd_r1, R2_0_r1_f0, ...]

        Only active residues contribute to the parameter vector.

        Args:
            p (list): Parameter list that is extended in place.
        """
        for r in self.dataset.rsds:
            if r.active:
                p.append(r.dd)
                for r2_0 in r.R2_0:
                    p.append(r2_0)

    def setLocalParams(self, p, covar=None):
        """Unpack optimised per-residue parameters from the solution vector.

        Updates each active residue's dd and R2_0 attributes.  If a covariance
        matrix is available, standard deviations are also computed.

        Args:
            p     (array-like): Optimised parameter vector from leastsq.
            covar (2-D array or None): Covariance matrix; None when the fit
                                       is singular or poorly constrained.
        """
        # Number of global parameters preceding the local ones
        if self.model in ('Matrix', 'London'):
            nGlob = 2   # kAB, kBA
        elif self.model == 'Meiboom':
            nGlob = 1   # kex

        j = nGlob  # index into the full parameter vector
        for i in range(0, len(self.dataset.rsds)):
            if self.dataset.rsds[i].active:
                self.dataset.rsds[i].dd = p[j]
                if covar is not None:
                    self.dataset.rsds[i].dd_std = sqrt(covar[j][j])
                j += 1
                for ii in range(len(self.dataset.rsds[i].R2_0)):
                    self.dataset.rsds[i].R2_0[ii] = p[j]
                    if covar is not None:
                        self.dataset.rsds[i].R2_0_std[ii] = sqrt(covar[j][j])
                    j += 1

    # -------------------------------------------------------------------------
    # R2 calculation models
    # -------------------------------------------------------------------------

    def R2_clc_matrix(self, kAB, kBA, domega, R2_0, vCPMG, tcp):
        """Compute R2eff using the exact Liouvillian matrix-exponential method.

        Constructs the 2×2 Liouvillian superoperator for two-site chemical
        exchange and propagates the magnetisation vector through one complete
        CPMG echo train of n repeats, each consisting of four delta intervals
        (delta = 1 / (4 * nu_CPMG))::

            B = exp(A*delta) · exp(A*·delta) · exp(A*·delta) · exp(A*delta)
            C = B^n
            R2eff = -1/(4*n*delta) * ln(C[0,0] / M0[0])

        Args:
            kAB    (float): A→B exchange rate (s^-1).
            kBA    (float): B→A exchange rate (s^-1).
            domega (float): Chemical shift difference in rad/s.
            R2_0   (float): Intrinsic transverse relaxation rate (s^-1).
            vCPMG  (float): CPMG pulsing frequency (Hz).
            tcp    (float): Total CPMG pulse delay (s); used to derive n.

        Returns:
            float: Calculated effective R2 relaxation rate (s^-1).
        """
        delta = 1.0 / (4.0 * vCPMG)   # time for one quarter-echo interval
        fn = tcp * vCPMG               # number of CPMG repeats (should be integer)

        n = int(fn)
        if fn - float(n) > 0.001:
            print('Error: n is not an integer')
            exit()

        # Build the Liouvillian matrix A for state A (real) and state B (complex)
        Alist = [[-R2_0 - kAB, kBA],
                 [kAB, -R2_0 - kBA + 1.0j * domega]]

        A = array(Alist)
        eA = expm(A * delta)           # propagator for forward interval

        Astar = A.conj()               # complex conjugate (180° pulse effect)
        eAstar = expm(Astar * delta)   # propagator for refocused interval

        # One full echo: delta → delta* → delta* → delta (CPMG refocusing)
        B = dot(dot(dot(eA, eAstar), eAstar), eA)
        # Repeat n times for the full CPMG train
        C = matrix_power(B, n)

        # Equilibrium populations
        pB = kAB / (kAB + kBA)
        M0 = matrix([1.0 - pB, pB]).T  # initial magnetisation vector

        MT = dot(C, M0)  # magnetisation after the full CPMG sequence

        # Recover R2eff from the decay of the ground-state magnetisation
        R2_clc = (-1.0 / (4.0 * n * delta)) * log(MT[0, 0] / M0[0, 0])

        return real(R2_clc)

    def R2_clc_matrix2(self, kAB, kBA, domega, R2_0, vCPMG, tcp):
        """Compute R2eff via eigenvalue decomposition of the Liouvillian.

        Alternative matrix formulation used inside errFunc during optimisation.
        The propagator is built from the analytic eigenvalues (l1, l2) and
        eigenvectors of A, which can be more numerically stable for certain
        parameter regimes.

        Args:
            kAB    (float): A→B exchange rate (s^-1).
            kBA    (float): B→A exchange rate (s^-1).
            domega (float): Chemical shift difference in rad/s.
            R2_0   (float): Intrinsic transverse relaxation rate (s^-1).
            vCPMG  (float): CPMG pulsing frequency (Hz).
            tcp    (float): Total CPMG pulse delay (s); used to derive n.

        Returns:
            float: Calculated effective R2 relaxation rate (s^-1).
        """
        delta = 1.0 / (4.0 * vCPMG)
        fn = tcp * vCPMG

        n = int(fn)
        if fn - float(n) > 0.001:
            print('Error: n is not an integer')
            exit()

        # Analytic square-root factor appearing in both eigenvalues
        AA = sqrt(kBA ** 2 + (2.0 * kAB - 2.0 * 1.0j * domega) * kBA
                  + kAB ** 2 + 2.0 * 1.0j * domega * kAB - domega ** 2)

        # Eigenvalues of A * delta
        l1 = -(delta * AA + delta * kBA + delta * kAB + (2.0 * R2_0 - 1.0j * domega) * delta) / 2.0
        l2 =  (delta * AA - delta * kBA - delta * kAB + (1.0j * domega - 2.0 * R2_0) * delta) / 2.0

        L = array([l1, l2])

        # Eigenvectors of A and their inverse
        v1 = (2 * kAB) / (-AA + kBA - kAB - 1.0j * domega)
        v2 = (2 * kAB) / ( AA + kBA - kAB - 1.0j * domega)
        X  = array([[1.0, 1.0], [v1, v2]])    # columns are eigenvectors
        Xi = array([
                [ v2 / (v2 - v1), -1.0 / (v2 - v1)],
                [-v1 / (v2 - v1),  1.0 / (v2 - v1)]
               ])                              # inverse of X

        # Complex conjugates for the refocused (Astar) propagator
        Lstar  = L.conj()
        Xstar  = X.conj()
        Xistar = Xi.conj()

        # Forward propagator:  exp(A*delta) = X · diag(exp(L)) · X^-1
        eA = dot(dot(X, diag(exp(L))), Xi)

        # Refocused propagator squared:  exp(Astar*2*delta)
        eAstar2 = dot(dot(Xstar, diag(exp(Lstar) ** 2)), Xistar)

        # One full echo and n repeats
        B = dot(dot(eA, eAstar2), eA)
        C = matrix_power(B, n)

        pB = kAB / (kAB + kBA)
        M0 = matrix([1.0 - pB, pB]).T
        MT = dot(C, M0)

        R2_clc = (-1.0 / (4.0 * n * delta)) * log(MT[0, 0] / M0[0, 0])

        return real(R2_clc)

    def R2_clc_london(self, kAB, kBA, domega, R2_0, vCPMG):
        """Compute R2eff using the London slow-exchange approximation.

        Closed-form analytical expression valid when kex << domega::

            Rex = (kex/2) - (1/tcp) * ln(sqrt(D+*cosh^2(eta+) - D-*cos^2(eta-))
                                         + sqrt(D+*sinh^2(eta+) + D-*sin^2(eta-)))

        where tcp = 1/(2*nu_CPMG).  Adds Rex to the intrinsic R2_0.

        Args:
            kAB    (float): A→B exchange rate (s^-1).
            kBA    (float): B→A exchange rate (s^-1).
            domega (float): Chemical shift difference in rad/s.
            R2_0   (float): Intrinsic transverse relaxation rate (s^-1).
            vCPMG  (float): CPMG pulsing frequency (Hz).

        Returns:
            float: Calculated effective R2 relaxation rate (s^-1).
        """
        # Auxiliary variables following London et al. notation
        zet = 2.0 * domega * (kAB - kBA)
        psi = (kAB - kBA) ** 2.0 - domega ** 2 + 4.0 * kAB * kBA
        ss  = sqrt(psi ** 2.0 + zet ** 2.0)

        tcp = 1.0 / (2.0 * vCPMG)  # effective inter-pulse delay for this nu_CPMG

        # eta± and D± factors of the London formula
        np = tcp / sqrt(8.0) * sqrt( psi + ss)
        nm = tcp / sqrt(8.0) * sqrt(-psi + ss)
        Dp = 1.0 / 2.0 * (1.0 + (psi + 2.0 * domega ** 2) / ss)
        Dm = 1.0 / 2.0 * (-1.0 + (psi + 2.0 * domega ** 2) / ss)

        Rex = ((kAB + kBA) / 2.0
               - (1.0 / tcp) * log(sqrt(Dp * cosh(np) ** 2 - Dm * cos(nm) ** 2)
                                   + sqrt(Dp * sinh(np) ** 2 + Dm * sin(nm) ** 2)))

        R2_clc = R2_0 + Rex
        return R2_clc

    def R2_clc_meiboom(self, kex, phi, R2_0, vCPMG, field):
        """Compute R2eff using the Meiboom fast-exchange approximation.

        Closed-form formula valid when kex >> domega::

            Rex = (phi_omega / kex) * [1 - (4*nu/kex) * tanh(kex / (4*nu))]

        where phi_omega = phi * (2*pi*field)^2 converts phi (in ppm^2) to
        rad^2/s^2 at the given field strength.

        Args:
            kex   (float): Total exchange rate kAB + kBA (s^-1).
            phi   (float): Effective exchange parameter
                           phi = pB*(1-pB)*delta_omega_ppm^2 (ppm^2).
            R2_0  (float): Intrinsic transverse relaxation rate (s^-1).
            vCPMG (float): CPMG pulsing frequency (Hz).
            field (float): Static field strength in MHz (proton Larmor frequency).

        Returns:
            float: Calculated effective R2 relaxation rate (s^-1).
        """
        # Convert phi from ppm^2 to rad^2/s^2 at this magnetic field
        phi_omega = phi * (2.0 * pi * field) ** 2
        Rex = (phi_omega / kex) * (1.0 - (4.0 * vCPMG / kex) * tanh(kex / (4.0 * vCPMG)))
        R2_clc = R2_0 + Rex
        return R2_clc

    def _cache_arrays(self):
        """Pre-cache numpy arrays and a batch plan for all dispersion curves.

        Converts list-based v, R2exp, R2stddev attributes to contiguous numpy
        arrays stored as v_arr, R2exp_arr, R2stddev_arr on each DispersionCurve.

        Also builds the errFunc "batch plan": active dispersion curves that
        share the same CPMG frequency grid and tcp are grouped so that all
        residues in a group can be evaluated in a single vectorised call.
        Per-group frequency-only invariants (delta, integer echo count n, and
        stacked R2exp / R2stddev) are precomputed here so they are not
        recomputed on every errFunc evaluation.

        Called once before fitting or grid-search begins.
        """
        for res in self.dataset.rsds:
            for dsp in res.dspCurves:
                dsp.v_arr        = np.asarray(dsp.v,       dtype=float)
                dsp.R2exp_arr    = np.asarray(dsp.R2exp,   dtype=float)
                dsp.R2stddev_arr = np.asarray(dsp.R2stddev, dtype=float)

        self._build_batch_plan()

    def _build_batch_plan(self):
        """Group active dispersion curves for vectorised errFunc evaluation.

        Builds ``self._batch_groups``: a list of groups, one per distinct
        (tcp, CPMG-frequency-grid) combination among the active residues.
        Each group is a dict holding

            * ``members``   : list of (residue_position, R2_0 field-index)
                              locating where each curve's parameters live in
                              the flat parameter vector unpacked by errFunc.
            * ``fields``    : per-member static field strength (MHz).
            * ``v_arr``     : shared CPMG frequency grid (Hz).
            * ``delta_arr`` : 1/(4 v)         — echo quarter-interval, shape (V,).
            * ``n_arr``     : round(tcp v)     — integer echo counts, shape (V,).
            * ``tcp``       : shared pulse delay (s).
            * ``R2exp``     : stacked experimental R2, shape (C, V).
            * ``R2stddev``  : stacked experimental sigma, shape (C, V).
            * ``slices``    : per-member (start, stop) offsets into the flat
                              residual vector, preserving errFunc ordering.

        The residual layout matches the original residue-major / field-order
        loop exactly, so the residual vector handed to leastsq is unchanged.
        """
        groups = {}          # key -> group dict
        order  = []          # group keys in first-seen order
        idx    = 0           # running offset into the residual vector

        for i_res, res in enumerate(self.dataset.rsds):
            if not res.active:
                continue
            for j, dsp in enumerate(res.dspCurves):
                v_arr = dsp.v_arr
                n     = len(v_arr)
                # Group curves sharing an identical (tcp, frequency grid)
                key = (dsp.tcp, v_arr.tobytes())
                if key not in groups:
                    fn_arr = dsp.tcp * v_arr
                    n_arr  = np.round(fn_arr).astype(int)
                    if np.any(np.abs(fn_arr - n_arr) > 0.001):
                        raise ValueError('n is not an integer for some vCPMG values')
                    groups[key] = {
                        'members':   [],
                        'fields':    [],
                        'v_arr':     v_arr,
                        'delta_arr': 1.0 / (4.0 * v_arr),
                        'n_arr':     n_arr,
                        'tcp':       dsp.tcp,
                        'R2exp':     [],
                        'R2stddev':  [],
                        'slices':    [],
                    }
                    order.append(key)
                g = groups[key]
                g['members'].append((i_res, j))
                g['fields'].append(dsp.field)
                g['R2exp'].append(dsp.R2exp_arr)
                g['R2stddev'].append(dsp.R2stddev_arr)
                g['slices'].append((idx, idx + n))
                idx += n

        # Finalise: convert per-group lists to stacked arrays
        self._batch_groups = []
        for key in order:
            g = groups[key]
            g['fields']   = np.asarray(g['fields'], dtype=float)
            g['R2exp']    = np.vstack(g['R2exp'])
            g['R2stddev'] = np.vstack(g['R2stddev'])
            self._batch_groups.append(g)
        self._n_residuals = idx

    # -------------------------------------------------------------------------
    # Vectorised R2 calculation models (accept v_array, return array)
    # -------------------------------------------------------------------------

    def R2_clc_matrix2_batch(self, kAB, kBA, domega_arr, R2_0_arr,
                             v_array, delta_arr, n_arr, tcp):
        """Vectorised R2eff over many residues *and* all CPMG frequencies.

        Evaluates a batch of curves that share the same CPMG frequency grid
        (and tcp).  Per-residue quantities (domega, R2_0)
        vary along the curve axis C; the eigen-decomposition and propagator are
        evaluated for the full (C, V) grid in one set of array operations.

        Falls back to the per-curve implementation for any curve whose
        eigenvectors are degenerate (kAB ≈ 0 or repeated eigenvalues), so the
        numerical result is identical to the scalar path in every regime.

        Args:
            kAB, kBA   (float)     : exchange rate constants (s^-1).
            domega_arr (np.ndarray): chemical shift differences, shape (C,) rad/s.
            R2_0_arr   (np.ndarray): intrinsic R2 per curve, shape (C,) s^-1.
            v_array    (np.ndarray): shared CPMG frequencies, shape (V,) Hz.
            delta_arr  (np.ndarray): 1/(4 v), shape (V,).
            n_arr      (np.ndarray): integer echo counts round(tcp v), shape (V,).
            tcp        (float)     : CPMG pulse delay (s), for the scalar fallback.

        Returns:
            np.ndarray: real R2eff values, shape (C, V).
        """
        domega = np.asarray(domega_arr, dtype=float)      # (C,)
        R2_0   = np.asarray(R2_0_arr,   dtype=float)       # (C,)
        C = domega.shape[0]

        # --- Per-curve scalar quantities (independent of vCPMG) ---
        AA = np.sqrt(kBA ** 2 + (2.0 * kAB - 2.0j * domega) * kBA
                     + kAB ** 2 + 2.0j * domega * kAB - domega ** 2)   # (C,)

        L1_base = -(AA + kBA + kAB + (2.0 * R2_0 - 1.0j * domega)) / 2.0  # (C,)
        L2_base =  (AA - kBA - kAB + (1.0j * domega - 2.0 * R2_0)) / 2.0  # (C,)

        denom_v1 = -AA + kBA - kAB - 1.0j * domega        # (C,)
        denom_v2 =  AA + kBA - kAB - 1.0j * domega
        v1 = (2.0 * kAB) / denom_v1
        v2 = (2.0 * kAB) / denom_v2
        diff_v = v2 - v1

        # Degenerate eigenvectors: hand those curves to the scalar path.
        bad = ((np.abs(denom_v1) < 1e-10) | (np.abs(denom_v2) < 1e-10)
               | (np.abs(diff_v) < 1e-10))
        if np.any(bad):
            out = np.empty((C, v_array.shape[0]))
            good = ~bad
            if np.any(good):
                out[good] = self.R2_clc_matrix2_batch(
                    kAB, kBA, domega[good], R2_0[good],
                    v_array, delta_arr, n_arr, tcp)
            for ci in np.nonzero(bad)[0]:
                out[ci] = [self.R2_clc_matrix2(kAB, kBA, domega[ci], R2_0[ci],
                                               v, tcp) for v in v_array]
            return out

        # Eigenvector matrices X and their inverse Xi, per curve — shape (C,2,2)
        X = np.empty((C, 2, 2), dtype=complex)
        X[:, 0, 0] = 1.0
        X[:, 0, 1] = 1.0
        X[:, 1, 0] = v1
        X[:, 1, 1] = v2
        Xi = np.empty((C, 2, 2), dtype=complex)
        Xi[:, 0, 0] =  v2 / diff_v
        Xi[:, 0, 1] = -1.0 / diff_v
        Xi[:, 1, 0] = -v1 / diff_v
        Xi[:, 1, 1] =  1.0 / diff_v
        Xstar  = X.conj()
        Xistar = Xi.conj()

        # --- (C, V) diagonal exponentials ---
        # Forward propagator eigenvalues: exp(L_base * delta)
        el1 = np.exp(L1_base[:, None] * delta_arr[None, :])   # (C, V)
        el2 = np.exp(L2_base[:, None] * delta_arr[None, :])
        # Conjugate propagator, doubled interval: exp(L_base* * 2 delta)
        el1s2 = np.exp(L1_base.conj()[:, None] * 2.0 * delta_arr[None, :])
        el2s2 = np.exp(L2_base.conj()[:, None] * 2.0 * delta_arr[None, :])

        # eA[c,v] = X[c] diag(el1,el2) Xi[c]; build by scaling X's columns.
        M_el  = np.stack([el1, el2], axis=2)                  # (C, V, 2)
        eA    = (X[:, None, :, :] * M_el[:, :, None, :]) @ Xi[:, None, :, :]
        M_es2 = np.stack([el1s2, el2s2], axis=2)              # (C, V, 2)
        eAstar2 = (Xstar[:, None, :, :] * M_es2[:, :, None, :]) @ Xistar[:, None, :, :]

        # One full CPMG echo, shape (C, V, 2, 2)
        B = eA @ eAstar2 @ eA

        # B^n via the analytic 2x2 matrix-function formula (diagonalisable 2x2)
        a = B[:, :, 0, 0]; b = B[:, :, 0, 1]
        c = B[:, :, 1, 0]; d = B[:, :, 1, 1]
        half_tr = (a + d) * 0.5
        disc    = np.sqrt(((a - d) * 0.5) ** 2 + b * c)
        lb1 = half_tr + disc
        lb2 = half_tr - disc
        denom = lb1 - lb2
        if np.any(np.abs(denom) < 1e-14):
            # Repeated eigenvalue somewhere — fall back per curve for safety
            return np.array([[self.R2_clc_matrix2(kAB, kBA, domega[ci], R2_0[ci],
                                                  v, tcp) for v in v_array]
                             for ci in range(C)])
        lb1_n = lb1 ** n_arr[None, :]
        lb2_n = lb2 ** n_arr[None, :]
        f1 = (lb1_n - lb2_n) / denom
        f2 = (lb1 * lb2_n - lb2 * lb1_n) / denom
        Cm = f1[:, :, None, None] * B
        Cm[:, :, 0, 0] += f2
        Cm[:, :, 1, 1] += f2

        pB   = kAB / (kAB + kBA)
        M0_0 = 1.0 - pB
        M0   = np.array([M0_0, pB], dtype=complex)
        MT0  = Cm[:, :, 0, :] @ M0                            # (C, V)

        R2 = (-1.0 / (4.0 * n_arr[None, :] * delta_arr[None, :])) * np.log(MT0 / M0_0)
        return np.real(R2)

    def R2_clc_london_vec(self, kAB, kBA, domega, R2_0, v_array):
        """Vectorised R2eff using the London slow-exchange approximation.

        Args:
            kAB     (float)     : A→B exchange rate (s^-1).
            kBA     (float)     : B→A exchange rate (s^-1).
            domega  (float)     : Chemical shift difference in rad/s.
            R2_0    (float)     : Intrinsic transverse relaxation rate (s^-1).
            v_array (np.ndarray): CPMG pulsing frequencies in Hz, shape (n_v,).

        Returns:
            np.ndarray: Real R2eff values, shape (n_v,).
        """
        zet = 2.0 * domega * (kAB - kBA)
        psi = (kAB - kBA) ** 2.0 - domega ** 2 + 4.0 * kAB * kBA
        ss  = sqrt(psi ** 2.0 + zet ** 2.0)

        tcp = 1.0 / (2.0 * v_array)   # effective inter-pulse delay, shape (n_v,)

        # Guard against negative argument under sqrt due to floating-point error
        eta_p = tcp / sqrt(8.0) * np.sqrt(np.maximum(0.0,  psi + ss))
        eta_m = tcp / sqrt(8.0) * np.sqrt(np.maximum(0.0, -psi + ss))

        Dp = 0.5 * (1.0 + (psi + 2.0 * domega ** 2) / ss)
        Dm = 0.5 * (-1.0 + (psi + 2.0 * domega ** 2) / ss)

        Rex = ((kAB + kBA) / 2.0
               - (1.0 / tcp) * np.log(
                   np.sqrt(Dp * np.cosh(eta_p) ** 2 - Dm * np.cos(eta_m) ** 2)
                   + np.sqrt(Dp * np.sinh(eta_p) ** 2 + Dm * np.sin(eta_m) ** 2)))
        return R2_0 + Rex

    def R2_clc_meiboom_vec(self, kex, phi, R2_0, v_array, field):
        """Vectorised R2eff using the Meiboom fast-exchange approximation.

        Args:
            kex     (float)     : Total exchange rate kAB + kBA (s^-1).
            phi     (float)     : phi = pB*(1-pB)*delta_omega_ppm^2 (ppm^2).
            R2_0    (float)     : Intrinsic transverse relaxation rate (s^-1).
            v_array (np.ndarray): CPMG pulsing frequencies in Hz, shape (n_v,).
            field   (float)     : Static field strength in MHz.

        Returns:
            np.ndarray: Real R2eff values, shape (n_v,).
        """
        phi_omega = phi * (2.0 * pi * field) ** 2
        Rex = (phi_omega / kex) * (1.0 - (4.0 * v_array / kex) * np.tanh(kex / (4.0 * v_array)))
        return R2_0 + Rex

    def _calc_r2(self, r, dspCurve, R2_0, v_array):
        """Compute R2 calculated values for a residue at given CPMG frequencies.

        Dispatches to the appropriate model (Matrix / London / Meiboom) and
        returns an array of R2eff values at all requested frequencies.

        For Matrix and London models, uses r.dd_copy (original csd) when
        available, falling back to r.dd.  r.dd_copy is set by toMeiboom() for
        active residues in fast-exchange runs; inactive residues retain r.dd
        as the original csd.

        Args:
            r        (Residue)         : Residue object with dd, R2_0, dd_copy.
            dspCurve (DispersionCurve) : Curve supplying field and tcp values.
            R2_0     (float)           : Intrinsic R2 for this field.
            v_array  (np.ndarray)      : CPMG frequencies at which to evaluate.

        Returns:
            np.ndarray: Calculated R2eff values (same length as v_array).
        """
        ty = zeros(v_array.shape)

        if self.model == 'Matrix':
            # Use the original csd (dd_copy) when toMeiboom() has already
            # replaced r.dd with phi; otherwise fall back to r.dd
            dd = r.dd_copy if r.dd_copy is not None else r.dd
            domega = dd * 2 * pi * dspCurve.field  # rad/s
            for ii in range(len(v_array)):
                ty[ii] = self.R2_clc_matrix(self.kAB, self.kBA, domega, R2_0, v_array[ii], dspCurve.tcp)

        elif self.model == 'London':
            # Same dd_copy guard as Matrix model above
            dd = r.dd_copy if r.dd_copy is not None else r.dd
            domega = dd * 2 * pi * dspCurve.field  # rad/s
            for ii in range(len(v_array)):
                ty[ii] = self.R2_clc_london(self.kAB, self.kBA, domega, R2_0, v_array[ii])

        elif self.model == 'Meiboom':
            # r.dd holds phi = pB*(1-pB)*dd^2 for the Meiboom model
            for ii in range(len(v_array)):
                ty[ii] = self.R2_clc_meiboom(self.kex, r.dd, R2_0, v_array[ii], dspCurve.field)

        return ty

    # -------------------------------------------------------------------------
    # Error function (called by leastsq)
    # -------------------------------------------------------------------------

    def errFunc(self, p):
        """Compute weighted residuals for all active residues and fields.

        This function is passed directly to scipy.optimize.leastsq and must
        return a 1-D array of residuals.  Each residual is::

            (R2exp - R2calc) / R2stddev

        Accumulates chi² and degrees-of-freedom for progress reporting.

        Args:
            p (array-like): Current parameter vector (global + local params).

        Returns:
            np.ndarray: 1-D array of weighted residuals.
        """
        # Lazy-initialise numpy array cache so errFunc works even when called
        # directly without a prior fit() / initGuessAll() call.
        if self.dataset.rsds and not hasattr(self.dataset.rsds[0].dspCurves[0], 'v_arr'):
            self._cache_arrays()

        # --- Unpack global exchange-rate parameters ---
        j = 0
        if self.model == 'Matrix':
            kAB = p[0]
            kBA = p[1]
            j = 2
        elif self.model == 'London':
            kAB = p[0]
            kBA = p[1]
            j = 2
        elif self.model == 'Meiboom':
            kex = p[0]
            j = 1
        else:
            print('Error: wrong model')
            exit()

        nres = len(self.dataset.rsds)  # total number of residues
        dd   = []   # per-residue csd or phi values (length nres)
        R2_0 = []   # per-residue R2_0 lists (length nres × nFields)

        # --- Unpack per-residue local parameters ---
        for i in range(0, nres):
            if self.dataset.rsds[i].active:
                dd.append(p[j])   # dd (csd or phi) for this residue
                j += 1
                R2_0.append([])
                for ii in range(0, len(self.dataset.rsds[i].R2_0)):
                    R2_0[i].append(p[j])   # R2_0 at each field
                    j += 1
            else:
                # Inactive residue: placeholder values (not used in residuals)
                dd.append(0.0)
                R2_0.append([])

        # --- Compute weighted residuals ---
        if self.model == 'Matrix':
            # Batched path: evaluate all residues sharing a CPMG grid at once.
            if not hasattr(self, '_batch_groups'):
                self._build_batch_plan()
            R2_residuals = np.empty(self._n_residuals)
            two_pi = 2.0 * pi
            for g in self._batch_groups:
                members = g['members']
                # Per-curve chemical shift (rad/s) and intrinsic R2 for this group
                domega_arr = np.fromiter(
                    (dd[i_res] * two_pi * fld
                     for (i_res, _), fld in zip(members, g['fields'])),
                    dtype=float, count=len(members))
                R2_0_arr = np.fromiter(
                    (R2_0[i_res][j] for (i_res, j) in members),
                    dtype=float, count=len(members))

                R2_clc = self.R2_clc_matrix2_batch(
                    kAB, kBA, domega_arr, R2_0_arr,
                    g['v_arr'], g['delta_arr'], g['n_arr'], g['tcp'])   # (C, V)

                res_block = (g['R2exp'] - R2_clc) / g['R2stddev']  # (C, V)
                for row, (start, stop) in enumerate(g['slices']):
                    R2_residuals[start:stop] = res_block[row]
        else:
            # Per-curve path for the analytic London / Meiboom models
            n_total = sum(len(res.dspCurves[j].R2exp)
                          for res in self.dataset.rsds if res.active
                          for j in range(len(res.dspCurves)))
            R2_residuals = np.empty(n_total)
            idx = 0
            for i in range(0, nres):
                res = self.dataset.rsds[i]
                if res.active:
                    for j in range(0, len(res.dspCurves)):
                        dsp = res.dspCurves[j]
                        n   = len(dsp.v_arr)

                        if self.model == 'London':
                            domega = dd[i] * 2 * pi * dsp.field
                            R2_clc = self.R2_clc_london_vec(
                                kAB, kBA, domega, R2_0[i][j], dsp.v_arr)
                        elif self.model == 'Meiboom':
                            R2_clc = self.R2_clc_meiboom_vec(
                                kex, dd[i], R2_0[i][j], dsp.v_arr, dsp.field)

                        R2_residuals[idx:idx + n] = (dsp.R2exp_arr - R2_clc) / dsp.R2stddev_arr
                        idx += n

        # --- Chi² and degrees of freedom ---
        self.chi2 = float(np.dot(R2_residuals, R2_residuals))

        self.dof  = len(R2_residuals) - len(p)  # degrees of freedom
        self.npar = len(p)
        self.nvar = len(R2_residuals)

        if self.verbose:
            # Guard against dof ≤ 0 (can happen with too few data points)
            red_chi2 = self.chi2 / self.dof if self.dof > 0 else float('nan')
            if self.model in ('Matrix', 'London'):
                print('kAB=%8.3f kBA=%8.3f  red_chi2=%12.3f  dof=%d  nv=%d  np=%d'
                      % (kAB, kBA, red_chi2, self.dof, len(R2_residuals), len(p)))
            elif self.model == 'Meiboom':
                print('kex=%8.3f            red_chi2=%12.3f  dof=%d  nv=%d  np=%d'
                      % (kex, red_chi2, self.dof, len(R2_residuals), len(p)))

        return R2_residuals

    # -------------------------------------------------------------------------
    # Fitting routines
    # -------------------------------------------------------------------------

    def fit(self):
        """Run nonlinear least-squares fitting for all active residues.

        Builds the initial parameter vector, calls scipy.optimize.leastsq,
        and updates the global exchange-rate parameters (kAB / kBA / kex) and
        per-residue local parameters (dd, R2_0) in place.

        If the covariance matrix is None (singular fit), standard deviations
        are set to zero and a warning is printed.
        """
        # Cache numpy arrays so errFunc can use pre-allocated v_arr, R2exp_arr, R2stddev_arr
        self._cache_arrays()

        # --- Build initial parameter vector p0 ---
        p0 = []

        if self.model in ('Matrix', 'London'):
            p0.append(self.kAB)  # global: kAB
            p0.append(self.kBA)  # global: kBA
        elif self.model == 'Meiboom':
            p0.append(self.kex)  # global: kex

        self.getLocalParams(p0)  # append per-residue dd and R2_0 values

        # --- Perform optimisation ---
        print(p0)
        out = leastsq(self.errFunc, x0=p0, full_output=1)

        p1    = out[0]  # solution vector
        covar = out[1]  # covariance matrix (None when fit is singular)

        if covar is None:
            print('Warning: covariance matrix is None (fit may be singular). '
                  'Standard deviations set to 0.')

        # --- Update global parameters from solution ---
        if self.model in ('Matrix', 'London'):
            self.kAB     = p1[0]
            self.kBA     = p1[1]
            self.kAB_std = sqrt(covar[0][0]) if covar is not None else 0.0
            self.kBA_std = sqrt(covar[1][1]) if covar is not None else 0.0
        elif self.model == 'Meiboom':
            self.kex     = p1[0]
            self.kex_std = sqrt(covar[0][0]) if covar is not None else 0.0

        self.setLocalParams(p1, covar)  # unpack per-residue parameters
        print()

    # -------------------------------------------------------------------------
    # Individual fitting and AIC-based model comparison
    # -------------------------------------------------------------------------

    def _build_p0_current(self):
        """Build the parameter vector from the current global + local state.

        Layout matches errFunc: [global params, dd, R2_0..., dd, R2_0..., ...]
        for every currently-active residue.

        Returns:
            list: initial-parameter vector p0.
        """
        p0 = []
        if self.model == 'Meiboom':
            p0.append(self.kex)
        else:
            p0.append(self.kAB)
            p0.append(self.kBA)
        self.getLocalParams(p0)
        return p0

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

            p0  = self._build_p0_current()
            out = leastsq(self.errFunc, x0=p0, full_output=1)
            p1  = out[0]
            fvec = out[2]['fvec']              # residuals at the solution
            chi2 = float(np.dot(fvec, fvec))
            npar = len(p1)
            nvar = len(fvec)

            # Decode fitted parameters for this residue
            r = rsds[i]
            rec = {'label': r.label, 'chi2': chi2, 'nvar': nvar, 'npar': npar}
            if self.model == 'Meiboom':
                kex = p1[0]
                rec['kex'] = kex
                rec['phi'] = p1[1]
                rec['R2_0'] = list(p1[2:])
            else:
                kAB, kBA = p1[0], p1[1]
                kex = kAB + kBA
                rec['kex'] = kex
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

        # Akaike weights from the (small-sample) AICc
        aicc_min = min(aicc_g, aicc_i)
        wg = exp(-0.5 * (aicc_g - aicc_min))
        wi = exp(-0.5 * (aicc_i - aicc_min))
        wsum = wg + wi
        wg, wi = wg / wsum, wi / wsum

        preferred = 'global' if aicc_g <= aicc_i else 'individual'

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
            'per_residue': indiv,
        }

    def aicReport(self, cmp):
        """Format the AIC model-comparison result as a text block.

        Args:
            cmp (dict): the dict returned by compareModelsAIC().

        Returns:
            str: multi-line report suitable for the log file / stdout.
        """
        g = cmp['global']
        i = cmp['individual']
        buf  = '\n========\n'
        buf += 'Model comparison: global vs individual fit (AIC)\n'
        buf += '========\n'
        buf += 'Model         %10s %6s %6s %12s %12s\n' % (
            'chi2', 'k', 'n', 'AIC', 'AICc')
        buf += 'global      %12.3f %6d %6d %12.3f %12.3f\n' % (
            g['chi2'], g['k'], g['n'], g['aic'], g['aicc'])
        buf += 'individual  %12.3f %6d %6d %12.3f %12.3f\n' % (
            i['chi2'], i['k'], i['n'], i['aic'], i['aicc'])
        buf += '---\n'
        buf += 'delta AIC  (individual - global): %12.3f\n' % cmp['delta_aic']
        buf += 'delta AICc (individual - global): %12.3f\n' % cmp['delta_aicc']
        buf += 'Akaike weights (AICc): global=%.4f  individual=%.4f\n' % (
            cmp['weight_global'], cmp['weight_individual'])
        buf += 'Preferred model: %s\n' % cmp['preferred']
        buf += '---\n'
        buf += 'Per-residue individual fit:\n'
        if self.model == 'Meiboom':
            buf += '%6s %10s %12s %10s\n' % ('resId', 'kex', 'phi', 'chi2')
            for d in cmp['per_residue']:
                buf += '%6s %10.3f %12.5f %10.3f\n' % (
                    d['label'], d['kex'], d['phi'], d['chi2'])
        else:
            buf += '%6s %10s %8s %10s %10s\n' % (
                'resId', 'kex', 'pB', 'csd', 'chi2')
            for d in cmp['per_residue']:
                buf += '%6s %10.3f %8.4f %10.4f %10.3f\n' % (
                    d['label'], d['kex'], d['pB'], d['dd'], d['chi2'])
        buf += '========\n'
        return buf

    def _build_grid_p0(self, kex_val, pB_val, csd_val):
        """Build a parameter vector for a grid search point without modifying state.

        Constructs the full parameter vector [global_params, dd_r0, R2_0_r0_f0, ...]
        directly from the given grid values, bypassing init_dd() / getLocalParams()
        so that multiple threads can call this concurrently.

        Args:
            kex_val (float): Total exchange rate for this grid point (s^-1).
            pB_val  (float): Minor-state population for this grid point.
            csd_val (float): Chemical shift difference (ppm) for this grid point.

        Returns:
            np.ndarray: Parameter vector ready to pass to errFunc.
        """
        p0 = []
        if self.model in ('Matrix', 'London'):
            kAB = pB_val * kex_val
            kBA = kex_val - kAB
            p0.append(kAB)
            p0.append(kBA)
            dd = csd_val
        else:  # Meiboom
            p0.append(kex_val)
            dd = pB_val * (1.0 - pB_val) * csd_val ** 2

        for r in self.dataset.rsds:
            if r.active:
                p0.append(dd)
                for r2_0 in r.R2_0:
                    p0.append(r2_0)
        return np.array(p0)

    def initGuessAll(self):
        """Grid-search for the best initial parameter values before fitting.

        Evaluates chi² on a coarse grid of (kex, pB, csd) combinations and
        selects the combination that gives the lowest chi².  The grid covers::

            kex  : [100, 500, 1000, 2000, 4000, 10000] s^-1
            pB   : [0.05, 0.1, 0.2, 0.4]
            csd  : [0.1, 0.5, 1.0, 2.0, 4.0, 10.0] ppm

        Sets self.kAB / self.kBA (or self.kex) and per-residue dd to the
        best-guess values, and initialises R2_0 from the data.
        """
        # Coarse parameter grids for the global search
        kex = [100.0, 500.0, 1000.0, 2000.0, 4000.0, 10000.0]
        pB  = [0.05, 0.1, 0.2, 0.4]
        csd = [0.1, 0.5, 1.0, 2.0, 4.0, 10.0]

        # Cache numpy arrays once before the grid-search loop
        self._cache_arrays()

        # Initialise R2_0 from the last (highest-nu) experimental data point
        self.init_R2_0()

        # --- Exhaustive grid search ---
        # Uses _build_grid_p0 to avoid mutating instance state per grid point.
        minChi2 = 9.0e99
        ii_kex = ii_pB = ii_csd = 0

        for i_kex, kv in enumerate(kex):
            for i_pB, pv in enumerate(pB):
                for i_csd, cv in enumerate(csd):
                    p0 = self._build_grid_p0(kv, pv, cv)
                    residuals = self.errFunc(p0)
                    chi2 = float(np.dot(residuals, residuals))
                    if chi2 < minChi2:
                        minChi2 = chi2
                        ii_kex, ii_pB, ii_csd = i_kex, i_pB, i_csd

        # --- Apply the best-guess values ---
        print('Initial values:')
        if self.model in ('Matrix', 'London'):
            self.kAB = pB[ii_pB] * kex[ii_kex]
            self.kBA = kex[ii_kex] - self.kAB
            self.init_dd(csd[ii_csd])
            print('kAB: %8.3f kBA: %8.3f CSD: %8.3f' % (self.kAB, self.kBA, csd[ii_csd]))
        elif self.model == 'Meiboom':
            self.kex = kex[ii_kex]
            self.init_dd(pB[ii_pB] * (1.0 - pB[ii_pB]) * csd[ii_csd] ** 2)
            print('kex: %8.3f phi: %8.3f'
                  % (self.kex, pB[ii_pB] * (1.0 - pB[ii_pB]) * csd[ii_csd] ** 2))

    # -------------------------------------------------------------------------
    # Reporting
    # -------------------------------------------------------------------------

    def _kex_pB(self):
        """Derive (kex, kex_std, pB, pB_std) from kAB/kBA and their stds."""
        kex     = self.kAB + self.kBA
        kex_std = sqrt(self.kAB_std ** 2 + self.kBA_std ** 2)
        pB      = self.kAB / kex
        pB_std  = sqrt(self.kAB_std ** 2 * (1.0 / kex - self.kAB / kex ** 2) ** 2
                       + self.kBA_std ** 2 * (self.kAB / kex ** 2) ** 2)
        return kex, kex_std, pB, pB_std

    def getLogBuffer(self, config_path=None):
        """Build a formatted text report of the fitting results.

        Includes the program header, run metadata (user, hostname, timestamp),
        the contents of the input config file, all datasets with experimental
        and calculated R2 values, fit statistics, and per-residue parameters
        with standard deviations.

        Returns:
            str: Multi-line text suitable for writing directly to a .log file.
        """
        buf = '********\n'
        buf += self.programName
        buf += '********\n\n'

        # Run metadata: user, hostname, and current time
        buf += '********\n'
        user     = getlogin()
        hostname = uname()[1]
        buf += 'User: %s@%s\n' % (user, hostname)
        buf += '%s\n' % ctime()
        buf += '********\n\n'

        # Echo the input config file for full reproducibility
        buf += '********\n\n'
        buf += 'Input file:\n'
        if config_path is None:
            config_path = argv[1]  # fallback to the CLI arg when not supplied
        with open(config_path, 'r') as file1:
            buf += file1.read()
        buf += '********\n\n'

        # Dataset summary: for each field, list experimental and calculated R2
        buf += 'Dataset(s):\n'
        for fi in range(0, len(self.dataset.fields)):
            buf += 'Field, tcp [MHz, s]:\n'
            buf += '%5.2f\n' % self.dataset.fields[fi]
            buf += '%f\n'    % self.dataset.tcps[fi]
            buf += '#v_cpmg(Hz)   R2(1/s)      Esd(R2)      R2_calc\n'

            for r in self.dataset.rsds:
                if r.active:
                    buf += '# %s\n' % r.label
                    for i in range(0, len(r.dspCurves)):
                        dspCurve = r.dspCurves[i]
                        R2_0 = r.R2_0[i]

                        va     = array(dspCurve.v)
                        R2_calc = self._calc_r2(r, dspCurve, R2_0, va)

                        # Only print the row for the field being reported
                        if dspCurve.field == self.dataset.fields[fi]:
                            for j in range(0, len(dspCurve.v)):
                                buf += '%8.3f %12.3f %12.3f %12.3f\n' % (
                                    dspCurve.v[j], dspCurve.R2exp[j],
                                    dspCurve.R2stddev[j], R2_calc[j])

        # Fit statistics
        buf += '\n\n\n---\n'
        buf += 'Exchange regime: %s\n' % self.dataset.exchange
        buf += 'Model: %s\n' % self.model
        buf += '---\n'
        # Guard against dof ≤ 0 to prevent division by zero
        red_chi2 = self.chi2 / self.dof if self.dof > 0 else float('nan')
        buf += 'npar=%d nvar=%d ndof=%d chi2=%8.3f chi2/dof=%8.3f\n' % (
            self.npar, self.nvar, self.dof, self.chi2, red_chi2)
        buf += '---\n'

        # Per-exchange-regime parameter summary
        if self.dataset.exchange == 'slow':
            # Derive kex and pB from the individual rate constants
            kex, kex_std, pB, pB_std = self._kex_pB()

            buf += 'kAB: %8.3f  +-  %8.3f\n' % (self.kAB, self.kAB_std)
            buf += 'kBA: %8.3f  +-  %8.3f\n' % (self.kBA, self.kBA_std)
            buf += '---\n'
            buf += 'kex: %8.3f  +-  %8.3f\n' % (kex, kex_std)
            buf += ' pB: %8.3f  +-  %8.3f\n' % (pB, pB_std)
            buf += '---\n'
            buf += 'resId:         csd [ppm]   '
            for field in self.dataset.fields:
                buf += '       R2_0 (%5.1f MHz) ' % field
            buf += '\n---\n'

        if self.dataset.exchange == 'fast':
            buf += 'kex: %8.3f  +-  %8.3f\n' % (self.kex, self.kex_std)
            buf += '---\n'
            buf += 'resId:         phi [ppm^2] '
            for field in self.dataset.fields:
                buf += '       R2_0 (%5.1f MHz) ' % field
            buf += '\n---\n'

        # Per-residue parameters (dd/phi and R2_0 at each field)
        for r in self.dataset.rsds:
            if r.active:
                buf += '%5s: %8.3f +- %8.3f' % (r.label, r.dd, r.dd_std)
                for cf in self.dataset.fields:
                    flag = False
                    for i in range(0, len(r.dspCurves)):
                        if r.dspCurves[i].field == cf:
                            buf += '    %8.3f' % r.R2_0[i]
                            buf += ' +- %8.3f' % r.R2_0_std[i]
                            flag = True
                    if not flag:
                        buf += '    --------'
                        buf += ' +- --------'
                buf += '\n'
            else:
                buf += '%5s:    excluded from the fit\n' % r.label
        buf += '---\nEnd.\n'

        return buf

    def pdf(self, pdfFileName):
        """Generate a multi-page PDF with one dispersion-curve plot per residue.

        Each page shows experimental R2 values (with error bars) and the
        corresponding calculated curve (solid line) for each magnetic field,
        colour-coded consistently.

        Args:
            pdfFileName (str): Path of the output PDF file to create.
        """
        pdf = PdfPages(pdfFileName)

        for r in self.dataset.rsds:
            if r.active:
                figure(figsize=(8, 6))

                # Fixed colour set, indexed by curve position
                colorsSet = ['b', 'g', 'r', 'c', 'm', 'y']

                max_x = []
                min_y = []
                max_y = []

                for i in range(0, len(r.dspCurves)):
                    dspCurve = r.dspCurves[i]
                    max_x.append(max(dspCurve.v))

                    # Assign colour by position (fall back to black when exhausted)
                    c = colorsSet[i] if i < len(colorsSet) else 'k'

                    # Plot experimental data as dots with error bars
                    errorbar(dspCurve.v, dspCurve.R2exp,
                             yerr=dspCurve.R2stddev,
                             fmt='%so' % c,
                             label=' exp %5.1f' % dspCurve.field)

                    # Evaluate calculated curve on a sorted frequency array
                    tx = sort(array(dspCurve.v))
                    R2_0 = r.R2_0[i]
                    ty = self._calc_r2(r, dspCurve, R2_0, tx)
                    plot(tx, ty, '%s-' % c, label='calc %5.1f' % dspCurve.field)

                    # Collect y-range for axis scaling
                    min_yv = [v - e for v, e in zip(dspCurve.R2exp, dspCurve.R2stddev)]
                    max_yv = [v + e for v, e in zip(dspCurve.R2exp, dspCurve.R2stddev)]
                    min_y.append(min(min_yv))
                    max_y.append(max(max_yv))

                title(r.label)
                xlabel('v [Hz]')
                ylabel('R2 [1/s]')

                # Shrink plot area by 20% on the right to make room for the legend
                box = subplot(111).get_position()
                subplot(111).set_position([box.x0, box.y0, box.width * 0.8, box.height])
                legend(loc='center left', bbox_to_anchor=(1, 0.5))

                grid(True)
                ylim(0.95 * min(min_y), 1.05 * max(max_y))
                xlim(0.0, 1.05 * max(max_x))
                pdf.savefig()
                close()

        pdf.close()

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
                r.dd_std = sqrt(
                    (dd_std_orig * 2.0 * pB * (1.0 - pB) * dd_orig) ** 2
                    + (pB_std * (1.0 - 2.0 * pB) * dd_orig ** 2) ** 2
                )

    # -------------------------------------------------------------------------
    # JSON result export
    # -------------------------------------------------------------------------

    def reportAllValues(self):
        """Export all per-residue fit results as a nested list of dictionaries.

        For each residue and each dispersion curve (one per field), builds a
        dictionary containing the calculated R2 curve, exchange parameters,
        and original experimental data.  The structure is::

            [
              [  # residue 0
                { field 0 data },
                { field 1 data },
                ...
              ],
              ...  # residue 1, 2, ...
            ]

        Returns:
            list: Nested list suitable for JSON serialisation via json.dumps().
        """
        all_data = []

        for r in self.dataset.rsds:
            all_data.append([])
            for i, dspCurve in enumerate(r.dspCurves):
                tx   = array(dspCurve.v)
                R2_0 = r.R2_0[i]
                ty   = self._calc_r2(r, dspCurve, R2_0, tx)  # calculated R2

                if self.dataset.exchange == 'slow':
                    kAB        = self.kAB
                    kBA        = self.kBA
                    kex        = kAB + kBA
                    pB         = kAB / kex
                    chemShift  = r.dd       # chemical shift difference (ppm)
                    chemShift_std = r.dd_std

                elif self.dataset.exchange == 'fast':
                    kex     = self.kex
                    phi     = r.dd          # effective Meiboom phi (ppm^2)
                    phi_std = r.dd_std

                # Append the per-field result dictionary
                if self.dataset.exchange == 'slow':
                    all_data[-1].append({
                        'tx':                 list(tx),
                        'r2_clc':             list(ty),
                        'kAB':                kAB,
                        'kBA':                kBA,
                        'kex':                kex,
                        'pB':                 pB,
                        'chemShiftDiff':      chemShift,
                        'chemShiftDiffStdDev': chemShift_std,
                        'r20':                R2_0,
                        'x1':                 dspCurve.v,
                        'y_val':              dspCurve.R2exp,
                        'y_val_err':          dspCurve.R2stddev,
                        'field':              dspCurve.field,
                        'label':              r.label
                    })
                elif self.dataset.exchange == 'fast':
                    all_data[-1].append({
                        'tx':        list(tx),
                        'r2_clc':    list(ty),
                        'kex':       kex,
                        'phi':       phi,
                        'phiStdDev': phi_std,
                        'r20':       R2_0,
                        'x1':        dspCurve.v,
                        'y_val':     dspCurve.R2exp,
                        'y_val_err': dspCurve.R2stddev,
                        'field':     dspCurve.field,
                        'label':     r.label
                    })

        return all_data
