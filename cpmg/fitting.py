"""
fitting.py - Nonlinear least-squares fitting engine for the two-state model.

Functions here take the CPMG_model instance as their first argument ``m``
and read and write its state in place, exactly as the former methods did.
The CPMG_model facade in model_2state.py delegates to these functions.

Attributes an ``m`` must expose: kAB, kBA, kex, kAB_std, kBA_std, kex_std,
dataset, model, verbose, chi2, dof, npar, nvar, covar, plus the R2 facade
methods R2_clc_matrix2_batch, R2_clc_london_vec and R2_clc_meiboom_vec.
"""

import sys
import numpy as np
from numpy import array, sqrt, pi, dot, real, log, sin, sinh, cos, cosh, exp, zeros, sort, tanh, isfinite
from scipy.optimize import leastsq

from .models import _echo_count


def _warn_if_not_converged(out):
    """Print a stderr warning when leastsq did not report convergence."""
    if out[4] not in (1, 2, 3, 4):
        msg = ' '.join(str(out[3]).split())
        print('Warning: optimizer did not converge (ier=%d): %s' % (out[4], msg),
              file=sys.stderr)


def init_dd(m, dd_AB):
    """Set the chemical shift difference (or phi) for all residues.

    Args:
        dd_AB (float): Chemical shift difference in ppm (Matrix/London)
                       or phi = pB*(1-pB)*dd^2 (Meiboom).
    """
    for r in m.dataset.rsds:
        r.dd = dd_AB
        r.dd_std = 0.0


def init_R2_0(m):
    """Initialise R2_0 for every residue from the last (highest nu) data point.

    The last CPMG frequency typically gives the smallest Rex contribution,
    so R2exp at that point is a good proxy for the intrinsic R2_0.
    """
    for r in m.dataset.rsds:
        # highest nu_CPMG (last point) → smallest Rex → best R2_0 proxy
        r.R2_0 = [dsp.R2exp[-1] for dsp in r.dspCurves]
        r.R2_0_std = [0.0] * len(r.dspCurves)


def get_local_params(m, p):
    """Append per-residue parameters (dd, R2_0 values) to parameter list p.

    The parameter vector layout for active residues is::

        [global params, dd_r0, R2_0_r0_f0, R2_0_r0_f1, ...,
                        dd_r1, R2_0_r1_f0, ...]

    Only active residues contribute to the parameter vector.

    Args:
        p (list): Parameter list that is extended in place.
    """
    for r in m.dataset.rsds:
        if r.active:
            p.append(r.dd)
            for r2_0 in r.R2_0:
                p.append(r2_0)


def set_local_params(m, p, covar=None):
    """Unpack optimised per-residue parameters from the solution vector.

    Updates each active residue's dd and R2_0 attributes.  If a covariance
    matrix is available, standard deviations are also computed.

    Args:
        p     (array-like): Optimised parameter vector from leastsq.
        covar (2-D array or None): Covariance matrix; None when the fit
                                   is singular or poorly constrained.
    """
    # Number of global parameters preceding the local ones
    if m.model in ('Matrix', 'London'):
        nGlob = 2   # kAB, kBA
    elif m.model == 'Meiboom':
        nGlob = 1   # kex

    j = nGlob  # index into the full parameter vector
    for i in range(0, len(m.dataset.rsds)):
        if m.dataset.rsds[i].active:
            m.dataset.rsds[i].dd = p[j]
            m.dataset.rsds[i]._dd_idx = j
            if covar is not None:
                m.dataset.rsds[i].dd_std = sqrt(covar[j][j])
            j += 1
            for ii in range(len(m.dataset.rsds[i].R2_0)):
                m.dataset.rsds[i].R2_0[ii] = p[j]
                if covar is not None:
                    m.dataset.rsds[i].R2_0_std[ii] = sqrt(covar[j][j])
                j += 1


def cache_arrays(m):
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
    for res in m.dataset.rsds:
        for dsp in res.dspCurves:
            dsp.v_arr        = np.asarray(dsp.v,       dtype=float)
            dsp.R2exp_arr    = np.asarray(dsp.R2exp,   dtype=float)
            dsp.R2stddev_arr = np.asarray(dsp.R2stddev, dtype=float)

    build_batch_plan(m)


def build_batch_plan(m):
    """Group active dispersion curves for vectorised errFunc evaluation.

    Builds ``m._batch_groups``: a list of groups, one per distinct
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

    for i_res, res in enumerate(m.dataset.rsds):
        if not res.active:
            continue
        for j, dsp in enumerate(res.dspCurves):
            v_arr = dsp.v_arr
            n     = len(v_arr)
            # Group curves sharing an identical (tcp, frequency grid)
            key = (dsp.tcp, v_arr.tobytes())
            if key not in groups:
                n_arr  = np.array([_echo_count(dsp.tcp, v) for v in v_arr], dtype=int)
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
    m._batch_groups = []
    for key in order:
        g = groups[key]
        g['fields']   = np.asarray(g['fields'], dtype=float)
        g['R2exp']    = np.vstack(g['R2exp'])
        g['R2stddev'] = np.vstack(g['R2stddev'])
        m._batch_groups.append(g)
    m._n_residuals = idx


def err_func(m, p):
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
    if m.dataset.rsds and not hasattr(m.dataset.rsds[0].dspCurves[0], 'v_arr'):
        cache_arrays(m)

    # --- Unpack global exchange-rate parameters ---
    j = 0
    if m.model == 'Matrix':
        kAB = p[0]
        kBA = p[1]
        j = 2
    elif m.model == 'London':
        kAB = p[0]
        kBA = p[1]
        j = 2
    elif m.model == 'Meiboom':
        kex = p[0]
        j = 1
    else:
        raise ValueError('Error: wrong model %r (expected Matrix, London or Meiboom)'
                         % m.model)

    nres = len(m.dataset.rsds)  # total number of residues
    dd   = []   # per-residue csd or phi values (length nres)
    R2_0 = []   # per-residue R2_0 lists (length nres × nFields)

    # --- Unpack per-residue local parameters ---
    for i in range(0, nres):
        if m.dataset.rsds[i].active:
            dd.append(p[j])   # dd (csd or phi) for this residue
            j += 1
            R2_0.append([])
            for ii in range(0, len(m.dataset.rsds[i].R2_0)):
                R2_0[i].append(p[j])   # R2_0 at each field
                j += 1
        else:
            # Inactive residue: placeholder values (not used in residuals)
            dd.append(0.0)
            R2_0.append([])

    # --- Compute weighted residuals ---
    if m.model == 'Matrix':
        # Batched path: evaluate all residues sharing a CPMG grid at once.
        if not hasattr(m, '_batch_groups'):
            build_batch_plan(m)
        R2_residuals = np.empty(m._n_residuals)
        two_pi = 2.0 * pi
        for g in m._batch_groups:
            members = g['members']
            # Per-curve chemical shift (rad/s) and intrinsic R2 for this group
            domega_arr = np.fromiter(
                (dd[i_res] * two_pi * fld
                 for (i_res, _), fld in zip(members, g['fields'])),
                dtype=float, count=len(members))
            R2_0_arr = np.fromiter(
                (R2_0[i_res][j] for (i_res, j) in members),
                dtype=float, count=len(members))

            R2_clc = m.R2_clc_matrix2_batch(
                kAB, kBA, domega_arr, R2_0_arr,
                g['v_arr'], g['delta_arr'], g['n_arr'], g['tcp'])   # (C, V)

            res_block = (g['R2exp'] - R2_clc) / g['R2stddev']  # (C, V)
            for row, (start, stop) in enumerate(g['slices']):
                R2_residuals[start:stop] = res_block[row]
    else:
        # Per-curve path for the analytic London / Meiboom models
        n_total = sum(len(res.dspCurves[j].R2exp)
                      for res in m.dataset.rsds if res.active
                      for j in range(len(res.dspCurves)))
        R2_residuals = np.empty(n_total)
        idx = 0
        for i in range(0, nres):
            res = m.dataset.rsds[i]
            if res.active:
                for j in range(0, len(res.dspCurves)):
                    dsp = res.dspCurves[j]
                    n   = len(dsp.v_arr)

                    if m.model == 'London':
                        domega = dd[i] * 2 * pi * dsp.field
                        R2_clc = m.R2_clc_london_vec(
                            kAB, kBA, domega, R2_0[i][j], dsp.v_arr)
                    elif m.model == 'Meiboom':
                        R2_clc = m.R2_clc_meiboom_vec(
                            kex, dd[i], R2_0[i][j], dsp.v_arr, dsp.field)

                    R2_residuals[idx:idx + n] = (dsp.R2exp_arr - R2_clc) / dsp.R2stddev_arr
                    idx += n

    # --- Chi² and degrees of freedom ---
    m.chi2 = float(np.dot(R2_residuals, R2_residuals))

    m.dof  = len(R2_residuals) - len(p)  # degrees of freedom
    m.npar = len(p)
    m.nvar = len(R2_residuals)

    if m.verbose:
        # Guard against dof ≤ 0 (can happen with too few data points)
        red_chi2 = m.chi2 / m.dof if m.dof > 0 else float('nan')
        if m.model in ('Matrix', 'London'):
            print('kAB=%8.3f kBA=%8.3f  red_chi2=%12.3f  dof=%d  nv=%d  np=%d'
                  % (kAB, kBA, red_chi2, m.dof, len(R2_residuals), len(p)))
        elif m.model == 'Meiboom':
            print('kex=%8.3f            red_chi2=%12.3f  dof=%d  nv=%d  np=%d'
                  % (kex, red_chi2, m.dof, len(R2_residuals), len(p)))

    return R2_residuals


def fit(m):
    """Run nonlinear least-squares fitting for all active residues.

    Builds the initial parameter vector, calls scipy.optimize.leastsq,
    and updates the global exchange-rate parameters (kAB / kBA / kex) and
    per-residue local parameters (dd, R2_0) in place.

    If the covariance matrix is None (singular fit), standard deviations
    are set to zero and a warning is printed.
    """
    # Cache numpy arrays so errFunc can use pre-allocated v_arr, R2exp_arr, R2stddev_arr
    cache_arrays(m)

    # Invalidate the leave-one-out cache: a new global solution changes it.
    m._loo_cache = None

    # --- Build initial parameter vector p0 ---
    p0 = []

    if m.model in ('Matrix', 'London'):
        p0.append(m.kAB)  # global: kAB
        p0.append(m.kBA)  # global: kBA
    elif m.model == 'Meiboom':
        p0.append(m.kex)  # global: kex

    get_local_params(m, p0)  # append per-residue dd and R2_0 values

    # --- Perform optimisation ---
    if m.verbose:
        print(p0)
    out = leastsq(m.errFunc, x0=p0, full_output=1)

    p1    = out[0]  # solution vector
    covar = out[1]  # covariance matrix (None when fit is singular)
    m.covar = covar
    _warn_if_not_converged(out)

    if covar is None:
        print('Warning: covariance matrix is None (fit may be singular). '
              'Standard deviations set to 0.')

    # --- Update global parameters from solution ---
    if m.model in ('Matrix', 'London'):
        m.kAB     = p1[0]
        m.kBA     = p1[1]
        m.kAB_std = sqrt(covar[0][0]) if covar is not None else 0.0
        m.kBA_std = sqrt(covar[1][1]) if covar is not None else 0.0
    elif m.model == 'Meiboom':
        m.kex     = p1[0]
        m.kex_std = sqrt(covar[0][0]) if covar is not None else 0.0

    set_local_params(m, p1, covar)  # unpack per-residue parameters
    if m.verbose:
        print()


def build_p0_current(m):
    """Build the parameter vector from the current global + local state.

    Layout matches errFunc: [global params, dd, R2_0..., dd, R2_0..., ...]
    for every currently-active residue.

    Returns:
        list: initial-parameter vector p0.
    """
    p0 = []
    if m.model == 'Meiboom':
        p0.append(m.kex)
    else:
        p0.append(m.kAB)
        p0.append(m.kBA)
    get_local_params(m, p0)
    return p0


def build_grid_p0(m, kex_val, pB_val, csd_val):
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
    if m.model in ('Matrix', 'London'):
        kAB = pB_val * kex_val
        kBA = kex_val - kAB
        p0.append(kAB)
        p0.append(kBA)
        dd = csd_val
    else:  # Meiboom
        p0.append(kex_val)
        dd = pB_val * (1.0 - pB_val) * csd_val ** 2

    for r in m.dataset.rsds:
        if r.active:
            p0.append(dd)
            for r2_0 in r.R2_0:
                p0.append(r2_0)
    return np.array(p0)


def init_guess_all(m):
    """Grid-search for the best initial parameter values before fitting.

    Evaluates chi² on a coarse grid of (kex, pB, csd) combinations and
    selects the combination that gives the lowest chi².  The grid covers::

        kex  : [100, 500, 1000, 2000, 4000, 10000] s^-1
        pB   : [0.05, 0.1, 0.2, 0.4]
        csd  : [0.1, 0.5, 1.0, 2.0, 4.0, 10.0] ppm

    Sets m.kAB / m.kBA (or m.kex) and per-residue dd to the
    best-guess values, and initialises R2_0 from the data.
    """
    # Coarse parameter grids for the global search
    kex = [100.0, 500.0, 1000.0, 2000.0, 4000.0, 10000.0]
    pB  = [0.05, 0.1, 0.2, 0.4]
    csd = [0.1, 0.5, 1.0, 2.0, 4.0, 10.0]

    # Cache numpy arrays once before the grid-search loop
    cache_arrays(m)

    # Initialise R2_0 from the last (highest-nu) experimental data point
    init_R2_0(m)

    # --- Exhaustive grid search ---
    # Uses _build_grid_p0 to avoid mutating instance state per grid point.
    minChi2 = 9.0e99
    ii_kex = ii_pB = ii_csd = 0

    for i_kex, kv in enumerate(kex):
        for i_pB, pv in enumerate(pB):
            for i_csd, cv in enumerate(csd):
                p0 = build_grid_p0(m, kv, pv, cv)
                residuals = err_func(m, p0)
                chi2 = float(np.dot(residuals, residuals))
                if chi2 < minChi2:
                    minChi2 = chi2
                    ii_kex, ii_pB, ii_csd = i_kex, i_pB, i_csd

    # --- Apply the best-guess values ---
    print('Initial values:')
    if m.model in ('Matrix', 'London'):
        m.kAB = pB[ii_pB] * kex[ii_kex]
        m.kBA = kex[ii_kex] - m.kAB
        init_dd(m, csd[ii_csd])
        print('kAB: %8.3f kBA: %8.3f CSD: %8.3f' % (m.kAB, m.kBA, csd[ii_csd]))
    elif m.model == 'Meiboom':
        m.kex = kex[ii_kex]
        init_dd(m, pB[ii_pB] * (1.0 - pB[ii_pB]) * csd[ii_csd] ** 2)
        print('kex: %8.3f phi: %8.3f'
              % (m.kex, pB[ii_pB] * (1.0 - pB[ii_pB]) * csd[ii_csd] ** 2))
