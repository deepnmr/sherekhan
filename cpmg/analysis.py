"""
analysis.py - Model comparison, jackknife, leave-one-out and derived-parameter analysis.

Functions here take the CPMG_model instance as their first argument ``m``
and read and write its state in place, exactly as the former methods did.
The CPMG_model facade in model_2state.py delegates to these functions.

Attributes an ``m`` must expose: kAB, kBA, kex, kAB_std, kBA_std, kex_std,
dataset, model, verbose, covar, plus errFunc and the fitting helpers used by the
analysis (fitting.cache_arrays, fitting.build_p0_current, fitting.fit).
"""

import sys
import numpy as np
from numpy import array, sqrt, pi, dot, real, log, sin, sinh, cos, cosh, exp, zeros, sort, tanh, isfinite
from scipy.optimize import leastsq

from . import fitting
from .fitting import _warn_if_not_converged


def derived_std(kAB, kBA, dd, C3):
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


def fit_individual(m):
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
    rsds = m.dataset.rsds

    # --- Snapshot global state so we can restore it afterwards ---
    snap_active = [r.active for r in rsds]
    snap_glob   = (m.kAB, m.kBA, m.kex)
    snap_local  = [(r.dd, getattr(r, 'dd_std', 0.0),
                    list(r.R2_0), list(r.R2_0_std)) for r in rsds]
    snap_verbose = m.verbose
    m.verbose = False   # individual fits would otherwise spam chi² lines

    active_idx = [i for i, r in enumerate(rsds) if snap_active[i]]
    results = []

    for i in active_idx:
        # Activate only residue i
        for k, r in enumerate(rsds):
            r.active = (k == i)
        fitting.cache_arrays(m)   # rebuild batch plan for the single residue

        p0    = fitting.build_p0_current(m)
        out   = leastsq(m.errFunc, x0=p0, full_output=1)
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
        if m.model == 'Meiboom':
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
    m.kAB, m.kBA, m.kex = snap_glob
    m.verbose = snap_verbose
    fitting.cache_arrays(m)   # rebuild plan for the full active set

    return results


def aic(chi2, k, n):
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


def compare_models_aic(m):
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
    fitting.cache_arrays(m)
    p_glob = fitting.build_p0_current(m)
    res_g  = m.errFunc(p_glob)
    chi2_g = float(np.dot(res_g, res_g))
    n      = len(res_g)
    k_g    = len(p_glob)
    aic_g, aicc_g = aic(chi2_g, k_g, n)

    # --- Individual model score ---
    indiv = fit_individual(m)
    chi2_i = sum(d['chi2'] for d in indiv)
    k_i    = sum(d['npar'] for d in indiv)
    aic_i, aicc_i = aic(chi2_i, k_i, n)

    # --- Per-residue global-vs-individual preference ---
    # Attribute to each residue its own chi2 contribution under the shared
    # global solution, then score global vs individual for that residue.
    # Parameter-count convention: the individual model pays for its extra
    # exchange parameters (nExch), the global model borrows the shared
    # exchange rate for free (local params only).  So the AIC per residue
    # answers "is this residue's own exchange rate worth its cost?".
    nExch      = 1 if m.model == 'Meiboom' else 2
    kex_global = m.kex if m.model == 'Meiboom' else m.kAB + m.kBA
    pr_glob    = per_residue_global_chi2(m)  # {label: (chi2, n)}
    # Leave-one-out kex per residue: the unbiased reference for z(kex).
    # Scoring a residue against the full global kex is circular (that kex
    # was fitted using the residue itself), so use the kex fitted from the
    # other residues only.  Falls back to the full global kex when LOO is
    # undefined (a single active residue).
    loo_kex = leave_one_out_kex(m)  # {label: kex_without_residue} or {}

    for d in indiv:
        chi2_g_r, n_r = pr_glob[d['label']]
        k_i_r = d['npar']
        k_g_r = k_i_r - nExch          # global shares the exchange param(s)
        aic_g_r, aicc_g_r = aic(chi2_g_r, k_g_r, n_r)
        aic_i_r, aicc_i_r = aic(d['chi2'], k_i_r, n_r)

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
    refresh_global_stats(m)

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


def refresh_global_stats(m):
    """Recompute m.chi2 / dof / npar / nvar for the full global solution.

    Individual fitting, the per-residue comparison, and the jackknife all
    make intermediate leastsq / errFunc calls that overwrite these scalars
    with sub-fit values.  Call this after such analyses so getLogBuffer
    reports the true full-active global-fit statistics.  Per-residue and
    global parameters are already restored by those methods; this only
    re-evaluates errFunc on the full active set to refresh the scalars.
    """
    fitting.cache_arrays(m)
    m.errFunc(fitting.build_p0_current(m))


def leave_one_out_kex(m):
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
    rsds = m.dataset.rsds
    active_labels = tuple(r.label for r in rsds if r.active)
    if len(active_labels) < 2:
        return {}

    cache = getattr(m, '_loo_cache', None)
    if cache is not None and cache[0] == active_labels:
        return dict(cache[1])

    # Snapshot the full global solution for warm-start and restore.
    snap_active = [r.active for r in rsds]
    snap_local  = [(r.dd, getattr(r, 'dd_std', 0.0),
                    list(r.R2_0), list(r.R2_0_std)) for r in rsds]
    snap_glob   = (m.kAB, m.kBA, m.kex)
    snap_std    = (m.kAB_std, m.kBA_std, m.kex_std)
    snap_verbose = m.verbose

    def restore_full():
        for k, r in enumerate(rsds):
            r.active = snap_active[k]
            dd, dd_std, r2_0, r2_0_std = snap_local[k]
            r.dd = dd
            r.dd_std = dd_std
            r.R2_0 = list(r2_0)
            r.R2_0_std = list(r2_0_std)
        m.kAB, m.kBA, m.kex = snap_glob
        m.kAB_std, m.kBA_std, m.kex_std = snap_std

    m.verbose = False   # keep the N refits quiet
    out = {}
    for i, r in enumerate(rsds):
        if not snap_active[i]:
            continue
        restore_full()          # each drop warm-starts from the full solution
        r.active = False
        fitting.fit(m)              # refit global on the remaining residues
        out[r.label] = m.kex if m.model == 'Meiboom' else m.kAB + m.kBA

    restore_full()
    m.verbose = snap_verbose
    refresh_global_stats(m)   # restore full-active global-fit statistics
    m._loo_cache = (active_labels, dict(out))
    return dict(out)


def per_residue_global_chi2(m):
    """Per-residue chi2 contribution under the current global solution.

    Evaluates each active residue in isolation at the shared global
    exchange parameters (and that residue's global-fit dd / R2_0 values),
    giving the residue's own contribution to the global chi2.  No refitting
    is done; the global state is restored before returning.

    Returns:
        dict: {residue_label: (chi2, n_points)} for every active residue.
    """
    rsds = m.dataset.rsds
    snap_active  = [r.active for r in rsds]
    snap_verbose = m.verbose
    m.verbose = False

    out = {}
    for i, r in enumerate(rsds):
        if not snap_active[i]:
            continue
        for k, rr in enumerate(rsds):
            rr.active = (k == i)
        fitting.cache_arrays(m)
        # _build_p0_current uses the (unchanged) global exchange rate and
        # this residue's global-fit dd / R2_0, so errFunc returns exactly
        # the residue's contribution to the global chi2.
        p   = fitting.build_p0_current(m)
        res = m.errFunc(p)
        out[r.label] = (float(np.dot(res, res)), len(res))

    for k, rr in enumerate(rsds):
        rr.active = snap_active[k]
    m.verbose = snap_verbose
    fitting.cache_arrays(m)
    return out


def jackknife_global(m):
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
    rsds = m.dataset.rsds
    # Ordered active labels; leave-one-out kex comes from the shared helper
    # (also used by compareModelsAIC, so the N refits happen only once).
    active_labels = [r.label for r in rsds if r.active]
    N = len(active_labels)
    if N < 2:
        return None

    kex_full = m.kex if m.model == 'Meiboom' else m.kAB + m.kBA
    loo   = leave_one_out_kex(m)               # {label: kex_without_residue}
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


def kex_pB(m):
    """Derive (kex, kex_std, pB, pB_std) from kAB/kBA and their stds."""
    kex     = m.kAB + m.kBA
    if m.covar is not None:
        # dd entries do not enter kex or pB; dd=0 only zeroes the phi part
        kex_std, pB_std, _ = derived_std(m.kAB, m.kBA, 0.0, m.covar[:3, :3])
        return kex, kex_std, m.kAB / kex, pB_std
    kex_std = sqrt(m.kAB_std ** 2 + m.kBA_std ** 2)
    pB      = m.kAB / kex
    pB_std  = sqrt(m.kAB_std ** 2 * (1.0 / kex - m.kAB / kex ** 2) ** 2
                   + m.kBA_std ** 2 * (m.kAB / kex ** 2) ** 2)
    return kex, kex_std, pB, pB_std


def to_meiboom(m):
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
    # (with error propagation) — also stores kex / kex_std on m.
    m.kex, m.kex_std, pB, pB_std = kex_pB(m)

    for r in m.dataset.rsds:
        if r.active:
            # Save original csd before overwriting r.dd with phi
            dd_orig     = r.dd
            dd_std_orig = r.dd_std
            r.dd_copy   = dd_orig        # kept for plot/report reconstruction

            # Replace dd with phi = pB*(1-pB)*dd^2
            r.dd = pB * (1.0 - pB) * dd_orig ** 2
            # Propagate uncertainty through phi = pB*(1-pB)*dd^2
            if m.covar is not None:
                idx = [0, 1, r._dd_idx]
                _, _, r.dd_std = derived_std(m.kAB, m.kBA, dd_orig,
                                              m.covar[np.ix_(idx, idx)])
            else:
                r.dd_std = sqrt(
                    (dd_std_orig * 2.0 * pB * (1.0 - pB) * dd_orig) ** 2
                    + (pB_std * (1.0 - 2.0 * pB) * dd_orig ** 2) ** 2
                )
