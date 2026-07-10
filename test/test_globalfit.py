#!/usr/bin/env python3
"""Self-check for the global/individual comparison and jackknife validation.

Runs a real Meiboom global fit on the bundled fast-exchange synthetic data,
then checks the invariants of compareModelsAIC() and jackknifeGlobal().

Run directly: `python3 test/test_globalfit.py` (asserts, no framework).
Requires numpy/scipy/matplotlib (same deps as the fitting code).
"""
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from cpmg.model_2state import CPMG_model


def _build_fitted_model():
    """Load the bundled fast synthetic data and run a Meiboom global fit."""
    m = CPMG_model()
    m.verbose = False
    m.dataset.addField(os.path.join(HERE, 'fast', 'in', 'fast-synth-60.dat'))
    m.dataset.addField(os.path.join(HERE, 'fast', 'in', 'fast-synth-90.dat'))
    m.dataset.exchange = 'fast'
    m.model = 'Meiboom'
    m.initGuessAll()   # grid-search starting point (what a real run does)
    m.init_R2_0()
    m.fit()
    return m


def test_per_residue_partition():
    """Per-residue global/individual chi2 must sum to the totals."""
    m = _build_fitted_model()
    cmp = m.compareModelsAIC()

    g_sum = sum(d['chi2_global'] for d in cmp['per_residue'])
    i_sum = sum(d['chi2'] for d in cmp['per_residue'])
    assert math.isclose(g_sum, cmp['global']['chi2'], rel_tol=1e-6), \
        'per-residue global chi2 must sum to the global total'
    assert math.isclose(i_sum, cmp['individual']['chi2'], rel_tol=1e-6), \
        'per-residue individual chi2 must sum to the individual total'

    # A residue always fits at least as well with its own exchange rate.
    for d in cmp['per_residue']:
        assert d['chi2'] <= d['chi2_global'] + 1e-6, \
            'individual chi2 cannot exceed global chi2 for a residue'
        assert d['preferred'] in ('global', 'individual')
        assert d['k_indiv'] - d['k_global'] == 1  # Meiboom shares one param (kex)

    n_glob = sum(1 for d in cmp['per_residue'] if d['preferred'] == 'global')
    assert cmp['n_pref_global'] == n_glob
    assert cmp['n_pref_individual'] == len(cmp['per_residue']) - n_glob


def test_global_stats_restored():
    """compareModelsAIC and jackknife must not corrupt the global fit stats."""
    m = _build_fitted_model()
    n_points = sum(len(dsp.R2exp) for r in m.dataset.rsds if r.active
                   for dsp in r.dspCurves)
    chi2_ref, nvar_ref, npar_ref = m.chi2, m.nvar, m.npar

    m.compareModelsAIC()
    assert m.nvar == nvar_ref == n_points, 'nvar clobbered by compareModelsAIC'
    assert m.npar == npar_ref, 'npar clobbered by compareModelsAIC'
    assert math.isclose(m.chi2, chi2_ref, rel_tol=1e-9), 'chi2 clobbered'

    m.jackknifeGlobal()
    assert m.nvar == nvar_ref, 'nvar clobbered by jackknifeGlobal'
    assert m.npar == npar_ref, 'npar clobbered by jackknifeGlobal'
    assert math.isclose(m.chi2, chi2_ref, rel_tol=1e-9), 'chi2 clobbered by jackknife'


def test_zkex_uses_leave_one_out_reference():
    """z(kex) must reference the leave-one-out kex (unbiased), not the full
    global kex, and the jackknife must reuse those same refits (shared cache)."""
    m = _build_fitted_model()
    cmp = m.compareModelsAIC()

    # The z-reference for each residue is the kex fit from the OTHER residues,
    # which is exactly what the jackknife reports as kex_drop.
    jk = m.jackknifeGlobal()
    kdrop = {d['label']: d['kex_drop'] for d in jk['per_residue']}
    for d in cmp['per_residue']:
        assert 'kex_ref' in d
        assert math.isclose(d['kex_ref'], kdrop[d['label']], rel_tol=1e-9), \
            'z(kex) reference must equal the leave-one-out (jackknife) kex'
        # The reference must NOT be the full global kex (that is the biased one).
        assert not math.isclose(d['kex_ref'], d['kex_global'], rel_tol=1e-12), \
            'leave-one-out reference should differ from the full global kex'


def test_jackknife():
    """Jackknife returns one entry per residue with consistent deltas."""
    m = _build_fitted_model()
    n_active = sum(1 for r in m.dataset.rsds if r.active)
    jk = m.jackknifeGlobal()

    assert jk is not None and jk['n'] == n_active
    assert len(jk['per_residue']) == n_active
    assert jk['se'] >= 0.0
    assert math.isfinite(jk['kex_full']) and jk['kex_full'] > 0.0
    for d in jk['per_residue']:
        assert math.isclose(d['delta_kex'], jk['kex_full'] - d['kex_drop'], rel_tol=1e-9)

    # kex is unchanged by running the jackknife (state fully restored).
    assert math.isclose(m.kex, jk['kex_full'], rel_tol=1e-9)


def test_jackknife_needs_two_residues():
    """Jackknife is undefined with fewer than 2 active residues."""
    m = _build_fitted_model()
    active = [r for r in m.dataset.rsds if r.active]
    for r in active[1:]:
        r.active = False
    assert m.jackknifeGlobal() is None


if __name__ == '__main__':
    test_per_residue_partition()
    test_global_stats_restored()
    test_zkex_uses_leave_one_out_reference()
    test_jackknife()
    test_jackknife_needs_two_residues()
    print('test_globalfit: OK')
