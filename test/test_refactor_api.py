"""Refactor guard: the CPMG_model facade keeps its names and delegates exactly to cpmg.models."""
import math
import os
import subprocess
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

from cpmg import analysis, fitting, models, report
from cpmg.model_2state import CPMG_model, _derived_std

# Facade method -> module function, for the seven pure R2 formulas.
PAIRS = [
    ('R2_clc_matrix', 'r2_matrix'),
    ('R2_clc_matrix2', 'r2_matrix_eig'),
    ('R2_clc_london', 'r2_london'),
    ('R2_clc_meiboom', 'r2_meiboom'),
    ('R2_clc_matrix2_batch', 'r2_matrix_eig_batch'),
    ('R2_clc_london_vec', 'r2_london_vec'),
    ('R2_clc_meiboom_vec', 'r2_meiboom_vec'),
]

# Names that must stay reachable on CPMG_model (plan: Must NOT have + delegators).
PUBLIC_AND_PRIVATE = [
    'R2_clc_matrix', 'R2_clc_matrix2', 'R2_clc_london', 'R2_clc_meiboom',
    'R2_clc_matrix2_batch', 'R2_clc_london_vec', 'R2_clc_meiboom_vec',
    '_aic', '_build_p0_current', '_cache_arrays', 'errFunc', '_leave_one_out_kex',
    '_per_residue_global_chi2', '_refresh_global_stats', '_kex_pB', '_calc_r2',
    '_build_grid_p0', '_build_batch_plan', 'getLocalParams', 'setLocalParams',
    'fit', 'fit_individual', 'compareModelsAIC', 'jackknifeGlobal', 'toMeiboom',
    'reportAllValues', 'getLogBuffer', 'pdf', 'initGuessAll', 'init_R2_0', 'init_dd',
    'jackknifeReport', 'aicReport',
]

N_SETS = 20


def test_names_exist_on_facade():
    missing = [n for n in PUBLIC_AND_PRIVATE if not hasattr(CPMG_model, n)]
    assert not missing, missing


def test_module_functions_exist():
    for _, fn in PAIRS:
        assert callable(getattr(models, fn)), fn


def test_fitting_facade_delegates():
    for name, fn in [('init_dd', 'init_dd'), ('init_R2_0', 'init_R2_0'),
                     ('getLocalParams', 'get_local_params'), ('setLocalParams', 'set_local_params'),
                     ('_cache_arrays', 'cache_arrays'), ('_build_batch_plan', 'build_batch_plan'),
                     ('errFunc', 'err_func'), ('fit', 'fit'), ('_build_p0_current', 'build_p0_current'),
                     ('_build_grid_p0', 'build_grid_p0'), ('initGuessAll', 'init_guess_all')]:
        assert hasattr(CPMG_model, name), name
        assert callable(getattr(fitting, fn)), fn


def test_report_facade_delegates():
    for name in ('_calc_r2', 'getLogBuffer', 'reportAllValues', 'aicReport', 'jackknifeReport', 'pdf'):
        assert hasattr(CPMG_model, name), name
    assert callable(report.write_pdf) and callable(report.log_buffer)


def test_analysis_facade_delegates():
    for name, fn in [('fit_individual', 'fit_individual'), ('_aic', 'aic'),
                     ('compareModelsAIC', 'compare_models_aic'),
                     ('_refresh_global_stats', 'refresh_global_stats'),
                     ('_leave_one_out_kex', 'leave_one_out_kex'),
                     ('_per_residue_global_chi2', 'per_residue_global_chi2'),
                     ('jackknifeGlobal', 'jackknife_global'), ('_kex_pB', 'kex_pB'),
                     ('toMeiboom', 'to_meiboom')]:
        assert hasattr(CPMG_model, name), name
        assert callable(getattr(analysis, fn)), fn


def test_aic_and_derived_std_exact():
    # Values recorded from the pre-step tree (bf260f1).
    assert CPMG_model._aic(10.5, 3, 50) == (16.5, 17.02173913043478)
    assert CPMG_model._aic(7., 2, 3) == (11.0, float('inf'))
    C = np.array([[100., -180., .3], [-180., 400., -.6], [.3, -.6, .01]])
    kex_std, pB_std, phi_std = _derived_std(100., 900., 2., C)
    assert abs(kex_std - 11.832159566199) < 1e-10
    assert abs(pB_std - 0.010835128056) < 1e-10
    assert abs(phi_std - 0.057083237470) < 1e-10


def test_derived_std_bit_exact():
    # Values recorded at commit c856545 (float.hex). Bit equality holds on one
    # numpy/libm build; across builds the last ulps of sqrt/products differ, so
    # compare at 4 ulp-scale relative tolerance (1e-14).
    cases = [
        (100., 900., 2., np.array([[100., -180., .3], [-180., 400., -.6], [.3, -.6, .01]]),
         ['0x1.7aa10d193c22dp+3', '0x1.630ba45343435p-7', '0x1.d3a039c29a407p-5']),
        (50., 950., 1.5, np.diag([4., 9., 0.01]),
         ['0x1.cd82b446159f3p+1', '0x1.f39f94597965dp-10', '0x1.e3c4630d8510fp-7']),
        (300., 700., -0.4, np.array([[30., 5., .1], [5., 20., -.2], [.1, -.2, .04]]),
         ['0x1.efbdeb14f4edap+2', '0x1.f1621dea44742p-9', '0x1.12eadac49beadp-5']),
    ]
    for kAB, kBA, dd, C, want in cases:
        got = [float(x) for x in _derived_std(kAB, kBA, dd, C)]
        exp = [float.fromhex(h) for h in want]
        assert all(math.isclose(g, e, rel_tol=1e-14, abs_tol=0.0) for g, e in zip(got, exp)), (got, exp)


def test_aic_report_takes_model_name():
    per = [{'label': 'A1', 'kex': 100., 'pB': .05, 'dd': 2., 'phi': .5, 'chi2': 1., 'chi2_global': 2.,
            'delta_aicc': -1., 'z_kex': 0.5, 'preferred': 'individual'}]
    cmp = {'global': {'chi2': 3., 'k': 3, 'n': 10, 'aic': 9., 'aicc': 10.},
           'individual': {'chi2': 1., 'k': 4, 'n': 10, 'aic': 7., 'aicc': 8.},
           'delta_aic': -2., 'delta_aicc': -2., 'weight_global': .3, 'weight_individual': .7,
           'preferred': 'individual', 'per_residue': per,
           'n_pref_global': 0, 'n_pref_individual': 1}
    m = CPMG_model()
    m.model = 'Matrix'
    assert m.aicReport(cmp) == report.aic_report('Matrix', cmp)
    assert report.aic_report('Meiboom', cmp) != report.aic_report('Matrix', cmp)


def test_import_does_not_load_matplotlib():
    out = subprocess.run(
        [sys.executable, '-c',
         "import sys, cpmg.model_2state; print('matplotlib' in sys.modules)"],
        cwd=os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'),
        capture_output=True, text=True, check=True)
    assert out.stdout.strip() == 'False', out.stdout


def _same(a, b):
    return np.array_equal(np.asarray(a), np.asarray(b), equal_nan=True)


def _random_sets(seed):
    """Deterministic parameter sets: normal range plus singular-prone corners.

    tcp*nu is an integer for every nu, so the echo-count check always passes.
    """
    rng = np.random.default_rng(seed)
    sets = []
    for _ in range(N_SETS):
        kAB = rng.uniform(0., 2000.)
        kBA = rng.uniform(0., 2000.)
        dw = rng.uniform(0., 5000.)
        R20 = rng.uniform(1., 30.)
        nu = 25. * rng.integers(1, 41)
        sets.append((kAB, kBA, dw, R20, nu, .04))
    return sets


def test_facade_equals_module_exactly():
    m = CPMG_model()
    v = np.array([50., 100., 200., 400., 600., 800., 1000.])
    delta = 1.0 / (4.0 * v)
    n_arr = np.round(.04 * v).astype(int)
    sets = _random_sets(20)
    sets += [(500., 500., 1000., 10., 100., .04), (0., 1000., 1000., 10., 100., .04),
             (20., 1000., 0., 10., 100., .04), (1000., 0., 1000., 10., 100., .04)]
    assert len(sets) >= N_SETS
    for kAB, kBA, dw, R20, nu, tcp in sets:
        assert _same(m.R2_clc_matrix(kAB, kBA, dw, R20, nu, tcp),
                     models.r2_matrix(kAB, kBA, dw, R20, nu, tcp))
        assert _same(m.R2_clc_matrix2(kAB, kBA, dw, R20, nu, tcp),
                     models.r2_matrix_eig(kAB, kBA, dw, R20, nu, tcp))
        dwa = np.array([dw, dw * 0.5])
        r2a = np.array([R20, R20 + 1.0])
        assert _same(m.R2_clc_matrix2_batch(kAB, kBA, dwa, r2a, v, delta, n_arr, tcp),
                     models.r2_matrix_eig_batch(kAB, kBA, dwa, r2a, v, delta, n_arr, tcp))
        assert _same(m.R2_clc_london(kAB, kBA, dw, R20, nu),
                     models.r2_london(kAB, kBA, dw, R20, nu))
        assert _same(m.R2_clc_london_vec(kAB, kBA, dw, R20, v),
                     models.r2_london_vec(kAB, kBA, dw, R20, v))
        assert _same(m.R2_clc_meiboom(kAB, dw / 1000., R20, nu, 100.),
                     models.r2_meiboom(kAB, dw / 1000., R20, nu, 100.))
        assert _same(m.R2_clc_meiboom_vec(kAB, dw / 1000., R20, v, 100.),
                     models.r2_meiboom_vec(kAB, dw / 1000., R20, v, 100.))


def test_edge_values_unchanged():
    assert abs(models.r2_matrix_eig(500., 500., 1000., 10., 100., .04) - 172.27349206218813) < 1e-9
    assert models.r2_matrix(0., 1000., 1000., 10., 100., .04) == 10.
    assert abs(models.r2_london(50000., 950000., 1000., 10., 25.) - 10.047495213750) < 1e-6
    assert models.r2_meiboom(0., 1., 10., 100., 100.) == 10.


if __name__ == '__main__':
    test_names_exist_on_facade()
    test_module_functions_exist()
    test_fitting_facade_delegates()
    test_report_facade_delegates()
    test_analysis_facade_delegates()
    test_aic_and_derived_std_exact()
    test_derived_std_bit_exact()
    test_aic_report_takes_model_name()
    test_import_does_not_load_matplotlib()
    test_facade_equals_module_exactly()
    test_edge_values_unchanged()
    print('test_refactor_api: OK')
