"""Refactor guard: the CPMG_model facade keeps its names and delegates exactly to cpmg.models."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

from cpmg import models
from cpmg.model_2state import CPMG_model

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


def test_names_exist_on_facade():
    missing = [n for n in PUBLIC_AND_PRIVATE if not hasattr(CPMG_model, n)]
    assert not missing, missing


def test_module_functions_exist():
    for _, fn in PAIRS:
        assert callable(getattr(models, fn)), fn


def _cases():
    """Parameter sets: normal values plus the todo-8 edge cases."""
    m = CPMG_model()
    out = []
    for kAB, kBA, dw, r20, nu, tcp in [
        (20., 1000., 3 * 2 * np.pi * 60.1, 20., 100., .04),
        (200., 800., 1000., 15., 300., .04),
        (500., 500., 1000., 10., 100., .04),
        (0., 1000., 1000., 10., 100., .04),
        (20., 1000., 0., 10., 100., .04),
        (50., 950., 2000., 12., 50., .04),
    ]:
        out.append(('matrix', (kAB, kBA, dw, r20, nu, tcp)))
    for kAB, kBA, dw, r20, nu in [
        (50000., 950000., 1000., 10., 25.),
        (500., 500., 1000., 10., 100.),
        (20., 1000., 3 * 2 * np.pi * 60.1, 20., 100.),
    ]:
        out.append(('london', (kAB, kBA, dw, r20, nu)))
    for kex, phi, r20, nu, field in [
        (0., 1., 10., 100., 100.),
        (1e-4, 1., 10., 100., 100.),
        (1e-6, 1., 10., 100., 100.),
        (100., 0.5, 12., 200., 60.),
    ]:
        out.append(('meiboom', (kex, phi, r20, nu, field)))
    assert m is not None
    return out


def test_facade_equals_module_exactly():
    m = CPMG_model()
    v = np.array([50., 100., 200., 400., 600., 800., 1000.])
    for kind, args in _cases():
        if kind == 'matrix':
            kAB, kBA, dw, r20, nu, tcp = args
            assert m.R2_clc_matrix(*args) == models.r2_matrix(*args)
            assert m.R2_clc_matrix2(*args) == models.r2_matrix_eig(*args)
            dwa = np.array([dw, dw * 0.5])
            r2a = np.array([r20, r20 + 1.0])
            got = m.R2_clc_matrix2_batch(kAB, kBA, dwa, r2a, v, 1.0 / (4.0 * v),
                                         np.round(tcp * v).astype(int), tcp)
            want = models.r2_matrix_eig_batch(kAB, kBA, dwa, r2a, v, 1.0 / (4.0 * v),
                                              np.round(tcp * v).astype(int), tcp)
            assert np.array_equal(got, want)
        elif kind == 'london':
            assert m.R2_clc_london(*args) == models.r2_london(*args)
            kAB, kBA, dw, r20, nu = args
            assert np.array_equal(m.R2_clc_london_vec(kAB, kBA, dw, r20, v),
                                  models.r2_london_vec(kAB, kBA, dw, r20, v))
        else:
            assert m.R2_clc_meiboom(*args) == models.r2_meiboom(*args)
            kex, phi, r20, nu, field = args
            assert np.array_equal(m.R2_clc_meiboom_vec(kex, phi, r20, v, field),
                                  models.r2_meiboom_vec(kex, phi, r20, v, field))


def test_edge_values_unchanged():
    assert abs(models.r2_matrix_eig(500., 500., 1000., 10., 100., .04) - 172.27349206218813) < 1e-9
    assert models.r2_matrix(0., 1000., 1000., 10., 100., .04) == 10.
    assert abs(models.r2_london(50000., 950000., 1000., 10., 25.) - 10.047495213750) < 1e-6
    assert models.r2_meiboom(0., 1., 10., 100., 100.) == 10.


if __name__ == '__main__':
    test_names_exist_on_facade()
    test_module_functions_exist()
    test_facade_equals_module_exactly()
    test_edge_values_unchanged()
    print('test_refactor_api: OK')
