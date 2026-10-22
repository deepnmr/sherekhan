#!/usr/bin/env python3
"""Self-checks for cpmg.model_2state hygiene (todo 5).

Run directly: `python3 test/test_models.py` (asserts, no framework).
"""
import contextlib
import io
import os
import re
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from cpmg import __version__
from cpmg.model_2state import CPMG_model


def _write_dat(field, tcp, residues):
    """Write a .dat file; residues is a list of (label, [(v, R2, sd), ...])."""
    lines = ['%f' % field, '%f' % tcp, '#nu_cpmg(Hz) R2(1/s) Esd(R2)']
    for label, rows in residues:
        lines.append('# %s' % label)
        lines += ['%f %f %f' % row for row in rows]
    fd, path = tempfile.mkstemp(suffix='.dat')
    with os.fdopen(fd, 'w') as f:
        f.write('\n'.join(lines) + '\n')
    return path


def test_echo_count_rounding_is_tolerant_and_strict():
    m = CPMG_model()
    base = m.R2_clc_matrix(20., 1000., 3 * 2 * np.pi * 60.1, 20., 100., .04)
    near = m.R2_clc_matrix(20., 1000., 3 * 2 * np.pi * 60.1, 20., 100., 0.039999999999999994)
    assert abs(base - near) < 1e-12, (base, near)
    near2 = m.R2_clc_matrix2(20., 1000., 3 * 2 * np.pi * 60.1, 20., 100., 0.039999999999999994)
    base2 = m.R2_clc_matrix2(20., 1000., 3 * 2 * np.pi * 60.1, 20., 100., .04)
    assert abs(base2 - near2) < 1e-12
    for fn in (m.R2_clc_matrix, m.R2_clc_matrix2):
        try:
            fn(20., 1000., 100., 20., 100., 1e-6)
        except ValueError as e:
            assert 'not a positive integer' in str(e)
        else:
            raise AssertionError('tcp*nu=1e-4 must raise ValueError')


def test_unknown_model_raises_value_error():
    m = CPMG_model()
    m.model = 'Foo'
    m.dataset = m.dataset  # empty dataset: errFunc still checks the model first
    try:
        m.errFunc([1.0])
    except ValueError as e:
        assert 'wrong model' in str(e)
    else:
        raise AssertionError('unknown model must raise ValueError')


def test_matrix_scalar_eig_and_batch_agree():
    v = np.arange(50, 1001, 50, dtype=float)
    m = CPMG_model()
    tcp = .04
    n_arr = np.round(tcp * v).astype(int)
    delta = 1.0 / (4.0 * v)
    for kAB, kBA, dw, r20 in ((20., 1000., 3 * 2 * np.pi * 60.1, 20.),
                              (200., 800., 1000., 15.)):
        scal = np.array([m.R2_clc_matrix(kAB, kBA, dw, r20, vi, tcp) for vi in v])
        eig = np.array([m.R2_clc_matrix2(kAB, kBA, dw, r20, vi, tcp) for vi in v])
        bat = m.R2_clc_matrix2_batch(kAB, kBA, np.array([dw]), np.array([r20]),
                                     v, delta, n_arr, tcp)[0]
        assert np.max(np.abs(scal - eig)) < 1e-10, np.max(np.abs(scal - eig))
        assert np.max(np.abs(scal - bat)) < 1e-10, np.max(np.abs(scal - bat))


def test_jackknife_influential_matches_report_flags():
    rng = np.random.RandomState(0)
    v = [100., 200., 300., 400., 500., 600., 800., 1000.]
    tcp, field = .04, 60.12
    m_true = CPMG_model()
    residues = []
    for k, (label, csd, r20) in enumerate([('A1', 1.5, 10.), ('A2', 2.0, 12.),
                                           ('A3', 2.5, 11.), ('A4', 1.8, 9.)]):
        dw = csd * 2 * np.pi * field
        rows = []
        for vi in v:
            y = m_true.R2_clc_matrix(100., 900., dw, r20, vi, tcp)
            rows.append((vi, y + rng.normal(0, 0.3), 0.3))
        residues.append((label, rows))
    path = _write_dat(field, tcp, residues)
    try:
        m = CPMG_model()
        m.dataset.addField(path)
        m.dataset.exchange = 'slow'
        m.model = 'Matrix'
        m.kAB, m.kBA = 90., 810.
        m.init_dd(2.0)
        m.init_R2_0()
        m.fit()
        jk = m.jackknifeGlobal()
        assert jk is not None and jk['se'] > 0
        text = m.jackknifeReport(jk)
        flagged = set()
        for line in text.splitlines():
            parts = line.split()
            if len(parts) == 5 and parts[-1] == '*' and parts[0] in ('A1', 'A2', 'A3', 'A4'):
                flagged.add(parts[0])
        expected = {d['label'] for d in jk['per_residue'] if d['influential']}
        assert all(isinstance(d['influential'], bool) for d in jk['per_residue'])
        assert flagged == expected, (flagged, expected)
    finally:
        os.remove(path)


def test_reportAllValues_undefined_exchange_matrix():
    rng = np.random.RandomState(1)
    v = [100., 200., 300., 400., 500., 600., 700., 800., 900., 1000.]
    tcp, field = .04, 60.12
    dw = 1 * 2 * np.pi * field
    m_true = CPMG_model()
    rows = [(vi, m_true.R2_clc_matrix(50., 950., dw, 10., vi, tcp) + rng.normal(0, 0.2), 0.2)
            for vi in v]
    path = _write_dat(field, tcp, [('K1', rows)])
    try:
        m = CPMG_model()
        m.dataset.addField(path)
        m.dataset.exchange = 'undefined'
        m.model = 'Matrix'
        m.kAB, m.kBA = 50., 950.
        m.init_dd(1.0)
        m.init_R2_0()
        m.fit()
        out = m.reportAllValues()
        assert len(out) == 1 and len(out[0]) == 1
        rec = out[0][0]
        assert {'kAB', 'kBA', 'kex', 'pB', 'chemShiftDiff', 'r2_clc'} <= set(rec)
        assert len(rec['r2_clc']) == len(v) and np.all(np.isfinite(rec['r2_clc']))
        log = m.getLogBuffer(None)
        assert 'Input file: (not provided)' in log
        assert 'kAB:' in log
    finally:
        os.remove(path)


def test_banner_version_and_citation():
    m = CPMG_model()
    assert ('ver. %s\n' % __version__) in m.programName
    assert 'Bioinformatics 29(14):1819-1820' in m.programName
    assert 'Please cite: xxx' not in m.programName


def test_no_process_exit_from_library():
    m = CPMG_model()
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        try:
            m.R2_clc_meiboom(100., 1., 10., 100., 60.)  # valid call: no exit
        except SystemExit:
            raise AssertionError('library must not call exit()')
    assert re.search(r'exit\(\)', open(os.path.join(os.path.dirname(__file__), '..',
                                                    'cpmg', 'model_2state.py')).read()) is None


def test_matrix_singular_cases_return_r2_0():
    m = CPMG_model()
    assert m.R2_clc_matrix(0., 1000., 1000., 10., 100., .04) == 10.
    assert m.R2_clc_matrix2(0., 1000., 1000., 10., 100., .04) == 10.
    assert m.R2_clc_matrix(0., 0., 1000., 10., 100., .04) == 10.
    assert m.R2_clc_matrix(20., 1000., 0., 10., 100., .04) == 10.


def test_matrix_nonfinite_eigen_point_falls_back_to_expm():
    m = CPMG_model()
    expected = 172.27349206218813
    assert abs(m.R2_clc_matrix(500., 500., 1000., 10., 100., .04) - expected) < 1e-9
    assert abs(m.R2_clc_matrix2(500., 500., 1000., 10., 100., .04) - expected) < 1e-9
    bat = m.R2_clc_matrix2_batch(500., 500., np.array([1000.]), np.array([10.]),
                                 np.array([100.]), np.array([1. / 400]), np.array([4]), .04)
    assert abs(bat[0, 0] - expected) < 1e-9, bat


def test_matrix_batch_singular_rows_return_r2_0_without_warnings():
    import warnings
    m = CPMG_model()
    v = np.arange(50, 1001, 50, dtype=float)
    tcp = .04
    n_arr = np.round(tcp * v).astype(int)
    delta = 1.0 / (4.0 * v)
    dw_n, r20_n = 3 * 2 * np.pi * 60.1, 20.
    domega = np.array([dw_n, 0.0])          # row 1 has no chemical shift difference
    R2_0 = np.array([r20_n, 7.5])
    cases = ((0., 1000.), (20., -20.))      # kAB == 0 and kAB + kBA == 0
    with warnings.catch_warnings():
        warnings.simplefilter('error', RuntimeWarning)
        for kAB, kBA in cases:
            out = m.R2_clc_matrix2_batch(kAB, kBA, domega, R2_0, v, delta, n_arr, tcp)
            assert out.shape == (2, len(v))
            assert np.all(out[0] == r20_n), out[0]
            assert np.all(out[1] == 7.5), out[1]
        # Normal row with a dw=0 row alongside: normal row matches the scalar path
        out = m.R2_clc_matrix2_batch(20., 1000., domega, R2_0, v, delta, n_arr, tcp)
    scal = np.array([m.R2_clc_matrix(20., 1000., dw_n, r20_n, vi, tcp) for vi in v])
    assert np.max(np.abs(out[0] - scal)) < 1e-10, np.max(np.abs(out[0] - scal))
    assert np.all(out[1] == 7.5), out[1]


def test_london_nonfinite_points_use_cr72_spectral_radius():
    m = CPMG_model()
    assert abs(m.R2_clc_london(50000., 950000., 1000., 10., 25.) - 10.047495213750) < 1e-6
    assert abs(m.R2_clc_london(500., 500., 1000., 10., 100.) - 180.553770725781) < 1e-6
    v = np.arange(50, 1001, 50, dtype=float)
    dw = 3 * 2 * np.pi * 60.1
    scal = np.array([m.R2_clc_london(20., 1000., dw, 20., x) for x in v])
    vec = m.R2_clc_london_vec(20., 1000., dw, 20., v)
    assert np.max(np.abs(scal - vec)) < 1e-12


def test_meiboom_zero_and_small_kex_limits():
    m = CPMG_model()
    assert m.R2_clc_meiboom(0., 1., 10., 100., 100.) == 10.
    assert abs(m.R2_clc_meiboom(1e-4, 1., 10., 100., 100.) - 10.000082246703343) < 1e-10
    assert abs(m.R2_clc_meiboom(1e-6, 1., 10., 100., 100.) - 10.000000822467033) < 1e-12
    vec = m.R2_clc_meiboom_vec(1e-4, 1., 10., np.array([100.]), 100.)
    assert abs(vec[0] - 10.000082246703343) < 1e-10


def test_meiboom_fast_limit_agrees_with_matrix():
    m = CPMG_model()
    kex, dw, r20, nu, field = 1e6, 1000., 10., 100., 100.
    pB = .05
    csd = dw / (2 * np.pi * field)
    mat = m.R2_clc_matrix(pB * kex, (1 - pB) * kex, dw, r20, nu, .04)
    mei = m.R2_clc_meiboom(kex, pB * (1 - pB) * csd ** 2, r20, nu, field)
    assert abs(mat - mei) < 1e-5, (mat, mei)


if __name__ == '__main__':
    tests = [v for k, v in sorted(globals().items()) if k.startswith('test_') and callable(v)]
    for t in tests:
        t()
        print('PASS', t.__name__)
    print('All %d model checks passed.' % len(tests))
