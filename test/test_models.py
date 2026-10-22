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


if __name__ == '__main__':
    tests = [v for k, v in sorted(globals().items()) if k.startswith('test_') and callable(v)]
    for t in tests:
        t()
        print('PASS', t.__name__)
    print('All %d model checks passed.' % len(tests))
