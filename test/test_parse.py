#!/usr/bin/env python3
"""Self-check for cpmg.CPMGDataSet.addField parsing and calcAlpha edge cases.

Run directly: `python3 test/test_parse.py` (asserts, no framework).
"""
import contextlib
import io
import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from cpmg.cpmg import CPMGDataSet, DispersionCurve, Residue

VALID = """60.12
0.040000
#nu_cpmg(Hz) R2(1/s) Esd(R2)
# K1f
50.000 22.415 0.448
100.000 22.000 0.440
"""

# Malformed: header present, but a data line appears where the first
# residue label ("# <label>") is expected — previously an unbound NameError.
MALFORMED = """60.12
0.040000
#nu_cpmg(Hz) R2(1/s) Esd(R2)
50.000 22.415 0.448
"""


def _write(text, newline=None):
    fd, path = tempfile.mkstemp(suffix='.dat')
    with os.fdopen(fd, 'w', newline=newline) as f:
        f.write(text)
    return path


def _parse_error(text, newline=None):
    """Parse text as a .dat file and return the ValueError message (file:line)."""
    path = _write(text, newline)
    try:
        try:
            CPMGDataSet().addField(path)
        except ValueError as e:
            assert os.path.basename(path) in str(e), str(e)
            return str(e)
        raise AssertionError('expected ValueError for input %r' % text)
    finally:
        os.remove(path)


HEADER = '60.12\n0.040000\n#h\n'


def test_valid_parse():
    path = _write(VALID)
    try:
        ds = CPMGDataSet()
        ds.addField(path)
        assert [r.label for r in ds.rsds] == ['K1f']
        assert ds.rsds[0].dspCurves[0].v == [50.0, 100.0]
        assert ds.rsds[0].dspCurves[0].R2exp == [22.415, 22.0]
        assert ds.rsds[0].dspCurves[0].R2stddev == [0.448, 0.440]
        assert ds.fields == [60.12] and ds.tcps == [0.04]
    finally:
        os.remove(path)


def test_malformed_raises_valueerror():
    path = _write(MALFORMED)
    try:
        ds = CPMGDataSet()
        try:
            ds.addField(path)
        except ValueError as e:
            assert 'residue label' in str(e)
            assert ':4:' in str(e)
        else:
            assert False, 'expected ValueError on missing residue label'
    finally:
        os.remove(path)


def test_invalid_inputs_name_file_and_line():
    assert ':5:' in _parse_error(HEADER + '# A1\n50 abc 0.4\n')
    assert ':5:' in _parse_error(HEADER + '# A1\n50 22 0.4 9\n')
    assert ':5:' in _parse_error(HEADER + '# A1\n50 22\n')
    assert ':5:' in _parse_error(HEADER + '# A1\n50 22 0\n')
    assert ':5:' in _parse_error(HEADER + '# A1\n50 22 nan\n')
    assert ':4:' in _parse_error(HEADER + '# A1 B2\n50 22 0.4\n')
    assert ':1:' in _parse_error('abc\n0.04\n#h\n# A1\n50 22 0.4\n')
    assert ':2:' in _parse_error('60.1\n0\n#h\n# A1\n50 22 0.4\n')
    assert ':1:' in _parse_error('-60.1\n0.04\n#h\n# A1\n50 22 0.4\n')
    assert ':4:' in _parse_error('60.1\n0.04\n#h\n50 22 0.4\n')


def test_short_files_report_missing_header_line():
    assert ':1: missing field (MHz) line' in _parse_error('')
    assert ':2: missing tcp (s) line' in _parse_error('60.12\n')
    assert ':3: missing column header line' in _parse_error('60.12\n0.04\n')


def test_label_forms_and_crlf():
    path = _write('60.1\r\n0.04\r\n#h\r\n#A1\r\n50 22 0.4\r\n\r\n# B2\r\n50 21 0.4\r\n',
                  newline='')
    try:
        ds = CPMGDataSet()
        ds.addField(path)
        assert [r.label for r in ds.rsds] == ['A1', 'B2']
        assert ds.rsds[0].dspCurves[0].R2exp == [22.0]
    finally:
        os.remove(path)


def test_label_without_data_is_dropped():
    path = _write(HEADER + '# A1\n# B2\n50 22 0.4\n')
    try:
        ds = CPMGDataSet()
        ds.addField(path)
        curves = {r.label: len(r.dspCurves) for r in ds.rsds}
        assert curves == {'A1': 0, 'B2': 1}, curves
    finally:
        os.remove(path)


def test_single_field_alpha_is_silent_on_stdout():
    ds = CPMGDataSet()
    ds.addField(_write(VALID))
    ds.fields = ds.fields[:1]
    out = io.StringIO()
    err = io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        ds.calcAlpha()
        ds.selectModelAlpha({})
    assert out.getvalue() == ''
    assert 'Two or more' in err.getvalue()


def test_alpha_matches_field_not_index():
    # Fields 60, 80, 100 MHz; residue R has curves at 60 and 100 only, so
    # its dspCurves list is shorter than the field list.
    ds = CPMGDataSet()
    low = DispersionCurve()
    low.field, low.v, low.R2exp, low.R2stddev = 60.0, [1.0, 2.0], [0.0, 1.0], [0.1, 0.1]
    high = DispersionCurve()
    high.field, high.v, high.R2exp, high.R2stddev = 100.0, [1.0, 2.0], [0.0, 2.0], [0.1, 0.1]
    res = Residue()
    res.label = 'R'
    res.dspCurves = [low, high]
    ds.rsds = [res]
    ds.fields = [60.0, 80.0, 100.0]
    ds.calcAlpha()
    # alpha = ((100+60)/(100-60)) * ((2-1)/(2+1)) = 4 * (1/3) = 4/3
    assert abs(res.alpha - 4.0 / 3.0) < 1e-12, res.alpha

    res.dspCurves = [low]
    ds.calcAlpha()
    assert res.alpha == 9999.0


if __name__ == '__main__':
    test_valid_parse()
    test_malformed_raises_valueerror()
    test_invalid_inputs_name_file_and_line()
    test_short_files_report_missing_header_line()
    test_label_forms_and_crlf()
    test_label_without_data_is_dropped()
    test_single_field_alpha_is_silent_on_stdout()
    test_alpha_matches_field_not_index()
    print('test_parse: OK')
