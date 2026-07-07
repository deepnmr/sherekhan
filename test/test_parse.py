#!/usr/bin/env python3
"""Self-check for cpmg.CPMGDataSet.addField parsing edge cases.

Run directly: `python3 test/test_parse.py` (asserts, no framework).
"""
import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from cpmg.cpmg import CPMGDataSet

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


def _write(text):
    fd, path = tempfile.mkstemp(suffix='.dat')
    with os.fdopen(fd, 'w') as f:
        f.write(text)
    return path


def test_valid_parse():
    path = _write(VALID)
    try:
        ds = CPMGDataSet()
        ds.addField(path)
        assert [r.label for r in ds.rsds] == ['K1f']
        assert ds.rsds[0].dspCurves[0].v == [50.0, 100.0]
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
        else:
            assert False, 'expected ValueError on missing residue label'
    finally:
        os.remove(path)


if __name__ == '__main__':
    test_valid_parse()
    test_malformed_raises_valueerror()
    print('test_parse: OK')
