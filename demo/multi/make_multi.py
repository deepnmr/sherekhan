#!/usr/bin/env python3
"""Generate structured mixed CPMG datasets: N individuals, and K global groups.

6 residues, 2 fields, exact Matrix forward model, pB=0.1, 1% noise.

Family A -- one shared global kex (1000) + N residues each with their own kex:
  m1: 1 individual, m2: 2, m3: 3, m4: 4

Family B -- residues cluster into K internally-consistent global groups:
  g2: 2 groups (600 | 1800), g3: 3 groups (400 | 1000 | 2500)

Each scenario writes f60.dat, f90.dat, <name>.conf, truth.json into its dir.
truth.json records per-residue kex, a 'group' id, and (Family A) the individuals.
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, REPO)
from cpmg.model_2state import CPMG_model

np.random.seed(2025)

FIELDS = [60.12, 90.23]
R20_FC = [0.0, 2.0]
TCP = 0.04
VCPMG = np.array([50., 100., 200., 300., 400., 500., 600., 700., 800., 900., 1000.])
NOISE = 0.01
R2_0 = 20.0
NAMES = ['R1', 'R2', 'R3', 'R4', 'R5', 'R6']
CSD = [2.0, 2.2, 2.4, 2.6, 2.8, 3.0]   # ppm, per residue
PB = 0.1

SHARED = 1000
# Family A: kex_list (shared residues use kex=1000, others are individuals).
FAMILY_A = {
    'm1': [SHARED, SHARED, SHARED, SHARED, SHARED, 2500],
    'm2': [SHARED, SHARED, SHARED, SHARED, 400,    2500],
    'm3': [SHARED, SHARED, SHARED, 400,    1600,   2500],
    'm4': [SHARED, SHARED, 400,    700,    1600,   2500],
}
# Family B: (kex_list, group-id list).  Residues sharing a group id share a kex.
FAMILY_B = {
    'g2': ([600, 600, 600, 1800, 1800, 1800], ['A', 'A', 'A', 'B', 'B', 'B']),
    'g3': ([400, 400, 1000, 1000, 2500, 2500], ['A', 'A', 'B', 'B', 'C', 'C']),
}

m = CPMG_model()


def write_scenario(name, kex_list, meta, outdir):
    os.makedirs(outdir, exist_ok=True)
    datasets = []
    for fi, field in enumerate(FIELDS):
        buff = '%5.2f\n%f\n#nu_cpmg(Hz)  R2(1/s)  Esd(R2)\n' % (field, TCP)
        for rname, kex, csd in zip(NAMES, kex_list, CSD):
            kAB, kBA = PB * kex, (1.0 - PB) * kex
            buff += '# %s\n' % rname
            for v in VCPMG:
                domega = csd * 2 * np.pi * field
                ideal = m.R2_clc_matrix(kAB, kBA, domega, R2_0 + R20_FC[fi], v, TCP)
                std = ideal * NOISE
                buff += '%8.3f %8.3f %8.3f\n' % (v, std * np.random.randn() + ideal, std)
        path = os.path.join(outdir, 'f%d.dat' % int(field))
        open(path, 'w').write(buff)
        datasets.append(os.path.basename(path))

    conf = {
        'Project Name': name, 'comments': 'scenario %s' % name,
        'datasets': datasets, 'exchange': 'slow', 'model': 'Matrix',
        'init': {'mode': 'values', 'kex': 1000.0, 'pB': PB, 'csd': 2.5},
        'compare_aic': True, 'jackknife': True,
        'residues': [{'name': n, 'flag': 'on'} for n in NAMES],
    }
    json.dump(conf, open(os.path.join(outdir, '%s.conf' % name), 'w'), indent=4)
    truth = {'name': name, 'kex': dict(zip(NAMES, kex_list)), **meta}
    json.dump(truth, open(os.path.join(outdir, 'truth.json'), 'w'), indent=4)
    print('generated %-4s -> %s' % (name, os.path.relpath(outdir, HERE)))


if __name__ == '__main__':
    for name, kex in FAMILY_A.items():
        individuals = [NAMES[i] for i, k in enumerate(kex) if k != SHARED]
        groups = ['shared' if k == SHARED else 'ind' for k in kex]
        write_scenario(name, kex, {'family': 'A', 'individuals': individuals,
                                   'group': dict(zip(NAMES, groups))}, os.path.join(HERE, name))
    for name, (kex, gids) in FAMILY_B.items():
        write_scenario(name, kex, {'family': 'B', 'group': dict(zip(NAMES, gids)),
                                   'n_groups': len(set(gids))}, os.path.join(HERE, name))
    print('done.')
