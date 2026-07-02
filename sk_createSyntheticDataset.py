#!/usr/bin/env python3
"""
sk_createSyntheticDataset.py — Generate synthetic CPMG dispersion data.

Usage
-----
    python sk_createSyntheticDataset.py

Produces two .dat files (synth-60.dat and synth-90.dat) containing
synthetic fast-exchange CPMG dispersion curves for five fictitious residues.
Gaussian noise proportional to the ideal R2eff value is added to mimic
experimental uncertainty.

The generated files can be used to test and validate the ShereKhan fitting
pipeline (sk_prepare.py → sk_run.py).

Output file format (example synth-60.dat)
------------------------------------------
    60.12                            # field (MHz)
    0.040000                         # tcp (s)
    #nu_cpmg(Hz)   R2(1/s)   Esd(R2) # header (ignored by reader)
    # K1f                            # residue label
    50.000   22.415   0.448          # nu_CPMG, R2exp, R2std
    ...

Authors
-------
Prof. Dr. Donghan Lee
Korea Basic Science Institute (KBSI)

Original code:
Adam Mazur, Bjoern Hammesfahr, Christian Griesinger,
Donghan Lee, Martin Kollmar
Max-Planck-Institute for Biophysical Chemistry, Goettingen, 2012

Copyright (c) 2025-2026 Prof. Dr. Donghan Lee, Korea Basic Science Institute (KBSI)
"""

###
# Copyright (c) 2025-2026 Prof. Dr. Donghan Lee, Korea Basic Science Institute (KBSI)
# Original code: Adam Mazur et al., MPI Goettingen, 2012
###

import numpy as np
from sys import stderr
from cpmg.model_2state import CPMG_model
from numpy.random import randn


def main():
    """Generate synthetic fast-exchange CPMG .dat files (one per field)."""
    # Use the model object only to call the R2_clc_matrix function
    m2 = CPMG_model()

    # -----------------------------------------------------------------------
    # Simulation parameters
    # -----------------------------------------------------------------------
    fields   = [60.12, 90.23]  # [MHz] — proton Larmor frequencies (1H fields)
    r20_fc   = [0.0, 2.0]      # [s^-1] — per-field R2_0 offset added to each residue's R2_0
    tcp      = 0.04             # [s]   — CPMG pulse delay (total echo time)

    # Exchange rate constants for the two-state model (fast exchange regime)
    kAB = 20.0     # [s^-1] — A→B rate (minor state population ~ kAB/(kAB+kBA) = 2%)
    kBA = 1000.0   # [s^-1] — B→A rate

    # Five residues with increasing chemical shift differences (ppm)
    res = []
    res.append({'name': 'K1f', 'csd': 1.0, 'R2_0': 20.0})  # fast exchange
    res.append({'name': 'L2f', 'csd': 1.2, 'R2_0': 20.0})
    res.append({'name': 'M3f', 'csd': 1.4, 'R2_0': 20.0})
    res.append({'name': 'N4f', 'csd': 1.6, 'R2_0': 20.0})
    res.append({'name': 'O5f', 'csd': 1.8, 'R2_0': 20.0})

    # CPMG pulsing frequencies (Hz) — 11 points from 50 to 1000 Hz
    vCPMG = np.array([50.0, 100.0, 200.0, 300.0, 400.0,
                       500.0, 600.0, 700.0, 800.0, 900.0, 1000.0])

    noiseRatio = 0.02  # Gaussian noise amplitude as a fraction of the ideal R2eff

    # -----------------------------------------------------------------------
    # Generate and write one .dat file per magnetic field
    # -----------------------------------------------------------------------
    for i in range(0, len(fields)):
        field = fields[i]

        # File header: field strength and tcp
        buff  = '%5.2f\n' % field
        buff += '%f\n'    % tcp
        buff += '#nu_cpmg(Hz)        R2(1/s)      Esd(R2)\n'

        for r in res:
            buff += '# %s\n' % r['name']  # residue label line

            for v in vCPMG:
                # Convert chemical shift difference from ppm to rad/s at this field
                domega = r['csd'] * 2 * np.pi * field

                # Intrinsic R2 with field-dependent offset
                R2_0 = r['R2_0'] + r20_fc[i]

                # Ideal (noise-free) R2eff from the exact matrix-exponential model
                r2eff_ideal = m2.R2_clc_matrix(kAB, kBA, domega, R2_0, v, tcp)

                # Noise standard deviation proportional to the ideal value
                r2effStd = r2eff_ideal * noiseRatio

                # Add Gaussian noise to simulate experimental scatter
                r2eff = r2effStd * randn() + r2eff_ideal

                buff += '%8.3f %8.3f %8.3f\n' % (v, r2eff, r2effStd)

        # Write the .dat file; filename encodes the field in MHz (integer)
        with open('synth-%d.dat' % int(field), 'w') as file1:
            file1.write(buff)

    stderr.write('Done.\n')


if __name__ == '__main__':
    main()
