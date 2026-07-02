#!/usr/bin/env python3
"""
sk_prepare.py — Generate a ShereKhan JSON config from raw .dat files.

Usage
-----
    python sk_prepare.py <file1.dat> [file2.dat ...] > config.json

Reads one or more CPMG dispersion .dat files (one per magnetic field),
computes per-residue alpha values to classify the exchange regime, and
prints a ready-to-use JSON configuration to stdout.  The config is
intended to be reviewed and passed to sk_run.py.

The auto-selected model follows this logic:
  * slow exchange majority → London
  * fast exchange majority → Meiboom
  * ambiguous / single field → Matrix

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

from sys import argv, stderr
from json import dumps, encoder
from cpmg.model_2state import CPMG_model

# Note: encoder.FLOAT_REPR is a Python 2 hook and has no effect in Python 3.
# Retained here for historical compatibility.
encoder.FLOAT_REPR = lambda o: format(o, '.4f')


def main():
    """Build a ShereKhan JSON config from .dat files passed on the command line."""
    if len(argv) < 2:
        stderr.write('Usage: sk_prepare.py <file1.dat> [file2.dat ...] > config.json\n')
        raise SystemExit(1)

    # -----------------------------------------------------------------------
    # Initialise model and output config skeleton
    # -----------------------------------------------------------------------
    m2 = CPMG_model()
    m2.verbose = True

    # Start the config dictionary with project defaults
    newDict = {
        'Project Name': 'noname',
        'comments': '',
        'init': {'mode': 'guess'}   # default to grid-search initialisation
    }
    expList = []  # will hold the paths of the loaded .dat files

    # -----------------------------------------------------------------------
    # Load all .dat files provided on the command line (one per field)
    # -----------------------------------------------------------------------
    for i in range(1, len(argv)):
        m2.dataset.addField(argv[i])  # parse file, extend dataset
        expList.append(argv[i])       # record the path for the config

    newDict['datasets'] = expList

    # -----------------------------------------------------------------------
    # Classify exchange regime and select model
    # -----------------------------------------------------------------------
    m2.dataset.calcAlpha()            # compute per-residue alpha values
    m2.dataset.selectModelAlpha(newDict)  # write 'exchange' and 'model' into newDict

    # Build the per-residue list (name, alpha, active flag)
    rdlist = m2.dataset.getResidues()
    newDict['residues'] = rdlist

    # -----------------------------------------------------------------------
    # Print the JSON config to stdout
    # -----------------------------------------------------------------------
    buf = dumps(newDict, sort_keys=True, indent=4)
    print(buf)

    stderr.write('Done.\n')


if __name__ == '__main__':
    main()
