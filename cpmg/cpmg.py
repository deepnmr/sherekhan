###
# Copyright (c) 2025-2026 Prof. Dr. Donghan Lee, Korea Basic Science Institute (KBSI)
# Original code: Adam Mazur et al., MPI Goettingen, 2012
###
"""
cpmg.py — Data structures for NMR CPMG dispersion experiments.

Classes
-------
DispersionCurve : R2 dispersion data at a single magnetic field strength.
Residue         : A protein residue observed at one or more field strengths.
CPMGDataSet     : Collection of all residues and their dispersion curves.
"""

import math
import sys


def _number(token, fileName, lineno, what):
    """Parse one finite number from a .dat file, raising ValueError with file:line."""
    try:
        value = float(token)
    except ValueError:
        value = None
    if value is None or not math.isfinite(value):
        raise ValueError('%s:%d: %s must be a finite number, got %r'
                         % (fileName, lineno, what, token))
    return value


class DispersionCurve:
    """R2 relaxation-dispersion data measured at a single magnetic field.

    Attributes
    ----------
    field    : float
        Larmor frequency in MHz of the observed nucleus (e.g. 15N at 600 MHz
        spectrometer frequency -> 60.8 MHz).
    tcp      : float
        Total CPMG relaxation delay in seconds.
    v        : list of float
        CPMG pulsing frequencies (nu_CPMG) in Hz.
    R2exp    : list of float
        Experimentally measured effective R2 relaxation rates in s^-1.
    R2stddev : list of float
        Standard deviations of the measured R2 values in s^-1.
    exchange : str
        Exchange regime label (e.g. 'fast', 'slow', 'undef').
    """

    def __init__(self):
        self.field = 0.0    # MHz  — Larmor frequency of the observed nucleus
        self.tcp = 0.0      # s    — total CPMG relaxation delay
        self.v = []         # Hz   — CPMG pulsing frequencies
        self.R2exp = []     # 1/s  — measured effective R2
        self.R2stddev = []  # 1/s  — experimental uncertainties

        self.exchange = 'undef'  # exchange regime, set later by model selection


class Residue:
    """A single protein residue with CPMG dispersion data at one or more fields.

    Attributes
    ----------
    label     : str
        Residue identifier string (e.g. 'K1f', 'A23').
    dspCurves : list of DispersionCurve
        One dispersion curve per magnetic field strength.
    active    : bool
        Whether this residue is included in the current fit.
    dd_copy   : float or None
        Copy of the original chemical shift difference (ppm) before
        toMeiboom() replaces r.dd with the phi parameter.  Set by
        toMeiboom() for active fast-exchange residues; None otherwise.
    """

    def __init__(self):
        self.label = ''
        self.dspCurves = []      # DispersionCurve objects, one per field
        self.active = True       # include this residue in the fit by default
        self.dd_copy = None      # set by toMeiboom() for fast+Matrix/London runs

        # Model-dependent attributes (set during fitting):
        #   dd       — chemical shift difference (ppm) for Matrix/London,
        #              or phi = pB*(1-pB)*dd^2 for Meiboom
        #   dd_std   — standard deviation of dd / phi
        #   R2_0     — intrinsic R2 (list, one per field)
        #   R2_0_std — standard deviations of R2_0 values


class CPMGDataSet:
    """Container for all residues and their CPMG dispersion curves.

    Data are loaded field-by-field from plain-text .dat files using
    addField().  Alpha values used for exchange-regime classification are
    computed by calcAlpha(), and model selection is performed by
    selectModelAlpha().

    Attributes
    ----------
    rsds     : list of Residue
        All residues found across all loaded data files.
    fields   : list of float
        Larmor frequencies (MHz) of the observed nucleus, in the order loaded.
    tcps     : list of float
        Total CPMG relaxation delays (s) corresponding to each field.
    exchange : str
        Exchange regime: 'fast', 'slow', or 'undefined'.
    """

    def __init__(self):
        self.rsds = []    # list of Residue objects
        self.fields = []  # field strengths in MHz
        self.tcps = []    # CPMG pulse delays in seconds

    def addField(self, fileName):
        """Parse a .dat file and add its dispersion curves to the dataset.

        The expected file format is::

            <field_MHz>
            <tcp_s>
            #nu_cpmg(Hz)  R2(1/s)  Esd(R2)   ← header line (ignored)
            # <ResidueLabel>                   ← residue label line
            <v1>  <R2_1>  <std_1>
            <v2>  <R2_2>  <std_2>
            ...
            # <NextResidueLabel>
            ...

        Blank lines are skipped.  Any other malformed line raises a
        ValueError naming the file and line number.  Residues that already
        exist (by label) are extended with a new DispersionCurve; new residues
        are created automatically.  Labels with no data points are discarded.

        Args:
            fileName (str): Path to the .dat input file.

        Raises:
            FileNotFoundError: if the file does not exist.
            ValueError: if the file content is malformed (message: file:line).
        """
        with open(fileName, 'r') as inFile:
            lines = inFile.read().splitlines()

        # Header lines 1-3 must exist; a missing one is reported at the line it should be on
        if len(lines) < 1:
            raise ValueError('%s:1: missing field (MHz) line' % fileName)
        if len(lines) < 2:
            raise ValueError('%s:2: missing tcp (s) line' % fileName)
        if len(lines) < 3:
            raise ValueError('%s:3: missing column header line' % fileName)

        # Line 1: Larmor frequency of the observed nucleus (MHz)
        currentField = _number(lines[0], fileName, 1, 'field (MHz)')
        if currentField <= 0.0:
            raise ValueError('%s:1: field (MHz) must be > 0, got %r' % (fileName, lines[0]))

        # Line 2: total CPMG relaxation delay (s)
        currentTcp = _number(lines[1], fileName, 2, 'tcp (s)')
        if currentTcp <= 0.0:
            raise ValueError('%s:2: tcp (s) must be > 0, got %r' % (fileName, lines[1]))

        # Line 3 is the column-header comment — skipped

        # Blocks of (label, DispersionCurve) in file order
        blocks = []
        label = None
        dc = None
        for lineno, line in enumerate(lines[3:], start=4):
            stripped = line.strip()
            if not stripped:
                continue

            if stripped.startswith('#'):
                tokens = stripped[1:].split()
                if len(tokens) != 1:
                    raise ValueError('%s:%d: residue label line must be "# <label>", got %r'
                                     % (fileName, lineno, line))
                if dc is not None:
                    blocks.append((label, dc))
                label = tokens[0]
                dc = DispersionCurve()
                dc.field = currentField
                dc.tcp = currentTcp
                continue

            if dc is None:
                raise ValueError('%s:%d: expected a residue label line ("# <label>") but got: %r'
                                 % (fileName, lineno, line))
            tokens = line.split()
            if len(tokens) != 3:
                raise ValueError('%s:%d: data row needs 3 numbers (nu_cpmg, R2, sigma), got %d: %r'
                                 % (fileName, lineno, len(tokens), line))
            v = _number(tokens[0], fileName, lineno, 'nu_cpmg (Hz)')
            r2 = _number(tokens[1], fileName, lineno, 'R2 (1/s)')
            sd = _number(tokens[2], fileName, lineno, 'sigma (1/s)')
            if v <= 0.0:
                raise ValueError('%s:%d: nu_cpmg (Hz) must be > 0, got %r' % (fileName, lineno, tokens[0]))
            if sd <= 0.0:
                raise ValueError('%s:%d: sigma (1/s) must be > 0, got %r' % (fileName, lineno, tokens[2]))
            dc.v.append(v)
            dc.R2exp.append(r2)
            dc.R2stddev.append(sd)
        if dc is not None:
            blocks.append((label, dc))

        for resLabel, block in blocks:
            # Find an existing Residue with this label, or create a new one
            res = None
            for r in self.rsds:
                if r.label == resLabel:
                    res = r
                    break
            if res is None:
                res = Residue()
                res.label = resLabel
                self.rsds.append(res)
            res.dspCurves.append(block)

        self.fields.append(currentField)
        self.tcps.append(currentTcp)

        # Discard any dispersion curves that ended up with no data points
        # (can happen when a residue label appears but has no following data)
        for r in self.rsds:
            r.dspCurves = [dc for dc in r.dspCurves if len(dc.v) > 0]
        # Then drop residues left without any curve (label-only residues)
        self.rsds = [r for r in self.rsds if r.dspCurves]

    # -------------------------------------------------------------------------
    # Model selection helpers
    # -------------------------------------------------------------------------

    def calcAlpha(self):
        """Calculate the alpha field-dependence parameter for every residue.

        Alpha is defined in JACS 122(12):2871 (2000), eq. 20 as::

            alpha = ((B0_high + B0_low) / (B0_high - B0_low))
                    * ((Rex_high - Rex_low) / (Rex_high + Rex_low))

        where Rex = R2_max - R2_min at a given field.

        * alpha ~ 1–2  →  fast exchange
        * alpha ~ 0–1  →  slow exchange

        The low and high field curves of each residue are located by their
        field value.  Residues lacking a curve at either extreme field, with
        flat dispersion curves, or with identical fields are assigned
        alpha = 9999.
        """
        if len(self.fields) <= 1:
            # Alpha calculation requires at least two field strengths
            sys.stderr.write('Two or more fields strengths are necessary for alpha calculation\n')
            for r in self.rsds:
                r.alpha = 9999.0
            return

        B0_low = min(self.fields)
        B0_high = max(self.fields)

        for r in self.rsds:
            low = next((dc for dc in r.dspCurves if dc.field == B0_low), None)
            high = next((dc for dc in r.dspCurves if dc.field == B0_high), None)
            if low is None or high is None:
                r.alpha = 9999.0
                continue

            # Rex estimate = max(R2) - min(R2) at each field
            Rex_low  = max(low.R2exp)  - min(low.R2exp)
            Rex_high = max(high.R2exp) - min(high.R2exp)

            # Guard: avoid division by zero for flat curves or identical fields
            field_diff = B0_high - B0_low
            rex_sum = Rex_high + Rex_low
            if field_diff == 0.0 or rex_sum == 0.0:
                r.alpha = 9999.0
                continue

            # JACS 122:2871 (2000) eq. 20
            r.alpha = ((B0_high + B0_low) / field_diff) * ((Rex_high - Rex_low) / rex_sum)

    def selectModelAlpha(self, dict1):
        """Classify the exchange regime and choose a fitting model from alpha values.

        Counts residues in fast exchange (1 < alpha < 2) and slow exchange
        (0 < alpha < 1).  The majority regime determines the global exchange
        flag.  Active/inactive residue flags are set accordingly, and the
        recommended model is written into *dict1*.

        Model selection logic:
          * slow exchange  →  London  (analytical slow-exchange formula)
          * fast exchange  →  Meiboom (analytical fast-exchange formula)
          * undefined      →  Matrix  (exact numerical solution, safest fallback)

        Args:
            dict1 (dict): Config dictionary that receives 'exchange' and
                          'model' keys (modified in place).
        """
        nSlow = 0  # residues in slow exchange (0 < alpha < 1)
        nFast = 0  # residues in fast exchange (1 < alpha < 2)

        if len(self.fields) <= 1:
            sys.stderr.write('Two or more fields strengths are necessary for alpha calculation\n')
            sys.stderr.write('Setting all alpha values to 9999\n')
            dict1['exchange'] = 'undefined'
            dict1['model'] = 'Matrix'
            return

        # Tally residues by exchange regime based on alpha thresholds
        for r in self.rsds:
            if r.alpha > 1.0 and r.alpha < 2.0:
                nFast += 1
            elif r.alpha < 1.0 and r.alpha > 0.0:
                nSlow += 1

        # Determine the dominant exchange regime
        if nSlow > nFast:
            self.exchange = 'slow'
        elif nSlow < nFast:
            self.exchange = 'fast'
        else:
            self.exchange = 'undefined'

        # Activate only residues that match the dominant regime
        if self.exchange == 'slow':
            for r in self.rsds:
                r.active = r.alpha < 1.0 and r.alpha > 0.0

        if self.exchange == 'fast':
            for r in self.rsds:
                r.active = r.alpha > 1.0 and r.alpha < 2.0

        if self.exchange == 'undefined':
            # Cannot classify — disable all residues; user must decide manually
            for r in self.rsds:
                r.active = False

        # Write results into the output config dictionary
        dict1['exchange'] = self.exchange
        if self.exchange == 'slow':
            dict1['model'] = 'London'
        elif self.exchange == 'fast':
            dict1['model'] = 'Meiboom'
        else:
            dict1['model'] = 'Matrix'

    def getResidues(self):
        """Build a list of residue metadata dictionaries for JSON output.

        Each dictionary contains:
          * 'name'  — residue label string
          * 'alpha' — computed alpha value
          * 'flag'  — 'on' if the residue is active, 'off' otherwise

        Returns:
            list of dict: One entry per residue in self.rsds.
        """
        rdlist = []
        for r in self.rsds:
            # With only one field, alpha cannot be computed so activate all residues
            if len(self.fields) <= 1:
                r.active = True
            rd = {'name': r.label, 'alpha': r.alpha}
            rd['flag'] = 'on' if r.active else 'off'
            rdlist.append(rd)
        return rdlist
