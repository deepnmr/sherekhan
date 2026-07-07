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


class DispersionCurve:
    """R2 relaxation-dispersion data measured at a single magnetic field.

    Attributes
    ----------
    field    : float
        Static magnetic field strength in MHz (proton Larmor frequency).
    tcp      : float
        CPMG pulse delay (time between 180° pulses) in seconds.
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
        self.field = 0.0    # MHz  — static field (proton Larmor frequency)
        self.tcp = 0.0      # s    — CPMG inter-pulse delay
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
        Magnetic field strengths (MHz) in the order they were loaded.
    tcps     : list of float
        CPMG pulse delays (s) corresponding to each field.
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

        Residues that already exist (by label) are extended with a new
        DispersionCurve; new residues are created automatically.
        Empty dispersion curves (no data points) are discarded at the end.

        Args:
            fileName (str): Path to the .dat input file.
        """
        with open(fileName, 'r') as inFile:

            # First line: static field strength in MHz
            line = inFile.readline()
            currentField = float(line)
            self.fields.append(currentField)

            # Second line: CPMG pulse delay (tcp) in seconds
            line = inFile.readline()
            currentTcp = float(line)
            self.tcps.append(currentTcp)

            # Third line is the column-header comment — skip it
            line = inFile.readline()

            # Read residue blocks until end of file
            line = inFile.readline()

            resLabel = None  # last-seen residue label; None until the first label line
            while line:
                sline = line.split()
                if len(sline) == 2:
                    # Lines with exactly 2 tokens are residue label lines: "# <label>"
                    resLabel = sline[1]

                # A data block must be preceded by a label line; fail clearly if not
                if resLabel is None:
                    raise ValueError(
                        '%s: expected a residue label line ("# <label>") but got: %r'
                        % (fileName, line.rstrip()))

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

                # Read the dispersion curve data points that follow the label line
                line = inFile.readline()
                sline = line.split()
                dc = DispersionCurve()
                dc.field = currentField
                dc.tcp = currentTcp
                while len(sline) == 3:
                    # Each data line: nu_CPMG(Hz)  R2exp(1/s)  R2stddev(1/s)
                    dc.v.append(float(sline[0]))
                    dc.R2exp.append(float(sline[1]))
                    dc.R2stddev.append(float(sline[2]))
                    line = inFile.readline()
                    if line:
                        sline = line.split()
                    else:
                        sline = []
                res.dspCurves.append(dc)

        # Discard any dispersion curves that ended up with no data points
        # (can happen when a residue label appears but has no following data)
        for r in self.rsds:
            r.dspCurves = [dc for dc in r.dspCurves if len(dc.v) > 0]

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

        Residues with flat dispersion curves, identical fields, or missing
        curves at the required field indices are assigned alpha = 9999.
        """
        if len(self.fields) <= 1:
            # Alpha calculation requires at least two field strengths
            print('Two or more fields strengths are necessary for alpha calculation')
            for r in self.rsds:
                r.alpha = 9999.0
            return

        # Identify the lowest and highest field indices
        # First occurrence of the lowest / highest field (matches strict-< / strict-> scan)
        B0_low_index = self.fields.index(min(self.fields))
        B0_high_index = self.fields.index(max(self.fields))

        for r in self.rsds:
            # Guard: residue must have dispersion curves at both field indices
            if len(r.dspCurves) <= max(B0_low_index, B0_high_index):
                r.alpha = 9999.0
                continue

            B0_low = r.dspCurves[B0_low_index].field
            B0_high = r.dspCurves[B0_high_index].field

            # Rex estimate = max(R2) - min(R2) at each field
            Rex_low  = max(r.dspCurves[B0_low_index].R2exp)  - min(r.dspCurves[B0_low_index].R2exp)
            Rex_high = max(r.dspCurves[B0_high_index].R2exp) - min(r.dspCurves[B0_high_index].R2exp)

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
            print('Two or more fields strengths are necessary for alpha calculation')
            print('Setting all alpha values to 9999')
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
