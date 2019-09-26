"""
Implementation for radiation diagnostics into ESMValTool
"""

from .c3s_511_basic import Basic_Diagnostic_SP
from .libs.MD_old.ESMValMD import ESMValMD
import os
import matplotlib.pyplot as plt
from .plots.basicplot import Plot2D


class fapar_Diagnostic(Basic_Diagnostic_SP):
    """
    class to implement radiation diagnostics, like e.g. global means,
    global differences, RMSD etc.
    """

    def run_diagnostic(self):
        """
        run parent diagnostic and the radiation specific diagnostic
        """
        super(fapar_Diagnostic, self).run_diagnostic()
        