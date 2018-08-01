"""
Implementation for XCO2 diagnostics into ESMValTool - This is when you do not add any of your diagnostics
"""

from c3s_511_basic import Basic_Diagnostic

class XCO2_Diagnostic(Basic_Diagnostic):
    """    
    class to implement XCO2 diagnostics, like e.g. global means,
    global differences, RMSD etc.
    """
