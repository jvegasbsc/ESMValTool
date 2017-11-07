"""
;;#############################################################################
;; ipcc_fig_9_42a.py
;; Author: Manuel Schlund (DLR, Germany)
;; ESMVal project
;;#############################################################################
;;
;; Description
;;      Calculate and plot the equilibrium climate sensitivity (ECS) vs. the
;;      global mean surface temperature (GMST) for several CMIP5 models
;;      (see IPCC AR5 WG1 ch. 9, fig. 9.42a)
;;
;; Required diag_script_info attributes (diagnostics specific)
;;      none
;;    
;; Optional diag_script_info attributes (diagnostic specific)
;;      none
;;
;; Required variable_info attributes (variable specific)
;;      none
;;
;; Required variable attributes (defined in namelist)
;;      none
;;
;; Caveats
;;
;; Modification history
;;      20171107-A_schl_m8: written
;;
;;#############################################################################
"""

from esmval_lib import ESMValProject
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

def main(project_info):
    """
    This is the main routine

    Parameters
    ----------
    project_info : dict
        Dictionary with project information
    """
    print("Hello, here is the dummy routine from the direct python interface!")
    print("Finished test successfully! :)")
