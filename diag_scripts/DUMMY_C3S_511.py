"""
;;#############################################################################
;; fAPAR Diagnostics
;; Author: Benjamin Mueller (LMU Munich, GER)
;; QA4ECV project
;;#############################################################################
;; Description
;;    Produces various general diagnostic plots and statistics for the
;;    reference data sets of the QA4ECV Project
;;
;; Required diag_script_info attributes (diagnostics specific)
;;    none
;;
;; Optional diag_script_info attributes (diagnostic specific)
;;    none
;;
;; Required variable_info attributes (variable specific)
;;    none
;;
;; Optional variable_info attributes (variable specific)
;;    none
;;
;; Caveats
;;
;; Modification history
;;    20161128-A_laue_ax: added call to write_references
;;    20160818-A_muel_bn: Routines written.
;;
;;#############################################################################
"""

# Basic Python packages
import sys
from copy import copy

# Add subfolder of the diagnostics to the path
sys.path.append('./diag_scripts/aux/C3S_511/')

from esmval_lib import ESMValProject

# Import full diagnostic routine
from c3s_511_basic import Basic_Diagnostic

def main(project_info):
    print('>>>>>>>> DUMMY_C3S_511.py is running! <<<<<<<<<<<<')

# A_laue_ax+
    E = ESMValProject(project_info)

    verbosity = E.get_verbosity()
    diag_script = E.get_diag_script_name()

# TODO Enter C3S_511 references
#    E.write_references(diag_script,              # diag script name
#                       ["A_muel_bn"],            # authors
#                       [""],                     # contributors
#                       [""],                     # diag_references
#                       ["E_qa4ecv_albedo"],      # obs_references
#                       ["P_qa4ecv"],             # proj_references
#                       project_info,
#                       verbosity,
#                       False)
# A_laue_ax-

    Diag = Basic_Diagnostic()
    
    Diag.set_info(proj_info=E)
    Diag.read_data()
    Diag.run_diagnostic()

   
    print('>>>>>>>> ENDED SUCESSFULLY!! <<<<<<<<<<<<')
    print('')
