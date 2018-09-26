"""
;;#############################################################################
;; File_packer_C3S_511
;; Author: Benjamin Mueller (LMU Munich, GER)
;; C3S_511 project
;;#############################################################################
;; Description
;;    Merges climo-files into one for subsequent sub-time-steps.
;;    This is limiting the required memory for preprocessing.
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
;;   Sub-time-steps may not overlap!
;;
;; Modification history
;;    20180921-A_muel_bn: Routines written.
;;
;;#############################################################################
"""

# Basic Python packages
import os
from cdo import Cdo as CDO 
from difflib import SequenceMatcher

from esmval_lib import ESMValProject

def string_matcher(string1,string2):
    match = SequenceMatcher(None, string1, string2).find_longest_match(0, len(string1), 0, len(string2))
    return string1[match.a: match.a + match.size]

def main(project_info):
    print('>>>>>>>> File_packer_C3S_511.py is running! <<<<<<<<<<<<')

    E = ESMValProject(project_info)
    
    
    climo_files = E.get_all_clim_models().keys()
    climo_paths = [os.path.dirname(cf) for cf in climo_files]
    
    matching_path_string = climo_paths[0]
    if len(climo_paths)>1:
        i=1
        while i<len(climo_paths):
            matching_path_string = string_matcher(climo_paths[i],matching_path_string)
            i+=1
    if not os.path.isdir(matching_path_string):
        raise ValueError("The requested model files do not origin in the same directory!")
        
    matching_file_string = climo_files[0]
    if len(climo_files)>1:
        i=1
        while i<len(climo_files):
            matching_file_string = string_matcher(climo_files[i],matching_file_string)
            i+=1
    
    years = [cf.replace(matching_file_string,'').split(".")[0].split("-") for cf in climo_files]
    s_years = min([y[0] for y in years])
    e_years = max([y[1] for y in years])
    
    output_file = matching_file_string + s_years + '-' + e_years + '.nc'
    
    if os.path.isfile(output_file):
        if E.get_configfile()=="overwrite":
            pass
        else:
            raise ValueError("The requested climo file is already existing! Please investigate: " + output_file)
    
    cdo = CDO()
    cdo.mergetime(input=" ".join(climo_files),output=output_file)
    
   
    print('>>>>>>>> ENDED SUCESSFULLY!! <<<<<<<<<<<<')
    print('')