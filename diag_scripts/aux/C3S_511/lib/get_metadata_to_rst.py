import sys  # WBD
import os   #WBD
sys.path.append("../../../diag_scripts/lib/python/")   #WBD
from METAdata import METAdata
import csv
path_out = "/media/bmueller/Work/ESMVAL_res/work/reports/sphinx/source"

def do_report(plot_list, report_title):
    """
    - plot_list    is a list of *full path* files (PNG, ...)
    - report_title is a string documenting the type of diagnostic 
    """
    
    # Create the output directory
    # THIS IS TO BE DONE WITH FUNCTION ensure_directory
    # WILL HAVE TO BE CHANGED DUE TO SPHINX
    if not (os.path.exists(path_out)):
      os.makedirs(path_out)    


    output_file = "report_" + report_title.split()[0].lower() + ".rst"  
    file = open(path_out + "/" + output_file, "w") 

    my_title ="DIAGNOSTICS OF " + report_title
    file.write(my_title + "\n")
    file.write("=" * len(my_title) + "\n\n")
    
    for f in plot_list:
#        MD = METAdata()
        
        caption = "test"#MD.read(f).get_dict()['ESMValTool']['caption']
    
        file.write(".. figure:: " + "../../../../diag_scripts/aux/C3S_511/" + f + "\n" 
                 "   :align:   center"   + "\n" 
                 "   :width:   95%" + "\n" 
                 "   " + caption[0] + "\n"
                  )
    file.close() 

    # Go to the directory and create PDF

#def do_smm_report(csv_expert, csv_definitions):
#    """
#    """
#csv_expert = "example_csvs/example_smm_expert.csv"
#with open(csv_expert, 'rb') as csvfile:
#    s = csv.reader(csvfile, delimiter = ",")
#    for row in s:
#      print(row)    
#
#output_file = "smm_report.rst"
#
#file = open(path_out + "/" + output_file, "w")
#my_title ="SYSTEM MATURITY MATRIX"
#file.write(my_title + "\n")
#file.write("=" * len(my_title) + "\n")
#table_length = 40
#file.write("+" + ("-" * (table_length - 2)) + "+\n")
#file.write("|" + " " + "*Namelists*" + (" "  * (table_length - 2 - 1 - len("*Namelists*"))) + "|" + "\n")
#file.write("+" + ("-" * (table_length - 2)) + "+\n")
#
#
#file.close()
#    
#
#plot_list = ["../../../diag_scripts/aux/C3S_511/example_images/albedo_QA4ECV_CMIP5_MRI-ESM1_historical_r1i1p1_regionalized_smean_ts.png" ,
#             "../../../diag_scripts/aux/C3S_511/example_images/albedo_QA4ECV_CMIP5_MRI-ESM1_historical_r1i1p1_4plots_gmd.png"]
#
#do_report(plot_list, "MEAN AND VARIABILITY")
