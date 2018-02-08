import sys  # WBD
import os   #WBD
import shutil

#sys.path.append("../../../diag_scripts/lib/python/")   #WBD
sys.path.append(os.path.abspath("./diag_scripts/lib/python"))
sys.path.append(os.path.abspath("./diag_scripts"))

from METAdata import METAdata
import csv
#path_out = "/media/bmueller/Work/ESMVAL_res/work/reports/sphinx/source"
work_dir = "/athome/laue_ax/sphinx/out"
path_out = work_dir + os.sep + "reporting"
src_dir = path_out + os.sep + "source12345"
bld_dir = path_out + os.sep + "build12345"

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

    if not (os.path.exists(src_dir)):
      os.makedirs(src_dir)    
    if not (os.path.exists(bld_dir)):
      os.makedirs(bld_dir)    

#    output_file = "report_" + report_title.split()[0].lower() + ".rst"  
    output_file = src_dir + os.sep + "report.rst"  
    file = open(output_file, "w") 

    my_title ="DIAGNOSTICS OF " + report_title
    file.write(my_title + "\n")
    file.write("=" * len(my_title) + "\n\n")
    
    for f in plot_list:
#        MD = METAdata()

        filename = f.rpartition(os.sep)[-1]
        filepath = f.rpartition(os.sep)[0]

        shutil.copy(f, src_dir)

        caption = "test"#MD.read(f).get_dict()['ESMValTool']['caption']
    
        file.write(".. figure:: " + filename + "\n" 
#        file.write(".. figure:: " + "../../../../diag_scripts/aux/C3S_511/" + f + "\n" 
                 "   :align:   center"   + "\n" 
                 "   :width:   95%" + "\n\n" 
                 "   " + caption[0] + "\n"
                  )
    file.close() 

    # Go to the directory and create PDF
    shutil.copy("doc/reporting/source/conf.py", src_dir)
    shutil.copy("doc/reporting/source/index.rst", src_dir)

    os.environ['SOURCEDIR'] = src_dir
    os.environ['BUILDDIR'] = bld_dir

    os.chdir("doc/reporting")

    os.system("make latexpdf")

    os.rename(bld_dir + os.sep + "latex" + os.sep + "ESMValToolC3S_511Report.pdf", path_out + os.sep + "report_" + report_title.split()[0].lower() + ".pdf")

    print("created reporting pdf!... or not.")

flist = ["/athome/laue_ax/sphinx/ESMValTool-private/diag_scripts/aux/C3S_511/example_images/albedo_QA4ECV_all_models_regionalized_smean_ts.png"]

do_report(flist, "mean and variability test")


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
