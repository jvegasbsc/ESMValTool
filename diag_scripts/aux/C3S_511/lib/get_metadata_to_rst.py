import sys
import os
import shutil

# TO DO: remove once this function is called within the ESMValTool
#        (included automatically then)
sys.path.append(os.path.abspath("./diag_scripts/lib/python"))
sys.path.append(os.path.abspath("./diag_scripts"))

from METAdata import METAdata
#import csv

# TO DO: get work directory from ESMValTool namelist
work_dir = "/athome/laue_ax/sphinx/out"

def do_report(plot_list, report_title):
    """
    - plot_list    is a list of plot file names  (.png, ...) including *full path*
    - report_title is a string used as title for the report 
    """

    # define output and temporary directories for Sphinx (source, build);
    # add process id to temporary directory names to allow for
    # execution of multiple instances in parallel

    pid = str(os.getpgid(0))

    path_out = work_dir + os.sep + "reporting"    # the final pdf will be put here
    src_dir = path_out + os.sep + "source_" + pid # Sphinx source code directory
    bld_dir = path_out + os.sep + "build_" + pid  # Sphinx build directory

    # create output and temporary directories
    if not (os.path.exists(path_out)):
      os.makedirs(path_out)
    if not (os.path.exists(src_dir)):
      os.makedirs(src_dir)    
    if not (os.path.exists(bld_dir)):
      os.makedirs(bld_dir)    

    # define filename of Sphinx source code file (.rst format)
    output_file = src_dir + os.sep + "report.rst"  
    file = open(output_file, "w") 

    # title (headline)
    my_title = report_title
    file.write(my_title + "\n")
    file.write("=" * len(my_title) + "\n\n")
    
    MD = METAdata()

    # add all plots in plot_list to Sphinx source code file;
    # the figure captions are extracted from the exif file headers
    # of the .png files (if present)

    for f in plot_list:

        filename = f.rpartition(os.sep)[-1]
        filepath = f.rpartition(os.sep)[0]

        shutil.copy(f, src_dir)

        try:
            caption = MD.read(f).get_dict()['ESMValTool']['caption']
        except:
            caption = ""
    
        file.write(".. figure:: " + filename + "\n" 
                 "   :align:   center" + "\n" 
                 "   :width:   95%" + "\n\n" 
                 "   " + caption[0] + "\n"
                  )
    file.close() 

    # copy Sphinx configuration and index file to temporary source directory
    shutil.copy("doc/reporting/source/conf.py", src_dir)
    shutil.copy("doc/reporting/source/index.rst", src_dir)

    # set environment variables for Sphinx (source directory and build directory)
    os.environ['SOURCEDIR'] = src_dir
    os.environ['BUILDDIR'] = bld_dir

    # run Sphinx to create a pdf
    os.chdir("doc/reporting")
    os.system("make latexpdf")

    # move pdf to the output directory and rename to report_xxx.pdf
    pdfname = path_out + os.sep + "report_" + report_title.split()[0].lower() + ".pdf"
    os.rename(bld_dir + os.sep + "latex" + os.sep + "ESMValToolC3S_511Report.pdf", pdfname)

    # clean up temporary directories
    if os.path.exists(src_dir):
         # remove if exists
         shutil.rmtree(src_dir)
    if os.path.exists(bld_dir):
         # remove if exists
         shutil.rmtree(bld_dir)

    print("created " + pdfname + "!")


# dummy code to test do_report

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
