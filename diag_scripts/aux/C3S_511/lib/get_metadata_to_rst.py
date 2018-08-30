import sys
import os
import shutil
import csv
import numpy as np
from collections import OrderedDict
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import subprocess
import re

# TO DO: remove once this function is called within the ESMValTool
#        (included automatically then)
sys.path.append(os.path.abspath("./diag_scripts/lib/python"))
sys.path.append(os.path.abspath("./diag_scripts"))

from METAdata import METAdata
#import csv

def do_report(report_data, report_title, work_dir, signature="", latex_opts=False):
    """
    - report_data  a dictionary of a list of plot file names  (.png, ...) including *full path*
                   OR a dictionary containing strings
                   keys are in ["listtext","freetext","plots"]
    - report_title is a string used as title for the report 
    
    Updated: February 18th, 2018 (B. Mueller)
    Updated: July 31st, 2018 (B. Mueller)
    """

    # define output and temporary directories for Sphinx (source, build);
    # add process id to temporary directory names to allow for
    # execution of multiple instances in parallel

    pid = str(os.getpgid(0))

    path_out = work_dir + os.sep + "reporting"    # the final pdf will be put here
    src_dir = path_out + os.sep + "source_" + report_title.split()[0].lower() + "_" + signature + "_" + pid # Sphinx source code directory
    bld_dir = path_out + os.sep + "build_" + report_title.split()[0].lower() + "_" + signature + "_" + pid  # Sphinx build directory

    # create output and temporary directories
    if not (os.path.exists(path_out)):
      os.makedirs(path_out)
    if not (os.path.exists(src_dir)):
      os.makedirs(src_dir)    
    if not (os.path.exists(bld_dir)):
      os.makedirs(bld_dir)    

    # define filename of Sphinx source code file (.rst format)
    output_file = src_dir + os.sep + "report.rst"  
    outfile = open(output_file, "w") 

    # title (headline)
#    my_title = report_title
#    outfile.write(my_title + "\n")
#    outfile.write("=" * len(my_title) + "\n\n")
    
    #search for text and plots instances in report_data
    
    if isinstance(report_data, dict):
        
        if "listtext" in report_data.keys():

            # if the input data are a dictionary, we create a bullet point list
            # from all key value pairs
            
            this_title = "Short Information"
            outfile.write(this_title + "\n")
            outfile.write("=" * len(this_title) + "\n\n")

            if isinstance(report_data["listtext"], dict):
        
            # if the input data are a list of filenames (plots), we simply put all
            # figures with their corresponding caption read from the plot meta data
            # into the report
        
                for key in report_data["listtext"]:
                    if isinstance(report_data["listtext"][key], dict):
                        outfile.write("* " + key + "\n\n")
                        for key2 in report_data["listtext"][key]:
                            if isinstance(report_data["listtext"][key][key2], dict):
                                outfile.write("  * " + key2 + "\n\n")
                                for key3 in report_data["listtext"][key][key2]:
                                    outfile.write("    * " + key3 + ": " + report_data["listtext"][key][key2][key3] + "\n")
                                outfile.write("\n")
                            else:
                                outfile.write("  * " + key2 + ": " + report_data["listtext"][key][key2] + "\n")
                        outfile.write("\n")
                    else:
                        outfile.write("* " + key + ": " + report_data["listtext"][key] + "\n")
                        
                outfile.write(".. raw:: latex \n\n")
                outfile.write("   \clearpage \n") #does not react on this
                        
            else:
                print("Wrong format in text entry, nothing can be written!") # TODO: ERROR function 
                        
        else:
            print("No writable list found! There was no 'listtext' in the dictionary!")
            
        if "freetext" in report_data.keys():
            
            this_title = "Description"
                        
            if isinstance(report_data["freetext"], (str,unicode)):
                
                if os.path.isfile(report_data["freetext"]):
                    with open(report_data["freetext"],"r") as freetext:
                        with open("./diag_scripts/aux/C3S_511/lib/predef/empty.txt","r") as empty:
                            text = freetext.read()
                            if len(set(text) - set(empty.read())):
                                outfile.write(this_title + "\n")
                                outfile.write("=" * len(this_title) + "\n\n")
                                outfile.write(text)
                                outfile.write("\n\n")
                            else: 
                                print("There is still the empty description from empty.txt!")
                else:
                    outfile.write(this_title + "\n")
                    outfile.write("=" * len(this_title) + "\n\n")
                    outfile.write(report_data["freetext"] + "\n\n")
                
            else:
                print("Wrong format in text entry, nothing can be written!") # TODO: ERROR function 
                
            outfile.write(".. raw:: latex \n\n")
            outfile.write("   \clearpage \n") #does not react on this
                
        else:
            print("No writable text found! There was no 'freetext' in the dictionary!")
            
        if "plots" in report_data.keys():
            
            this_title = "Figure(s)"
            outfile.write(this_title + "\n")
            outfile.write("=" * len(this_title) + "\n\n")

            if isinstance(report_data["plots"], list):
                MD = METAdata()
        
                # add all plots in plot_list to Sphinx source code file;
                # the figure captions are extracted from the exif file headers
                # of the .png files (if present)
        
                for f in report_data["plots"]:
        
                    filename = f.rpartition(os.sep)[-1]
                    filepath = f.rpartition(os.sep)[0]
        
                    shutil.copy(f, src_dir)
        
                    try:
                        caption = MD.read(f).get_dict()['ESMValTool']['caption']
                    except:
                        caption = ["Error: empty caption in " + filename]
            
                    outfile.write(".. figure:: " + filename + "\n" 
                                  "   :align:   center" + "\n" 
                                  "   :width:   95%" + "\n\n" 
                                  "   " + caption[0] + "\n"
                                 )
                    outfile.write(".. raw:: latex \n\n")
                    outfile.write("   \FloatBarrier \n") #does not react on this
        
            else:
                print("Wrong format in plots entry, nothing can be written!") # TODO: ERROR function 
                
        else:
            print("No plottable links found! There was no 'plots' in the dictionary!")
            
        if "freetext" not in report_data.keys() and "listtext" not in report_data.keys() and "plots" not in report_data.keys():
            print("Nothing to write reports from!! There was no 'plots' nor 'text' in the dictionary!") # TODO: ERROR function 
            return

    outfile.close() 

    # copy Sphinx configuration and index file to temporary source directory
    shutil.copy("doc/reporting/source/conf.py", src_dir)
    shutil.copy("doc/reporting/source/index.rst", src_dir)

    # set environment variables for Sphinx (source directory and build directory)
    os.environ['SOURCEDIR'] = src_dir
    os.environ['BUILDDIR'] = bld_dir

    
    if latex_opts is not None:
        # run Sphinx to create a pdf
        oldpath = os.getcwd()
        os.chdir("doc/reporting")
        if not latex_opts:
            with open(os.devnull, 'wb') as devnull:
                subprocess.call("make latexpdf", shell=True,
                                stdout=devnull, stderr=subprocess.STDOUT)
        else:
            subprocess.call("make latexpdf", shell=True)
         
        os.chdir(oldpath)

        # move pdf to the output directory and rename to report_xxx.pdf
        pdfname = path_out + os.sep + "report_" + report_title.split()[0].lower() + "_" + signature + ".pdf"
        os.rename(bld_dir + os.sep + "latex" + os.sep + "ESMValToolC3S_511Report.pdf", pdfname)

        # clean up temporary directories
        if os.path.exists(src_dir):
             # remove if exists
             shutil.rmtree(src_dir)
             pass
        if os.path.exists(bld_dir):
             # remove if exists
             shutil.rmtree(bld_dir)
             pass
    
        print("Successfully created " + pdfname + "!")
    
    
def do_smm_table(csv_expert, csv_definitions):
    """
    Authors:  F. Massonnet and B. Hassler
    Creation: February 8th, 2018
    Updated: February 15th, 2018 (B. Mueller)

    Inputs:
        - csv_expert : A CSV file giving the grades assigned to each
                       cell of the System Maturity Matrix (1 to 6)
                       The first line is supposed to be a header with
                       the names of the different categories. This row is not read,
                       but it useful to keep it if the CSV themselves have to be opened
                       by the expert.
                       Fields in the CSV file must be separated by commas. 
                       White spaces can be present around the commas. 
                       The strings should be quoted (" ... "), especially if there are 
                       commas within them.
                       Non relevant information (i.e., where no grade is given) shall be
                       reported by not filling the field.
                       The CSV file should not have an empty line at the end!
                       Example:
                       >> cat input_file.csv
                             "Software Readiness" , "Metadata", "User Documentation", "Uncertainty Characterisation", "Public access feedback and update","Usage"
                              1                   , 4         , 3                   , 2                             , 1                                 , 6
                              3                   , 2         , 1                   , 4                             , 6                                 , 5
                              1                   , 3         , 1                   , 2                             , 4                                 ,
                              6                   ,           , 3                   , 2                             , 6                                 ,

        - csv_definitions : A CSV file (independent of the product to be assessed) giving
                            the standard names for criteria in the System Maturity Matrix.
                            Same conventions as for the csv_expert file apply.
                            The sign "\n" can be used if a linebreak should appear in the 
                            name. 
                          
                            Example: 
                            >> cat definition_file.csv
                               "Software\nReadiness", "Metadata",  "User\nDocumentation", "Uncertainty\nCharacterisation", "Public access,\nfeedback,\nand update", "Usage"
                               "Coding\nStandards"  , "Standards", "Formal description\nof scientific\nmethodology", "Standards", "Public\nAccess/Archive", "Research"
                               "Software\nDocumentation", "Collection\nlevel", "Formal validation\nreport", "Validation", "Version", "Decision\nsupport\nsystem"
                               "Numerical\nReproducibility\nand portability", "File level", "Formal product\nuser guide", "Uncertainty\nquantification", "User\nfeedback",
                               "Security", , "Formal description\nof operations\nconcept", "Automated quality\nmonitoring", "Updates to record" ,


    Output: a .rst file including the System Maturity Matrix
 
    """

#    # WILL HAVE TO BE DELETED
#    # -----------------------
#    path_out = work_dir + "/plot_scripts/python/system_maturity_matrix"
#
#    path_report_out = work_dir + "/reports/sphinx/source"    

    # ----------------------
 
    max_grade = 6 # Maximum possible grade. Grades are integers to be given between 
                  # 1 and max_grade (1, ..., max_grade)
 
    # Create a list. Each item of the list will be itself a list of strings, corresponding to each
    # word to appear in the System Maturity Matrix (header and words in the matrix itsef)
    definitions = list()
    with open(csv_definitions, 'rb') as csvfile:
        s = csv.reader(csvfile, delimiter = ",", skipinitialspace = True)
        for row in s:
            definitions.append(row)

    # Get the dimensions of the System Maturity Matrix, including header.
    # The x-dimension goes horizontally (along columns) while the y-dimension goes
    # vertically (along rows)
    
    ny = len(definitions)
    nx = len(definitions[0])

    # Check if, possibly, one of the rows of the CSV has not the same number of items
    # TO BE IMPROVED WITH NEW BACK END
    if sum([1.0 * (len(definitions[i])==nx) for i in range(len(definitions))] ) != ny:
        sys.exit("(do_smm_report) STOP: uneven number of columns in definition file")
       

    # The grades to be used as color in the System maturity matrix
    grades = np.empty((ny, nx)) 
    grades[:] = np.nan
    
    with open(csv_expert, 'rb') as csvfile:
        s = csv.reader(csvfile, delimiter = ",", skipinitialspace = True)
        counter_y = 0 # counter to iterate through rows
        for row in s:
            # Check if row length matches definition file
            if counter_y > 0: # Ignore header
                counter_x = 0
                for item in row:
                    if item.strip() != '':
                        # TO BE IMPROVED WITH NEW BACK END
                        try:
                            grades[counter_y, counter_x] = int(item)
                        except IndexError:
                            sys.exit("(do_smm_table) number of columns of input " + 
                                     " csv file for system maturity matrix differs" +
                                     " from definition file")
                    counter_x += 1
            counter_y += 1
    
    # Create the figure
    fig = plt.figure(figsize = (9, 5))
    # X-Y mesh to plot the color array
    # The Y dimension is written from ny to zero as to write items in visually
    # descending order.
    x_grid, y_grid = np.meshgrid(np.arange(nx + 1), np.arange(ny, -1, -1))

    # The color map
    cmap = plt.get_cmap('RdYlGn', max_grade + 2)

    # Plot the cells in color
    plt.pcolor(x_grid, y_grid, grades, cmap = cmap, vmin = -0.5, vmax = max_grade + 1.5)
    # Create colorbar at the bottom
    cb = plt.colorbar(boundaries = np.arange(0.5, max_grade + 1), 
                      ticks = np.arange(1, nx + 1), orientation = "horizontal", pad = 0.05)
    # Remove ticks
    cb.set_ticks([])
    # Put colorbar as the most bottom layer
    cb.ax.zorder = -1
 
    # Write legend inside colorbar
    [plt.text(0.5 * nx / max_grade + 1.0 * (i - 1) * nx / max_grade, -0.65, str(i), fontsize = 14, 
     fontweight = "bold", ha = "center", va = "center") for i in np.arange(1, max_grade + 1)]


    # Finish polishing the figure
    plt.title("System Maturity Matrix", fontsize = 18)
    # Grid lines
    [plt.plot((0, nx), (y, y), color = 'k') for y in range(ny)]
    [plt.plot((x, x), (0, ny), color = 'k') for x in range(nx)]

    for y in range(ny):
        if y == ny - 1: # if header
          fontweight = "bold"
          fontsize   = 10
        else:
          fontweight = "normal"
          fontsize   = 10
        for x in range(nx):
            # Read in the "go to line" in the csv and convert it to "go to line" instruction
            # When \n stands in a CSV, python reads \\n
            instring = "\n".join(definitions[ny - y - 1][x].split("\\n")).strip()
            plt.text(x + 0.5, y + 0.5, instring, ha = 'center', va = 'center', 
                     fontweight = fontweight, fontsize = fontsize)

    plt.xticks([])
    plt.yticks([])
    plt.xlim(0, nx)
    plt.ylim(0, ny)
    plt.tight_layout()
    
    return fig

#    # Create output path for figure
#    if not (os.path.exists(path_out)):
#      os.makedirs(path_out)
#    plt.savefig(path_out + "/" + "system_maturity_matrix.png", dpi = 400)
#    plt.close("fig")
#    # Create *.rst report
#    if not (os.path.exists(path_report_out)):
#      os.makedirs(path_report_out)
#
#    output_file = "report_smm.rst"
#    file = open(path_report_out + "/" + output_file, "w")
#
#    my_title ="SYSTEM MATURITY MATRIX"
#    file.write(my_title + "\n")
#    file.write("=" * len(my_title) + "\n\n")
#
#    file.write(".. figure:: " + "../../../../plot_scripts/python/system_maturity_matrix/system_maturity_matrix.png" + "\n"
#               "   :align:   center"   + "\n"
#               "   :width:   95%" + "\n"
#              )
#    file.close()    


def do_eval_table(varname, eval_expert, eval_data, ecv_length):
    """
    Author  : Arun Rana
    Creation: July 26th, 2018
    Updated : August 08th, 2018
    
    Inputs:
        - varname: the variable name to read the right evaluation/validation table parts
        
        - eval_expert: A CSV file giving the Evaluation/validation criteria value for the dataset
                       or product assessed.
                       The first line is supposed to be a header with
                       the names of the different criteria. This provides
                       the criteria for.
                       Fields in the CSV file must be separated by commas. 
                       White spaces can be present around the commas. 
                       The strings should be quoted (" ... "), especially if there are 
                       commas within them.
                       Non relevant information (i.e., where no validation is given) shall be
                       reported by not filling the field.
                       The sign "\n" can be used if a linebreak should appear in the 
                       name.
                       The CSV file should not have an empty line at the end!

                       Example:
                       >> cat eval_expert.csv
                            "ECV", "Mean\n climatology", "Trend", "Variability"
                             tas,	20,	42,	0.2

        - eval_data: A dictionary of values (read from in the product). If the expert input
			  is smaller than the  length of the dataset for ECV, that cell would get the number for "green".
			  If the expert input is bigger than the  length of the dataset for ECV,
			  that cell would get the number for "red".
                          It is a constant value refering to the total length of data that is read in and is pointed in overview report
			  and saved for GCOS requirement as well. (tim_range_spec)
    Outputs:
        - A *.png file that includes a one-row table with green or red or white shading depending
          on whether the evaluation standards are met or not

    TODO:
       
    """

#    # WILL HAVE TO BE DELETED
#    # -----------------------
#    path_out = work_dir + "/plot_scripts/python/system_maturity_matrix"
#
#    path_report_out = work_dir + "/reports/sphinx/source"    

    # ----------------------
 
    max_grade_eval = 3 # Maximum possible grade. Grades are integers to be given between 
                  # 1 and max_grade_eval (1, ..., max_grade_eval)
 
    # Create a list. Items of the list will be itself a list of strings, corresponding to each
    # word to appear in the Evaluation Matrix (header of the matrix itself - parameter for evaluation)
    definitions = list()
    with open(eval_data, 'rb') as csvfile:
        s = csv.reader(csvfile, delimiter = ",", skipinitialspace = True)
        for row in s:
            definitions.append(row)
    
    # Get the dimensions of the Validation Matrix, which is just header in this case.
    # The x-dimension goes horizontally (along columns) while the y-dimension goes
    # vertically (along rows)
    ny = len(definitions)
    nx = len(definitions[0])
   
    # The grades to be used as color in the Validation matrix
    grades = np.empty((ny, nx)) 
    grades[:] = np.nan
       
    # perform the comparison with validation data and prepare the grading matrix for plotting
    # here the grades matrix is to be filled with number assigned with colors below it give it right feel
    
    time_range = list()    
    with open(eval_expert, 'rb') as csvfile:
	s = csv.reader(csvfile, delimiter = ",", skipinitialspace = True)
        for row in s:
            time_range.append(row)
    
    time_range = time_range[1][:]  #excluding the header for comparison with range (was [1][1:])  
        
    # making sure that the conversion to integers does not crash the namelist
    try:
        time_range = [int(i) for i in time_range] #conversion to integer for comparison to ECV data length
    except:
        time_range = [np.nan for i in time_range]
    #print time_range - done during coding
    
    for i in range(len(time_range)):
        if time_range[i] <= ecv_length:    # it's ok if the data are just as long as the limit given by the expert user
           grades[0][i] = 3
        elif time_range[i] > ecv_length:
           grades[0][i] = 1
        else:
           grades[0][i] = 2
        #print grades.item (i)
     
    # Create the figure
    fig = plt.figure(figsize = (6, 2))
    # X-Y mesh to plot the color array
    # The Y dimension is written from ny to zero as to write items in visually
    # descending order.
    x_grid, y_grid = np.meshgrid(np.arange(nx + 1), np.arange(ny, -1, -1))

    # The color map
    cmap = plt.get_cmap('RdYlGn', max_grade_eval )

    # Plot the cells in color
    plt.pcolor(x_grid, y_grid, grades, cmap = cmap, vmin = -0.5, vmax = max_grade_eval + 1.5)
    

    # Finish polishing the figure
    plt.title("ESM Evaluation Table and Grading", fontsize = 18)
    # Grid lines
    [plt.plot((0, nx), (y, y), color = 'k') for y in range(ny)]
    [plt.plot((x, x), (0, ny), color = 'k') for x in range(nx)]

    for y in range(ny):
        if y == ny - 1: # if header
          fontweight = "bold"
          fontsize   = 15
        else:
          fontweight = "normal"
          fontsize   = 15
        for x in range(nx):
            # Read in the "go to line" in the csv and convert it to "go to line" instruction
            # When \n stands in a CSV, python reads \\n
            instring = "\n".join(definitions[ny - y - 1][x].split("\\n")).strip()
            plt.text(x + 0.5, y + 0.5, instring, ha = 'center', va = 'center', 
                     fontweight = fontweight, fontsize = fontsize)

    plt.xticks([])
    plt.yticks([])
    plt.xlim(0, nx)
    plt.ylim(0, ny)
    plt.tight_layout()
    
    return fig


def do_gcos_table(varname, gcos_expert, gcos_reference):
    """
    Author  :  F. Massonnet
    Creation: February 9th, 2018
    Updated: February 15th, 2018 (B. Mueller)
             March 1st, 2018 (B. Mueller) (no csv for expert)

    Inputs:
        - varname: the variable name to read the right gcos table parts
        
        - gcos_expert: A CSV file giving the GCOS criteria value for the dataset
                       or product assessed.
                       The first line is supposed to be a header with
                       the names of the different criteria. This row is not read,
                       but it useful to keep it if the CSV themselves have to be opened
                       by the expert.
                       Fields in the CSV file must be separated by commas. 
                       White spaces can be present around the commas. 
                       The strings should be quoted (" ... "), especially if there are 
                       commas within them.
                       Non relevant information (i.e., where no grade is given) shall be
                       reported by not filling the field.
                       The sign "\n" can be used if a linebreak should appear in the 
                       name.
                       The CSV file should not have an empty line at the end!

                       Example:
                       >> cat gcos_expert.csv
                            "Accuracy", "Stability", "Frequency",         "Resolution"
                            0.3       , 0.02       , "Hourly\nto weekly", 0.9

        - gcos_reference: A dictionary the reference values (GCOS standards for that product)
                          keys must of course match the expert file!
                          Example:
                          {"Accuracy":{"value":0.3,"unit":"m2"},
                           "Stability":{"value":0.5,"unit":"m2 0.1 y-1"},
                           "Frequency":{"value":28,"unit":"days"},
                           "Resolution":{"value":0.5,"unit":"degrees"}}
    Outputs:
        - A *.png file that includes a two-row table: first row = header,
          second row = numbers from expert with green or red shading depending
          on whether the GCOS standards are met or not

    TODO:
        - is a GCOS standard met if the expert value is less in an absolute sense?
        - how about frequency? reported as "hourly", "decadal", not numbers
    """

#    # WILL HAVE TO BE DELETED
#    # -----------------------
#    path_out = work_dir + "/plot_scripts/python/gcos"
#
#    path_report_out = work_dir + "/reports/sphinx/source"

    # ----------------------

    # Create a list. Each item of the list will be itself a list of strings, corresponding either to the 
    # headers or to the GCOS reference values
    
    print gcos_expert
    
    contents = list()
    with open(gcos_reference, 'rb') as csvfile:
        s = csv.reader(csvfile, delimiter = ";", skipinitialspace = True)
        for row in s:
            contents.append(row)

    variable = varname
    #variable = 'Water vapour'
    
    # convert the list into an array for easier checks for entries
    gcos_ref_array = np.asarray(contents)
	
    gcos_idx = []
    for i in range(len(gcos_ref_array[:,2])):        # look into the third column of the table to find the CMOR ECV name
        if re.search(r'\b' + variable.strip().upper() + r'\b', gcos_ref_array[i,2].upper()):
            gcos_idx.append(i)

    gcos_contents = OrderedDict()
    gcos_contents["Frequency"] = [""]
    gcos_contents["Resolution"] = [""]
    gcos_contents["Accuracy"] = [""]
    gcos_contents["Stability"] = [""]
    
    if len(gcos_idx) >= 1:
        gcos_contents["Frequency"] = [gcos_ref_array[gcos_idx[0],3]]
        gcos_contents["Resolution"] = [gcos_ref_array[gcos_idx[0],4]]
        gcos_contents["Accuracy"] = [gcos_ref_array[gcos_idx[0],5]]
        gcos_contents["Stability"] = [gcos_ref_array[gcos_idx[0],6]]

				
    # Get the horizontal dimension of the GCOS table
    nx = len(gcos_contents)
    
#    # Read in the product table and store the data. ignore the header!
#    counter_y = 0 # Counter along rows of the CSV file
#    with open(gcos_expert, 'rb') as csvfile:
#        # Check the match up between headers
#        s = csv.reader(csvfile, delimiter = ",", skipinitialspace = True)
#        for row in s:
#            if counter_y == 0:
#                if not [l.strip().lower() for l in row] == [l.strip().lower() for l in contents[0]]:
#                    sys.exit("(do_gcos_report) NO MATCH-UP BETWEEN HEADER NAMES IN REFERENCE AND PRODUCT CSV FILES. MAKE SURE NAMES AND ORDER MATCH UP.")
#            else:
#                contents.append(row)
#            counter_y += 1
     
    if not isinstance(gcos_expert, dict):
        assert False, "wrong input type in gcos"
    is_equal = np.array_equal(np.sort(np.array(gcos_expert.keys())), np.sort(np.array(gcos_contents.keys())))
    if not is_equal:
        assert False, "wrong names in gcos"
    else:
        #data_contents=[]
        for key in gcos_contents:
            if gcos_expert[key]["unit"] is None:
                this_unit = ""
            elif gcos_expert[key]["unit"] in ["1","unkown","no-unit"]:
                this_unit = ""
            else:
                this_unit = " " + gcos_expert[key]["unit"]
            gcos_contents[key].append(str(gcos_expert[key]["value"]) + this_unit)
            ny = len(gcos_contents[key]) + 1
        #gcos_contents.append(data_contents)
#    print(gcos_contents)
    

    # Check if, possibly, one of the rows of the CSV has not the same number of items
    # TO BE IMPROVED WITH NEW BACK END
    #if sum([1.0 * (len(gcos_contents[i])==nx) for i in range(len(gcos_contents))] ) != ny:
    #    sys.exit("(do_gcos_report) STOP: uneven number of columns in reference file")

    # Create the figure
    fig = plt.figure(figsize = (12, 7))
    # X-Y mesh to plot the color array
    # The Y dimension is written from ny to zero as to write items in visually
    # descending order.
    x_grid, y_grid = np.meshgrid(np.arange(nx + 1), np.arange(ny, -1, -1))

    plt.title("GCOS requirements", fontsize = 36)
    # Grid lines
    [plt.plot((0, nx), (y, y), color = 'k') for y in range(ny)]
    [plt.plot((x, x), (0, ny), color = 'k') for x in range(nx)]

    
    for y in range(ny - 1):
        for (x, key) in enumerate(gcos_contents.keys()):
            #print(gcos_contents[key][y])
            instring = "\n".join(gcos_contents[key][y].split("\\n")).strip()
            #print(instring)
            plt.text(x + 0.5, ny - y - 1.5, instring, ha = 'center', va = 'center',
                     fontweight = "normal", fontsize = 20)
    for (x, key) in enumerate(gcos_contents.keys()):
        instring = "\n".join(key.split("\\n")).strip()
        plt.text(x + 0.5, ny - 0.5, instring, ha = 'center', va = 'center',
                 fontweight = "bold", fontsize = 20)
      
    # Add legend on the left
    plt.text(-0.3, ny - 1.5, "GCOS", rotation = 90, ha = "center", va = "center", fontsize=15)
    plt.text(-0.3, ny - 2.5, "ECV\n(averages)", rotation = 90, ha = "center", va = "center", fontsize=15)
    plt.xticks([])
    plt.yticks([])
    plt.xlim(0, nx)
    plt.ylim(0, ny)
    plt.tight_layout()

    return fig

#    # Create output path for figure
#    if not (os.path.exists(path_out)):
#      os.makedirs(path_out)
#    plt.savefig(path_out + "/" + "gcos_requirements.png", dpi = 400)
#    plt.close(fig)
#
#
#    # Create *.rst report
#    if not (os.path.exists(path_report_out)):
#      os.makedirs(path_report_out)
#
#    output_file = "report_gcos.rst"
#    file = open(path_report_out + "/" + output_file, "w")
#
#    my_title ="GCOS"
#    file.write(my_title + "\n")
#    file.write("=" * len(my_title) + "\n\n")
#
#    file.write(".. figure:: " + "../../../../plot_scripts/python/gcos/gcos_requirements.png" + "\n"
#               "   :align:   center"   + "\n"
#               "   :width:   95%" + "\n"
#              )
#    file.close()


# dummy code to test do_report

#flist = ["diag_scripts/aux/C3S_511/example_images/albedo_QA4ECV_all_models_regionalized_smean_ts.png"]
#testdict = {"resolution" : {"dx":"1 deg", "dy":"1 deg"}, "time":{"period":{"start":"2000-01-01", "end":"2000-12-31"}}, "attribute1":"none"}
#
#do_report(testdict, "dictionary test")
#do_report(flist, "mean and variability test")


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
