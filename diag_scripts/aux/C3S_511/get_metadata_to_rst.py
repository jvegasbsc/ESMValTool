import sys  # WBD
import os   #WBD
import matplotlib.pyplot as plt
sys.path.append("../../../diag_scripts/lib/python/")   #WBD
from METAdata import METAdata
import csv
import numpy as np

path_out = "../../../work/reports/sphinx/source"

def wrapper(string, splitter):
    subs = string.split(splitter)
    subs = [i + splitter for i in subs[:-1]] + [subs[-1]]
    if len(subs) > 1:
        return subs
    else:
      return subs

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
        MD = METAdata()
        
        caption = MD.read(f).get_dict()['ESMValTool']['caption']
    
        file.write(".. figure:: " + "../../../../diag_scripts/aux/C3S_511/" + f + "\n" 
                 "   :align:   center"   + "\n" 
                 "   :width:   95%" + "\n" 
                 "   " + caption[0] + "\n"
                  )
    file.close() 

    # Go to the directory and create PDF

def do_smm_report(csv_expert, csv_definitions):
    """
    Authors:  F. Massonnet and B. Hassler
    Creation: February 8th, 2018

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

    # WILL HAVE TO BE DELETED
    # -----------------------
    path_out = "../../../plot_scripts/python/system_maturity_matrix"

    path_report_out = "../../../work/reports/sphinx/source"    

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
                            sys.exit("(do_smm_report) number of columns of input " + 
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

    # Create output path for figure
    if not (os.path.exists(path_out)):
      os.makedirs(path_out)
    plt.savefig(path_out + "/" + "system_maturity_matrix.png", dpi = 400)
    plt.close("fig")
    # Create *.rst report
    if not (os.path.exists(path_report_out)):
      os.makedirs(path_report_out)

    output_file = "report_smm.rst"
    file = open(path_report_out + "/" + output_file, "w")

    my_title ="SYSTEM MATURITY MATRIX"
    file.write(my_title + "\n")
    file.write("=" * len(my_title) + "\n\n")

    file.write(".. figure:: " + "../../../../plot_scripts/python/system_maturity_matrix/system_maturity_matrix.png" + "\n"
               "   :align:   center"   + "\n"
               "   :width:   95%" + "\n"
              )
    file.close()    

def do_gcos_report(gcos_expert, gcos_reference):
    """
    Author  :  F. Massonnet
    Creation: February 9th, 2018

    Inputs:
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

        - gcos_reference: A CSV file with the reference values (GCOS standards for that product)
                          Order of columns must of course match the expert file!
                          Example:
                          >> cat gcos_reference.csv
                               "Accuracy", "Stability", "Frequency", "Resolution"
                               0.05       , 0.05       , 0.2        , 0.5
    Outputs:
        - A *.rst file that includes a two-row table: first row = header,
          second row = numbers from expert with green or red shading depending
          on whether the GCOS standards are met or not

    TO DO:
        - is a GCOS standard met if the expert value is less in an absolute sense?
        - how about frequency? reported as "hourly", "decadal", not numbers
    """

    # WILL HAVE TO BE DELETED
    # -----------------------
    path_out = "../../../plot_scripts/python/gcos"

    path_report_out = "../../../work/reports/sphinx/source"

    # ----------------------

    # Create a list. Each item of the list will be itself a list of strings, corresponding either to the 
    # headers or to the GCOS reference values
    contents = list()
    with open(gcos_reference, 'rb') as csvfile:
        s = csv.reader(csvfile, delimiter = ",", skipinitialspace = True)
        for row in s:
            contents.append(row)

    # Get the horizontal dimension of the GCOS table
    nx = len(contents[0])
    
    # Read in the product table and store the data. ignore the header!
    counter_y = 0 # Counter along rows of the CSV file
    with open(gcos_expert, 'rb') as csvfile:
        # Check the match up between headers
        s = csv.reader(csvfile, delimiter = ",", skipinitialspace = True)
        for row in s:
            if counter_y == 0:
                if not [l.strip().lower() for l in row] == [l.strip().lower() for l in contents[0]]:
                    sys.exit("(do_gcos_report) NO MATCH-UP BETWEEN HEADER NAMES IN REFERENCE AND PRODUCT CSV FILES. MAKE SURE NAMES AND ORDER MATCH UP.")
            else:
                contents.append(row)
            counter_y += 1

    ny = len(contents)


    # Check if, possibly, one of the rows of the CSV has not the same number of items
    # TO BE IMPROVED WITH NEW BACK END
    if sum([1.0 * (len(contents[i])==nx) for i in range(len(contents))] ) != ny:
        sys.exit("(do_gcos_report) STOP: uneven number of columns in reference file")

    # Create the figure
    fig = plt.figure(figsize = (5, 2))
    # X-Y mesh to plot the color array
    # The Y dimension is written from ny to zero as to write items in visually
    # descending order.
    x_grid, y_grid = np.meshgrid(np.arange(nx + 1), np.arange(ny, -1, -1))

    plt.title("GCOS requirements", fontsize = 18)
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
            instring = "\n".join(contents[ny - y - 1][x].split("\\n")).strip()
            plt.text(x + 0.5, y + 0.5, instring, ha = 'center', va = 'center',
                     fontweight = fontweight, fontsize = fontsize)

    # Add legend on the left
    plt.text(-0.1, ny - 1.5, "GCOS", rotation = 90, ha = "center", va = "center")
    plt.text(-0.1, ny - 2.5, "ECV", rotation = 90, ha = "center", va = "center")
    plt.xticks([])
    plt.yticks([])
    plt.xlim(0, nx)
    plt.ylim(0, ny)
    plt.tight_layout()

    # Create output path for figure
    if not (os.path.exists(path_out)):
      os.makedirs(path_out)
    plt.savefig(path_out + "/" + "gcos_requirements.png", dpi = 400)
    plt.close(fig)


    # Create *.rst report
    if not (os.path.exists(path_report_out)):
      os.makedirs(path_report_out)

    output_file = "report_gcos.rst"
    file = open(path_report_out + "/" + output_file, "w")

    my_title ="GCOS"
    file.write(my_title + "\n")
    file.write("=" * len(my_title) + "\n\n")

    file.write(".. figure:: " + "../../../../plot_scripts/python/gcos/gcos_requirements.png" + "\n"
               "   :align:   center"   + "\n"
               "   :width:   95%" + "\n"
              )
    file.close()


# ===========================
# END OF FUNCTIONS
# ===========================

# SANDBOX HERE

do_gcos_report("./example_csvs/example_gcos_expert.csv", "./example_csvs/example_gcos_reference.csv")

plot_list = ["../../../diag_scripts/aux/C3S_511/example_images/albedo_QA4ECV_CMIP5_MRI-ESM1_historical_r1i1p1_regionalized_smean_ts.png" ,
             "../../../diag_scripts/aux/C3S_511/example_images/albedo_QA4ECV_CMIP5_MRI-ESM1_historical_r1i1p1_4plots_gmd.png"]

do_report(plot_list, "MEAN AND VARIABILITY")

do_smm_report("example_csvs/example_smm_expert.csv", "example_csvs/smm_definitions.csv")
