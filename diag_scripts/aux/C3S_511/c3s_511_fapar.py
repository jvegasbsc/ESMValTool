"""
Implementation for FAPAR diagnostics into ESMValTool
"""

import iris
from c3s_511_basic import Basic_Diagnostic
from ESMValMD import ESMValMD
import sys
import os
import matplotlib.pyplot as plt
[sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.abspath(__file__)), dir)) for dir in ["lib", "plots"]]
from plot import Plot2D

class FAPAR_Diagnostic(Basic_Diagnostic):
    """    
    class to implement fapar diagnostics, like e.g. global means,
    global differences, RMSD etc.
    """

    def run_diagnostic(self):
        """
        run parent diagnostic and the FAPAR specific diagnostic
        """
        super(FAPAR_Diagnostic, self).run_diagnostic()
        self.__interquartile_range__()
        
 
    
    def __interquartile_range__(self):
        """
        method for calculating interquartile range (iqrange) and 
        relative interquartile range (riqrange)
        """
        this_function="interquartile range"
        
        percentiles = [25.,50.,75.]
        
        # calculate quartiles
        quartiles_cube = self.sp_data.collapsed("time", iris.analysis.PERCENTILE, percent=percentiles)

        # extract quartiles
        Q25_perc = quartiles_cube.extract(iris.Constraint(percentile_over_time=25.))
        Q75_perc = quartiles_cube.extract(iris.Constraint(percentile_over_time=75.))
        
        # calculate differences
        iqrange=Q75_perc - Q25_perc
        
        # make list of plots
        list_of_plots=[]
        
        # define filename
        filename = self.__plot_dir__ + os.sep + self.__basic_filename__ + "_IQR" + "." + self.__output_type__
        
        # add current plot to list
        list_of_plots.append(filename)
        
        # produce the plot object
        x=Plot2D(iqrange)
        
        # initialize the figure
        fig = plt.figure()
        
        # intialize the axis where we plot
        ax = [plt.subplot(1,1,1)]
        
        # call the plot
        x.plot(ax=ax, title=" ".join([self.__dataset_id__[indx] for indx in [0,2,1,3]]) + " (" + self.__time_period__ + ")")
        
        # save figure
        fig.savefig(filename)
        
        # close figure (memory reasons)
        plt.close(fig.number)
                        
        # add metadata
        ESMValMD("meta",
                 filename,
                 self.__basetags__ + ['DM_global', 'C3S_mean_var'],
                 # caption below EDIT!
                 str("IQR" + ' values of ' + self.__varname__ + ' for the data set "' + "_".join(self.__dataset_id__) + '" (' + self.__time_period__ + ')'),
                 # identifier EDIT!
                 '#C3S' + "iqr" + self.__varname__,
                 self.__infile__,
                 self.diagname,
                 self.authors)
        
        # extract median
        Q50_perc = quartiles_cube.extract(iris.Constraint(percentile_over_time=50.))
        
        # calculate realtive IQR
        riqrange=iqrange/Q50_perc
        
        # define filename
        filename = self.__plot_dir__ + os.sep + self.__basic_filename__ + "_rIQR" + "." + self.__output_type__
        
        # add current plot to list
        list_of_plots.append(filename)
        
        # produce the plot object
        x=Plot2D(riqrange)
        
        # initialize the figure
        fig = plt.figure()
        
        # intialize the axis where we plot
        ax = [plt.subplot(1,1,1)]
        
        # call the plot
        x.plot(ax=ax, title=" ".join([self.__dataset_id__[indx] for indx in [0,2,1,3]]) + " (" + self.__time_period__ + ")")
        
        # save figure
        fig.savefig(filename)
        
        # close figure (memory reasons)
        plt.close(fig.number)
                        
        # add metadata               
        ESMValMD("meta",
                 filename,
                 self.__basetags__ + ['DM_global', 'C3S_mean_var'],
                 # caption below EDIT!
                 str("rIQR" + ' values of ' + self.__varname__ + ' for the data set "' + "_".join(self.__dataset_id__) + '" (' + self.__time_period__ + ')'),
                 # identifier EDIT!
                 '#C3S' + "riqr" + self.__varname__,
                 self.__infile__,
                 self.diagname,
                 self.authors)
        
        # produce report
                # produce report
        expected_input, found = \
            self.__file_anouncement__(subdir="c3s_511/ECV_specific_input",
                                      expfile="_ECV_specific.txt",
                                      protofile="empty.txt",
                                      function=this_function)
            
        if found:    
            self.__do_report__(content={"plots":list_of_plots,"freetext":expected_input}, filename=this_function.upper())
        else:
            self.__do_report__(content={"plots":list_of_plots}, filename=this_function.upper())
        