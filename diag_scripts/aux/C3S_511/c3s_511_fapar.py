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
        #super(FAPAR_Diagnostic, self).run_diagnostic()
        self.__my_diag__()
        
 
    
    def __my_diag__(self):

        percentiles = [25.,75.]
        
        interquart = self.sp_data.collapsed("time", iris.analysis.PERCENTILE, percent=percentiles)
        
        print interquart
        
        
        
        
        
        
        Q25_perc = interquart.extract(iris.Constraint(percentile_over_time=25.))
        Q75_perc = interquart.extract(iris.Constraint(percentile_over_time=75.))
        
        iqrange=Q25_perc #TODO: MINUS
        
        list_of_plots=[]
        filename = self.__plot_dir__ + os.sep + self.__basic_filename__ + "_IQR" + "." + self.__output_type__
        list_of_plots.append(filename)
        
        x=Plot2D(iqrange)
        fig = plt.figure()
        ax = [plt.subplot(1,1,1)]
        
        x.plot(ax=ax, title=" ".join([self.__dataset_id__[indx] for indx in [0,2,1,3]]) + " (" + self.__time_period__ + ")",vminmax=vminmax)
        fig.savefig(filename)
        plt.close(fig.number)
                        
                        
        ESMValMD("meta",
                 filename,
                 self.__basetags__ + ['DM_global', 'C3S_mean_var'],
                 str("IQR" + ' values of ' + self.__varname__ + ' for the data set "' + "_".join(self.__dataset_id__) + '" (' + self.__time_period__ + ')'),
                 '#C3S' + "igr" + self.__varname__,
                 self.__infile__,
                 self.diagname,
                 self.authors)