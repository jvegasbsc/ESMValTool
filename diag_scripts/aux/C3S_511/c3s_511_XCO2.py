"""
Implementation for XCO2 diagnostics into ESMValTool - with changes that needs to be included when you add your diagnostics (import iris, for plotting stuff etc.)
"""
import iris
#these are for plots
from c3s_511_basic import Basic_Diagnostic
from ESMValMD import ESMValMD
import sys
import os
import matplotlib.pyplot as plt
[sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.abspath(__file__)), dir)) for dir in ["lib", "plots"]]


from plot import Plot2D
from c3s_511_basic import Basic_Diagnostic

class XCO2_Diagnostic(Basic_Diagnostic):
    """    
    class to implement XCO2 diagnostics, like e.g. global means,
    global differences, RMSD etc.
    """


#If you want to add your own diagnostic then you add that here below as it is done in basic file in run_diagnostic and you define you function here

    def run_diagnostic(self):
#        self.sp_data = self.__spatiotemp_subsets__()["Germany_2001-2005"]
#        self.__do_overview__()
#        self.__do_mean_var__()
#        self.__do_trends__()
#        self.__do_extremes__()
#        self.__do_sectors__()
#        self.__do_maturity_matrix__()
#        self.__do_gcos_requirements__()
#        self.__do_esmvalidation__()
	super(XCO2_Diagnostic, self).run_diagnostic()
	self._my_diag_()   #these diagnostics name should be different for all the ones you want to run esle if it is same then it would only perform the last one below
	self._interquartile_()

# Example below of what you can add
    def _my_diag_(self):
	print("I did my diag!")

# Lets calculate 25th and 75 th percentile and calculate interquartile difference

    def _interquartile_(self):
	percentiles=[25., 75.]
	interquart=self.sp_data.collapsed("time", iris.analysis.PERCENTILE, percent=percentiles)
	#print interquart

	#calculating the interquartile range here
	Q25_perc=interquart.extract(iris.Constraint(percentile_over_time=25.))
	Q75_perc=interquart.extract(iris.Constraint(percentile_over_time=75.))
	iqrange=Q25_perc-Q75_perc
	#print iqrange
	
        list_of_plots=[]
        # define filename
        filename = self.__plot_dir__ + os.sep + self.__basic_filename__ + "_IQR" + "." + self.__output_type__
        list_of_plots.append(filename)
        
	#produce the plot object
        x=Plot2D(iqrange)
        # initialize the fuigure
        fig = plt.figure()
	#intialize the axis wehere we plot
        ax = [plt.subplot(1,1,1)]
        
        # call the plot
        x.plot(ax=ax, title=" ".join([self.__dataset_id__[indx] for indx in [0,2,1,3]]) + " (" + self.__time_period__ + ")")
        # save figure
        fig.savefig(filename)
        # close figure
        plt.close(fig.number)
                        
                        
        ESMValMD("meta",
                 filename,
                 self.__basetags__ + ['DM_global', 'C3S_mean_var'],
                 # caption below EDIT!
                 str("IQR" + ' values of ' + self.__varname__ + ' for the data set "' + "_".join(self.__dataset_id__) + '" (' + self.__time_period__ + ')'),
                 # identifier EDIT!
                 '#C3S' + "igr" + self.__varname__,
                 self.__infile__,
                 self.diagname,
                 self.authors)
	#produce report
        self.__do_report__(content={"plots":list_of_plots}, filename="quartiles")


# TODO: inuput

