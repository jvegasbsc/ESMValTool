# -*- coding: utf-8 -*-
"""
Created on Fri Feb 09 09:29:45 2018

@author: arunranain
"""

import csv
import numpy as np
from get_metadata_to_rst import do_report


    def __do_esm_validation__(self):
        
        this_function = "ESM validation"
		
        # read in the ESM evaluation grading csv file
        esm_grad_input = os.path.dirname(os.path.realpath(__file__)) + "/lib/predef/esmgrad_expert.csv"

        # plotting routines
        filename = self.__plot_dir__ + os.sep + self.__basic_filename__ + "_" + "".join(this_function.split()) + "." + self.__output_type__
        fig = do_grad_table(self.__varname__, esm_grad_input, os.path.dirname(os.path.realpath(__file__)) + "/example_csvs/example_grad_reference.csv")
        fig.savefig(filename)
        plt.close(fig)
        
        caption = str(this_function + ' for the variable ' + self.__varname__ + ' in the data set "' + "_".join(self.__dataset_id__) + '" (' + self.__time_period__ + ')')

        ESMValMD("meta",
                 filename,
                 self.__basetags__ + ['C3S_Eval'],
                 caption,
                 '#C3S' + 'Eval' + self.__varname__,
                 self.__infile__,
                 self.diagname,
                 self.authors)
        
        # produce report
        self.__do_report__(content={"plots":[filename]}, filename="".join(this_function.upper().split()))
        
        return
