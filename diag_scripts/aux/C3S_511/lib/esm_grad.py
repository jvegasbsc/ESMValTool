# -*- coding: utf-8 -*-
"""
Created on Fri Feb 09 09:29:45 2018

@author: hass_bg
"""

import csv
import numpy as np
from get_metadata_to_rst import do_report


    def __do_esm_validation__(self):
        
        this_function = "ESM validation"
		
        # read in the ESM evaluation grading csv file
        esm_grad_input = os.path.dirname(os.path.realpath(__file__)) + "/lib/predef/esmgrad_expert.csv"
		
        # Create a list. Each item of the list will be itself a list of strings, corresponding either to the 
        # headers or to the ESM evaluation entries for the different ECVs
        contents = list()
        with open(esm_grad_input, 'rb') as csvfile:
            s = csv.reader(csvfile, delimiter = ";", skipinitialspace = True)
            for row in s:
                contents.append(row)
		
        # convert the list into an array
        esmgrad_data = np.asarray(contents)
		
        # check the number of entries for the ECV in question, and write all of the available entries in an ordered dictionary
        # for easy output in the reports
        esmgrad_dict=collections.OrderedDict()
		
        for num_entries in range(0, len(np.nonzero(esmgrad_data == self.__varname__)[0])):
            insert_dict=collections.OrderedDict()
            for column in range(1, len(esmgrad_data[0,:])): 
                insert_dict.update({esmgrad_data[0, column]: esmgrad_data[np.nonzero(esmgrad_data == self.__varname__)[0][num_entries], column]}) 
            esmgrad_dict.update({'R' + str(num_entries + 1):insert_dict})

		              
        # produce report
        self.__do_report__(content={"text":esmgrad_dict}, filename = this_function.upper())
        
        return
