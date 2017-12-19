diag_script_info = True

;; Used by: tas_ts_ipcc_1205.ncl
diag_script_info@styleset = "CMIP5"         ; "CMIP5", "DEFAULT"
diag_script_info@ncdf = "tas_ts_line.nc" ;enable to output to netCDF; either use "default" or give a full file name

diag_script_info@scenarios = (/"rcp26", "rcp45", "rcp60", "rcp85"/) ;  scenarios
diag_script_info@syears = (/"1850", "2006", "2101"/) ; start year in time periods
diag_script_info@eyears = (/"2005", "2100", "2300"/) ; end year in time periods
diag_script_info@begin_ref_year = 1986
diag_script_info@end_ref_year = 2005
diag_script_info@label = (/"historical", "RCP2.6", "RCP4.5", "RCP6.0", "RCP8.5"/)
diag_script_info@title = "Global surface temperature change (~S~o~N~ C)"
diag_script_info@yaxis = "(~S~o~N~ C)"
diag_script_info@ymin = -2
diag_script_info@ymax = 12
diag_script_info@model_nr = False
diag_script_info@figure_nr = "12.05"

diag_script_info@diagscript = "ts_line_collins_ipcc_ch12"

;; Or read rgb-file for colormap
 ;delete(diag_script_info@colormap) ;; if defined above it is necessary to delete
                                   ;; the colormap here
colormap_file = "diag_scripts/lib/ncl/rgb/ipcc_color_tseries.rgb"
diag_script_info@colormap = colormap_file

