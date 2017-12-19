diag_script_info = True

;; Used by: radTOA_ts_ipcc_1215.ncl
diag_script_info@styleset = "CMIP5"         ; "CMIP5", "DEFAULT"
diag_script_info@ncdf = "rtmt_ts_line.nc"  ; enable to output to netCDF; either use "default" or give a full file name
diag_script_info@scenarios = (/"rcp26", "rcp45", "rcp60", "rcp85"/)  ;  scenarios
diag_script_info@syears = (/"1850", "2006"/)  ; start year in time periods
diag_script_info@eyears = (/"2005", "2100"/)  ; end year in time periods
diag_script_info@begin_ref_year = 1900
diag_script_info@end_ref_year = 1950
diag_script_info@label = (/"historical", "RCP2.6", "RCP4.5", "RCP6.0", "RCP8.5"/)
diag_script_info@title = "a) net flux TOA"
diag_script_info@yaxis = "(Wm~S~-2~N~)"
diag_script_info@ymin = -3
diag_script_info@ymax = 4
diag_script_info@model_nr = False
diag_script_info@figure_nr = "12.15"

diag_script_info@diagscript = "ts_line_collins_ipcc_ch12"

;; Or read rgb-file for colormap
colormap_file = "diag_scripts/lib/ncl/rgb/ipcc_color_tseries.rgb"
diag_script_info@colormap = colormap_file
