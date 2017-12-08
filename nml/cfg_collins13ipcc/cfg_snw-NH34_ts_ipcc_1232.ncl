diag_script_info = True

;; Used by: snw_area_seas_ipcc1232.ncl
diag_script_info@styleset = "CMIP5"         ; "CMIP5", "DEFAULT"
;; enable to output to netCDF, either use "default" or give a full file name
diag_script_info@ncdf = "snw_area_seas_ipcc1232.nc" 

diag_script_info@scenarios = (/"rcp26", "rcp45", "rcp60", "rcp85"/) ;  scenarios
diag_script_info@syears = (/"1980", "2006"/) ; start year in time periods
diag_script_info@eyears = (/"2005", "2100"/) ; end year in time periods
diag_script_info@label = (/"historical", "RCP2.6", "RCP4.5", "RCP6.0", "RCP8.5"/)
diag_script_info@title = "Relative snow cover extent change"
diag_script_info@yaxis = "(-)"
diag_script_info@ts_minlat = (/0/)
diag_script_info@ts_maxlat = (/90/)
diag_script_info@ts_minlon = (/0/)
diag_script_info@ts_maxlon = (/360/)
diag_script_info@begin_ref_year = 1986
diag_script_info@end_ref_year = 2005
diag_script_info@months = (/"3", "4"/)
diag_script_info@figure_nr = "12.32"

diag_script_info@ymin = -0.4
diag_script_info@ymax = 0.3

diag_script_info@model_nr = False
;; diagscript name which uses file, to find output netcdf when plotting
diag_script_info@diagscript = "snw_area_seas_ipcc1232"

;; Or read rgb-file for colormap
;delete(diag_script_info@colormap) ;; if defined above it is necessary to delete
                                   ;; the colormap here
colormap_file = "diag_scripts/lib/ncl/rgb/ipcc_color_tseries.rgb"
diag_script_info@colormap = colormap_file

