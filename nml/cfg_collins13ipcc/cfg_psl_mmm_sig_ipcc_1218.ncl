diag_script_info = True

;; Used by: radTOA_mmm_sig_ipcc_1216.ncl
;; map projection, e.g., Mollweide, Mercator
diag_script_info@projection = "Robinson"   
diag_script_info@styleset = "CMIP5"         ; "CMIP5", "DEFAULT"
diag_script_info@scenarios = (/"rcp26", "rcp45", "rcp85"/) ;  scenarios
diag_script_info@periods = (/"2081"/) ; start year in time periods
diag_script_info@seasons = (/0, 2/)  ;season index, DJF:0, MAM:1, JJA:2, SON:3
diag_script_info@time_avg = (/"seasonalclim"/)
diag_script_info@grid = "1x1"
;; max allowed number of rows on a panel page (vertical)
diag_script_info@max_vert     =    2
;; max allowed number of columns on a panel page (horizontal)
diag_script_info@max_hori     =    3        
diag_script_info@title = "Seasonal mean sea level pressure change (2081-2100)"
diag_script_info@label = (/"RCP2.6 - DJF","RCP4.5 - DJF","RCP8.5 - DJF", "RCP2.6 - JJA","RCP4.5 - JJA","RCP8.5 - JJA"/)
;; a default colormap, e.g., WhiteBlueGreenYellowRed, rainbow
;diag_script_info@colormap = "WhiteBlueGreenYellowRed"
diag_script_info@sig = True
diag_script_info@not_sig = True
diag_script_info@model_nr = True
diag_script_info@figure_nr = "12.18"
  
;enable to output to netCDF; either use "default" or give a full file name
diag_script_info@ncdf = "psl_mmm_sig_ipcc_1218.nc" 
diag_script_info@diagscript = "map_diff_mmm_stippandhatch_ipcc12"

;; Or read rgb-file for colormap
;delete(diag_script_info@colormap) ;; if defined above it is necessary to delete
                                   ;; the colormap here
colormap_file = "diag_scripts/lib/ncl/rgb/ipcc_colors_blu2red_centered.rgb"
diag_script_info@colormap = read_colormap_file(colormap_file)

;; Set the contour levels explicitly
diag_script_info@explicit_cn_levels = fspan(-4, 4, 17)
