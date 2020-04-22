diag_script_info = True

;; Used by: pr_mmm_sig_ipcc_1216.ncl
;; map projection, e.g., Mollweide, Mercator
diag_script_info@projection = "Robinson"   
diag_script_info@styleset = "CMIP5"         ; "CMIP5", "DEFAULT"
diag_script_info@scenarios = (/"rcp85"/) ;  scenarios
diag_script_info@periods = (/"2046", "2081", "2181"/) ; start year in time periods
diag_script_info@time_avg = (/"seasonalclim"/)
diag_script_info@seasons = (/0, 1, 2, 3/)  ;season index, DJF:0, MAM:1, JJA:2, SON:3
diag_script_info@percent = 1 ; 0 or 1, difference in percent=1, default = 0
diag_script_info@grid = "1x1"
;diag_script_info@iavmode = "each"
;; max allowed number of rows on a panel page (vertical)
diag_script_info@max_vert     =    4
;; max allowed number of columns on a panel page (horizontal)
diag_script_info@max_hori     =    3        
diag_script_info@title = "Seasonal mean percentage precipitation change (RCP8.5)"
diag_script_info@label = (/"2046-2065 - DJF","2081-2100 - DJF","2181-2200 - DJF", "2046-2065 - MAM","2081-2100 - MAM","2181-2200 - MAM", "2046-2065 - JJA","2081-2100 - JJA","2181-2200 - JJA", "2046-2065 - SON","2081-2100 - SON","2181-2200 - SON"/)
diag_script_info@units = (/"(%)"/)
;; a default colormap, e.g., WhiteBlueGreenYellowRed, rainbow
;diag_script_info@colormap = "WhiteBlueGreenYellowRed"
diag_script_info@sig = True
diag_script_info@not_sig = True
diag_script_info@model_nr = True
diag_script_info@figure_nr = "12.22"
  
;enable to output to netCDF; either use "default" or give a full file name
diag_script_info@ncdf = "pr_mmm_sig_ipcc_1222_seasonalclim.nc" 
diag_script_info@diagscript = "map_diff_mmm_stippandhatch_ipcc12"

colormap_file = "diag_scripts/lib/ncl/rgb/ipcc12_precipitation.rgb"
diag_script_info@colormap = read_colormap_file(colormap_file)

;; Set the contour levels explicitly
diag_script_info@explicit_cn_levels = fspan(-50, 50, 11)
