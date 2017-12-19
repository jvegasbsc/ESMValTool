diag_script_info = True

;; Used by: map_diff_mmm_stippandhatch_ipcc12.ncl
;; map projection, e.g., Mollweide, Mercator
diag_script_info@projection = "Robinson"   
diag_script_info@styleset = "CMIP5"         ; "CMIP5", "DEFAULT"
diag_script_info@scenarios = (/"rcp85"/) ;  scenarios
diag_script_info@periods = (/"2081"/) ; start year in time periods
diag_script_info@time_avg = (/"annualclim"/)
;diag_script_info@iavmode = "each"
diag_script_info@grid = "1x1"
;; max allowed number of rows on a panel page (vertical), number of scenarios
diag_script_info@max_vert     =    1
;; max allowed number of columns on a panel page (horizontal), number of periods
diag_script_info@max_hori     =    1       
diag_script_info@title = "Annual mean surface salinity change"
diag_script_info@label = (/"RCP8.5: 2081-2100"/)
diag_script_info@units = (/"(psu)"/)
diag_script_info@sig = True
diag_script_info@not_sig = True
diag_script_info@model_nr = True
diag_script_info@figure_nr = "12.34"
;; a default colormap, e.g., WhiteBlueGreenYellowRed, rainbow
;diag_script_info@colormap = "WhiteBlueGreenYellowRed"
  
;enable to output to netCDF; either use "default" or give a full file name
diag_script_info@ncdf = "sos_mmm_sig_ipcc_1211.nc"
diag_script_info@diagscript = "map_diff_mmm_stippandhatch_ipcc12"

;; Or read rgb-file for colormap
;delete(diag_script_info@colormap) ;; if defined above it is necessary to delete
                                   ;; the colormap here
colormap_file = "diag_scripts/lib/ncl/rgb/ipcc_colors_blu2red_centered.rgb"
diag_script_info@colormap = read_colormap_file(colormap_file)
diag_script_info@span = True
;; Set the contour levels explicitly
diag_script_info@explicit_cn_levels = (/-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2/)
