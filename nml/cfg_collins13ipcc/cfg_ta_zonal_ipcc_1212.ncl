diag_script_info = True

;; Used by: ta_zonal_ipcc_1212.ncl
;; map projection, e.g., Mollweide, Mercator
;diag_script_info@projection = "Robinson"   ;
diag_script_info@styleset = "CMIP5"         ; "CMIP5", "DEFAULT"
diag_script_info@scenarios = (/"rcp26", "rcp45", "rcp85"/) ;  scenarios
diag_script_info@label = (/"RCP2.6", "RCP4.5", "RCP8.5"/)
diag_script_info@periods = (/"2081"/)
diag_script_info@title = "Annual mean atmospheric temperature change (2081-2100)"
diag_script_info@units = (/"(~F35~J~F~C)"/)
diag_script_info@time_avg = (/"annualclim"/)
diag_script_info@grid = "1x1"
;; max allowed number of rows on a panel page (vertical)
diag_script_info@max_vert     =    1
;; max allowed number of columns on a panel page (horizontal)     
diag_script_info@max_hori     =    3        
diag_script_info@sig = True
diag_script_info@not_sig = True
diag_script_info@model_nr = True
diag_script_info@figure_nr = "12.12"

;; a default colormap, e.g., WhiteBlueGreenYellowRed, rainbow
;diag_script_info@colormap = "WhiteBlueGreenYellowRed"

;enable to output to netCDF; either use "default" or give a full file name
diag_script_info@ncdf = "ta_zonal_ipcc_1212.nc" 
diag_script_info@diagscript = "zonal_cont_diff_mmm_stippandhatch_ipcc12"

;; Or read rgb-file for colormap
;delete(diag_script_info@colormap) ;; if defined above it is necessary to delete
                                   ;; the colormap here
colormap_file = "diag_scripts/lib/ncl/rgb/ipcc_temperature_with_grey.rgb"
diag_script_info@colormap = read_colormap_file(colormap_file)
diag_script_info@span = False
;; Set the contour levels explicitly
diag_script_info@explicit_cn_levels = (/-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 3, 4, 5, 7, 9, 11/)
