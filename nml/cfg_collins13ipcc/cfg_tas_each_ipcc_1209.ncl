diag_script_info = True

;; Used by: map_diff_tas_each_model_1209.ncl
;; map projection, e.g., Mollweide, Mercator
diag_script_info@projection = "Robinson"   
diag_script_info@styleset = "CMIP5"         ; "CMIP5", "DEFAULT"
;; max allowed number of rows on a panel page (vertical)
diag_script_info@max_vert     =    9
;; max allowed number of columns on a panel page (horizontal)     
diag_script_info@max_hori     =    5 ; 5 is good for all models        
diag_script_info@title = "Annual mean surface air temperature change (RCP4.5: 2081-2100)"
diag_script_info@time_avg = "annualclim"
diag_script_info@experiment = "rcp45"
diag_script_info@figure_nr = "12.09"
;; a default colormap, e.g., WhiteBlueGreenYellowRed, rainbow
;diag_script_info@colormap = "WhiteBlueGreenYellowRed"

;enable to output to netCDF; either use "default" or give a full file name
diag_script_info@ncdf = "_tas_each_model_1209.nc" 

;; Or read rgb-file for colormap
;delete(diag_script_info@colormap) ;; if defined above it is necessary to delete
                                   ;; the colormap here
colormap_file = "diag_scripts/lib/ncl/rgb/ipcc_temperature_with_grey.rgb"
diag_script_info@colormap = read_colormap_file(colormap_file)

;; Set the contour levels explicitly
diag_script_info@explicit_cn_levels = (/-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 3, 4, 5, 7, 9, 11/)
