diag_script_info = True

;; Used by: map_diff_mmm_stippandhatch_ipcc12.ncl
;; map projection, e.g., Mollweide, Mercator
diag_script_info@projection = "Robinson"   
diag_script_info@styleset = "CMIP5"         ; "CMIP5", "DEFAULT"
diag_script_info@scenarios = (/"rcp26", "rcp45", "rcp60", "rcp85"/) ;  scenarios
diag_script_info@periods = (/"2046", "2081", "2181"/) ; start year in time periods
diag_script_info@time_avg = (/"annualclim"/)
;diag_script_info@iavmode = "each"
diag_script_info@grid = "1x1"
;; max allowed number of rows on a panel page (vertical), number of scenarios
diag_script_info@max_vert     =    4
;; max allowed number of columns on a panel page (horizontal), number of periods
diag_script_info@max_hori     =    3        
diag_script_info@title = "Annual mean surface air temperature change"
diag_script_info@label = (/"RCP2.6: 2046-2065", "RCP2.6: 2081-2100", "RCP2.6: 2181-2200", "RCP4.5: 2046-2065", "RCP4.5: 2081-2100", "RCP4.5: 2181-2200", "RCP6.0: 2046-2065", "RCP6.0: 2081-2100", "RCP6.0: 2181-2200", "RCP8.5: 2046-2065", "RCP8.5: 2081-2100", "RCP8.5: 2181-2200"/)
diag_script_info@units = (/"(~F35~J~F~C)"/)
diag_script_info@figure_nr = "12.11"
;; a default colormap, e.g., WhiteBlueGreenYellowRed, rainbow
;diag_script_info@colormap = "WhiteBlueGreenYellowRed"
diag_script_info@sig = True
diag_script_info@not_sig = True
diag_script_info@model_nr = True
  
;enable to output to netCDF; either use "default" or give a full file name
diag_script_info@ncdf = "tas_mmm_sig_ipcc_1211.nc"
diag_script_info@diagscript = "map_diff_mmm_stippandhatch_ipcc12"

;; Or read rgb-file for colormap
;delete(diag_script_info@colormap) ;; if defined above it is necessary to delete
                                   ;; the colormap here
colormap_file = "diag_scripts/lib/ncl/rgb/ipcc_temperature_with_grey.rgb"
diag_script_info@colormap = read_colormap_file(colormap_file)
diag_script_info@span = False
;; Set the contour levels explicitly
diag_script_info@explicit_cn_levels = (/-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,3,4,5,7,9,11/)
