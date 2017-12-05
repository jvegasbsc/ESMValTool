diag_script_info = True

;; Used by: map_diff_mmm_stippandhatch_ipcc12.ncl
;; map projection, e.g., Mollweide, Mercator
diag_script_info@projection = "Robinson"   
diag_script_info@styleset = "CMIP5"         ; "CMIP5", "DEFAULT"
diag_script_info@scenarios = (/"rcp26", "rcp45", "rcp60", "rcp85"/) ;  scenarios
diag_script_info@periods = (/"2081", "2181"/) ; start year in time periods
diag_script_info@time_avg = (/"annualclim"/)
diag_script_info@grid = "1x1"
;; max allowed number of rows on a panel page (vertical), number of scenarios
diag_script_info@max_vert     =    2
;; max allowed number of columns on a panel page (horizontal), number of periods
diag_script_info@max_hori     =    2        
diag_script_info@title = "Temperature scaled by global T (~F35~J~F~C per ~F35~J~F~C)"
diag_script_info@label = (/"2081-2100", "2181-2200"/)
diag_script_info@units = (/"(~F35~J~F~C per ~F35~J~F~C)"/)
diag_script_info@sig = True
diag_script_info@figure_nr = "12.10"
;; a default colormap, e.g., WhiteBlueGreenYellowRed, rainbow
;diag_script_info@colormap = "WhiteBlueGreenYellowRed"
  
;enable to output to netCDF; either use "default" or give a full file name
diag_script_info@ncdf = "tas_mmm_sig_ipcc_1210.nc"
diag_script_info@diagscript = "calc_map_diff_scaleT_mmm_stipp_ipcc12"
;; Or read rgb-file for colormap
;delete(diag_script_info@colormap) ;; if defined above it is necessary to delete
                                   ;; the colormap here
colormap_file = "diag_scripts/lib/ncl/rgb/ipcc_tx_grey2darkred.rgb"
diag_script_info@colormap = read_colormap_file(colormap_file)
diag_script_info@span = False
;; Set the contour levels explicitly
diag_script_info@explicit_cn_levels = (/0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2/)
