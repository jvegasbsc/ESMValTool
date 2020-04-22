diag_script_info = True

;; Used by: IAV_for_stippandhatch_ipcc12.ncl
diag_script_info@styleset = "CMIP5"         ; "CMIP5", "DEFAULT"
diag_script_info@periodlength = 20.         ; in years
diag_script_info@time_avg = (/"seasonalclim"/)
diag_script_info@seasons = (/0, 1, 2, 3/)
diag_script_info@overwrite = True
diag_script_info@grid = "1x1"
;enable to output to netCDF; either use "default" or give a full file name
diag_script_info@ncdf = "default"

