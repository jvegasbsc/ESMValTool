diag_script_info = True

;; Define region for the spatial average
diag_script_info@region = "Global"

;; Plot legend in a separate file
diag_script_info@legend_outside = False

;; Plot style (EMAC, CMIP5, DEFAULT)
diag_script_info@styleset = "CMIP5"

;; Sort model alphabetically
diag_script_info@sort = True

;; Draw seasonal cycle plot, including diff y-axis setting
diag_script_info@seasonal_cycle_plot = False

;; Draw errorbar plot
diag_script_info@errorbar_plot = True

;; Draw Mean (x-axis), IAV (y-axis) plot
diag_script_info@mean_IAV_plot = False

;; For regridding; default: 2deg - others not supported yet
diag_script_info@grid = "2deg"

;; GCP obs external read-in
diag_script_info@obsname = "GCP"
path = getenv("ESMValTool_OBSPATH") + "/Tier2/GCP/"
;; Currently only one file will be read in
if variables(0).eq."NBP" then
    diag_script_info@ref_varname = "nbp"
else
    diag_script_info@ref_varname = "fgco2"
end if
diag_script_info@obsfile = path + \
       "OBS_GCP_reanaly_1_T0M_" + \
       diag_script_info@ref_varname + \
       "_195901-201512.nc"
