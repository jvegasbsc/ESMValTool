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
if (any(variables(0).eq.(/"tas-degC", "pr-mmday", "tos"/))) then
    diag_script_info@seasonal_cycle_plot = False
else
    diag_script_info@seasonal_cycle_plot = True
    if (any(variables(0).eq.(/"GPP", "NBP", "FGCO2"/))) then
        diag_script_info@anav_month = True
    end if
end if

;; Draw errorbar plot
if (any(variables(0).eq.(/"tas-degC", "pr-mmday", "tos"/))) then
    diag_script_info@errorbar_plot = False
else
    diag_script_info@errorbar_plot = True
end if

;; Draw Mean (x-axis), IAV (y-axis) plot
diag_script_info@mean_IAV_plot = False

;; For regridding; default: 2deg - others not supported yet
diag_script_info@grid = "2deg"

;; Reference model and evolution plot
;; Evolution plot only with 2 or more non-obs models
if (any(variables(0).eq.(/"tas-degC", "pr-mmday"/))) then
    diag_script_info@evolution_plot_ref_model = "CRU"
    diag_script_info@evolution_plot = True
    diag_script_info@evolution_plot_anomaly = True
    diag_script_info@evolution_plot_volcanoes = True
end if
if (variables(0).eq."tos") then
    diag_script_info@evolution_plot_ref_model = "HadISST"
    diag_script_info@evolution_plot = True
    diag_script_info@evolution_plot_anomaly = True
    diag_script_info@evolution_plot_volcanoes = True
end if
if (any(variables(0).eq.(/"NBP", "FGCO2"/))) then
    diag_script_info@evolution_plot = True
    diag_script_info@evolution_plot_volcanoes = False
    diag_script_info@evolution_plot_ref_model = "GCP"
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
end if
