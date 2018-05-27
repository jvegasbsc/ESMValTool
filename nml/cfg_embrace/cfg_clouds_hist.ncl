diag_script_info = True

diag_script_info@projection = "CylindricalEquidistant"   ; map projection, e.g., Mollweide, Mercator

diag_script_info@showdiff = True
diag_script_info@extralegend = True
diag_script_info@calcmm = False
diag_script_info@embrace_setup = True      ; True  = 2 plots per line
                                           ; False = 4 plots per line

diag_script_info@filename_add = "hist"

; --- time averaging ---
diag_script_info@timemean = "annualclim"   ; options: 1) annualclim (1 plots: annual mean)
                                           ;          2) seasonalclim (4 plots: DJF, MAM, JJA, SON)

; target grid for regridding: use grid of reference model
; (valid options are: "ref", "coarsest", "finest", "1x1")
diag_script_info@grid = "1x1"

if (variables(0).eq."lwp") then
    diag_script_info@ref_model = (/"UWisc"/)
end if
if (any(variables(0).eq.(/"rlutcs", "rlut", "rsutcs", "rsut","SW_CRE","LW_CRE","NET_CRE"/))) then
    diag_script_info@ref_model = (/"CERES-EBAF"/)
end if
if (any(variables(0).eq.(/"clt"/))) then
    diag_script_info@ref_model = (/"MODIS"/)
end if
if (any(variables(0).eq.(/"clivi"/))) then
    diag_script_info@ref_model = (/"MODIS-CFMIP"/)
end if
if (any(variables(0).eq.(/"pr","pr-mmday"/))) then
    diag_script_info@ref_model = (/"GPCP-SG"/)
end if
if (any(variables(0).eq.(/"tas-degC","tas"/))) then
    diag_script_info@ref_model = (/"ERA-Interim"/)
end if
; -------------------------------------------------------------

