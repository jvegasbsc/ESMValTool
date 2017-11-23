diag_script_info = True

;; Call the collect and plot script after last variable
if (variables(0).eq."tas") then
    diag_script_info@alt_obs = "NCEP"
    diag_script_info@collect = False
end if
if (variables(0).eq."rlut") then
    diag_script_info@collect = False
end if
if (variables(0).eq."pr") then
    diag_script_info@alt_obs = "CMAP"
    diag_script_info@collect = False
end if
if (variables(0).eq."SW_CRE") then
    diag_script_info@collect = True
end if
