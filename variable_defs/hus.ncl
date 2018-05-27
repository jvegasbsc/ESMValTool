;
;  Requires: none
;
variable_info = True
variable_info@derived = False

; Countour levels for zonalmean_profile.ncl
variable_info@zonal_ref_Levels = ispan(0, 15, 1)
;variable_info@zonal_diff_Levels = ispan(-5, 5, 1)
;variable_info@zonal_diff_Levels = fspan(-0.5, 0.5, 11)
variable_info@zonal_diff_Levels = (/-1, -0.5, -0.1, -0.05, -0.01, -0.005, -0.001, 0, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1/)

; Minimum pressure for zonalmean_profile.ncl
;variable_info@zonal_ymin = 300.
variable_info@zonal_ymin = 5.

;; Plot units
variable_info@plot_units = "g/kg"
