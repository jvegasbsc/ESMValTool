;
; Requires: none
;
variable_info = True
variable_info@derived = False

; Countour levels for zonalmean_profile.ncl
variable_info@zonal_ref_Levels = fspan(0, 50, 11)
variable_info@zonal_diff_Levels = (/-20, -10, -5, -3, -2, -1, 0, 1, 2, 3, 5, 10, 20/)

; Minimum pressure for zonalmean_profile.ncl
variable_info@zonal_ymin = 5.

;; Plot units
variable_info@plot_units = "%"

