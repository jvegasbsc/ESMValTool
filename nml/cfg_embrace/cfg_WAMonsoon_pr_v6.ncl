; Settings used by most WAMonsoon plot_types
diag_script_info = True
diag_script_info@season = "JJAS"

diag_script_info@exact_panel_position_file = "./diag_scripts/aux/WAMonsoon/exact_panel_positions_5_by_2.ncl"

; Settings specific for the precip basic plot_type
diag_script_info@latrange_basic = (/-10, 40/)
diag_script_info@lonrange_basic = (/-30, 50/)
;diag_script_info@latrange_basic = (/-25, 35/)
;diag_script_info@lonrange_basic = (/-30, 50/)

diag_script_info@diff_colormap = "diag_scripts/aux/WAMonsoon/cmap_difference_inv.rgb"
;diag_script_info@diff_colormap = "ncl_default"

diag_script_info@cn_levels_mean_basic = ispan(1, 10, 1)
diag_script_info@cn_levels_mean_basic_sprintf = "%1.0f"
diag_script_info@cn_levels_mean_diff_basic_sprintf = "%1.0f"

diag_script_info@cn_levels_mean_diff_basic = (/-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5/)

diag_script_info@cn_levels_stddev_basic = (/0.0, 0.3, 0.6, 0.9, 1.2, 1.5, 1.8, 2.1, 2.4, 2.7, 3.0/)
diag_script_info@cn_levels_stddev_diff_basic = (/-2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5/)


diag_script_info@my_region = "WA"
diag_script_info@plot_stddev = False
diag_script_info@plot_diff_reference = False
diag_script_info@three_panel_plots = False

;diag_script_info@invert_diff_palette = True

