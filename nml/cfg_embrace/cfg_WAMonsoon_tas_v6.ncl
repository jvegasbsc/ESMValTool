; Settings used by most WAMonsoon plot_types
diag_script_info = True
diag_script_info@season = "JJAS"

diag_script_info@exact_panel_position_file = "./diag_scripts/aux/WAMonsoon/exact_panel_positions_5_by_2.ncl"

; Settings specific for the precip basic plot_type
diag_script_info@latrange_basic = (/-10, 40/)
diag_script_info@lonrange_basic = (/-30, 50/)
;diag_script_info@latrange_basic = (/-25, 35/)
;diag_script_info@lonrange_basic = (/-30, 50/)

diag_script_info@diff_colormap = "diag_scripts/aux/WAMonsoon/cmap_difference_tas.rgb"

diag_script_info@cn_levels_mean_basic = ispan(16, 40, 2)
diag_script_info@cn_levels_mean_diff_basic = ispan(-7, 7, 1)
diag_script_info@cn_levels_mean_basic_sprintf = "%3.0f"
diag_script_info@cn_levels_mean_diff_basic_sprintf = "%1.0f"

diag_script_info@cn_levels_stddev_basic = (/0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0/)
diag_script_info@cn_levels_stddev_diff_basic = (/-2.0, -1.7, -1.4, -1.1, -0.8, -0.5, -0.2, 0.1, 0.4, 0.7, 1.0, 1.3, 1.6, 1.9, 2.2, 2.5/)

diag_script_info@my_region = "WA"
diag_script_info@plot_stddev = False
diag_script_info@plot_diff_reference = False
diag_script_info@three_panel_plots = False

