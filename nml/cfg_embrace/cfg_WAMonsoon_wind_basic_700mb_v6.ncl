diag_script_info = True

diag_script_info@exact_panel_position_file = "./diag_scripts/aux/WAMonsoon/exact_panel_positions_5_by_2.ncl"

; Style
diag_script_info@styleset = "CMIP5"

diag_script_info@season = "JJAS"
diag_script_info@latrange_basic = (/-10, 40/)
diag_script_info@lonrange_basic = (/-30, 50/)

diag_script_info@cn_levels_mean_basic = (/-8, -6, -4, -2, 0, 2, 4, 6, 8/)
diag_script_info@cn_levels_mean_diff_basic = (/-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5/)
diag_script_info@diff_colormap = "diag_scripts/aux/WAMonsoon/cmap_difference_wind.rgb"
diag_script_info@cn_levels_mean_basic_sprintf= "%3.1f"
diag_script_info@cn_levels_mean_diff_basic_sprintf= "%3.0f"

diag_script_info@cn_levels_stddev_basic = (/0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0/)
diag_script_info@cn_levels_stddev_diff_basic = (/-3.0, -2.6, -2.2, -1.8, -1.4, -1.0, -0.6, -0.2, 0.2, 0.6, 1.0, 1.4, 1.8, 2.2, 2.6, 3.0/)

diag_script_info@plottype_lbres = False

diag_script_info@my_region = "WA"
diag_script_info@use_for_contour = "zonal"  ; 'speed', 'zonal' (ua) or 'meridional' (va)
diag_script_info@plot_stddev = False
diag_script_info@plot_diff_reference = False
diag_script_info@three_panel_plots = False

