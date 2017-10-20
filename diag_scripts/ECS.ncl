;;#############################################################################
;; CALCULATES ECS (Equilibrium climate sensitivity)
;; Author: Ricarda Rettinger (DLR, Germany)
;; CRESCENDO project
;;#############################################################################
;;
;; Description
;;    Flatens a 2D/3D array into a 1D time-series and calculates ECS  
;;    following Andrews et al. (2012).
;;    After doing a linear fit on the piControl runs (variables = tas/rtmt) 
;;    it is substracted from the corresponding abr4xCO2 experiment.
;;    The result is plotted and a netcdf file is generated where the ECS value  
;;    of every model is listed.
;;
;; Required diag_script_info attributes (diagnostics specific)
;;    none
;;    
;; Optional diag_script_info attributes (diagnostic specific)
;;    ts_minlat, ts_maxlat, ts_minlon, ts_maxlon: possible to define extent 
;;    default is global: (-90, 90, 0, 360)
;;
;; Required variable_info attributes (variable specific)
;;    none
;;
;; Required variable attributes (defined in namelist)
;;    none
;;
;; Caveats
;;
;; Modification history
;;    20171013-A_laue_ax: * bugfix minlat, maxlat, minlon, maxlon
;;                        * added meta data to plots for reporting
;;    20170320-A_rett_ri: written.
;;
;;#############################################################################
load "./interface_data/ncl.interface"
load "./interface_scripts/auxiliary.ncl"
load "./interface_scripts/data_handling.ncl"
load "./interface_scripts/messaging.ncl"

;load "./diag_scripts/lib/ncl/scaling.ncl"
load "./diag_scripts/lib/ncl/statistics.ncl"
load "./diag_scripts/lib/ncl/latlon.ncl"
load "./diag_scripts/lib/ncl/style.ncl"

load "./plot_scripts/ncl/xy_line.ncl"

load "./diag_scripts/lib/ncl/meta_data.ncl"

begin
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering " + diag_script, verbosity, 4)
    info_output("++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)
    info_output(diag_script + " (var: " + variables(0) + ")", verbosity, 1)
    info_output("++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)

    ;; References
;     write_references(diag_script, "A_rett_ri", "A_laue_ax", \
;     "D_0001",  "", "P_crescendo")
    dim_MOD = dimsizes(models@name)
    dim_VAR = dimsizes(variables)
    var0 = variables(0)
    var1 = variables(1)
    field_type0 = field_types(0)
    field_type1 = field_types(1)

    ;; Get time range
    start_year = min(stringtoint(models@start_year))
    end_year   = max(stringtoint(models@end_year))
    ntime = 150

    ;; Get area extent
    if (isatt(diag_script_info, "minlat")) then
        minlat = diag_script_info@minlat
    else
        minlat = -90
    end if

    if (isatt(diag_script_info, "maxlat")) then
        maxlat = diag_script_info@maxlat
    else
        maxlat = 90
    end if

    if (isatt(diag_script_info, "minlon")) then
        minlon = diag_script_info@minlon
    else
        minlon = 0
    end if

    if (isatt(diag_script_info, "maxlon")) then
        maxlon = diag_script_info@maxlon
    else
        maxlon = 360
    end if
 
end


begin
    vardeffile = "interface_data/" + var0 + "_info.tmp"
    loadscript(vardeffile)
    vardeffile2 = "interface_data/" + var1 + "_info.tmp"
    loadscript(vardeffile2)
end


begin
    ;; Create empty model arrays for rtmt and tas
    rtmt_arr = new((/dim_MOD, ntime/), float)
    rtmt_arr!0 = "model"
    rtmt_arr!1 = "time"

    tas_arr = new((/dim_MOD, ntime/), float)
    tas_arr!0 = "model"
    tas_arr!1 = "time"

    ;; Loop over models
    do imod = 0, dim_MOD - 1

        ;; Read data
        data_rtmt = read_data(imod, var1, field_type1)           
        d_rtmt = extract_data(imod, data_rtmt, -1, 0, 0)

        data_tas = read_data(imod, var0, field_type0)           
        d_tas = extract_data(imod, data_tas, -1, 0, 0) 

        ;; Calculate time average
        d_rtmt_yrl = time_operations(d_rtmt, -1, -1, "average", "yearly", True)
        d_tas_yrl = time_operations(d_tas, -1, -1, "average", "yearly", True)
        
        ;; Calculate area average
        d_rtmt_yrl_area = area_operations(d_rtmt_yrl, minlat, \
                          maxlat, \
                          minlon, \
                          maxlon, "average", \
                          True)
        d_tas_yrl_area = area_operations(d_tas_yrl, minlat, \
                         maxlat, \
                         minlon, \
                         maxlon, "average", \
                         True)

        ;; Extend shorter time length with NAs
        period = stringtoint(models@end_year(imod))-(stringtoint(models@start_year(imod)))
        if (period .lt. 149) then
            t = 149 - period
            zero = new((/t/),float)
            d_rtmt_yrl_area := array_append_record(d_rtmt_yrl_area,zero, 0)
            d_tas_yrl_area := array_append_record(d_tas_yrl_area,zero, 0)

            delete(t)
            delete(zero)
        end if

        ;; Write model to model array
        rtmt_arr(imod, :) = d_rtmt_yrl_area
        tas_arr(imod, :) = d_tas_yrl_area
 
        delete(d_rtmt_yrl_area)
        delete(d_tas_yrl_area)
        delete(d_tas_yrl)
        delete(d_rtmt_yrl)
        delete(d_tas)
        delete(d_rtmt)
        delete(data_rtmt)
        delete(data_tas)
        delete(period)
             
    end do
    delete(imod)

    ;; get Modelnames
    modelnames = project_style(diag_script_info, "annots")
    
    ;; Create substrings to get index of models
    ind_4xCO2 = str_match_ind(modelnames, "abrupt4xCO2")
    ind_piCon = str_match_ind(modelnames, "piControl")
 
    ;; regression of piControl (tas/ rtmt)
    piCon_tas  = new((/dimsizes(ind_piCon), ntime/), float)
    piCon_rtmt = new((/dimsizes(ind_piCon), ntime/), float)

    do n = 0, dimsizes(ind_piCon) -1

        ;; select only the piControl Experiment
        tas_var = tas_arr(ind_piCon(n), :)
        tas_nm  = dimsizes(tas_var)

        rtmt_var = rtmt_arr(ind_piCon(n), :)
        rtmt_nm  = dimsizes(rtmt_var)

        if (rtmt_nm - tas_nm .ne. 0) then
            error_msg( "f", diag_script, "", \
                      "time steps of tas and rtmt differ")
        end if

        years_ind = ispan(0, ntime-1, 1)

        ;; least squares regression of rtmt/tas (as in Adrews: x = tas)
        var_tas = regline( years_ind, tas_var)
        reg_tas = var_tas * (years_ind - var_tas@xave) + var_tas@yave

        var_rtmt = regline( years_ind ,rtmt_var)
        reg_rtmt = var_rtmt * (years_ind - var_rtmt@xave) + var_rtmt@yave

        piCon_tas(n, :) = reg_tas
        piCon_rtmt(n, :) = reg_rtmt
    
        delete(tas_var)
        delete(tas_nm)
        delete(rtmt_var)
        delete(rtmt_nm) 
        delete(var_tas)
        delete(var_rtmt)
        delete(reg_tas)
        delete(reg_rtmt)
        delete(years_ind)
        
    end do
    delete(n)

    ;; Define coordinates
    piCon_tas!0 = "model"
    piCon_tas!1 = "time"
    piCon_tas&time = ispan(0, ntime-1, 1)
    piCon_tas&model = str_match(modelnames, "piControl")
  
    piCon_rtmt!0 = "model"
    piCon_rtmt!1 = "time"
    piCon_rtmt&time = ispan(0, ntime-1, 1)
    piCon_rtmt&model = str_match(modelnames, "piControl")

;    piCon_reg!0 = "model"
;    piCon_reg!1 = "time"
;    piCon_reg&time = ispan(0, ntime-1, 1)
;    piCon_reg&model = str_match(modelnames, "piControl")

    ;; Create abrupt4xCO2 arrays
    abr4xCO2_tas = new((/dimsizes(ind_4xCO2), ntime/), float)
    abr4xCO2_rtmt = new((/dimsizes(ind_4xCO2), ntime/), float)

    do n = 0, dimsizes(ind_4xCO2) -1
       
       ;; select only the abr4xCO2 - Experiment
       tas_var = tas_arr(ind_4xCO2(n), :)
       tas_nm  = dimsizes(tas_var)

       rtmt_var = rtmt_arr(ind_4xCO2(n), :)
       rtmt_nm  = dimsizes(rtmt_var)
       
       ;; Control if rtmt and tas have same number of years
       if (rtmt_nm - tas_nm .ne. 0) then
           error_msg( "f", diag_script, "", \
                      "time steps of tas and rtmt differ")
       end if
       
       ;; extract years
       abr4xCO2_rtmt(n, :) = rtmt_var
       abr4xCO2_tas(n, :) = tas_var
  
       delete(rtmt_var)
       delete(rtmt_nm)
       delete(tas_var)
       delete(tas_nm)

    end do

    ;; Define coordinates
    abr4xCO2_rtmt!0 = "model"
    abr4xCO2_rtmt!1 = "time"
    abr4xCO2_rtmt&time = ispan(0, ntime-1, 1)
    abr4xCO2_rtmt&model = str_match(modelnames, "abrupt4xCO2")

    abr4xCO2_tas!0 = "model"
    abr4xCO2_tas!1 = "time"
    abr4xCO2_tas&time = ispan(0, ntime-1, 1)
    abr4xCO2_tas&model = str_match(modelnames, "abrupt4xCO2")
    
    if (dimsizes(abr4xCO2_rtmt&model) .ne. dimsizes(piCon_rtmt&model)) then
        error_msg("f", diag_script, "", \
                  "different Models in use!")
    end if

    delete(n)
end

begin
    ;; Substract piControl from abrupt4xCO2
    subs_rtmt = abr4xCO2_rtmt - piCon_rtmt
    subs_tas  = abr4xCO2_tas  - piCon_tas
     
    ;; Defining Dimensions and Coordinates
    subs_rtmt!0 = "model"
    subs_rtmt!1 = "time"
    subs_rtmt&time = piCon_tas&time
    subs_rtmt&model = str_match(modelnames, "piControl")

    subs_tas!0 = "model"
    subs_tas!1 = "time"
    subs_tas&time = piCon_tas&time
    subs_tas&model = str_match(modelnames, "piControl")


    ;; Time averaging
    if (isatt(diag_script_info, "time_avg")) then
        time_avg = diag_script_info@time_avg
        if (time_avg.ne."yearly") then
            error_msg("f", diag_script, "", \
                      "time averaging option " + time_avg + " not yet " + \
                      "implemented")
        end if
    else
        time_avg = "yearly"
    end if

    out_path = get_output_dir(get_work_dir(), basename(diag_script))
    out_dir = out_path + "ECS.txt"

    mod_names       = new((/dimsizes(subs_tas&model)/), string)
    rtmt_regression = new((/dimsizes(subs_tas&model)/), float)
    alpha_net       = new((/dimsizes(subs_tas&model)/), float)
    ECS_2           = new((/dimsizes(subs_tas&model)/), float)

; -------------------------------------------------------------------------------
; collect information for adding meta data to plots (for reporting)

climofiles = new(2, string)

alltags = array_append_record(tags, (/"PT_scatter", "DM_global"/), 0)
id = diag_script + "_ECS"
contrib_authors = (/"A_rett_ri", "A_laue_ax"/) 
; -------------------------------------------------------------------------------

do n = 0, dimsizes(subs_tas&model) -1

    mod_name = str_split(subs_tas&model(n), "_")
    model_name = mod_name(0)

    ;; Create data frame for scatterplot
    data = new((/2, dimsizes(subs_tas&time) + 1/), typeof(subs_tas))
    fillval = (/(/-10/)/)
    data (0, :) = array_append_record(tofloat(fillval),subs_rtmt(n, :), 0)
    
    ;; regress rtmt/tas and append to data
    rc     = regline_stats(subs_tas(n, :), subs_rtmt(n, :))
    YAx    = rc@yintercept
    YAx!0  = "time"
    regr   = rc * (subs_tas(n, :) - rc@xave) + rc@yave
    regr!0 = "time"
    YAx@_FillValue = regr@_FillValue
    data (1, :) = array_append_record(tofloat(YAx), regr, 0)
     
    ;; Add 0 to Delta T
    deltaT = array_append_record(tofloat(0), subs_tas(n, :), 0)

    ;; Define workstation
    outfile = var0 + "_" +  "_" +  model_name
    wks = get_wks(0, diag_script, outfile)

    ;; Set resources
    res                   = True                     ; plot mods desired
    res@gsnDraw           = False                    ; don't draw
    res@gsnFrame          = False                    ; don't advance frame
    res@xyMarkLineModes   = (/"Markers", "Lines"/)   ; choose which have markers
    res@xyMarkers         =  9                       ; choose type of marker  
    res@xyMarkerColor     = (/"blue", "black"/)      ; Marker color
    res@xyMarkerSizeF     = 0.01                     ; Marker size (default 0.01)
    res@tiMainString      = model_name
    res@xyDashPatterns    = 1                        ; solid Line

    res@tiYAxisString     = "N" + " [" + rtmt_arr@units + "]"
    res@tiXAxisString     = "Delta"+" T [" + tas_arr@units + "]"
  
    res@gsnYRefLine       = 0                        ; create horizontal line
    res@gsnYRefLineColor  = "black"
    res@gsnYRefLineDashPattern = 1

    res@trXMinF           = 0
    res@trYMinF           = -2
    res@trXMaxF           = 7
    res@trYMaxF           = 10

    plot = gsn_csm_xy (wks, deltaT, data, res)

    txres                 = True                            
    txres@txFontHeightF   = 0.018                       ; Set the font height

    label = " r = -" + decimalPlaces(rc@r, 2, True) + \
            "~C~ Y-Intercept = " + \ 
            decimalPlaces(rc@yintercept, 2, True) +\
            "    ECS = " + \
            decimalPlaces((rc@yintercept/(2 *(-rc))), 2, True)
    gsn_text(wks, plot, label, 2.2, -1, txres)

    draw(plot)
    frame(wks)
 
    info_output(" gv " + wks@fullname + " &", verbosity, 1)

    ; -------------------------------------------------------------------------------
    ; add meta data to plot (for reporting)
    caption = "Equilibrium climate sensitivity for model " + model_name + "."
    climofiles(0) = interface_get_inpaths(n) + "/" \
        + interface_get_infile(variables(0), field_types(0), n)
    climofiles(1) = interface_get_inpaths(n) + "/" \
        + interface_get_infile(variables(1), field_types(1), n)

    ESMValMD(wks@fullname, alltags, caption, id, variables, \
             model_name, climofiles, diag_script, contrib_authors)
    ; -------------------------------------------------------------------------------

    ;; Write Strings for Table 
    mod_names(n) = tochar(model_name)
    rtmt_regression(n) = rc@yintercept
    alpha_net(n) = rc
    ECS_2(n) = rtmt_regression(n) / (2 * (-rc))

    Fill_Val =  regr@_FillValue
    delete(YAx)
    delete(regr)
    delete(data)
    delete(rc)
    delete(model_name)
    delete(mod_name)
    delete(plot)
    delete(res)
    delete(deltaT)

end do

    ;; Export Data to .nc
    out_d = out_path + "ECS.nc"
    output_nc = ECS_2    
    output_nc!0 =  "model"
    output_nc&model =   mod_names 
    output_nc@var = "ecs"
    output_nc@ncdf = out_path
    output_nc@diag_script = diag_script
    output_nc@long_name = "Equilibrium Climate Sensitivity"
    output_nc@standard_name = "ECS"
    output_nc@units = "K"
    output_nc@_FillValue = Fill_Val
    ncdf_outfile = ncdf_write(output_nc, out_d)
    info_output(" ncdump " + out_d + " &", verbosity, 1)

    delete(out_d)
    delete(output_nc)
    delete(Fill_Val)
    
   ;; References
   write_references(diag_script, "A_rett_ri", "A_laue_ax", "D_0001", \
                     "", "P_crescendo")
   info_output("<<<<<<<< Leaving " + diag_script, verbosity, 4)

end

