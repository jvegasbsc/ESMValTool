;;#############################################################################
;; CALCULATES TCR (Transient climate response)
;; Author: Ricarda Rettinger (DLR, Germany)
;; CRESCENDO project
;;#############################################################################
;;
;; Description
;;    To produce IPCC AR5 figure 9.42b (Flato et al., 2013) the Equilibrium Climate
;;    Sensitivity (ECS) is plotted against the Transcient Climate Response (TCR).
;;    Therefore, run ./diag_scripts/ECS.ncl before using this script.
;;    For TCR the 1pctCO2 experiment is used and it is calculated
;;    taking the difference between the start of the experiment and a 20-year
;;    period centred on the time of CO2 doubling.
;;
;; References
;;    Flato et al. (2013). Evaluation of Climate Models. In: Climate Change 2013:
;;    The Physical Science Basis. Contribution of Working Group I to the Fifth
;;    Assessment Report of the Intergovernmental Panel on Climate Change. Climate
;;    Change 2013, 5, 741-866.
;;
;; Required variable_info attributes (variable specific)
;;    none
;;
;; Required diag_script_info attributes (diagnostics specific)
;;   none
;;
;; Optional diag_script_info attributes (diagnostic specific)
;;   calcmm:     - False (default) : do not calculate multimodel mean
;;               - True : calculate multimodel mean
;;   ecs_file:   path relative to ESMValTool work_dir + file name of
;;               netCDF containing precalculated ECS values;
;;               default: "ECS/ECS.nc"
;;   styleset:   "CMIP5" (if not set, a color table and symbols for
;;               plotting will be created)
;;
;; Caveats
;;
;; Modification history
;;    20171013-A_laue_ax: * cleaned up header
;;                        * added meta data to plots for reporting
;;    20161010-A_rett_ri: written.
;;
;;#############################################################################
load "./interface_data/ncl.interface"
load "./interface_scripts/auxiliary.ncl"
;load "./interface_scripts/data_handling.ncl"
;load "./interface_scripts/messaging.ncl"

load "./diag_scripts/lib/ncl/scaling.ncl"
load "./diag_scripts/lib/ncl/statistics.ncl"
load "./diag_scripts/lib/ncl/latlon.ncl"
load "./diag_scripts/lib/ncl/style.ncl"

;load "./plot_scripts/ncl/xy_line.ncl"
load "./plot_scripts/ncl/scatterplot.ncl"

load "./diag_scripts/lib/ncl/meta_data.ncl"

begin

    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering " + diag_script, verbosity, 4)
    info_output("++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)
    info_output(diag_script + " (var: " + variables(0) + ")", verbosity, 1)
    info_output("++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)

    ;; References
     write_references(diag_script,      \  ; script name
                      "A_rett_ri",      \  ; authors
                      "A_laue_ax",      \  ; contributors
                      "D_flato13ipcc",  \  ; diag_references
                      "",               \  ; obs_references
                      (/"P_crescendo"/))   ; proj_references

    ;; basename and filetype
    diag_script_base = basename(diag_script)
    file_type = getenv("ESMValTool_output_file_type")
    if (ismissing(file_type)) then
        file_type = "ps"
    end if

    ;; output path
    out_path = get_output_dir(get_work_dir(), basename(diag_script))

    dim_MOD = dimsizes(models@name)
    dim_VAR = dimsizes(variables)
    var0 = variables(0)
    field_type0 = field_types(0)

    ;; Model indices
    do idx_mod = 0, dim_MOD - 1
        data_tas = read_data(idx_mod, var0, field_type0)
    end do
    delete(idx_mod)

    ;; Calculate MMM -> append empty field
    flag_mean = False
    mmm_ind = -1
    i = 1 ;; needed for data extraction
    if (isatt(diag_script_info, "calcmm")) then
        if (diag_script_info@calcmm) then
            i := 2
            flag_mean = True
            modelsatname = array_append_record(models@name, (/"multi-model-mean"/), 0)
            mmm_ind := dim_MOD
            dim_MOD := dim_MOD + 1
            delete(models@name)
            models@name = modelsatname
            delete(modelsatname)
        end if
    end if
    modelnames = models@name

end

begin
    vardeffile = "interface_data/" + var0 + "_info.tmp"
    loadscript(vardeffile)
end

begin

    ;; calculate TCR

    ;; empty index array
    TCR = new((/dim_MOD/),"float")

    ;; Extract model corresponding data without MMM!
    do imod = 0, dim_MOD - i
        var_tas = extract_data_extended(imod, var0, data_tas, -1, 0, 0)

        if (models@experiment(0) .eq. "historical") then
            avg_area = area_operations(var_tas, -90, 90, 0, 360, "average", True)
            TCR(imod) = avg(avg_area) -273.15
            TCR@var = "GMSAT"
            TCR@long_name = "Global mean surface temperature"
            TCR@standard_name = "GMSAT"
        else
            avg_area = area_operations(var_tas, -90, 90, 0, 360, "average", True)
            avg_tas_mon = month_to_annual(avg_area,1)
            TCR_i = (avg(avg_tas_mon(60:79)))-avg_tas_mon(0)
            TCR(imod) = TCR_i; - 273.15
            TCR@var = "TCR"
            TCR@long_name = "Transcient Climate Response"
            TCR@standard_name = "TCR"
            delete(avg_tas_mon)
            delete(TCR_i)
        end if
        delete(avg_area)
        delete(var_tas)
    end do
    delete(imod)

    ;; calculate MMM if calcmm = True
    if (isatt(diag_script_info, "calcmm")) then
        TCR(dim_MOD-1) = avg(TCR)
    end if

    ;; Write .nc file
    outfile = out_path + TCR@var + ".nc"
    TCR!0 = "model"
    TCR&model = models@name
    TCR@units = "°C"
   ; = nam0
    TCR@ncdf = out_path
    TCR@diag_script = diag_script

    ncdf_outfile = ncdf_write(TCR, outfile)
    info_output("ncdump " + outfile + " &", verbosity, 1)

    ; ========================================================================
    ; ==================== read precalculated ECS from file ==================
    ; ========================================================================

    work_dir = getenv("ESMValTool_wrk_dir")
    if (isatt(diag_script_info, "ecs_file")) then
        ecsfile = work_dir + "/" + diag_script_info@ecs_file
    else
        ecsfile = work_dir + "/ECS/ECS.nc"
    end if

    if (isfilepresent_esmval(ecsfile)) then
        info_output("Read in: " + ecsfile, verbosity, 1)
        ecs_tmp = addfile(ecsfile, "r")
        ecs_att = ecs_tmp->ecs ;; get ecs values
        ecs_tmp_model = tostring(ecs_tmp->model) ; get modelnames

    else
        error_msg("f", diag_script, "", "file with precalculated ECS data (" + \
            ecsfile + ") not found. Please run ./diag_scripts/ECS.ncl.")
    end if

    ; now copy the precalculated ECS values read from the external file to the
    ; array ECS
    dim_MOD = dimsizes(models@name)
    ECS = new(dim_MOD, float)
    ECS = ecs_att@_FillValue
    ECS@_FillValue = ecs_att@_FillValue
    ECS@units = ecs_att@units
    ECS@standard_name = ecs_att@standard_name

    exit_flag = False
    wECS_idx = new((/dim_MOD/),float)
    wECS_idx@_FillValue = ecs_att@_FillValue

    do i = 0, dim_MOD - 1
         if (i.eq.mmm_ind) then
             wECS_idx(i) =  1
             continue
         end if
        idx = ind(models@name(i).eq.ecs_tmp_model); &model
        if (ismissing(idx)) then
            print("no pre-calcuated ECS value for model " + models@name(i) + " available")
           ;  exit_flag = True
        else
            ECS(i) = ecs_att(idx)
            wECS_idx(i) =  1
            print("pre-calcuated ECS value for model " + models@name(i) + " = " + ECS(i))
        end if
    end do
    delete(i)

    if (exit_flag) then
        error_msg("f", diag_script, "", "precalculated ECS data (" + \
            ecsfile + ") not available for all models.")
    end if

    if (flag_mean) then
        ECS(mmm_ind) = avg(ECS)
    end if

    ; remove models from ECS/TCR if ECS value is missing
    org_wECS_idx = wECS_idx
    wECS_idx := ind(wECS_idx.eq.1)
    ECS := ECS(wECS_idx)
    TCR := TCR(wECS_idx)
    dim_MOD := dimsizes(ECS)

    ; update mmm_ind if needed

    if (flag_mean) then
        i = mmm_ind - 1
        mmm_ind = toint(sum(org_wECS_idx(0:i)))
    end if

    delete(org_wECS_idx)

    ;############################
    ;            Plot
    ;############################

    ; output directory
    plot_dir = getenv("ESMValTool_plot_dir")
    outputdir = get_output_dir(plot_dir, diag_script_base)

    if (isatt(diag_script_info, "legend_outside")) then
        legend_outside = diag_script_info@legend_outside
    else
        legend_outside = False
    end if

    ;; Define workstation
    if (models@experiment(0) .eq. "historical") then
        outfile = outputdir + "Figure_9.42a"
    else
        outfile = outputdir + "Figure_9.42b"
    end if
    wks = gsn_open_wks(file_type, outfile)
    wks@legendfile = outfile + "_legend"
    wks@fullname = outfile + "." + file_type

    if (isatt(diag_script_info, "styleset")) then
        colortab = project_style(diag_script_info, "colors")
        markertab = project_style(diag_script_info, "markers")

    else
        colortab = (/"(/0.00, 0.00, 0.59/)", "(/0.00, 0.39, 1.00/)", \
                     "(/0.20, 1.00, 1.00/)", "(/0.20, 0.88, 0.00/)", \
                     "(/1.00, 0.88, 0.00/)", "(/1.00, 0.59, 0.00/)", \
                     "(/1.00, 0.20, 0.00/)", "(/0.59, 0.00, 0.00/)", \
                     "(/0.78, 0.00, 0.78/)", "(/0.59, 0.00, 0.59/)", \
                     "(/0.90, 0.90, 0.90/)", "(/0.70, 0.70, 0.70/)", \
                     "(/0.50, 0.50, 0.50/)", "(/0.30, 0.30, 0.30/)"/)
        markertab = (/16, 4, 5, 0/)
    end if

    colors  = new(dim_MOD, string)
    markers = new(dim_MOD, integer)

    ; create new marker: filled star
    mstring = "z"
    fontnum = 35
    size    = 1.5
    angle   = 0.0
    new_index = NhlNewMarker(wks, mstring, fontnum, 0.0, 0.0, 1.0, size, angle)

    if (isatt(diag_script_info, "styleset")) then
        colors = colortab(wECS_idx)
        markers = markertab(wECS_idx)
        if (flag_mean) then
            colors(mmm_ind)  = "(/0.00, 0.00, 0.00/)"  ; black
            markers(mmm_ind) = new_index
        end if
     else
        i = 0
        idx = 0
        nmod = dim_MOD
        do while (i.lt.nmod)
            do n = 0, (nmod - 1) / dimsizes(colortab)
                colors(i)  = colortab(idx)
                markers(i) = markertab(n)
                i = i + 1
                if (i.ge.nmod) then
                    break
                end if
            end do
            idx = idx + 1
            if (idx.ge.dimsizes(colortab)) then
                idx = dimsizes(colortab) - 1
            end if
        end do
    end if

    ;; collect data
    data_arr = new((/2, dim_MOD /), float)
    data_arr!0 = "statistic"
    data_arr!1 = "model"
    data_arr&statistic = (/"ECS",TCR@var/)
    data_arr&model = models@name(wECS_idx)
    data_arr(0, :) = (/ECS/)
    data_arr(1, :) = (/TCR/)
    data_arr@units = (/"~S~o~N~C","~S~o~N~C"/)
    data_arr@legend_outside = legend_outside
    data_arr@colors = colors
    data_arr@markers = markers
    data_arr@annots = models@name(wECS_idx)
    data_arr@res_trXMinF = min(ECS) - 0.5
    data_arr@res_trYMinF = min(TCR) - 0.5
    data_arr@res_trXMaxF = max(ECS) + 0.5
    data_arr@res_trYMaxF = max(TCR) + 0.5

    ; draw scatter plot
    data_arr@res_tiMainString = ""
    plot = scatterplot_markers(wks, data_arr, var0)

    ; add regression line (drawn in the background)
    ; (can only be calculated if there are at least 3 models)
    if (dim_MOD.ge.3) then
        rc = regline(ECS, TCR)
        ; draw only if regression coefficient is greater/less than 0.1/-0.1
        if ((rc .gt. 0.1).or.(rc .lt. -0.1)) then
            xmin = min(ECS)
            xmax = max(ECS)
            xdelta = xmax - xmin
            x0 = xmin - 0.5 * xdelta
            x1 = xmax + 0.5 * xdelta
            x = (/x0, x1/)
            y = (/x0 * rc + rc@yintercept, x1 * rc + rc@yintercept/)
            res = True
            res@gsLineColor = (/0.00, 0.00, 0.00/)
            res@gsLineThicknessF = 4.0
            res@tfPolyDrawOrder = "Draw"
            regl = gsn_add_polyline(wks, plot, x, y, res)
        end if
    end if

    draw(plot)
    frame(wks)

    ; ========================================================================

    ; -------------------------------------------------------------------------------
    ; add meta data to plot (for reporting)

    climofiles = new(dim_MOD + 1, string)

    alltags = array_append_record(tags, (/"PT_scatter", "DM_global"/), 0)
    id = diag_script + "_TCR"
    contrib_authors = (/"A_rett_ri", "A_laue_ax"/)

    caption = "Equilibrium climate sensitivity against transient climate " \
            + "response (TCR). The ECS and TCR information are based on " \
            + "Andrews et al. (2012) and Forster et al. (2013) and updated " \
            + "from the CMIP5 archive. Similar to Flato et al. 2013, fig. " \
            + "9.42b"

    idx = 0
    do i = 0, dim_MOD - 1
        if (i.eq.mmm_ind) then
            continue
        end if

        climofiles(idx) = interface_get_inpaths(i) + "/" \
            + interface_get_infile(variables(0), field_types(0), i)

        idx = idx + 1
    end do

    climofiles(idx) = ecsfile

    vars = array_append_record(variables, (/"ECS"/), 0)

    ESMValMD(wks@fullname, alltags, caption, id, vars, \
             models@name(wECS_idx), climofiles(0:idx), diag_script, contrib_authors)
    ; -------------------------------------------------------------------------------

    info_output(" Wrote " + outfile + "." + str_lower(file_type) + " & " , verbosity, 1)
    info_output(">>>>>>>> Leaving " + diag_script, verbosity, 4)

end
