;;#############################################################################
;; Multi-model mean change scled by global T change map with significance
;; Author: Ruth Lorenz (ETH, Switzerland)
;; CRESCENDO project
;;#############################################################################
;;
;; Description
;;    Annual mean or seasonal change sclaed by global T change per model
;;    multi-model mean with signficance in different periods
;;    (e.g. 2081-2100 and 2181-2200 with respect to 1986-2005)
;;
;; Required diag_script_info attributes (diagnostics specific)
;;    scenarios: list with scenarios to be included in the
;;               figure, e.g (/"rcp26","rcp45","rcp60","rcp85"/)
;;    periods: list with start years of periods to be included
;;    	       e.g. (/"2046","2081","2180"/)
;;    time_avg: list with seasons or annual to be included
;;              e.g (/"annualclim"/)
;;
;; Optional diag_script_info attributes (diagnostic specific)
;;    seasons: list with seasons index to be included if
;;             averagetime is "seasonalclim" (required for
;;             seasonalclim), DJF:0, MAM:1, JJA:2, SON:3
;;    iavmode: use IAV data from multi-model mean calculated by
;;             IAV_mmm_for_stippandhatch_ipcc12.ncl or include identical
;;             models to the ones plotted, calculated for each model on
;;             native grid with IAV_for_stippandhatch_ipcc12.ncl?
;;             default if not given is mmm, otherwise set iavmode = "each"
;;    grid:    grid for regridding (coarsest, finest, ref, 1x1, 2.5x2.5)
;;    percent: 0 or 1, difference in percent = 1, default = 0
;;    explicit_cn_levels: explicit levels for all contour plots
;;    max_vert: maximum number of plots in vertical
;;              default determined by number of scenarios
;;    max_hori: maximum number of plots in horizontal
;;              default determined by number of periods
;;    title: Figure title
;;    label: labels on top of each panel, e.g. scenario + rcp,
;;           loop to plot is seasons, scenarios, periods
;;    colormap: specify colormap for contour plots
;;    units: unit string in ncl formatting for legend title
;;
;; Required variable_info attributes (variable specific)
;;    att1: short description
;;    att2: short description
;;
;; Optional variable_info attributes (variable specific)
;;    att1: short description
;;    att2: short description
;;
;; Caveats
;;    Only handles one variable at a time
;;
;; Modification history
;;    20161027-A_lore_ru: written for ESMValTool
;;
;;#############################################################################

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

; A temporary file written by the invoking Python script
; Passes on a number of variables from Python to NCL
load "./interface_data/ncl.interface"

; Auxiliary NCL routines
load "./interface_scripts/auxiliary.ncl"
load "./interface_scripts/data_handling.ncl"
load "./interface_scripts/messaging.ncl"

load "./diag_scripts/lib/ncl/style.ncl"
load "./diag_scripts/lib/ncl/statistics.ncl"
load "./diag_scripts/lib/ncl/ensemble.ncl"

begin
    ;;##############################################################
    ;;# Fetch general parameters, set in namelist_collins13ipcc.xml#
    ;;# passed via environment variables by python code            #
    ;;##############################################################
    verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering " + diag_script, verbosity, 4)

    ;; 'diag_script' and 'variables' are fetched from the
    ;; above 'ncl.interface' file
    info_output("++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)
    info_output(diag_script + " (var: " + variables(0) + ")", verbosity, 1)
    info_output("++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)

    ;; 'models', 'variables' and 'field_types' are fetched from the
    ;; above 'ncl.interface' file
    dim_MOD = dimsizes(models@name)
    dim_VAR = dimsizes(variables)
    if (dim_VAR .eq. 2) then
        var0 = variables(0)
        var1 = variables(1)
    else
        var0 = variables(0)
        var1 = variables(0)
    end if
    field_type0 = field_types(0)

    ;;###########################################
    ;;# References Information                  #
    ;;###########################################
    write_references(diag_script, \  ; script name
                     "A_lore_ru", \  ; authors
                     "A_sedl_ja", \  ; contributors
                     "D_collins13ipcc",    \  ; diag_references
                     "",          \  ; obs_references
                     "P_cresce")   ; proj_references

    do idx_mod = 0, dim_MOD - 1
        data = read_data(idx_mod, var0, field_type0)
        temp = read_data(idx_mod, var1, field_type0)
    end do
end

; Load variable_info
begin
    vardeffile = "interface_data/" + var0 + "_info.tmp"
    loadscript(vardeffile)
end

begin
    ;; Get environment variables
    file_type = getenv("ESMValTool_output_file_type")
    if (ismissing(file_type)) then
        file_type = "PS"
    end if

    ;; Check required diag_script_info attributes
    req_atts = (/"scenarios", "periods", "time_avg"/)
    exit_if_missing_atts(diag_script_info, req_atts)

    ;; Basename of diag_script
    diag_script_base = basename(diag_script)

    ;; Output netcdf directory
    work_dir = getenv("ESMValTool_wrk_dir")
    netcdf_dir = get_output_dir(work_dir, diag_script_base)
end

begin
    ;;#############################################
    ;;# Get parameters from ./variable_defs/*.ncl #
    ;;# passed via the 'info' attribute           #
    ;;#############################################
    if(isvar("MyParam")) then
        delete(MyParam)
    end if
    if(isatt(variable_info, "long_name")) then
        MyParam = variable_info@long_name
        info_output(" MyParam = " + MyParam, verbosity, 5)
    else
        error_msg("f", diag_script, "", "missing variable attribute " + \
                  "'long_name'")
    end if

    ;;###########################################
    ;;# Get data and average time               #
    ;;###########################################
    ;; get data from first model
    imod = 0  ; NCL array indicies start from zero
    info_output("processing " + models@name(imod) + "_"\
                              + models@case_name(imod) + "_"\
                              + models@ensemble(imod), verbosity, 5)

    ;; See ./interface_scripts/data_handling.ncl
    A0 = extract_data(imod, data, -1, 0, 0)
    info_output(" ++++++++++++++ Variable summary of original data (A0): ", \
                verbosity, 5)
    if (verbosity .ge. 4) then
        printVarSummary(A0)
    end if

    ;; Check dimensions
    dims = getvardims(A0)
    if (dimsizes(dims) .lt. 2) then
        error_msg("f", diag_script, "", dimsizes(dims) + " dimensions, " + \
                  "need 2 or 3")
    end if
    idx = ind(dims .eq. "lat")
    if (ismissing(idx)) then
        error_msg("f", diag_script, "", "no lat dimension")
    end if
    idx = ind(dims .eq. "lon")
    if (ismissing(idx)) then
        error_msg("f", diag_script, "", "no lon dimension")
    end if

    ;; determine destination grid (regridding.ncl)
    index = ispan(0, dim_MOD - 1, 1)
    if (isatt(diag_script_info, "seasons")) then
        dim_seas = dimsizes(diag_script_info@seasons)
    else
        dim_seas = 1
    end if
    if (.not.(isatt(diag_script_info, "grid"))) then
        error_msg("w", diag_script, "", \
                  "Regridding method not defined. " + "Using 2.5x2.5")
        selected_grid = "2.5x2.5"
    else
        selected_grid = diag_script_info@grid
    end if
    info_output("Setting target grid:", verbosity, 2)
    if (selected_grid.eq."ref") then
        ; set reference model
        ref_ind = ind(models@name.eq.diag_script_info@ref_model(0))
        ref_model = diag_script_info@ref_model
        if (ismissing(ref_ind)) then
            if (diag_script_info@ref_model(0).eq."n.a.") then
                ref_ind = 0
            else
                error_msg("f", diag_script, "", \
                          "no adequate reference model provided")
            end if
        end if
        dest_grd = find_destination_grid(ref_ind, \
                                         var0, field_type0, "coarsest")
    else if (selected_grid.eq."1x1") then
        dest_grd     = new((/180, 360/), float)
        dest_grd!0   = "lat"
        dest_grd!1   = "lon"
        dest_grd&lat = flt2dble(fspan(-89.5,  89.5, 180))  ; CMOR standard
        dest_grd&lon = flt2dble(fspan(0.5, 359.5, 360))    ; CMOR standard
        dest_grd&lat@standard_name = "latitude"
        dest_grd&lon@standard_name = "longitude"
        dest_grd&lat@long_name     = "latitude"
        dest_grd&lon@long_name     = "longitude"
        dest_grd&lat@units         = "degrees_north"
        dest_grd&lon@units         = "degrees_east"
        dest_grd&lat@axis          = "Y"
        dest_grd&lon@axis          = "X"
    else if (selected_grid.eq."2.5x2.5") then
        dest_grd     = new((/72, 144/), float)
        dest_grd!0   = "lat"
        dest_grd!1   = "lon"
        dest_grd&lat = flt2dble(fspan(-88.75, 88.75, 72))  ; CMOR standard
        dest_grd&lon = flt2dble(fspan(1.25, 358.75, 144))    ; CMOR standard
        dest_grd&lat@standard_name = "latitude"
        dest_grd&lon@standard_name = "longitude"
        dest_grd&lat@long_name     = "latitude"
        dest_grd&lon@long_name     = "longitude"
        dest_grd&lat@units         = "degrees_north"
        dest_grd&lon@units         = "degrees_east"
        dest_grd&lat@axis          = "Y"
        dest_grd&lon@axis          = "X"
    else
        dest_grd = find_destination_grid(index, \
                                         var0, field_type0, selected_grid)
    end if
    end if
    end if
    annots = project_style(diag_script_info, "annots")
    Fill = default_fillvalue(typeof(A0))
    var_grd = new((/dim_MOD, dim_seas, dimsizes(dest_grd&lat), \
                   dimsizes(dest_grd&lon)/), typeof(A0), Fill)
    var_grd!0 = "models"
    var_grd!1 = "season"
    var_grd!2 = "lat"
    var_grd!3 = "lon"
    var_grd&models = annots
    if (isatt(diag_script_info, "seasons")) then
        var_grd&season = (/ispan(0, dim_seas - 1, 1)/)
    else
        var_grd&season = (/0/)
    end if
    var_grd&lat = dest_grd&lat
    var_grd&lon = dest_grd&lon
    if (isatt(diag_script_info, "iavmode")) then
        if (diag_script_info@iavmode .eq. "each") then
            iav = new((/dim_MOD, dim_seas, dimsizes(dest_grd&lat), \
                       dimsizes(dest_grd&lon)/), typeof(A0), Fill)
            copy_VarMeta(var_grd, iav)
        end if
    end if
    temp_avg = new((/dim_MOD/), typeof(A0), Fill)
    temp_avg!0 = "models"
    temp_avg&models = annots
    ;; Average over time (function in ./diag_scripts/lib/ncl/statistics.ncl)
    if ((diag_script_info@time_avg .eq. "seasonalclim") .and. \
        (dim_seas .ne. 4)) then
        tmp_seas = time_operations(A0, -1, -1, "average", \
                                   diag_script_info@time_avg, True)
        data1 = tmp_seas(diag_script_info@seasons, :, :)
    else
        data1 = time_operations(A0, -1, -1, "average", \
                                diag_script_info@time_avg, True)
    end if
    info_output(" ++++++++++++++ Variable summary of time averaged data " + \
                "(data1): ", verbosity, 5)
    if (verbosity .ge. 5) then
        printVarSummary(data1)
    end if
    rad = (4.0 * atan(1.0) / 180.0)
    do imod = 0, dim_MOD - 1
        if (imod .eq. 0) then
            A0_temp = extract_data(imod, temp, -1, 0, 0)
            temp_tavg = dim_avg_n_Wrap(A0_temp, 0)
            latw = cos(data1&lat * rad)
            temp_avg(imod) = wgt_areaave_Wrap(temp_tavg, latw, 1.0, 0)
            if (dim_seas .eq. 1) then
                var_grd(imod, 0, :, :) = linint2_Wrap(data1&lon, \
                                                      data1&lat, \
                                                      data1, True, \
                                                      dest_grd&lon, \
                                                      dest_grd&lat, 0)
            else
                var_grd(imod, :, :, :) = linint2_Wrap(data1&lon, \
                                                      data1&lat, \
                                                      data1, True, \
                                                      dest_grd&lon, \
                                                      dest_grd&lat, 0)
            end if
            delete(A0)
            delete(data1)
            delete(temp_tavg)
            delete(latw)
        else
            A1_temp = extract_data(imod, temp, -1, 0, 0)
            temp_tavg = dim_avg_n_Wrap(A1_temp, 0)
            latw = cos(temp_tavg&lat * rad)
            temp_avg(imod) = wgt_areaave_Wrap(temp_tavg, latw, 1.0, 0)
            A1 = extract_data(imod, data, -1, 0, 0)
            ;; Average over time,
            ;; function in ./diag_scripts/lib/ncl/statistics.ncl)
            if ((diag_script_info@time_avg .eq. "seasonalclim") .and. \
                (dim_seas .ne. 4)) then
                tmp_seas = time_operations(A1, -1, -1, "average", \
                                           diag_script_info@time_avg, True)
                data2 = tmp_seas(diag_script_info@seasons, :, :)
                delete(tmp_seas)
            else
                data2 = time_operations(A1, -1, -1, "average", \
                                        diag_script_info@time_avg, True)
            end if
            ;; Regrid to common grid
            if (dim_seas .eq. 1) then
                var_grd(imod, 0, :, :) = linint2_Wrap(data2&lon, \
                                                      data2&lat, \
                                                      data2, True, \
                                                      dest_grd&lon, \
                                                      dest_grd&lat, 0)
            else
                var_grd(imod, :, :, :) = linint2_Wrap(data2&lon, \
                                                      data2&lat, \
                                                      data2, True, \
                                                      dest_grd&lon, \
                                                      dest_grd&lat, 0)
            end if
            delete(A1)
            delete(A1_temp)
            delete(data2)
            delete(temp_tavg)
            delete(latw)
        end if
    end do
    delete(imod)
    if (verbosity .ge. 5) then
        printVarSummary(var_grd)
        printMinMax(var_grd, True)
    end if
    printVarSummary(temp_avg)
    printMinMax(temp_avg, True)
    ;;##############################################
    ;;# Calculate change from hist to periods      #
    ;;##############################################
    ;loop over rcps and periods
    dim_scen = dimsizes(diag_script_info@scenarios)
    dim_per = dimsizes(diag_script_info@periods)

    to_plot = new((/dim_per * dim_seas, dimsizes(dest_grd&lat), \
                   dimsizes(dest_grd&lon)/), typeof(var_grd), Fill)
    if (isatt(diag_script_info, "percent")) then
        percent = diag_script_info@percent
    else
        percent = 0
    end if
    to_plot_signif = new((/dim_per * dim_seas,  \
                          dimsizes(dest_grd&lat), dimsizes(dest_grd&lon)/), \
                         typeof(var_grd), Fill)
    if (var0 .eq. "mrsos" .or. var0 .eq. "mrro") then
        ; read mask to mask data where not defined (ocean and ice)
        in = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc", "r")
        msk_data = in->LSMASK
        msk = landsea_mask(msk_data, dest_grd&lat, dest_grd&lon)
    end if
    i = 0
    do seas = 0, dim_seas - 1  ; loop over seasons
        do per = 0, dim_per - 1  ; loop over periods
            diff_rcp_scal = new((/dim_scen, dimsizes(dest_grd&lat), \
                                dimsizes(dest_grd&lon)/), typeof(var_grd), \
                                Fill)
            do rcp = 0, dim_scen - 1  ; loop over rcps
                idx_rcp = ind((models@case_name .eq. \
                               diag_script_info@scenarios(rcp)) .and. \
                              (models@start_year .eq. \
                               diag_script_info@periods(per)))
                proj = var_grd(idx_rcp, seas, :, :)
                temp_proj = temp_avg(idx_rcp)
                ;; **cut low values ***
                proj = where(abs(proj) .gt. 1.e-14, proj, proj@_FillValue)
                ;; find historical runs from same models as in rcp
                idx_hist = find_model_index(models@name(idx_rcp), \
                                            models@name, "historical", \
                                            models@experiment, \
                                            models@ensemble(idx_rcp), \
                                            models@ensemble, 0)
                base = var_grd(idx_hist, seas, :, :)
                temp_base = temp_avg(idx_hist)
                ;; **cut low values ***
                base = where(abs(base) .gt. 1.e-14, base, base@_FillValue)
                ;; scale each model by global T change
                dim_mod = dimsizes(idx_rcp)
                if (dim_mod .lt. 2) then
                    error_msg("f", diag_script, "", \
                    "Only one model found, multi-model needs at least two.")
                end if
                var_diff_scal = new((/dim_mod, dimsizes(dest_grd&lat), \
                                    dimsizes(dest_grd&lon)/), typeof(var_grd),\
                                    Fill)
                do imod = 0, dim_mod - 1
                    if (percent .eq. 1) then
                        var_diff_scal(imod, :, :) = \
                            ((100 * (proj(imod, :, :) - base(imod, :, :))) / \
                             where(base(imod, :, :) .ne. 0.,\
                                   base(imod, :, :), base@_FillValue)) / \
                            (temp_proj(imod) - temp_base(imod))
                    else
                        var_diff_scal(imod, :, :) = (proj(imod, :, :) - \
                                                     base(imod, :, :)) / \
                                                    (temp_proj(imod) - \
                                                     temp_base(imod))
                    end if
                end do
                delete(base)
                delete(proj)
                delete(temp_proj)
                delete(temp_base)
                ;;  average over models
                diff_rcp_scal(rcp, :, :) = dim_avg_n_Wrap(var_diff_scal, 0)
                delete(idx_hist)
                delete(idx_rcp)
                delete(var_diff_scal)
            end do  ; rcp
            ;; average over rcps
            to_plot(i, :, :) = dim_avg_n_Wrap(diff_rcp_scal, 0)
            ;; calculate 95-percentile over all realizations
            diff_sort = dim_pqsort_n(diff_rcp_scal, 2, 0)
            ndiff = dim_num_n(.not.ismissing(diff_sort), 0)  ; ndiff(nlat,mlon)
            ind_95 = tointeger(floor(0.95 * ndiff) - 1)
            perc95 = new(dimsizes(ind_95), float)
            do ilat = 0, dimsizes(dest_grd&lat) - 1
                do ilon = 0, dimsizes(dest_grd&lon) - 1
                    index = ind_95(ilat, ilon)
                    perc95(ilat, ilon) = diff_sort(index(0), ilat, ilon)
                end do
            end do
            to_plot_signif(i, :, :) = where(abs(to_plot(i, :, :)) .gt. perc95,\
                                            1.0, 0)
            if (var0 .eq. "mrsos" .or. var0 .eq. "mrro") then
                to_plot(i, :, :) = mask(to_plot(i, :, :), msk .eq. 0, False)
                to_plot(i, :, :) = mask(to_plot(i, :, :), msk .eq. 4, False)
            end if
            i = i + 1
        end do  ; per
    end do  ; seas
    to_plot!0 = "panel"
    to_plot&panel = diag_script_info@label
    to_plot!1 = "lat"
    to_plot&lat = dest_grd&lat
    to_plot!2 = "lon"
    to_plot&lon = dest_grd&lon

    if (percent .eq. 1) then
        to_plot@units     = "% per ~F35~J~F~C"
    end if

    ;;###########################################
    ;;# Other Metadata: diag_script, var        #
    ;;###########################################
    ;; Add to to_plot, as attributes without prefix
    if (isatt(to_plot, "diag_script")) then  ; Add to existing entries
        tmp = to_plot@diag_script
        delete(to_plot@diag_script)
        to_plot@diag_script = array_append_record(tmp, (/diag_script/), 0)
        delete(tmp)
    else  ; Add as new attribute
        to_plot@diag_script = (/diag_script/)
    end if
    to_plot@var = var0  ; Overwrite existing entry
    if (isatt(variable_info, "long_name")) then
        to_plot@var_long_name = variable_info@long_name
    end if
    ;; check unit and adjust for plotting if necessary
    if (isatt(variable_info, "units")) then
        if ((var0 .eq. "mrro" .or. var0 .eq. "evspsbl") \
            .and. variable_info@units .eq. "kg m-2 s-1") then
            ;; convert [kg m-2 s-1] to [mm day-1]
            to_plot = to_plot * 60 * 60 * 24
            variable_info@units = "mm day-1"
            to_plot@var_units   = variable_info@units
        else if (var0 .eq. "mrsos" .and. variable_info@units .eq. "kg m-2") \
            then
            ;; convert [kg m-2] to [mm]
            variable_info@units = "mm"
            to_plot@var_units   = variable_info@units
        else if (var0 .eq. "hfls" .and. variable_info@units .eq. "W m-2") then
            ;; convert [W m-2] to [mm day-1]
            ;; ET[kg m**-2 s**-1] = LH / lambda
            ;;                    = LE / (2.5 * 10 ^ 6 [J kg**-1])
            ;; [kg m**-2 s**-1] * 24 * 60 * 60 = [mm / day]
            to_plot = to_plot / (2.5 * 10 ^ 6) * 24 * 60 * 60
            variable_info@units = "mm day-1"
            to_plot@var_units   = variable_info@units
        else if (variable_info@units .eq. "Pa") then
            ;; convert to hPa
            to_plot = to_plot / 100
            to_plot@var_units     = "hPa"
        else
            to_plot@var_units     = variable_info@units
        end if
        end if
        end if
        end if
    end if
    ;; swap downward/upward radiation fluxes according to IPCC Chapter12
    if ((var0 .eq. "rsut") .or. (var0 .eq. "rlut")) then
        to_plot = - to_plot
    end if
    copy_VarCoords(to_plot, to_plot_signif)

    to_plot_signif@diag_script = (/diag_script/)
    to_plot_signif@var = "signif"
    to_plot_signif@var_long_name = "significance: where average change " + \
                                   "larger than 95-percentile over all " + \
                                   "realizations"
    to_plot_signif@var_units     = 1

    if (verbosity .ge. 5) then
        printVarSummary(to_plot)
        printVarSummary(to_plot_signif)
    end if

    ;;###########################################
    ;;# Optional output to netCDF               #
    ;;###########################################
    ;Optional output (controlled by diag_script_info)
    if (isatt(diag_script_info, "ncdf")) then
        to_plot@ncdf = diag_script_info@ncdf
        ncdf_file = netcdf_dir + "/" + diag_script_info@ncdf
        ncdf_file@existing = "overwrite"
        ;; Function in ~/interface_scripts/auxiliary.ncl
        ncdf_outfile = ncdf_write(to_plot, ncdf_file)
        ncdf_file@existing = "append"
        ncdf_outfile = ncdf_write(to_plot_signif, ncdf_file)
    end if
end
