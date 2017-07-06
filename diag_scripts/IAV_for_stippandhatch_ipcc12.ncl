;;#############################################################################
;; INTERANNUAL VARIABILITY FOR STIPPLING SIGNIFICANCE IPCC CH12
;; Author: Ruth Lorenz (ETH Zurich, Switzerland)
;; CRESCENDO project
;;#############################################################################
;;
;; Description
;;    Calculate Interannual variability from piControl runs for plotting
;;    significance with stippling and hatching, save data for each model
;;    on native grid.
;;    Can either be calculated over full time period of piControl run
;;    (no periodlength given)
;;    or first averaged annually/seasonally over periodlength and then standard
;;    deviation calculated based on averaged periods.
;;
;; Required diag_script_info attributes (diagnostics specific)
;;    time_avg: time period to calculate IAV, e.g. annualclim, seasonalclim
;;    overwrite: True or False, if False files which already exist are not
;;               recalculated
;;
;; Optional diag_script_info attributes (diagnostic specific)
;;    periodlength: length of periods to calculate IAV accross, depends on
;;                  period lengths used in
;;                  map_diff_mmm_stippandhatch_ipcc12.ncl if not given whole
;;                  time period calculated at once
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
;;    List possible caveats or limitations of this diagnostic
;;    Features to-be-implemented shall also be mentioned here
;;
;; Modification history
;;    20161219-A_lore_ru: remove seasonal cycle before std if seasonal
;;    20161024-A_lore_ru: adapted to ESMValTool
;;    YYYYMMDD-A_sedl_ja: written for IPCC AR5.
;;
;;#############################################################################

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/calendar_decode2.ncl"

; A temporary file written by the invoking Python script
; Passes on a number of variables from Python to NCL
load "./interface_data/ncl.interface"

; Auxiliary NCL routines
load "./interface_scripts/auxiliary.ncl"
load "./interface_scripts/data_handling.ncl"
load "./interface_scripts/messaging.ncl"

load "./diag_scripts/lib/ncl/statistics.ncl"

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
    var0 = variables(0)
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
    end do
end

; Load variable_info
begin
    vardeffile = "interface_data/" + var0 + "_info.tmp"
    loadscript(vardeffile)
end
begin
    ;; Check required diag_script_info attributes
    req_atts = (/"time_avg", "overwrite"/)
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
    if (isvar("MyParam")) then
        delete(MyParam)
    end if
    if (isatt(variable_info, "long_name")) then
        MyParam = variable_info@long_name
        info_output(" MyParam = " + MyParam, verbosity, 5)
    else
        error_msg("f", diag_script, "", "missing variable attribute " + \
                  "'long_name'")
    end if

end

begin
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
    rank = dimsizes(dims)
    print("rank: " + rank)
    if (dimsizes(dims) .lt. 3) then
        error_msg("f", diag_script, "", dimsizes(dims) + " dimensions, " + \
                  "need 3 or 4")
    end if
    idx = ind(dims .eq. "lat" .or. dims .eq. "rlat" .or. dims .eq. "j")
    if (ismissing(idx)) then
        error_msg("f", diag_script, "", "no lat dimension")
    end if
    lat_name = dims(idx)

    idx = ind(dims .eq. "lon" .or. dims .eq. "rlon" .or. dims .eq. "i")
    if (ismissing(idx)) then
        error_msg("f", diag_script, "", "no lon dimension")
    end if
    lon_name = dims(idx)

    if (rank .eq. 4) then
        idx = ind(dims .eq. "lev" .or. dims .eq. "plev")
        if (ismissing(idx)) then
            error_msg("f", diag_script, "", "no level dimension")
        end if
        lev_name = dims(idx)
    end if

    if (isatt(diag_script_info, "seasons")) then
        dim_seas = dimsizes(diag_script_info@seasons)
    else
        dim_seas = 1
    end if
    Fill = default_fillvalue(typeof(A0))

    do imod = 0, dim_MOD - 1
        outfile = "IAV_" + models@experiment(imod) + "_" + models@name(imod) +\
                  "_" + models@ensemble(imod) + "_" + field_type0 + "_" + \
                  var0 + "_" + diag_script_info@time_avg + ".nc"
        file_exist = isfilepresent(outfile)
        if (file_exist .and. diag_script_info@overwrite .eq. False) then
            continue
        end if
        if (imod .ne. 0) then
            A0 = extract_data(imod, data, -1, 0, 0)
            dims = getvardims(A0)
            idx = ind(dims .eq. "lat" .or. dims .eq. "rlat" .or. dims .eq. "j")
            if (ismissing(idx)) then
                error_msg("f", diag_script, "", "no lat dimension")
            end if
            lat_name = dims(idx)
            idx = ind(dims .eq. "lon" .or. dims .eq. "rlon" .or. dims .eq. "i")
            if (ismissing(idx)) then
                error_msg("f", diag_script, "", "no lon dimension")
            end if
            lon_name = dims(idx)
        end if
        ;; Calculate IAV (function in ./diag_scripts/lib/ncl/statistics.ncl)
        ;; cut the first 100 yr
        time = A0&time
        start_yr = cd_calendar(time(0), 0)
        ind_end = dimsizes(time) - 1
        end_yr = cd_calendar(time(ind_end), 0)
        if ((end_yr(0, 0) - start_yr(0, 0) + 1) .lt. 500) then
            print("warning: Control run is less than 500 years.")
        end if
        new_start_yr = toint((start_yr(0, 0) + 100))
        end_yr_int = toint(end_yr(0, 0))
        if (isatt(diag_script_info, "periodlength")) then
            length_of_period = toint(diag_script_info@periodlength)
            nr_periods = toint(floor((dimsizes(time) - 12 * 100) / \
                                     (length_of_period * 12.)))
            yr_possible = (dimsizes(time) - 12 * 100) / (length_of_period * 12)
            rest = yr_possible - nr_periods
            start_yrs = ispan(new_start_yr, end_yr_int, length_of_period)
            if (rank .eq. 3) then
                data_tmp =  new((/nr_periods, dim_seas, \
                                 dimsizes(A0&$lat_name$), \
                                 dimsizes(A0&$lon_name$)/), typeof(A0), Fill)
            else if (rank .eq. 4) then
                data_tmp =  new((/nr_periods, dim_seas,  \
                                 dimsizes(A0&$lev_name$), \
                                 dimsizes(A0&$lat_name$), \
                                 dimsizes(A0&$lon_name$)/), typeof(A0), Fill)
            end if
            end if
            if (diag_script_info@time_avg .eq. "seasonalclim") then
                ;; remove seasonal cycle
                deseas_data = rmMonAnnCycTLL(A0)
                delete(A0)
                A0 = deseas_data
            end if
            do per = 0, nr_periods - 1
                if ((rest .gt. 0) .and. (per .eq. nr_periods - 1)) then
                    new_end_yr = end_yr_int
                else
                    new_end_yr = toint(start_yrs(per) + length_of_period - 1.0)
                end if
                ;; calculate seasonal/annual averages over periods
                if (rank .eq. 3) then
                    if (dim_seas .eq. 1) then
                        data_tmp(per, 0, :, :) = \
                            time_operations(A0, start_yrs(per), new_end_yr, \
                                            "average", \
                                            diag_script_info@time_avg, True)
                    else
                        data_tmp(per, :, :, :) = \
                            time_operations(A0, start_yrs(per), new_end_yr, \
                                            "average", \
                                            diag_script_info@time_avg, True)
                    end if
                else if (rank .eq. 4) then
                    if (dim_seas .eq. 1) then
                        data_tmp(per, 0, :, :, :) = \
                            time_operations(A0, tart_yrs(per), new_end_yr, \
                                            "average", \
                                            diag_script_info@time_avg, True)
                    else
                        data_tmp(per, :, :, :, :) = \
                            time_operations(A0, start_yrs(per), new_end_yr, \
                                            "average", \
                                            diag_script_info@time_avg, True)
                    end if
                end if
                end if
            end do  ; end per
            data_dtr = dtrend_quadratic_msg_n(data_tmp, False, False, 0)
            data1 = dim_stddev_n(data_dtr, 0)
            data1!0 = "season"
            data1&season = (/ispan(0, dim_seas - 1, 1)/)
            if (rank .eq. 3) then
                data1!1 = lat_name
                data1!2 = lon_name
            else if (rank .eq. 4) then
                data1!1 = lev_name
                data1!2 = lat_name
                data1!3 = lon_name
                data1&$lev_name$ = data_tmp&$lev_name$
            end if
            end if
            data1&$lat_name$ = data_tmp&$lat_name$
            data1&$lon_name$ = data_tmp&$lon_name$
            delete(start_yrs)
            delete(data_tmp)
            delete(data_dtr)
        else
            data1 = interannual_variability(A0, new_start_yr, end_yr_int, \
                                            diag_script_info@time_avg, \
                                            "quadratic")
        end if
        info_output(" ++++++++++++++ Variable summary of IAV annual data " + \
                    "(data1): ", verbosity, 5)
        if (verbosity .ge. 5) then
            printVarSummary(data1)
        end if
        delete(A0)
        if (diag_script_info@time_avg .eq. "seasonalclim") then
            delete(deseas_data)
        end if
        delete(start_yr)
        delete(end_yr)
        delete(time)
        ;;###########################################
        ;;# Output to netCDF                        #
        ;;###########################################
        data1@diag_script = (/diag_script/)
        data1@var = "iav"
        if (isatt(variable_info, "long_name")) then
            data1@var_long_name = "iav of " + variable_info@long_name
        end if
        if (isatt(variable_info, "units")) then
            data1@var_units     = variable_info@units
        end if
        ;; Function in ~/interface_scripts/auxiliary.ncl
        ncdf_outfile = ncdf_write(data1, netcdf_dir + outfile)
        delete(data1)
    end do  ; imod
end
