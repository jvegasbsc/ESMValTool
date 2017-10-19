;;#############################################################################
;; SNOW AREA CHANGE
;; Author: Ruth Lorenz (ETH Zurich, Switzerland)
;; CRESCENDO project
;;#############################################################################
;;
;; Description
;;    Calculate Snow area in a region (e.g Northern Hemisphere)
;;    and season (e.g. NH spring March & April)
;;    relative to reference period (e.g 1986-2005) as in Fig. 12.32 of IPCC AR5
;;
;; Required diag_script_info attributes (diagnostics specific)
;;    scenarios: which scenarios included in figure
;;    syears : start years in time periods
;;             (e.g. start of historical period and rcps)
;;    eyears : end years in time periods
;;             (end year of historical runs and rcps)
;;    begin_ref_year: start year of reference period (e.g. 1986)
;;    end_ref_year: end year of reference period (e.g 2005)
;;    months: months in year to calculate snow area from
;;
;; Optional diag_script_info attributes (diagnostic specific)
;;    ncdf : enable to output to netCDF; either use "default"
;;           or give a full file name
;;    ts_minlat: minimum latitude in case not global e.g(/0/)
;;    ts_maxlat: maximum latitude in case not global e.g.(/90/)
;;    ts_minlon: minimum longitude in case not global e.g.(/0/)
;;    ts_maxlon: maximum longitude in case not global e.g.(/360/)
;;    colormap: choose file to create colormap from
;;    yaxis: y-axis title
;;    ymin: minimim value on y-axis
;;    ymax: maximum value on y-axis
;;    title: title in plot
;;    label: labels to use in legend depending on scenarios
;;
;; Required variable_info attributes (variable specific)
;;    long_name: longer variable name
;;    units: variable units
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
;;    20170120-A_lore_ru: written.
;;
;;#############################################################################

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

; A temporary file written by the invoking Python script
; Passes on a number of variables from Python to NCL
load "./interface_data/ncl.interface"

; Auxiliary NCL routines
load "./interface_scripts/auxiliary.ncl"
load "interface_scripts/data_handling.ncl"
load "./interface_scripts/messaging.ncl"
load "diag_scripts/lib/ncl/latlon.ncl"
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
                     "P_crescendo")   ; proj_references

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
    ;; Get environment variables
    file_type = getenv("ESMValTool_output_file_type")
    if (ismissing(file_type)) then
        file_type = "PS"
    end if

    ;; Check required diag_script_info attributes
    req_atts = (/"scenarios", "syears", "eyears", "begin_ref_year", \
                "end_ref_year", "months"/)
    exit_if_missing_atts(diag_script_info, req_atts)

    ;; Basename of diag_script
    diag_script_base = basename(diag_script)

    ;; Output netcdf directory
    work_dir = getenv("ESMValTool_wrk_dir")
    netcdf_dir = get_output_dir(work_dir, diag_script_base)

    ;; Output plot directory
    plot_dir = getenv("ESMValTool_plot_dir")
    output_dir = get_output_dir(plot_dir, diag_script_base)
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
    ;;# Get data and average annual and globally#
    ;;###########################################
    ;; get data from first model for historical period
    imod = 0  ; NCL array indices start from zero
    info_output("processing " + models@name(imod) + "_" \
                              + models@case_name(imod) + "_" \
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

    ;; how many months included in analysis?
    months = diag_script_info@months
    dim_SEAS = dimsizes(months)

    ;; Extract data for time
    ;; (function in ./diag_scripts/lib/ncl/statistics.ncl)
    data3 = time_operations(A0, -1, -1, "extract", [months(0)], False)
    dsize = dimsizes(data3)
    data0 = new((/dim_SEAS, dsize(0), dsize(1), dsize(2)/), float)
    data0(0, :, :, :) = data3
    do m = 1, dim_SEAS - 1
        data0(m, :, :, :) = time_operations(A0, -1, -1, "extract", \
                                            [months(m)], False)
    end do
    data0!0 = "seas"
    data0&seas = months
    date    = cd_calendar(data3&time, 0)
    years = date(:, 0)
    nyears = dimsizes(years)
    nlat = dimsizes(A0&lat)
    nlon = dimsizes(A0&lon)
    data1 = dim_avg_n_Wrap(data0, 0)

    ;; if variable is "snw" we need to calculate the snow area using sftgif
    ;; and sftlf, files need to be provided as reformatted data
    ;; (trick ESMValTool into providing those files), if snw_area already
    ;; derived using variable_defs/snw_area.ncl nothing needs to be done
    if (var0 .eq. "snw") then
        inpath = interface_get_inpaths(imod)
        infile = inpath + "/" + models@project(imod) + "_fx_historical_" + \
            models@name(imod) + "_r0i0p0_F2Ms_sftgif_0-0.nc"

        f = addfile(infile, "r")
        ;; read sftgif (fraction of grid covered with land ice)
        sftgif = f->sftgif
        sftgif_0 = sftgif * (- 0.01)
        sftgif_1 = sftgif_0 + 1

        ;; read sftlf (land sea mask)
        infile = inpath + "/" + models@project(imod) + "_fx_historical_" + \
            models@name(imod) + "_r0i0p0_F2Ms_sftlf_0-0.nc"
        f = addfile(infile, "r")
        sftlf = f->sftlf
        sftifl = sftgif_1 * sftlf

        tmp2 = data1 / 2.
        tmp3 = tmp2 * 5.10067e14
        ;; where > 5
        tmp4 = where(tmp3 .gt. 5., tmp3, sftgif@_FillValue) / 100.
        data2 = new((/nyears, nlat, nlon/), typeof(A0), A0@_FillValue)
        do yr = 0, nyears - 1
            data2(yr, :, :) = tmp4(yr, :, :) * sftifl
        end do
        delete([/tmp2, tmp3, tmp4/])
        ;; Add coordinates
        data2!0 = "year"
        data2&year = date(:, 0)
        data2!1 = "lat"
        data2&lat = A0&lat
        data2!2 = "lon"
        data2&lon = A0&lon
    else
        data2 = data1
    end if
    info_output(" ++++++++++++++ Variable summary of time averaged data " + \
                "(data2): ", verbosity, 5)
    if (verbosity .ge. 5) then
        printVarSummary(data2)
    end if
    delete(data3)
    dim_file = dimsizes(data2)  ; file size of historical file
    ;; Average over globe/region
    ;; (function in ./diag_scripts/lib/ncl/latlon.ncl)
    if (isatt(diag_script_info, "ts_minlat")) then
        data1_glob = area_operations(data2, diag_script_info@ts_minlat, \
                                     diag_script_info@ts_maxlat, \
                                     diag_script_info@ts_minlon, \
                                     diag_script_info@ts_maxlon, "average", \
                                     True)
    else
        data1_glob = coswgt_areaave(data2)
    end if
    info_output(" ++++++++++++++ Variable summary of area averaged data " + \
                "(data1_glob): ", verbosity, 5)

    if (verbosity .ge. 5) then
        printVarSummary(data1_glob)
    end if

    ;; How many historical model runs?
    ind_hist = ind(models@case_name .eq. "historical")
    dim_MOD_hist = dimsizes(models@name(ind_hist))
    ref_data2 = new((/dim_MOD_hist/), float, data1_glob@_FillValue)
    data2_rel = new((/dim_MOD_hist, dim_file(0)/), float, \
                    data1_glob@_FillValue)

    ;; How many scenarios?
    dim_scen = dimsizes(diag_script_info@scenarios)
    ;; start and end years of timeperiods
    syears = diag_script_info@syears
    eyears = diag_script_info@eyears
    nperiods = dimsizes(syears)
    nyears_full = toint(eyears(nperiods - 1)) - toint(syears(0)) + 1
    nyearsp1 =  toint(eyears(nperiods - 1)) - toint(eyears(0))

    ;; How many rcp2.6 model runs?
    ind_26 = ind(models@case_name .eq. diag_script_info@scenarios(0))
    dim_MOD_26 = dimsizes(models@name(ind_26))
    ;; indicate array for scenarios
    data26_rel = new((/dim_MOD_26, nyearsp1/), float, \
                     data1_glob@_FillValue)

    ;; How many rcp4.5 model runs?
    ind_45 = ind(models@case_name .eq. diag_script_info@scenarios(1))
    dim_MOD_45 = dimsizes(models@name(ind_45))
    data45_rel = new((/dim_MOD_45, nyearsp1/), float, \
                     data1_glob@_FillValue)

    ;; How many rcp6.0 model runs?
    ind_60 = ind(models@case_name .eq. diag_script_info@scenarios(2))
    dim_MOD_60 = dimsizes(models@name(ind_60))
    data60_rel = new((/dim_MOD_60, nyearsp1/), float, \
                     data1_glob@_FillValue)

    ;; How many rcp8.5 model runs?
    ind_85 = ind(models@case_name .eq. diag_script_info@scenarios(3))
    dim_MOD_85 = dimsizes(models@name(ind_85))
    data85_rel = new((/dim_MOD_85, nyearsp1/), float, \
                     data1_glob@_FillValue)

    i = 0
    i26 = 0
    i45 = 0
    i60 = 0
    i85 = 0
    do imod = 0, dim_MOD_hist - 1
        if (imod .eq. 0) then
            ;; Add coordinates
            data2!0 = "year"
            data2&year = date(:, 0)
            data2!1 = "lat"
            data2&lat = A0&lat
            data2!2 = "lon"
            data2&lon = A0&lon
            delete(A0)
            delete(data0)
        else
            A0 = extract_data(imod, data, -1, 0, 0)
            ;; Extact over time
            ;; (function in ./diag_scripts/lib/ncl/statistics.ncl)
            data3 = time_operations(A0, -1, -1, "extract", [months(0)], False)
            dsize = dimsizes(data3)
            data0 = new((/dim_SEAS, dsize(0), dsize(1), dsize(2)/), float)
            data0(0, :, :, :) = data3
            do m = 1, dim_SEAS - 1
                data0(m, :, :, :) = time_operations(A0, -1, -1, "extract", \
                                                    [months(m)], False)
            end do
            data0!0 = "seas"
            data0&seas = months
            date    = cd_calendar(data3&time, 0)
            years = date(:, 0)
            nyears = dimsizes(years)
            nlat = dimsizes(A0&lat)
            nlon = dimsizes(A0&lon)
            data1 = dim_avg_n_Wrap(data0, 0)
            if (var0 .eq. "snw") then
                inpath = interface_get_inpaths(imod)
                infile = inpath + "/" + models@project(imod) + \
                    "_fx_historical_" + models@name(imod) + \
                    "_r0i0p0_F2Ms_sftgif_0-0.nc"
                f = addfile(infile, "r")
                ;; read sftgif (fraction of grid covered with land ice)
                sftgif = f->sftgif
                sftgif_0 = sftgif * (- 0.01)
                sftgif_1 = sftgif_0 + 1

                ;; read sftlf (land sea mask)
                infile = inpath + "/" + models@project(imod) + \
                    "_fx_historical_" + models@name(imod) + \
                    "_r0i0p0_F2Ms_sftlf_0-0.nc"
                f = addfile(infile, "r")
                sftlf = f->sftlf
                sftifl = sftgif_1 * sftlf

                tmp2 = data1 / 2.
                tmp3 = tmp2 * 5.10067e14
                ;; where > 5
                tmp4 = where(tmp3 .gt. 5., tmp3, sftgif@_FillValue) / 100.
                data2 = new((/nyears, nlat, nlon/), typeof(A0), A0@_FillValue)
                do yr = 0, nyears - 1
                    data2(yr, :, :) = tmp4(yr, :, :) * sftifl
                end do
                ;; Add coordinates
                data2!0 = "year"
                data2&year = date(:, 0)
                data2!1 = "lat"
                data2&lat = A0&lat
                data2!2 = "lon"
                data2&lon = A0&lon
                delete([/tmp2, tmp3, tmp4/])
            else
                data2 = data1
                ;; Add coordinates
                data2!0 = "year"
                data2&year = date(:, 0)
            end if
            delete([/A0, data0, data3/])
        end if
        ;; Average over globe (function in
        ;; ./diag_scripts/lib/ncl/statistics.ncl)
        if (isatt(diag_script_info, "ts_minlat")) then
            data2_glob = area_operations(data2, diag_script_info@ts_minlat, \
                                         diag_script_info@ts_maxlat, \
                                         diag_script_info@ts_minlon, \
                                         diag_script_info@ts_maxlon, \
                                         "average", True)
        else
            data2_glob = coswgt_areaave(data2)
        end if
        ind_start = ind(data2&year .eq. diag_script_info@begin_ref_year)
        ind_end = ind(data2&year .eq. diag_script_info@end_ref_year)
        ref_data2(imod) = dim_avg(data2_glob(ind_start : ind_end))
        dim_data2g = dimsizes(data2_glob)
        if (dim_data2g .ne. dim_file(0)) then
            tdiff = dim_file(0) - dim_data2g
            data2_rel(imod, 0 : tdiff - 1) = data2_glob@_FillValue
            data2_rel(imod, tdiff :) =  (data2_glob / ref_data2(imod)) - 1
        else
            data2_rel(imod, :) = (data2_glob / ref_data2(imod)) - 1
        end if
        ;; clean up
        delete([/data1, data2, data2_glob/])
        delete([/date, years, nyears/])

        do s = 0, dim_scen - 1
            ;; find all other runs from this model and calculate relative ts
            imodrcp = ind(models@name .eq. models@name(imod) .and. \
                          models@case_name .eq. diag_script_info@scenarios(s))
            if (.not.all(ismissing(imodrcp))) then
                if (dimsizes(imodrcp) .gt. 1) then
                    do ens = 0, dimsizes(imodrcp) - 1
                        A1 = extract_data(imodrcp(ens), data, -1, 0, 0)
                    end do
                else
                    A1 = extract_data(imodrcp, data, -1, 0, 0)
                end if
                ;; Extract over time
                ;; (function in ./diag_scripts/lib/ncl/statistics.ncl)
                data3 = time_operations(A1, -1, -1, "extract", [months(0)], \
                                        False)
                dsize = dimsizes(data3)
                data0 = new((/dim_SEAS, dsize(0), dsize(1), dsize(2)/), float)
                data0(0, :, :, :) = data3
                do m = 1, dim_SEAS - 1
                    data0(m, :, :, :) = time_operations(A1, -1, -1, "extract",\
                                                        [months(m)], False)
                end do
                data0!0 = "seas"
                data0&seas = months
                date    = cd_calendar(data3&time, 0)
                years = date(:, 0)
                nyears = dimsizes(years)
                nlat = dimsizes(A1&lat)
                nlon = dimsizes(A1&lon)
                data1 = dim_avg_n_Wrap(data0, 0)

                if (var0 .eq. "snw") then
                    tmp2 = data1 / 2.
                    tmp3 = tmp2 * 5.10067e14
                    ;; where > 5
                    tmp4 = where(tmp3 .gt. 5., tmp3, sftgif@_FillValue) / 100.
                    data2 = new((/nyears, nlat, nlon/), typeof(A1), \
                                A1@_FillValue)
                    do yr = 0, nyears - 1
                        data2(yr, :, :) = tmp4(yr, :, :) * sftifl
                    end do
                    ;; Add coordinates
                    data2!0 = "year"
                    data2&year = years
                    data2!1 = "lat"
                    data2&lat = A1&lat
                    data2!2 = "lon"
                    data2&lon = A1&lon
                else
                    data2 = data1
                end if
                delete(A1)
                if (isatt(diag_script_info, "ts_minlat")) then
                    data2_glob = area_operations(data2, \
                                                 diag_script_info@ts_minlat, \
                                                 diag_script_info@ts_maxlat, \
                                                 diag_script_info@ts_minlon, \
                                                 diag_script_info@ts_maxlon, \
                                                 "average", True)
                else
                    data2_glob = coswgt_areaave(data2)
                end if
                dim_data2g = dimsizes(data2_glob)
                if (s .eq. 0) then
                    if (dim_data2g .ne. nyearsp1) then
                        data26_rel(i26, 0 : dim_data2g - 1) = \
                            (data2_glob / ref_data2(imod)) - 1
                        data26_rel(i26, dim_data2g :) =  data2_glob@_FillValue
                    else
                        data26_rel(i26, :) = (data2_glob / ref_data2(imod)) - 1
                    end if
                    i26 = i26 + 1
                end if
                if (s .eq. 1) then
                    if (dim_data2g .ne. nyearsp1) then
                        data45_rel(i45, 0 : dim_data2g - 1) = \
                            (data2_glob / ref_data2(imod)) - 1
                        data45_rel(i45, dim_data2g:) =  data2_glob@_FillValue
                    else
                        data45_rel(i45, :) = (data2_glob / ref_data2(imod)) - 1
                    end if
                    i45 = i45 + 1
                end if
                if (s .eq. 2) then
                    if (dim_data2g.ne.nyearsp1) then
                        data60_rel(i60, 0 : dim_data2g - 1) = \
                            (data2_glob / ref_data2(imod)) - 1
                        data60_rel(i60, dim_data2g:) =  data2_glob@_FillValue
                    else
                        data60_rel(i60, :) = (data2_glob / ref_data2(imod)) - 1
                    end if
                    i60 = i60 + 1
                end if
                if (s .eq. 3) then
                    if (dim_data2g.ne.nyearsp1) then
                        data85_rel(i85, 0 : dim_data2g - 1) = \
                            (data2_glob / ref_data2(imod)) - 1
                        data85_rel(i85, dim_data2g :) = data2_glob@_FillValue
                    else
                        data85_rel(i85, :) = (data2_glob / ref_data2(imod)) - 1
                    end if
                    i85 = i85 + 1
                end if
                ;; clean up
                delete([/data1, data2, data3, data0, data2_glob/])
                delete([/date, years, nyears/])
                if (var0 .eq. "snw") then
                    delete([/tmp2, tmp3, tmp4/])
                end if
            end if
            delete(imodrcp)
        end do  ; scenario
        if (var0 .eq. "snw") then
            delete([/sftgif, sftgif_0, sftgif_1, sftlf, sftifl/])
        end if
        delete([/nlat, nlon/])
    end do  ; imod

    ;; average and standard deviation over all models
    datahist_avg = dim_avg_n_Wrap(data2_rel, 0)
    datahist_std = dim_stddev_n_Wrap(data2_rel, 0)

    data26_avg = dim_avg_n_Wrap(data26_rel, 0)
    data26_std = dim_stddev_n_Wrap(data26_rel, 0)

    data45_avg = dim_avg_n_Wrap(data45_rel, 0)
    data45_std = dim_stddev_n_Wrap(data45_rel, 0)

    data60_avg = dim_avg_n_Wrap(data60_rel, 0)
    data60_std = dim_stddev_n_Wrap(data60_rel, 0)

    data85_avg = dim_avg_n_Wrap(data85_rel, 0)
    data85_std = dim_stddev_n_Wrap(data85_rel, 0)

    ;; put all data from 1850-2100 into mean_val and std_val
    mean_val = new((/dim_scen + 1, nyears_full/), typeof(data2_rel),\
                   data2_rel@_FillValue)
    std_val = new((/dim_scen + 1, nyears_full/), typeof(data2_rel),\
                  data2_rel@_FillValue)
    p0 = toint(eyears(0)) - toint(syears(0))
    p1 = toint(eyears(0)) - toint(syears(0)) + 1
    mean_val(0, 0 : p0) = datahist_avg
    std_val(0, 0 : p0) = datahist_std

    mean_val(1, p1 :) = data26_avg
    std_val(1, p1 :) = data26_std
    mean_val(2, p1 :) = data45_avg
    std_val(2, p1 :) = data45_std
    mean_val(3, p1 :) = data60_avg
    std_val(3, p1 :) = data60_std
    mean_val(4, p1 :) = data85_avg
    std_val(4, p1 :) = data85_std

    ;; check number of models for timeperiods in scenarios
    ;; if there are less than 5 models set to missing
    runs = new((/1 + dim_scen, nperiods/), integer)
    runs(0, 0) = dimsizes(ind(models@case_name .eq. "historical" .and. \
                              models@end_year .eq. diag_script_info@eyears(0)))
    do scen = 0, dim_scen - 1
        do p = 0, nperiods - 1
            runs(scen + 1, p) = dimsizes(ind(models@case_name .eq. \
                                             diag_script_info@scenarios(scen) \
                                             .and. models@end_year .eq. \
                                             diag_script_info@eyears(p)))
            if (runs(scen + 1, p) .lt. 5) .and. (p .ne. 0) then
                p2 = toint(eyears(p - 1)) - toint(syears(0))
                p3 = toint(eyears(p)) - toint(syears(0))
                mean_val(scen + 1, p2 : p3) = mean_val@_FillValue
                std_val(scen + 1, p2 : p3) = std_val@_FillValue
            end if
        end do
    end do

    dims_data = dimsizes(mean_val)
    to_plot = new((/3 * dims_data(0), dims_data(1)/), float)

    do j = 0, dims_data(0) - 1
        ;;  if variable sw or lw TOA flux: change direction
        if (var0 .eq. "rlut" .or. var0 .eq. "rsut") then
            to_plot(j * 3, :) = - mean_val(j, :)
        else
            to_plot(j * 3, :) = mean_val(j, :)
        end if
        to_plot((j * 3) + 1, :) = to_plot(j * 3, :) + 1. * std_val(j, :)
        to_plot((j * 3) + 2, :) = to_plot(j * 3, :) - 1. * std_val(j, :)
    end do
    copy_VarMeta(mean_val, to_plot)
    to_plot!0 = "metric"
    to_plot&metric = (/"mean", "+1.std", "-1.std", \
                       "mean", "+1.std", "-1.std", \
                       "mean", "+1.std", "-1.std", \
                       "mean", "+1.std", "-1.std", \
                       "mean", "+1.std", "-1.std"/)
    to_plot!1 = "year"
    to_plot&year = ispan(toint(syears(0)), toint(eyears(nperiods - 1)), 1)
    ;;###########################################
    ;;# Other Metadata: diag_script, var        #
    ;;###########################################
    ;; Add to to_plot, as attributes without prefix
    if(isatt(to_plot, "diag_script")) then  ; Add to existing entries
        temp = to_plot@diag_script
        delete(to_plot@diag_script)
        to_plot@diag_script = array_append_record(temp, (/diag_script/), 0)
        delete(temp)
    else  ; Add as new attribute
        to_plot@diag_script = (/diag_script/)
    end if
    to_plot@var = var0  ; Overwrite existing entry
    if(isatt(variable_info, "long_name")) then
        to_plot@var_long_name = variable_info@long_name
    end if
    if(isatt(variable_info, "units")) then
        to_plot@var_units     = variable_info@units
    end if
    if (verbosity .ge. 5) then
        printVarSummary(to_plot)
    end if
    runs@var = "nr"
    runs!0 = "scenario"
    runs&scenario = diag_script_info@label
    runs!1 = "period"
    runs&period = diag_script_info@syears
    runs@diag_script = (/diag_script/)
    runs@var_long_name = "number of model runs per scenario and period"
    if (verbosity .ge. 5) then
        printVarSummary(runs)
    end if
    ;;###########################################
    ;;# Output to netCDF               #
    ;;###########################################
    ;; Output (controlled by diag_script_info)
    if(isatt(diag_script_info, "ncdf")) then
        to_plot@ncdf = netcdf_dir + "/" + diag_script_info@ncdf
    else
        to_plot@ncdf = netcdf_dir + "/" + diag_script + "_" + var0 + ".nc"
    end if
    ncdf_file = to_plot@ncdf
    ncdf_file@existing = "overwrite"
    ;; Function in ~/interface_scripts/auxiliary.ncl
    ncdf_outfile = ncdf_write(to_plot, ncdf_file)
    if(isatt(diag_script_info, "model_nr") .and. \
        diag_script_info@model_nr .eq. True) then
        ncdf_file@existing = "append"
        ncdf_outfile = ncdf_write(runs, ncdf_file)
    end if

end
