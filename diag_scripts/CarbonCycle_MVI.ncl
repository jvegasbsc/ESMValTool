;;#############################################################################
;; DIAGNOSTIC SCRIPT TO CREATE MVI-IAV-TREND PLOT
;; Author: Alessandro Anav (University of Exeter, UK)
;; EMBRACE project
;;#############################################################################
;;
;; Description
;;    Calculates the model variability index (MVI), interannual variability
;;    (IAV) and mean, and draws them in a 3D scatter plot.
;;
;; Required diag_script_info attributes (diagnostics specific)
;;    grid: target grid for the regridding
;;    region: region to be averaged
;;    ref_model: reference for the MVI calculation
;;
;; Optional diag_script_info attributes (diagnostic specific)
;;
;; Required variable_info attributes (variable specific)
;;
;; Optional variable_info attributes (variable specific)
;;
;; Caveats
;;    The derived variables pr-mmday and tas-degC needs to be used in order
;;    to get the correct linear trend. For land variables nbp, gpp and lai, it
;;    is recommended to use the masked version (NBP, GPP and LAI, respectively)
;;    in order to have a consistent land mask for all models.
;;
;; Modification history
;;    20170316-A_gier_be: Added cdo regridding, landmask for CRU obs
;;    20151130-A_laue_ax: moved call to write_references to beginning of code
;;    20151104-A_righ_ma: graphic part moved to external plot_scripts.
;;    20151102-A_righ_ma: replacing area functions with ESMValTool functions.
;;    201507??-A_wenz_sa: adapted to ESMValTool structure.
;;    201505??-A_anav_al: written.
;;
;;#############################################################################
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

load "./interface_data/ncl.interface"

load "./interface_scripts/auxiliary.ncl"
load "./interface_scripts/data_handling.ncl"
load "./interface_scripts/messaging.ncl"
load "./interface_scripts/constants.ncl"

load "./diag_scripts/lib/ncl/latlon.ncl"
load "./diag_scripts/lib/ncl/statistics.ncl"
load "./diag_scripts/lib/ncl/regridding.ncl"
load "./diag_scripts/lib/ncl/scaling.ncl"

load "./plot_scripts/ncl/scatterplot.ncl"

begin

    verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(diag_script, "", 2)
    info_output("++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)
    info_output(diag_script + " (var: " + variables(0) + ")", verbosity, 1)
    info_output("++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)

    dim_MOD = dimsizes(models@name)
    dim_VAR = dimsizes(variables)
    var0 = variables(0)
    field_type0 = field_types(0)

    ;; Read models data
    data = True
    do idx_mod = 0, dim_MOD - 1
        data = read_data(idx_mod, var0, field_type0)
    end do

end

begin
    vardeffile = "interface_data/" + var0 + "_info.tmp"
    loadscript(vardeffile)
end

begin
    write_references(diag_script,                  \  ; script name
                     (/"A_anav_al"/),              \  ; authors
                     (/"A_wenz_sa", "A_righ_ma"/), \  ; contributors
                     (/"D_anav13jclim"/),          \  ; diag_references
                     (/""/),                       \  ; obs_references
                     (/"P_embrace"/))                 ; proj_references

    ;; Basename of diag_script
    diag_script_base = basename(diag_script)

    ;; Define file type
    file_type = getenv("ESMValTool_output_file_type")
    if(ismissing(file_type)) then
        file_type = "ps"
    end if

    ;; Check required diag_script_info attributes
    req_atts = (/"grid", "region", "ref_model"/)
    exit_if_missing_atts(diag_script_info, req_atts)

    ;; Check field type
    plot_fields = (/"T2Ms"/)
    if (field_type0.ne.plot_fields) then
        error_msg("f", diag_script, "", "can't use field " + field_type + \
                  " for this diag_script")
    end if
    delete(plot_fields)

    ;; Set region - supports 2D array for lat and strings
    if isstring(diag_script_info@region) then
        region = select_region(diag_script_info@region)
    else if dimsizes(diag_script_info@region).eq.2 then
        region = (/diag_script_info@region(0), \
                   diag_script_info@region(1), 0., 360./)        
        region@name = flt2string(diag_script_info@region(0)) + "_-_" + \
                      flt2string(diag_script_info@region(1))
        else
            error_msg("f", diag_script, "", \
                  diag_script_info@region + " is an invalid region type")
        end if
    end if

    ;; Output plot directory
    plot_dir = getenv("ESMValTool_plot_dir")
    output_dir = get_output_dir(plot_dir, diag_script_base)

    ;; Set time range
    start_year = min(stringtoint(models@start_year))
    end_year = max(stringtoint(models@end_year))
    all_years = ispan(start_year, end_year, 1)

    ;; Computed mean over all years, remember that some models have 
    ;; longer timespan than what we analyze here
    ;; Anav mean and MVI over 1986-2005! -> put this in cfg as var?
    mean_minyear = 1986
    mean_maxyear = 2005
    minyear_ind = ind(all_years.eq.mean_minyear)
    maxyear_ind = ind(all_years.eq.mean_maxyear)

    ;; Check for irregular grid
    l_irreg = new(dim_MOD, logical)
    l_irreg = False
    do imod = 0, dim_MOD - 1
        cnt = interface_get_dictkeys(var0, imod) + "_dimNames"
        dims = data@$cnt$
        ndims = dimsizes(dims)
        if (ndims.lt.2) then
            continue
        end if
        if (ndims.eq.2 .and. dims(0).eq."time") then
            continue
        end if
        if (dims(ndims - 1).ne."lon" .and. dims(ndims - 2).ne."lat") then
            l_irreg(imod) = True
            continue
        end if
    end do
    delete([/cnt, dims, ndims/])


    ;; Regrid to destination grid
    grid_resolution        = diag_script_info@grid

    if (grid_resolution.eq."0.5deg")    then
        dst_grid = "./diag_scripts/aux/CDO/cdo_dst_grid_g005"
        dNLON = 720
        dNLAT = 360
    else if (grid_resolution.eq."1deg") then
        dst_grid = "./diag_scripts/aux/CDO/cdo_dst_grid_g010"
        dNLON = 360
        dNLAT = 180
    else if (grid_resolution.eq."2deg") then
        dst_grid = "./diag_scripts/aux/CDO/cdo_dst_grid_g020"
        dNLON = 180
        dNLAT = 90
    else
        error_msg("f", diag_script, "", "Wrong destination grid defined: " + \
                  diag_script_info@grid)
    end if
    end if
    end if

    ;; Gridded matrix for a given resolution where all yearly data are stored
    all_data_yearly = \
        new((/dim_MOD, end_year - start_year + 1, dNLAT, dNLON/), float)
    all_data_yearly!0 = "model"
    all_data_yearly!1 = "year"
    all_data_yearly!2 = "lat"
    all_data_yearly!3 = "lon"
    all_data_yearly&model = models@name

    ;; Gridded matrix for a given resolution where all MVI data are stored
    MVI_all_data = new((/dim_MOD, dNLAT, dNLON/), float)
    MVI_all_data!0 = "model"
    MVI_all_data!1 = "lat"
    MVI_all_data!2 = "lon"

    ;; Gridded matrix for a given resolution where all trend data are stored
    Trend_all_data = new((/dim_MOD, dNLAT, dNLON/), float)
    Trend_all_data!0 = "model"
    Trend_all_data!1 = "lat"
    Trend_all_data!2 = "lon"


    ;; Loop over models and regrid
    do imod = 0, dim_MOD - 1

        info_output("Processing " + models@name(imod), verbosity, 1)

        ;; Extract model corresponding data
        var  = extract_data_extended(imod, var0, data, -1, 0, 0)

        ;; Regridding - do always since regrid to 2x2 after?
        ;;if (l_irreg(imod)) then
            regr_file = basename(diag_script) + "_" + models@name(imod)+ "_" + \
                        var0 + "_" + start_year + "-" + end_year + "_grid_" + \
                       grid_resolution
            regr_path = wrk_dir + "/CDO/" + regr_file + ".nc"            
            if (isfilepresent(regr_path)) then
                var_regrid = ncdf_read(regr_path, var0)
            else
                var_regrid = cdo_remapdis(var0, var, imod, data, \
                             dst_grid, regr_file)
            end if
;;            if (var_regrid!0.eq."lon") then  ; reorder dimensions
                var_regrid := var_regrid(time|:, lat|:, lon|:)
;;            end if

            delete(var)
            var = var_regrid
            delete(var_regrid)
        ;;else
            
        ;;end if


        ;; CRU is only land-based - Need landmask and additionally cut off Antarctica
        if (isatt(diag_script_info, "ref_model")).and.diag_script_info@ref_model.eq."CRU" then
            f = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc", "r")
            a = f->LSMASK
            sftlf = byte2flt(landsea_mask(a, var&lat, var&lon))
            sftlf = where(sftlf.gt.1., 1., sftlf)
            sftlf = where(sftlf.eq.0, sftlf@_FillValue, sftlf)
        ;; Mask out Antarctica (below 60degS)
            sftlf({:-60}, :) = sftlf@_FillValue
        ;; Apply mask
            var = var * conform(var, sftlf, (/1, 2/))
            delete([/sftlf, a, f/])
        end if

        ;; Computing annual mean
        data_annual = time_operations(var, toint(models@start_year(imod)), \
                                      toint(models@end_year(imod)), "average",\
                                      "yearly", True)

        ;; Place array in model specific time range
        idx1 = start_year - toint(models@start_year(imod))
        idx2 = dimsizes(all_years) - \
            (end_year - toint(models@end_year(imod))) - 1
        
        all_data_yearly(imod, idx1:idx2, :, :) = data_annual
;;                regrid_3D_to_rectilinear_grid(data_annual, var&lon, var&lat, \
;;                                              grid_resolution, "", True)

        ;;If lat and lon are 2D arrays, assume irregular grid
;;        if (dimsizes(dimsizes(var&lat)).ge.2 .and. \
;;            dimsizes(dimsizes(var&lon)).ge.2) then
;;            filename = interface_get_fullpath(var0, field_type0, imod)
;;            all_data_yearly(imod, idx1:idx2, :, :) = \
;;                regrid_3D_to_rectilinear_grid(data_annual, var&lon, var&lat, \
;;                                              grid_resolution, filename, False)

;;        else  ; Regular grid
;;            all_data_yearly(imod, idx1:idx2, :, :) = \
;;                regrid_3D_to_rectilinear_grid(data_annual, var&lon, var&lat, \
;;                                              grid_resolution, "", True)
;;        end if
;;        all_data_yearly&model = models@name

        ;; Cleaning up
        delete(var)
        delete(data_annual)

    end do  ; end loop over models

    ;; Compute annual mean, MVI, trend
    MVI_all_data&model = models@name
    MVI_all_data&lat = all_data_yearly&lat
    MVI_all_data&lon = all_data_yearly&lon
    MVI_all_data@_FillValue = 1.e20

    Trend_all_data&model = models@name
    Trend_all_data&lat = all_data_yearly&lat
    Trend_all_data&lon = all_data_yearly&lon
    Trend_all_data@_FillValue = 1.e20

    ;; Determine reference model and compute MVI
    ref_ind = ind(models@name.eq.diag_script_info@ref_model)
    if (ismissing(ref_ind)) then
        error_msg("f", diag_script, "", "no adequate reference model provided")
    end if

    info_output("Reference model: " + models@name(ref_ind), verbosity, 1)
    ; Remember to restrict MVI to time period 1986-2005
    reference = all_data_yearly(ref_ind, minyear_ind:maxyear_ind, :, :)
    ref = dim_stddev_n(reference, 0)
    ref = where(ref.ne.0, ref, all_data_yearly@_FillValue)

    ;; Loop over models
    do imod = 0, dim_MOD - 1

        temp  = all_data_yearly(imod, minyear_ind:maxyear_ind, :, :)
        model = dim_stddev_n(temp, 0)

        ;; Prevents a division by 0
        model = where(model.ne.0, model, all_data_yearly@_FillValue)

        ;; Commpute MVI
        MVI_all_data(imod, :, :) = (model/ ref -  ref/model) ^ 2
        delete(model)
        delete(temp)

        model = all_data_yearly(imod, :, :, :)

         ;; Commpute Trend (for ALL years)
        rc = regCoef_n(1.0 * all_years, model, 0, 0)
        rc!0 = "lat"
        rc!1 = "lon"
        rc&lat = model&lat
        rc&lon = model&lon
        Trend_all_data(imod, :, :) = rc
        delete(model)
        delete(rc)

    end do

    ;; Put a cap. Where MVI > 1000 we assign a FillValue
    ;;info_output(" Total points: " + num(.not.ismissing(MVI_all_data)), verbosity, 1)
    ;;info_output(" Points with MVI>10: " +  num(MVI_all_data.gt.10) + \
    ;; "  MVI > 50: " + num(MVI_all_data.gt.50) + "  MVI > 100: " + \
    ;; num(MVI_all_data.gt.100) + "  MVI > 1000: " + num(MVI_all_data.gt.1000), verbosity, 1)
    ;; How is 50 a fill value? Often more than half the MVI lie above 10!
    ;; This doesn't make any sense....
    MVI_all_data = where(MVI_all_data.lt.10, MVI_all_data, 50)

    ;; Special case for LAI
    ;; Since models have different land covers, all sparse vegetated points
    ;; must be masked out
    if (var0.eq."LAI" .or. var0.eq."LAI") then
        all_data_yearly = where(all_data_yearly.ge.0.025, \
                                all_data_yearly, 1.e20)
    end if

    all_annual_mean = dim_avg_n(all_data_yearly(:, minyear_ind:maxyear_ind, :, :), 1)
    all_annual_mean!0 = "model"
    all_annual_mean!1 = "lat"
    all_annual_mean!2 = "lon"

    all_annual_mean&model = models@name
    all_annual_mean&lat = all_data_yearly&lat
    all_annual_mean&lon = all_data_yearly&lon
    all_annual_mean@units = all_data_yearly@units
    delete(all_data_yearly)

    Mean_var = new((/dim_MOD/), typeof(all_annual_mean))
    Mean_var!0 = "model"
    Mean_var&model = models@name

    Trend_var = new((/dim_MOD/), typeof(all_annual_mean))
    Trend_var!0 = "model"
    Trend_var&model = models@name

    MVI_var = new((/dim_MOD/), typeof(all_annual_mean))
    MVI_var!0 = "model"
    MVI_var&model = models@name

    ;; Compute global area
    area = map_area(all_annual_mean&lat, all_annual_mean&lon)

    ;; Set type of area operation (flux vs. non-flux variables)
    if (all_annual_mean@units.eq."kg m-2 s-1") then
        sp_opt = "sum"  ; for flux variable
    else
        sp_opt = "average"  ; for non-flux variable
    end if

    ;; Loop over models
    do imod = 0, dim_MOD - 1
        Mean_var(imod) = area_operations(all_annual_mean(imod, :, :), \
                                         region(0), region(1), region(2), \
                                         region(3), sp_opt, True)
        Trend_var(imod) = area_operations(Trend_all_data(imod, :, :), \
                                          region(0), region(1), region(2), \
                                          region(3), sp_opt, True)
        MVI_var(imod) = area_operations(MVI_all_data(imod, :, :), \
                                        region(0), region(1), region(2), \
                                        region(3), "average", True)
    end do

    ;; Diagnostic-specific conversions
    ;; This cannot be done via the @plot_units attribute, since it also
    ;; depends on the area integration/average
    MVI_var@units = "1"
    l_conv = False
    if (var0.eq."pr-mmday") then
        ;; [mm d-1] --> [mm y-1]
        Mean_var = Mean_var * 365.
        Mean_var@units = "mm y-1"
        ;; [mm d-1 y-1] --> [mm y-1 decade-1]
        Trend_var = Trend_var * 365. * 10.
        Trend_var@units = "mm y-1 decade-1"
        l_conv = True
    end if
    if (var0.eq."tas-degC") then
        ;; [degC] --> [degC decade-1]
        Trend_var = Trend_var * 10.
        Trend_var@units = "degC decade-1"
        l_conv = True
    end if
    if (var0.eq."tos") then
        Mean_var = Mean_var - 273.15
        Mean_var@units = "degC"
        Trend_var = Trend_var * 10.
        Trend_var@units = "degC decade-1"
        l_conv = True
    end if
    if (any(var0.eq.(/"NBP", "GPP"/))) then
        ;; [kg s-1] --> [Pg y-1]
        Mean_var = Mean_var * 3600. * 24. * 365 / 1.e12
        Mean_var@units = "PgC y-1"
        ;; [kg s-1 y-1] --> [Pg y-2]
        Trend_var = Trend_var * 3600. * 24. * 365 / 1.e12
        Trend_var@units = "PgC y-2"
        l_conv = True
    end if
    if (var0.eq."LAI") then
        Trend_var@units = "y-1"
        l_conv = True
    end if


    ;; Check conversions
    if (.not.l_conv) then
        error_msg("f", diag_script, "", "no unit conversion applied to " + \
                  "variable " + var0)
    end if

    ;; Optional output to NetCDF
    if (getenv("ESMValTool_write_netcdf").eq."True") then

        ;; Mean
        new_path = get_output_dir(get_work_dir(), basename(diag_script))
        new_path = new_path + basename(diag_script) + \
                   "_mean_" + var0 + "_" + region@name + ".nc"
        Mean_var@var = var0
        Mean_var@diag_script = diag_script
        Mean_var@ncdf = new_path
        ncdf_outfile = ncdf_write(Mean_var, new_path)

        ;; Trend
        new_path = get_output_dir(get_work_dir(), basename(diag_script))
        new_path = new_path + basename(diag_script) + \
                   "_trend_" + var0 +  "_" + region@name + ".nc"
        Trend_var@var = var0
        Trend_var@diag_script = diag_script
        Trend_var@ncdf = new_path
        ncdf_outfile = ncdf_write(Trend_var, new_path)

        ;; MVI
        new_path = get_output_dir(get_work_dir(), basename(diag_script))
        new_path = new_path + basename(diag_script) + \
                   "_MVI_" + var0 +  "_" + region@name + ".nc"
        MVI_var@var = var0
        MVI_var@diag_script = diag_script
        MVI_var@ncdf = new_path
        ncdf_outfile = ncdf_write(MVI_var, new_path)

    end if

    ;; Set plot output file
    outfile = output_dir + diag_script + "_" + var0 + "_" + \
        region@name
    wks = gsn_open_wks(file_type, outfile)
    gsn_define_colormap(wks, "BlAqGrYeOrReVi200")
    gsn_reverse_colormap(wks)

    ;; Pack data
    data_arr = new((/3, dimsizes(Mean_var)/), float)
    data_arr(0, :) = (/Mean_var/)
    data_arr(1, :) = (/Trend_var/)
    data_arr(2, :) = (/MVI_var/)
    data_arr!0 = "statistic"
    data_arr!1 = "models"
    data_arr&statistic = (/"mean", "linear trend", "MVI"/)
    data_arr&models = models@name
    data_arr@units = (/Mean_var@units, Trend_var@units, "1"/)

    ;; Set levels for the color coding
    tmp = MVI_var
    tmp(ref_ind) = tmp@_FillValue
    nlevs = 11
    tmp := nice_mnmxintvl(min(tmp), max(tmp), nlevs, True)
    levs = decimalPlaces(fspan(tmp(0), tmp(1), 11), 2, True)
    delete(tmp)

    ;; Draw a 3D scatterplot
    data_arr@res_tiMainString = var0 + " - " + region@name
    plot = scatterplot3D(wks, data_arr, var0, levs)
    draw(plot)
    frame(wks)
    delete(data_arr)

    leave_msg(diag_script, "", 2)

end
