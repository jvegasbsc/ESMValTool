;
;  Requires: ta:T3M,tas:T2Ms
;

variable_info = True
variable_info@derived = True
variable_info@long_name = "Lower tropospheric temperature"
variable_info@short_name = "TLT"
variable_info@units = "K"

;;###############################################################
load "interface_scripts/data_handling.ncl"
load "interface_scripts/auxiliary.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

undef("calculate")
function calculate(index [1] : integer,
                   variable [1] : string,
                   field_type [1] : string)
;;                 return_val [1] : logical
;; Arguments
;;    index    - index to current infile defined in the
;;               'interface_data/ncl.interface'-file
;;    variable - Current variable as string
;;    field_type  - string with field type classification
;;
;; Return value
;;    data_new - logical
;;
;; Description 
;;      Computes TLT from ta and tas following MSU satellite
;;      weights. Currently using NCL landsea mask, as land
;;      and ocean have different weights. tas will be regridded
;;      to ta coordinates if necessary. Missing Value treatment
;;      following MSU satellite convention for RSS data.
;;
;; References
;;    The default NCL mask is $NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc
;;    Weights from ftp://ftp.remss.com/msu/weighting_functions/
;;    accessed 2017-04-06
;;
;; Modification history
;;    20171025-A_laue_ax: changed procedure to calculate vertical weights
;;             added scaling of weight for surface temperature
;;    20171024-A_gier_be: Added checks for pressure levels
;;    20170517-A_gier_be: Added treatment for missing values
;;    20170410-A_gier_be: written.
;;
local tmp1, tmp2, dum, verbosity, land_weights, land_w, oce_weights, oce_w, \
      f, a, sftlf, l_weights, o_weights, tas_lw, tas_ow, \
      tmp_ta, tmp_tas, tmp_land, tmp_ocean, land_tlt, oce_tlt, xregrid, \
      yregrid, var_hreg, p_cut, l_scale, land_scale, ocean_scale, \
      o_scale, land_height, oce_height, pfine, o_weights_fine, l_weights_fine, \
      o_heights_fine, l_heights_fine, n, odz, ldz, ibot, itop, pbot, ptop, \
      omean, lmean
begin

    funcname = "calculate"
    scriptname = "variable_defs/TLT.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 4)

    ;; Read original variable
    tmp1 = read_data(index, "ta", "T3M")
    tmp2 = read_data(index, "tas", "T2Ms")
    tmp_ta = extract_data_extended(index, "ta", tmp1, -1, 0, 0)
    tmp_tas = extract_data_extended(index, "tas", tmp2, -1, 0, 0)
    delete([/tmp1, tmp2/])

    ;; Check if plev satisfies assumptions: dim 1
    ;; Give error if not - do not want to reorder 3D fields
    ;; CAVEAT: does not check for known aliases like p, plevs, lev_p
    if (tmp_ta!1).ne."plev" then
        error_msg("f", scriptname, funcname, "Input field does not " + \
                  "have plev as second coordinate.")
    end if

    ;; Check that plev is Monotonic -> required for weight calculation
    if isMonotonic(tmp_ta&plev).eq.0 then
        error_msg("f", scriptname, funcname, "Pressure levels are not " + \
                  "monotonic. Please adapt input field.")
    end if

    ;; Only works if tas, ta on same grid, regrid tas if not
    if ((dimsizes(tmp_ta&lat).ne.dimsizes(tmp_tas&lat)) .or. \
        (dimsizes(tmp_ta&lon).ne.dimsizes(tmp_tas&lon))) then
        info_output("Regridding tas", verbosity, 1)
        ; Need to regrid, determine dimensions to regrid
        xregrid = True ;lon
        yregrid = True ;lat
        if (dimsizes(tmp_tas&lon).eq.dimsizes(tmp_ta&lon)) then
            if (max(abs(tmp_tas&lon - tmp_ta&lon)).lt.1.0d-5) then
                xregrid = False
            end if
        end if
        if (dimsizes(tmp_tas&lat).eq.dimsizes(tmp_ta&lat)) then
            if (max(abs(tmp_tas&lat - tmp_ta&lat)).lt.1.0d-5) then
                yregrid = False
            end if
        end if

        dnames = getVarDimNames(tmp_tas)
        if (xregrid.and.yregrid) then
            var_hreg = linint2_Wrap(tmp_tas&lon, tmp_tas&lat, tmp_tas, True, \
                                      tmp_ta&lon, tmp_ta&lat, 0)
        else if (xregrid.and..not.yregrid) then
            dx = ind(dnames.eq."lon")
            var_hreg = linint1_n_Wrap(tmp_tas&lon, tmp_tas, True, \
                                                  tmp_ta&lon, 0, dx)
        else if (.not.xregrid.and.yregrid) then
            dx = ind(dnames.eq."lat")
            var_hreg = linint1_n_Wrap(tmp_tas&lat, tmp_tas, False, \
                                                   tmp_ta&lat, 0, dx)
        end if
        end if
        end if

        delete(tmp_tas)
        tmp_tas = var_hreg
        delete(var_hreg)
        delete(xregrid)
        delete(yregrid)
    end if

    ;; Requires landsea mask due to different weights
    ;; Read NCL standard mask
    f = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc", "r")
    a = f->LSMASK
    sftof = byte2flt(landsea_mask(a, tmp_ta&lat, tmp_ta&lon))
    sftof = where(sftof.eq.0., 1., 0.)
    ;; exactly opposite for landmask
    sftlf = where(sftof.eq.0., 1., 0.)

    ;; Apply mask
    tmp_land = tmp_ta * conform(tmp_ta, sftlf, (/2,3/))
    tas_land = tmp_tas * conform(tmp_tas, sftlf, (/1,2/))
    tmp_oce = tmp_ta * conform(tmp_ta, sftof, (/2,3/))
    tas_oce = tmp_tas * conform(tmp_tas, sftof, (/1,2/))    

    ;; Loading weight functions
    weight_folder = "variable_defs/variable_defs_aux/"
    ;; Pressure levels in column 4, weights col 6
    land_weights = readAsciiTable(weight_folder + \
       "std_atmosphere_wt_function_chan_tlt_land.txt", 6, "float", 7)
    land_w = land_weights(:, 5)
    land_w!0 = "pressure"
    land_w&pressure = land_weights(:, 3)
    land_height = land_weights(:,1) ;height in m

    oce_weights = readAsciiTable(weight_folder + \
       "std_atmosphere_wt_function_chan_tlt_ocean.txt", 6, "float", 7)
    oce_w = oce_weights(:, 5)
    oce_w!0 = "pressure"
    oce_w&pressure = oce_weights(:, 3)
    oce_height = oce_weights(:,1) ;height in m

    ;; surface weight is listed before header, not in table
    tas_lw = asciiread(weight_folder + \
       "std_atmosphere_wt_function_chan_tlt_land.txt", 1, "float")
    tas_ow = asciiread(weight_folder + \
       "std_atmosphere_wt_function_chan_tlt_ocean.txt", 1, "float")

    ;; calculate weights for levels in ta array
    ;; weight = average over level times layer thickness in km

    pfine = fspan(0., 101325., 1013)

    o_weights_fine = tofloat(ftcurv(oce_w&pressure(::-1), oce_w(::-1), pfine))
    l_weights_fine = tofloat(ftcurv(land_w&pressure(::-1), land_w(::-1), pfine))
    o_heights_fine = tofloat(ftcurv(oce_w&pressure(::-1), oce_height(::-1), pfine))
    l_heights_fine = tofloat(ftcurv(land_w&pressure(::-1), land_height(::-1), pfine))

    n = dimsizes(tmp_ta&plev)

    o_weights = new((/n/), float)
    l_weights = new((/n/), float)

    do i = 0, n - 1
        if (i.gt.0) then
            pbot = 0.5*tofloat((tmp_ta&plev(i) + tmp_ta&plev(i-1)))
        else
            pbot = 101325.
        end if
        if (i.lt.n-1) then
            ptop = 0.5*tofloat((tmp_ta&plev(i) + tmp_ta&plev(i+1)))
        else
            ptop = 0.
        end if

        ibot = ind_nearest_coord(pbot, pfine, 0)
        itop = ind_nearest_coord(ptop, pfine, 0)

        omean = avg(o_weights_fine(ibot:itop))
        odz = (o_heights_fine(itop) - o_heights_fine(ibot)) / 1000.0 ; in km
        o_weights(i) = omean * odz

        lmean = avg(l_weights_fine(ibot:itop))
        ldz = (l_heights_fine(itop) - l_heights_fine(ibot)) / 1000.0 ; in km
        l_weights(i) = lmean * ldz
    end do

    ;; Some land over mountains have missing values in pressure levels
    ;; Treat these like in TLT from RSS - "Land areas w/ surface height
    ;; exceeding threshold of 1500m (~85 000 Pa) were excluded"
    ;; Missing values change per time and height, so 4D weighting
    ;; array is required

    p_cut = ind_nearest_coord(85000, tmp_ta&plev, 0)

    ;; Variant 1: add up all weights and apply to lowest level available data
    ;; used because it resembles RSS obs the most, smoother
    ;; Only execute loop if missing values exist!
    if any(ismissing(tmp_ta)) then
        info_output("Treating missing values.", verbosity, 1)
        land_weight_arr = conform(tmp_ta, l_weights, 1)
        land_weight_arr@_FillValue = default_fillvalue(typeof(land_weight_arr))
        ocean_weight_arr = conform(tmp_ta, o_weights, 1)
        ocean_weight_arr@_FillValue = default_fillvalue(typeof(ocean_weight_arr))

        if isMonotonic(tmp_ta&plev).eq.(-1) then
            do iplev = 0, p_cut
                land_weight_arr(:, iplev + 1, :, :)  = where(ismissing( \
                                tmp_ta(:, iplev, :, :)), \
                                land_weight_arr(:, iplev, :, :) + \
                                land_weight_arr(:, iplev + 1, :, :), \
                                land_weight_arr(:, iplev + 1, :, :))
                land_weight_arr(:, iplev, :, :)  = where(ismissing( \
                                tmp_ta(:, iplev, :, :)), \
                                land_weight_arr@_FillValue, \ 
                                land_weight_arr(:, iplev, :, :))
                ocean_weight_arr(:, iplev + 1, :, :)  = where(ismissing( \
                                 tmp_ta(:, iplev, :, :)), \
                                 ocean_weight_arr(:, iplev, :, :) + \ 
                                 ocean_weight_arr(:, iplev + 1, :, :), \ 
                                 ocean_weight_arr(:, iplev + 1, :, :))
                ocean_weight_arr(:, iplev, :, :)  = where(ismissing( \
                                 tmp_ta(:, iplev, :, :)), \
                                 ocean_weight_arr@_FillValue, \ 
                                 ocean_weight_arr(:, iplev, :, :))
            end do
        else
            do nplev = 0, dimsizes(tmp_ta&plev)-1 - p_cut
                iplev = dimsizes(tmp_ta&plev) - 1 - nplev
                land_weight_arr(:, iplev - 1, :, :)  = where(ismissing( \
                                tmp_ta(:, iplev, :, :)), \
                                land_weight_arr(:, iplev, :, :) + \ 
                                land_weight_arr(:, iplev - 1, :, :), \ 
                                land_weight_arr(:, iplev - 1, :, :))
                land_weight_arr(:, iplev, :, :)  = where(ismissing( \
                                tmp_ta(:, iplev, :, :)), \
                                land_weight_arr@_FillValue, \ 
                                land_weight_arr(:, iplev, :, :))
                ocean_weight_arr(:, iplev - 1, :, :)  = where(ismissing( \
                                 tmp_ta(:, iplev, :, :)), \
                                 ocean_weight_arr(:, iplev, :, :) + \ 
                                 ocean_weight_arr(:, iplev - 1, :, :), \ 
                                 ocean_weight_arr(:, iplev - 1, :, :))
                ocean_weight_arr(:, iplev, :, :)  = where(ismissing( \
                                 tmp_ta(:, iplev, :, :)), \
                                 ocean_weight_arr@_FillValue, \ 
                                 ocean_weight_arr(:, iplev, :, :))
            end do
        end if

        l_scale = dim_sum_n(land_weight_arr, 1) + tas_lw
        land_scale = conform(land_weight_arr, l_scale, (/0, 2, 3/))
        land_weight_arr = land_weight_arr / land_scale
        o_scale = dim_sum_n(ocean_weight_arr, 1) + tas_ow
        ocean_scale = conform(ocean_weight_arr, o_scale, (/0, 2, 3/))
        ocean_weight_arr = ocean_weight_arr / ocean_scale

    else
        l_scale = sum(l_weights)+tas_lw
        o_scale = sum(o_weights)+tas_ow

        ;; If no missing values -> weights the same
        land_weight_arr = conform(tmp_land, l_weights, 1) / l_scale
        ocean_weight_arr = conform(tmp_oce, o_weights, 1) / o_scale

    end if

    ;; Apply weights and calc TLT for land, oce seperately
    tmp_land = tmp_land * land_weight_arr
    tmp_oce = tmp_oce * ocean_weight_arr

    delete(land_weight_arr)
    delete(ocean_weight_arr)

    copy_VarMeta(tmp_ta, tmp_land)
    copy_VarMeta(tmp_ta, tmp_oce)
    land_tlt = dim_sum_n_Wrap(tmp_land, 1)
    land_tlt = land_tlt + tas_land*tas_lw/l_scale
    oce_tlt = dim_sum_n_Wrap(tmp_oce, 1)
    oce_tlt = oce_tlt + tas_oce*tas_ow/o_scale  

    ;; Add TLT
    dum = land_tlt + oce_tlt
    copy_VarMeta(land_tlt, dum)

    delete(land_tlt)
    delete(oce_tlt)

    ;; Take out all points with missing values above the allowed cutoff
    dum = where(ismissing(tmp_ta(:, p_cut, :, :)), \
               dum@_FillValue, dum)
    
    delete(tmp_ta)

    data_new = True
    dum@long_name = variable_info@long_name
    dum@units = variable_info@units
    add_data_var(index, data_new, dum, variable)

    leave_msg(scriptname, funcname, 4)
    return(data_new)
end
