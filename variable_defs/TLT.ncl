;
;  Requires: ta:T3M,tas:T2Ms
;

variable_info = True
variable_info@derived = True
variable_info@long_name = "Lower tropospheric temperature"
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
;;      Computes TLT from ta array following MSU satellite
;;      weights. Land mask for given model used, otherwise NCL
;;      mask, as land and ocean have different weights.
;;
;; References
;;    The default NCL mask is $NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc
;;    Weights from ftp://ftp.remss.com/msu/weighting_functions/
;;    accessed 2017-04-06
;;
;; Modification history
;;    20170410-A_gier_be: written.
;;
local tmp1, tmp2, dum, verbosity, land_weights, land_w, oce_weights, oce_w, \
      f, a, sftlf, l_weights, o_weights, p_ini, p_fin, tas_lw, tas_ow, \
      tmp_ta, tmp_tas, tmp_land, tmp_ocean, land_tlt, oce_tlt
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

    ;; Only works if tas, ta on same grid, do not wanna regrid here
    if .not.((dimsizes(tmp_ta&lat).eq.dimsizes(tmp_tas&lat)) .or. \
        (dimsizes(tmp_ta&lon).eq.dimsizes(tmp_tas&lon))) then
        error_msg("f", scriptname, "", "ta and tas are not using" \
        + " the same grid, please regrid manually")
    end if

    ;; Loading weight functions
    ;; Pressure levels in column 4, weights col 6
    land_weights = readAsciiTable( \
       "variable_defs/variable_defs_aux/std_atmosphere_wt_function_chan_tlt_land.txt", \
       6, "float", 7)
    land_w = land_weights(:, 5)
    land_w!0 = "pressure"
    land_w&pressure = land_weights(:, 3)

    oce_weights = readAsciiTable( \
       "variable_defs/variable_defs_aux/std_atmosphere_wt_function_chan_tlt_ocean.txt", \
       6, "float", 7)
    oce_w = oce_weights(:, 5)
    oce_w!0 = "pressure"
    oce_w&pressure = oce_weights(:, 3)

    ;; surface weight is listed before header, not in table
    tas_lw = asciiread( \ 
       "variable_defs/variable_defs_aux/std_atmosphere_wt_function_chan_tlt_land.txt", \
       1, "float")
    tas_ow = asciiread( \ 
       "variable_defs/variable_defs_aux/std_atmosphere_wt_function_chan_tlt_ocean.txt", \
       1, "float")

    
    ;; Calculate weights for data specific pressure layers
    l_weights = new(dimsizes(tmp_ta&plev), "float")
    o_weights = new(dimsizes(tmp_ta&plev), "float")
    do nplev = 1, dimsizes(tmp_ta&plev)-1
        ;; Get indizes - pressure levels are the same for both
        ;; weight functions
        p_ini = ind_nearest_coord(tmp_ta&plev(nplev-1), oce_w&pressure, 0)
        p_fin = ind_nearest_coord(tmp_ta&plev(nplev), oce_w&pressure, 0)
        ;; avg of weights in the area. weights given in km-1
        ;; each step resembles 300m -> sum * 0.3
        ;; instead of avg * (0.3*n)
        ;; weights on file smaller spacing than plev from data
        ;; p_fin-1 to not count some levels double
        l_weights(nplev-1) = sum(land_w(p_ini:p_fin-1))*0.3
        o_weights(nplev-1) = sum(oce_w(p_ini:p_fin-1))*0.3
        if nplev.eq.(dimsizes(tmp_ta&plev)-1) then
            l_weights(nplev) = sum(land_w(p_fin:))*0.3
            o_weights(nplev) = sum(oce_w(p_fin:))*0.3
        end if
        delete([/p_ini, p_fin/])
    end do

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

    ;; Apply weights and calc TLT for land, oce seperately
    tmp_land = tmp_land * conform(tmp_land, l_weights, 1)
    tmp_oce = tmp_oce * conform(tmp_oce, o_weights, 1)
    copy_VarMeta(tmp_ta, tmp_land)
    copy_VarMeta(tmp_ta, tmp_oce)
    land_tlt = dim_sum_n_Wrap(tmp_land, 1)
    land_tlt = land_tlt + tas_land*tas_lw
    oce_tlt = dim_sum_n_Wrap(tmp_oce, 1)
    oce_tlt = oce_tlt + tas_oce*tas_ow

    ;; Add TLT
    dum = land_tlt + oce_tlt
    copy_VarMeta(land_tlt, dum)

    data_new = True
    dum@long_name = variable_info@long_name
    dum@units = variable_info@units
    add_data_var(index, data_new, dum, variable)

    leave_msg(scriptname, funcname, 4)
    return(data_new)
end
