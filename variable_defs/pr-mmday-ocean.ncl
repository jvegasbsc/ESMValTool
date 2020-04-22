;
; Requires: pr:T2*s
;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "interface_scripts/data_handling.ncl"
load "./interface_scripts/messaging.ncl"

variable_info = True
variable_info@derived = True
variable_info@long_name = "Precipitation Rate (ocean only)"
variable_info@units = "mm day-1"

; Contour levels
variable_info@map_ref_Levels = ispan(10, 100, 15) * 0.1
variable_info@map_diff_Levels = ispan(-30, 30, 5) * 0.1

undef("calculate")
function calculate(index [1] : integer,
                   variable [1] : string,
                   field_type [1] : string)
;;                 return_val [1] : logical
;; Arguments:
;;    index    - index to current infile defined in the 'interface_data/ncl.interface'-file
;;    variable - string with current variable
;;    field_type  - string with field type classification
;; Return value:
;;    data_new - logical
;; Description:
;;    Converts precipitation rate from units "kg m-2 s-1" to "mm day-1" and
;;    masks precipitation field over the ocean. If available, the land mask for the
;;    given model is used, otherwise the NCL mask is applied.
;; References
;;    The default NCL mask is $NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc
;;
;; Modification history
;;    20161102-A_laue_ax: written (based on variable_defs/CSOIL.ncl)
;;

local result, T, verbosity, funcname, scriptname, f, a, sftlf, prOcean

begin
    funcname = "calculate"
    scriptname = "variable_defs/pr-mmday-ocean.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 4)

    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering calculate (pr-mmday)", verbosity, 4)
    result = read_data(index, "pr", "T2*s")
    T = extract_data_extended(index, "pr", result, -1, 0, 0)
    T = T*24*3600
    T@units = variable_info@units

    f = addfile(interface_get_fullpath("pr", field_types(0), index), "r")
    if (isfilevar(f, "sftlf")) then
        sftlf = f->sftlf
        sftlf = sftlf / 100.  ; [%] --> [1]
    else
        f = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc", "r")
        a = f->LSMASK
        sftlf = byte2flt(landsea_mask(a, T&lat, T&lon))
        sftlf = where(sftlf.gt.1., 1., sftlf)
        info_output("Using NCL default mask", verbosity, 1)
    end if
    sftlf = where(sftlf.ge.1, sftlf@_FillValue, 0.0)

    ;; Apply mask
    prOcean = T
    delete(prOcean@_FillValue)
    printVarSummary(prOcean)
    prOcean = prOcean + conform(T, sftlf, (/1, 2/))
    prOcean = where(ismissing(prOcean), -999., prOcean)
;    prOcean@_FillValue = -999.

    data_new = True
    add_data_var(index, data_new, prOcean, "pr-mmday-ocean")

    leave_msg(scriptname, funcname, 4)
    return(data_new)
end

