;
;  Requires: ts:*2*s
;
variable_info = True
variable_info@derived = True 
variable_info@units = "K"
variable_info@long_name = "Skin temperature (ocean only)"

;;#############################################################################
load "interface_scripts/data_handling.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

undef("calculate")
function calculate(index[1]:integer,
                   variable[1]:string,
                   field_type[1]:string)
;;
;; Arguments
;;    index: index to current infile defined in the 
;;           'interface_data/ncl.interface'-file
;;    variable: a string with current variable
;;    field_type: a string with field type classification
;;
;; Return value
;;    data_new: logical
;;
;; Description
;;    Masks skin temperature over the ocean. If available, the land mask for the
;;    given model is used, otherwise the NCL mask is applied.
;;
;; Caveats
;;
;; References
;;    The default NCL mask is $NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc
;;
;; Modification history
;;    20161020-A_laue_ax: written (based on variable_defs/CSOIL.ncl)
;;
local funcname, scriptname, verbosity, result, tmp, f, a, sftlf, tsOcean
begin

    funcname = "calculate"
    scriptname = "variable_defs/ts-ocean.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 4)

    ;; Read original variable
    result = read_data(index, "ts", "*2*s")
    tmp = extract_data_extended(index, "ts" , result, -1, 0, 0)

    f = addfile(interface_get_fullpath("ts", field_types(0), index), "r")
    if (isfilevar(f, "sftlf")) then
        sftlf = f->sftlf
        sftlf = sftlf / 100.  ; [%] --> [1]
    else
        f = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc", "r")
        a = f->LSMASK
        sftlf = byte2flt(landsea_mask(a, tmp&lat, tmp&lon))
        sftlf = where(sftlf.gt.1., 1., sftlf)
        info_output("Using NCL default mask", verbosity, 1)
    end if
    sftlf = where(sftlf.ge.1, sftlf@_FillValue, 0.0)

    ;; Apply mask
    tsOcean = tmp
    tsOcean = tsOcean + conform(tmp, sftlf, (/1, 2/))
    tsOcean = where(ismissing(tsOcean), -999., tsOcean)
;    tsOcean@_FillValue = -999.

    data_new = True
    add_data_var(index, data_new, tsOcean, "ts-ocean")

    leave_msg(scriptname, funcname, 4)
    return(data_new)

end
