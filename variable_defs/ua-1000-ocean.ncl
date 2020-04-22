;
;  Requires: ua:*3*
;
load "interface_scripts/data_handling.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

variable_info = True
variable_info@derived = True
variable_info@long_name="1000mb wind (ocean only)"
variable_info@units="m s-1"

variable_info@map_ref_Levels = ispan(-70, 70, 10)
variable_info@map_diff_Levels = ispan(-20, 20, 5)

undef("calculate")
function calculate(index [1] : integer,
                   variable [1] : string,
                   field_type [1] : string)
;;                 return_val [1] : logical
;; Arguments:
;;    index    - index to current infile defined in the 'interface_data/ncl.interface'-file
;;    variable - Current variable as string
;;    field_type  - string with field type classification
;; Return value:
;;    data_new - logical
;; Description:
;;    Masks 1000 hPa u-wind cover over the ocean using the NCL mask.

local tmp, dum, i, verbosity
begin
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering ua-1000-ocean.ncl", verbosity, 4)
    tmp = read_data(index, "ua", "*3*")
    dum = extract_data(index, tmp, 1, 100000., 100000.)

    dum@long_name = variable_info@long_name
    dum@units = variable_info@units

    f = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc", "r")
    a = f->LSMASK
    sftlf = byte2flt(landsea_mask(a, dum&lat, dum&lon))
    sftlf = where(sftlf.gt.1., 1., sftlf)

    info_output("Using NCL default mask", verbosity, 1)

    sftlf = where(sftlf.ge.1, sftlf@_FillValue, 0.0)

    ;; Apply mask
    dum = dum + conform(dum, sftlf, (/1, 2/))
    dum = where(ismissing(dum), -999., dum)
;    dum@_FillValue = -999.

    data_new = True
    add_data_var(index, data_new, dum, variable)

    info_output(">>>>>>>> Leaving ua-1000-ocean.ncl", verbosity, 4)
    return(data_new)
end
