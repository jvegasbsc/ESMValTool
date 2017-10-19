;
;  Requires: snw:*2*s
;
variable_info = True
variable_info@derived = True
variable_info@long_name = "Surface Snow Area"
variable_info@units = "m-2"

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
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20170615-A_lore_ru: written.
;;
local funcname, scriptname, verbosity, data, tmp, f, sftlf, sftgif, SNW_A
begin
    funcname = "calculate"
    scriptname = "variable_defs/snw_area.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 4)

    ;; Read original variable
    data = read_data(index, "snw", "*2*s")
    tmp = extract_data_extended(index, "snw" , data, -1, 0, 0)

    inpath = models@dir(index)
    infile = inpath + "/" + models@experiment(index) + "/fx/sftgif/" + \
             models@name(index)+ "/r0i0p0/sftgif_fx_" + models@name(index) \
             + "_" + models@experiment(index) + "_r0i0p0.nc"
    print(infile)
    f = addfile(infile, "r")
    sftgif = f->sftgif

    infile2 = inpath + "/" + models@experiment(index) + "/fx/sftlf/" + \
             models@name(index)+ "/r0i0p0/sftlf_fx_" + models@name(index) \
             + "_" + models@experiment(index) + "_r0i0p0.nc"

    f = addfile(infile2, "r")
    sftlf = f->sftlf

    sftgif_0 = sftgif * (- 0.01)
    sftgif_1 = sftgif_0 + 1
    sftifl = sftgif_1 * sftlf
    sftifl_conf = conform(tmp, sftifl, (/1, 2/))

    tmp2 = tmp / 2.
    tmp3 = tmp2 * 5.10067e14
    ;; where > 5
    tmp4 = where(tmp3 .gt. 5., tmp3, sftgif@_FillValue) / 100.

    SNW_A = tmp4 * sftifl_conf
    copy_VarCoords(tmp, SNW_A)
    SNW_A@units = "m-2"

    data_new = True
    add_data_var(index, data_new, SNW_A, "snw_area")

    leave_msg(scriptname, funcname, 4)
    return(data_new)
end


