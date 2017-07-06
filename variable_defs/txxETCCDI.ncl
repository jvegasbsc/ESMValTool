;
;  Requires: tasmax:T2*s
;
variable_info = True
variable_info@derived = True
variable_info@long_name = "annual maximum temperature"

;;#############################################################################
load "interface_scripts/data_handling.ncl"
load "interface_scripts/constants.ncl"
load "diag_scripts/lib/ncl/latlon.ncl"
load "./diag_scripts/lib/ncl/statistics.ncl"
load "./interface_data/ncl.interface"

undef("calculate")
function calculate(index[1] : integer,
                   variable[1] : string,
                   field_type[1] : string)
begin
    funcname = "calculate"
    scriptname = "variable_defs/txxETCCDI.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 4)

    ;; Fetch data
    tmp1 = read_data(index, "txxETCCDI", "T2*s")
    data0 = extract_data(index, tmp1, -1, 0, 0)

    return(data0)

end

