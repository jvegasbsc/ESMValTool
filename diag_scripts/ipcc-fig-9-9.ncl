;;#############################################################################
;; DIAGNOSTIC SCRIPT for reproducing IPCC ch. 9 fig. 9.9
;; Author: Bettina Gier (DLR, Germany)
;; CRESCENDO project
;;#############################################################################
;;
;; Description
;;    Calculated centred pattern correlations for annual mean climatologies
;;    and plots them. Like IPCC ch. 9 fig 9.6
;;
;; Required diag_script_info attributes (diagnostics specific)
;;     diag_script_info@collect: Set to True for last variable for collecting
;;                               all variables and plotting
;;
;; Optional diag_script_info attributes (diagnostic specific)
;;
;; Required variable_info attributes (variable specific)
;;     none
;;
;; Optional variable_info attributes (variable specific)
;;     none
;;
;; Required variable attributes (defined in namelist)
;;    ref_model: name of reference data set (observations)
;;
;; Caveats
;;     Effect of different regridding methods not yet determined
;;
;; Modification history
;;    201704010-A_gier_be: written.
;;
;;#############################################################################

load "./interface_data/ncl.interface"

load "./interface_scripts/auxiliary.ncl"
load "./interface_scripts/data_handling.ncl"
load "./interface_scripts/messaging.ncl"

load "./diag_scripts/lib/ncl/style.ncl"
load "./diag_scripts/lib/ncl/latlon.ncl"
load "./diag_scripts/lib/ncl/statistics.ncl"

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

    ;; Read model data
    data = True
    do imod = 0, dim_MOD - 1
        data = read_data(imod, var0, field_type0)
    end do

end

begin
    vardeffile = "interface_data/" + var0 + "_info.tmp"
    loadscript(vardeffile)
end

begin
    write_references(diag_script,                  \
                     (/""/),                       \  ; authors
                     (/""/),                       \  ; contributors
                     (/""/),                       \  ; diag_references
                     (/""/),                       \  ; obs_references
                     (/"P_crescendo"/))               ; proj_references

    ;; Basename of diag_script
    diag_script_base = basename(diag_script)

    ;; Define file type
    file_type = getenv("ESMValTool_output_file_type")
    if (ismissing(file_type)) then
        file_type = "ps"
    end if

    ;; Check required diag_script_info attributes
    ;;req_atts = (/"dst_grid"/)
    ;;exit_if_missing_atts(diag_script_info, req_atts)

    ;; Check field type
    plot_fields = (/"T2Ms"/)
    if (field_type0.ne.plot_fields) then
        error_msg("f", diag_script, "", "can't use field " + field_type + \
                  " for this diag_script")
    end if
    delete(plot_fields)


    ;; Output directories
    plot_dir = getenv("ESMValTool_plot_dir")
    plot_output_dir = get_output_dir(plot_dir, diag_script_base)
    output_dir = get_output_dir(get_work_dir(), diag_script_base)
    xml_name = getenv("ESMValTool_xml_name")
    suffix = get_file_suffix(xml_name, 0)
    name = str_sub_str(xml_name, suffix, "")

end
