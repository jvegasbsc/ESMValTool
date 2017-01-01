;;#############################################################################
;; DIAGNOSTIC SCRIPT for reproducing IPCC ch. 9 fig. 9.6
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
;;    20170404-A_gier_be: written.
;;
;;#############################################################################

load "./interface_data/ncl.interface"

load "./interface_scripts/auxiliary.ncl"
load "./interface_scripts/data_handling.ncl"
load "./interface_scripts/messaging.ncl"

load "./diag_scripts/lib/ncl/style.ncl"
load "./diag_scripts/lib/ncl/latlon.ncl"
load "./diag_scripts/lib/ncl/statistics.ncl"

load "./plot_scripts/ncl/scatterplot.ncl"

begin    
    verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(diag_script, "", 2)
    info_output("++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)
    info_output(diag_script + " (var: " + variables(0) + ")", verbosity, 1)
    info_output("++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)

    ref_model = var_attr_ref(0)  ; produces an error if ref_model is not set in the namelist
    info_output("reference model = " + ref_model, verbosity, 1)
    idx_ref = get_ref_model_idx(models, var_attr_ref(0))

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

    ;; Set time range
    ;start_year = min(stringtoint(models@start_year))
    ;end_year = max(stringtoint(models@end_year))

    ;; Set destination grid
    grid_resolution = "5degx4deg"
    dst_grid = "./diag_scripts/aux/CDO/cdo_dst_grid_g050x040"

;; ------------------------------------------------------------------------------
;; ------------------- Compute correlation for one var --------------------------
;; ------------------------------------------------------------------------------

    ;; Reference model
    ref_ind = get_ref_model_idx(models, var_attr_ref(0))
    if (ismissing(ref_ind)) then
        error_msg("f", diag_script, "", "reference model missing")
    end if
    ref_model = models@name(ref_ind)
    mod_ind_woref = ind(ref_model.ne.models@name)

    ;; Make output array
    all_cor = new((/dim_MOD-1/), float)
    all_cor!0 = "models"
    all_cor&models = models@name(mod_ind_woref)
    ;all_cor&models@project = models@project(mod_ind_woref)
    ;; Pass on alt models
    if isatt(diag_script_info, "alt_obs").and. \
        any(models@name.eq.diag_script_info@alt_obs) then
        all_cor@alt_obs =  diag_script_info@alt_obs
    else
        all_cor@alt_obs = "none"
    end if

    ;; Loop over models
    do iloop = 0, dim_MOD - 1
        ;; Loop over models, with ref model processed first
        model_ind = array_append_record(ref_ind, mod_ind_woref, 0)
        imod = model_ind(iloop)
        info_output("Processing " + models@name(imod), verbosity, 0)

        ;; Extract model corresponding data
        var  = extract_data_extended(imod, var0, data, -1, 0, 0)

        ;; Regridding
        regr_file = basename(diag_script) + "_" + models@name(imod)+ "_" + \
                      var0 + "_" + models@start_year(imod) + "-" + \
                      models@end_year(imod) + "_grid_" + \
                      grid_resolution
        regr_path = wrk_dir + "/CDO/" + regr_file + ".nc"            
        if (isfilepresent(regr_path)) then
            var_regrid = ncdf_read(regr_path, var0)
        else
            var_regrid = cdo_remapdis(var0, var, imod, data, \
                          dst_grid, regr_file)                
        end if
        var_regrid := var_regrid(time|:, lat|:, lon|:)
        delete(var)
        var = var_regrid
        delete(var_regrid)

        ;; Calculate annual mean
        data_yearly = \
              time_operations(var, toint(models@start_year(imod)), \
                                   toint(models@end_year(imod)), \
                                    "average", "yearly", True)
                                    
        ;; Mean over the years
        mean_years = dim_avg_n_Wrap(data_yearly, 0)
        delete(data_yearly)

        ;; Compute centred pattern correlation (ref_model will be first so no error)
        if models@name(imod).eq.ref_model then
            data_ref = mean_years
        else
            ;; Since ref model processed first, move all models up one spot
            all_cor(iloop-1) = pattern_cor(data_ref, mean_years, 1.0, 0)
        end if
        delete(mean_years)
        delete(var)

    end do

    ;; Write data
    ncdf_dir = output_dir + name + "_" + var0 + ".nc"
    all_cor@corvar = var0
    all_cor@var = "cor"
    all_cor@diag_script = diag_script
    all_cor@ncdf = ncdf_dir
    ncdf_outfile = ncdf_write(all_cor, ncdf_dir)

    ;; Write temporary correlation collection list
    temp_dir = get_output_dir(getenv("ESMValTool_wrk_dir"), "temp") + \
        name + "_patterncor.nc"

    if (isfilepresent_esmval(temp_dir)) then
        ; Read correlation list
        temp_file = addfile(temp_dir, "r")
        temp_list = temp_file->temp_list
        temp_list_string = tostring(temp_list)
        ; Add new entry
        temp_list_new_string = array_append_record(temp_list_string, \
                                                   ncdf_dir, 0)
        temp_list_new = tochar(temp_list_new_string)
        ; Remove previous file
        system("rm -f " + temp_dir)
    else
        ncdf_char = tochar(ncdf_dir)
        temp_list_new = new((/1, dimsizes(ncdf_char)/), character)
        temp_list_new(0, :) = ncdf_char
    end if

    ;; Create new file and add list
    temp = addfile(temp_dir, "c")
    temp->temp_list = temp_list_new


;; ---------------------------------------------------------------------------
;; ------------------ Collecting Data ----------------------------------------
;; ---------------------------------------------------------------------------

    if diag_script_info@collect then
        if (isfilepresent_esmval(temp_dir)) then
            ;; Read correlation list
            files = addfile(temp_dir, "r")
            filename_list = files->temp_list
            file_name_list = tostring(filename_list)
            file_name_list@_FillValue = default_fillvalue(typeof(file_name_list))
            delete(files)
            delete(filename_list)

            ;; Remove duplicates
            do ii = 0, dimsizes(file_name_list)-1
                curr = file_name_list(ii)
                idx = ind(curr.eq.file_name_list)
                if (dimsizes(idx).gt.1) then
                    file_name_list(idx) = file_name_list@_FillValue
                    file_name_list(ii) = curr
                end if
                delete(curr)
                delete(idx)
            end do
            file_list = file_name_list(ind(.not.ismissing(file_name_list)))
            delete(file_name_list)
            var_collect = new(dimsizes(file_list), string)
            alt_obs = new(dimsizes(file_list), string)

            ;; Loop over files, append data
            do ii=0, dimsizes(file_list)-1
                data_temp = ncdf_read(file_list(ii), "cor")
                var_collect(ii) = data_temp@corvar
                alt_obs(ii) = data_temp@alt_obs

                ;; Make 2D array to store all data
                if (.not.isdefined("data_all")) then
                    data_all = new((/dimsizes(data_temp), dimsizes(file_list)/), float)
                    data_all(:, 0) = data_temp

                    data_all!0 = "models"
                    data_all&models = data_temp&models
                else
                    ;; If model coordinates are identical
                    if (dimsizes(data_temp&models).eq.dimsizes(data_all&models)) \
                        .and. all(data_temp&models.eq.data_all&models) then
                        data_all(:, ii) = (/data_temp/)
                    else
                        ;; Loop over models in new data entry
                        do imod_temp = 0, dimsizes(data_temp&models) -1
                            ;; If current model is not already part
                            ;; of the model coordinate
                            if (.not.any(data_temp&models(imod_temp) .eq. \
                                data_all&models)) then
                                ;; Append record for model(imod)
                                data_new = extend_var_at(data_all, 0, \
                                                        dimsizes(data_all&models))
                                data_new(dimsizes(data_all&models), ii) = \
                                                         (/data_temp(imod_temp)/)
                                data_new&models(dimsizes(data_all&models)) = \
                                                 (/data_temp&models(imod_temp)/)
                                delete(data_all)
                                data_all = data_new
                                delete(data_new)
                            else
                                ;; Loop over models of data
                                do imod = 0, dimsizes(data_all&models)-1
                                    ;; if neq data model is similar to current
                                    ;; entry, write data entry
                                    if (data_all&models(imod).eq. \
                                        data_temp&models(imod_temp)) then
                                        data_all(imod_temp, ii) = \
                                                 (/data_temp(imod_temp)/)
                                    end if
                                end do
                            end if
                        end do                                  
                    end if  
                end if
                delete(data_temp)
            end do
            data_all!1 = "vars"
            data_all&vars = var_collect
            delete(var_collect)
        else
            error_msg("f", diag_script, "", "cannot find correlation list")
        end if



        ;; Remove temporary correlation list
        system("rm -f " + temp_dir)

;; -------------------------------------------------------------------------
;; ----------------- Interim Functions (unused atm) ------------------------
;; -------------------------------------------------------------------------

undef("get_unique_entries")
function get_unique_entries (array)
;;
;; Arguments:
;;      array: 1D array
;;
;; Return value: 1D array of unique entries in array
;;
;; Modification history:
;;    20170406-A_gier_be: written.
local dummy_array, unique_new, new_array, nodupes
begin
    dummy_array = array
    do while (dimsizes(dummy_array).ne.0)
        if (.not.isdefined("unique")) then
            unique = dummy_array(0)
        else
            unique_new = array_append_record(unique, dummy_array(0), 0)
            delete(unique)
            unique = unique_new
            delete(unique_new)
        end if
        nodupes = ind(dummy_array.ne.dummy_array(0))
        ;; Missing value index are dim 1 and would give an error
        if (dimsizes(dummy_array).eq.dimsizes(ind(dummy_array.eq.dummy_array(0)))) then
            break
        end if
        new_array = dummy_array(nodupes)
        delete(nodupes)
        delete(dummy_array)
        dummy_array = new_array
        delete(new_array)
    end do
    return(unique)
end   

;; ---------------------------------------------------------------------------
;; ---------------------------- Plotting -------------------------------------
;; ---------------------------------------------------------------------------

        ;; Calculating necessary values
        ;; Number of Projects needed to determine span
        ;; For now just CMIP5
        projects = "CMIP5"
        nr_projects = dimsizes(projects)
        x_val = ispan(1, dimsizes(data_all(0, :)), 1)
        x_val = x_val * nr_projects
        
        ;; Mean and Median of Ensemble - without alt obs
        obs_ind = get1Dindex(data_all&models, alt_obs)
        ex_ind = obs_ind(ind(.not.ismissing(obs_ind)))
        mods_ind = ispan(0, dimsizes(data_all&models)-1, 1)
        ex_ind@_FillValue = default_fillvalue("integer")
        mods_ind@_FillValue = default_fillvalue("integer")
        mod_ind = get1Dindex_Collapse(mods_ind, ex_ind)
        delete(mods_ind)
        delete(ex_ind)
        means = dim_avg_n(data_all(mod_ind, :), 0)
        median = dim_median_n(data_all(mod_ind, :), 0)
       
        
        ;; Plotting preparation
        outfile = plot_output_dir + name + "_patterncor." + file_type
        wks = gsn_open_wks(file_type, outfile)
        wks@fullname = outfile    

        ;; Set one marker size for all plot markers
        marker_size = 0.02

        res = True
        res@gsnDraw          = False
        res@gsnFrame         = False
        ;;res@gsnMaximize      = True
        res@xyMarkLineMode = "Markers"
        res@xyMarkerSizeF = marker_size
        res@xyMarker = 4
        res@tiYAxisString = "Correlation"
        res@trYMinF = 0.5
        res@trYMaxF = 1
        res@trXMinF = 0
        res@trXMaxF = (dimsizes(data_all(0,:))+1)*nr_projects
        res@tmXBLabels = data_all&vars
        res@tmXBValues = x_val
        res@tmXBMode = "Explicit"
        res@tmYLMode = "Manual"
        res@tmYLTickSpacingF = 0.05
        res@tmYLLabelStride = 1
    
        ;; Resources for model lines
        res_lines                   = True                  ; polyline mods desired
        res_lines@gsLineDashPattern = 0.                    ; solid line
        res_lines@gsLineThicknessF  = 0.5                   ; line thickness
        res_lines@tfPolyDrawOrder = "PreDraw"

        ;; Resources for mean lines
        res_mlines                   = True                  ; polyline mods desired
        res_mlines@gsLineDashPattern = 0.                    ; solid line
        res_mlines@gsLineThicknessF  = 2                     ; line thicker
        res_mlines@tfPolyDrawOrder = "PreDraw"
        
        ;; Resources for obs data markers
        res_circ = True
        res_circ@gsMarkerIndex = 16
        res_circ@gsMarkerColor = "green"
        res_circ@gsMarkerSizeF = marker_size
        res_circ@gsMarkerOpacityF = 0.4

        ;; Resources for white markers below median
        res_circw = True
        res_circw@gsMarkerIndex = 16
        res_circw@gsMarkerColor = "white"
        res_circw@gsMarkerSizeF = marker_size
        res_circw@tfPolyDrawOrder = "PreDraw"
        
        ;; Resources for legend text
        res_text                    = True                  ; text mods desired
        res_text@txFontHeightF      = 0.02                  ; change text size
        res_text@txJust             = "CenterLeft"          ; text justification
        
           
        ;; gsn_csm_xy won't plot just one point -> double up on it
        if dimsizes(x_val).eq.1 then
            x_val_new = (/x_val, x_val/)
            median_new = (/median, median/)
            delete([/x_val, median/])
            x_val = x_val_new
            median = median_new
            delete([/x_val_new, median_new/])
        end if

        plot = gsn_csm_xy(wks, x_val, median, res)
        
        ;; add lines for individual models
        do ivar = 0, dimsizes(data_all(0, :))-1
            do iproj = 0, dimsizes(projects)-1 ;; Placeholder for now
                proj_center = ivar*nr_projects + iproj + 1 
                xx = (/proj_center-0.15, proj_center+0.15/)
                ;; Plot lines for mean
                xx_mean = (/xx(0)-0.1, xx(1)+0.1/)
                yy_mean = (/means(ivar), means(ivar)/)
                plot@$unique_string("dum")$ = gsn_add_polyline \
                     (wks, plot, xx_mean, yy_mean, res_mlines)
                do imod = 0, dimsizes(data_all(:,0))-1
                    ;; Don't plot obs as lines
                    if (.not.ismissing(data_all(imod, ivar))) then
                        if (data_all&models(imod).ne.alt_obs(ivar)) then
                            yy = (/data_all(imod, ivar), data_all(imod, ivar)/)
                            plot@$unique_string("dum")$ = gsn_add_polyline \
                                            (wks, plot, xx, yy, res_lines)                        
                        end if
                    end if
                end do
            end do
            plot@$unique_string("dum_ci")$ = gsn_add_polymarker \
                   (wks, plot, x_val(ivar), median(ivar), res_circw)
            if (alt_obs(ivar).ne."none") then
                ;; Plot obs as circles
                plot@$unique_string("dum_circ")$ = gsn_add_polymarker \
                  (wks, plot, x_val(ivar), data_all(obs_ind(ivar), ivar), res_circ)
            end if
        end do
        
        ;; Draw Legend manually
        plot@$unique_string("dum_l")$ = gsn_add_text(wks, plot, "CMIP5", \
                                                  x_val(0), 0.65, res_text)
        res_text@txFontColor = "green"
        plot@$unique_string("dum_l")$ = gsn_add_text(wks, plot, "OBS", \
                                                  x_val(0), 0.62, res_text)



        draw(plot)
        frame(wks)
    end if

    leave_msg(diag_script, "", 2)

end
