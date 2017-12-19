;;#############################################################################
;; DIAGNOSTIC SCRIPT for reproducing IPCC ch. 9 fig. 9.9
;; Author: Bettina Gier (DLR, Germany)
;; CRESCENDO project
;;#############################################################################
;;
;; Description
;;    Compute decadal trends in tropical prw and TLT. Like IPCC ch. 9 fig 9.9
;;
;; Required diag_script_info attributes (diagnostics specific)
;;     diag_script_info@region: Region to be averaged over
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
;;     none
;;
;; Caveats
;;     Effect of different regridding methods not yet determined
;;
;; Modification history
;;    20171023-A_laue_ax: replaced annual means with low-pass filter
;;                        bug fix: using prw average instead of sum
;;    20171020-A_laue_ax: added tags for reporting
;;    20170410-A_gier_be: written.
;;
;;#############################################################################
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

load "./interface_data/ncl.interface"

load "./interface_scripts/auxiliary.ncl"
load "./interface_scripts/data_handling.ncl"
load "./interface_scripts/messaging.ncl"

load "./diag_scripts/lib/ncl/style.ncl"
load "./diag_scripts/lib/ncl/latlon.ncl"
load "./diag_scripts/lib/ncl/statistics.ncl"

load "./diag_scripts/lib/ncl/meta_data.ncl"

begin    
    verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(diag_script, "", 2)
    info_output("++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)
    info_output(diag_script + " (var: " + variables(0) + ")", verbosity, 1)
    info_output("++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)

    dim_MOD = dimsizes(models@name)
    dim_VAR = dimsizes(variables)
    var0 = variables(0)
    var1 = variables(1)
    field_type0 = field_types(0)
    field_type1 = field_types(1)

    climofiles = new(dim_MOD * dim_VAR, string)

    i = 0
    do jj = 0, dim_VAR - 1
        do ii = 0, dim_MOD - 1
            climofiles(i) = \
                interface_get_inpaths(ii) + "/" \
                + interface_get_infile(variables(jj), field_types(jj), ii)
            i = i + 1
        end do
    end do

    alltags = array_append_record(tags, (/"PT_scatter", "ST_trend", "DM_trop"/), 0)
end

begin
    vardeffile = "interface_data/" + var0 + "_info.tmp"
    loadscript(vardeffile)
end

begin
    authors = "A_gier_be"

    write_references(diag_script,                  \
                     authors,                      \  ; authors
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
    req_atts = (/"region"/)
    exit_if_missing_atts(diag_script_info, req_atts)

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

    ;; Set region
    region = select_region(diag_script_info@region)

    ;; GLobal array for collecting all data
    all_data_A = new((/dim_MOD/), float)
    all_data_A!0 = "model"
    all_data_A&model = models@name

    all_data_B = new((/dim_MOD/), float)
    all_data_B!0 = "model"
    all_data_B&model = models@name

    ;; Loop over models
    do imod = 0, dim_MOD - 1

        info_output("Processing " + all_data_A&model(imod), verbosity, 1)

        ;; First variable (TLT)
        if (var_attr_exclude(0).ne.model_attr_id(imod)) then

            ;; Read data
            dataA = read_data(imod, var0, field_type0)
            var = extract_data_extended(imod, var0, dataA, -1, 0, 0)
            delete(dataA)

            ;; Compute and apply NCL standard ocean mask
            f = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc", "r")
            a = f->LSMASK
            sftof = byte2flt(landsea_mask(a, var&lat, var&lon))
            sftof = where(sftof.eq.0., 1., 0.)
            sftof = where(sftof.eq.0, sftof@_FillValue, sftof)
            varA = var
            varA = varA * conform(var, sftof, (/1, 2/))
            delete(var)
;            filename = models@name(imod) + ".nc"
;            system ("rm " + filename)
;            debugfile = addfile(filename, "c")
;            debugfile->mask = sftof
;            debugfile->var = varA
;            delete(filename)
;            delete(debugfile)
            delete(sftof)

            yyyymmdd = cd_calendar (varA&time, 2)
            if isatt(yyyymmdd, "calendar") then
                ; The "proleptic_gregorian" calendar is not supported by the NCL
                ; function "yyyymmdd_to_yyyyfrac". The error by changing the calendar
                ; to standard is negligible.
                if (yyyymmdd@calendar.eq."proleptic_gregorian") then
                    yyyymmdd@calendar = "standard"
                end if
            end if
            regtime = yyyymmdd_to_yyyyfrac(yyyymmdd, 0)
            delete(yyyymmdd)

            ;; Computing area-weighted average
            avg_A = area_operations(varA, region(0), region(1), region(2), \
                                    region(3), "average", True)
            delete(varA)

            ; apply low-pass filter with a cutoff period of 12 months

            ihp   = 0                             ; low_pass
            sigma = 1.0                           ; Lanczos sigma
            nWgt  = 25                            ; loose 12 months each end                            
            fca   = 1./12.                        ; 12 months
            wgt   = filwgts_lanczos(nWgt, ihp, fca, -999., sigma)

            filtered_A = wgt_runave(avg_A, wgt, 0)

            delete(avg_A)
            delete(wgt)

            rc = regline(regtime, filtered_A)
            ;; Assign to global array
            all_data_A(imod) = tofloat(rc)
            delete(filtered_A)
            delete(rc)
            delete(regtime)

        end if

        ;; Second variable (prw)
        if (var_attr_exclude(1).ne.model_attr_id(imod)) then

            ;; Read data
            dataB = read_data(imod, var1, field_type1)
            var = extract_data_extended(imod, var1, dataB, -1, 0, 0)
            delete(dataB)

            ;; Compute and apply NCL standard ocean mask
            f = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc", "r")
            a = f->LSMASK
            sftof = byte2flt(landsea_mask(a, var&lat, var&lon))
            sftof = where(sftof.eq.0., 1., 0.)
            sftof = where(sftof.eq.0, sftof@_FillValue, sftof)
            varB = var
            varB = varB * conform(varB, sftof, (/1, 2/))
            delete(var)
            delete(sftof)

            yyyymmdd = cd_calendar (varB&time, 2)
            if isatt(yyyymmdd, "calendar") then
                ; The "proleptic_gregorian" calendar is not supported by the NCL
                ; function "yyyymmdd_to_yyyyfrac". The error by changing the calendar
                ; to standard is negligible.
                if (yyyymmdd@calendar.eq."proleptic_gregorian") then
                    yyyymmdd@calendar = "standard"
                end if
            end if
            regtime = yyyymmdd_to_yyyyfrac(yyyymmdd, 0)
            delete(yyyymmdd)

            ;; Computing area-weighted average (prw)
            avg_B = area_operations(varB, region(0), region(1), region(2), \
                                    region(3), "average", True)
            delete(varB)

            ; apply low-pass filter with a cutoff period of 12 months

            ihp   = 0                             ; low_pass
            sigma = 1.0                           ; Lanczos sigma
            nWgt  = 25                            ; loose 12 months each end                            
            fca   = 1./12.                        ; 12 months
            wgt   = filwgts_lanczos(nWgt, ihp, fca, -999., sigma)

            filtered_B = wgt_runave(avg_B, wgt, 0)

            delete(avg_B)
            delete(wgt)

            ;; Computing anomaly in percent
            anomaly_B = (filtered_B - avg(filtered_B)) / avg(filtered_B) * 100.
            rc = regline(regtime, anomaly_B)
            ;; Assign to global array
            all_data_B(imod) = tofloat(rc)
            delete(rc)
            delete(anomaly_B)
            delete(regtime)
            delete(filtered_B)
        end if

    end do

    ;; UAH uses the same prw as RSS, thus is excluded in prw read-in
    if any(models@name.eq."UAH") then
        idx_UAH = ind(var_attr_exclude(1).eq.model_attr_id)
        idx_RSS = ind(models@name.eq."RSS")
        all_data_B(idx_UAH) = all_data_B(idx_RSS)
    end if


;; Plotting routine from Hongmei, only changed where necessary

;*************************************************
; Scatter plot of lower tropospheric temperature vs precipitable water
; Author: H. Li, MPI-M
; Date: May 2013
;************************************************

;----------------------------------------------------------------------
; Data is hard-coded here. It was originally calculated by running
; a compute-intensive script.
;----------------------------------------------------------------------

;; models = (/"NorESM1-ME","NorESM1-M","MRI-CGCM3","MPI-ESM-MR","MPI-ESM-LR", \
;;            "MIROC-ESM-CHEM", "MIROC-ESM","MIROC5","IPSL-CM5B-LR", \
;;            "IPSL-CM5A-MR","IPSL-CM5A-LR","INM-CM4", "HadGEM2-ES", \
;;            "HadGEM2-CC","GISS-E2-R","GISS-E2-H","GFDL-ESM2M","GFDL-ESM2G", \
;;            "GFDL-CM3","FIO-ESM","FGOALS-g2","EC-EARTH","CSIRO-Mk3.6.0", \
;;            "CNRM-CM5","CMCC-CMS", "CMCC-CESM","CESM1(CAM5)","CESM1(BGC)", \
;;            "CCSM4","CanESM2","BCC-CSM1.1(m)", "BCC-CSM1.1","ACCESS1.3", \
;;            "ACCESS1.0","RSS V3.3-SSM/I","UAH V5.4-SSM/I","MERRA", \ 
;;            "JRA-25","ERA-Interim"/)
;;  data = asciiread("figure_9_9.rev2_HL_ncl.txt",(/39,2/),"float")
;;  tr_tlt = data(38:0,0)
;;  tr_pr = data(38:0,1)
  model_names = models@name
  tr_tlt = all_data_A*10
  tr_pr = all_data_B*10
  nmodel = dimsizes(model_names)

;;  print(tr_tlt)
;;  print(tr_pr)

  markertab = project_style(diag_script_info, "markers")
  colortab  = project_style(diag_script_info, "colors")

  markers = new(dim_MOD, integer)
  colors  = new(dim_MOD, string)

  colors = colortab
  markers = markertab


;----------------------------------------------------------------------
; Open a PS file and define the color table
;----------------------------------------------------------------------
  
  outfile = plot_output_dir + diag_script_base + "_" + diag_script_info@region + \
            "." + file_type
  ;wks = gsn_open_wks ("pdf","Fig9_09")
  wks = gsn_open_wks(file_type, outfile)
  ;gsn_define_colormap(wks,"IPCC5")
  plot = new(2,graphic)

;----------------------------------------------------------------------
; Set some plot resources
;----------------------------------------------------------------------
  res                    = True

  res@gsnDraw            = False       ; Don't draw plot
  res@gsnFrame           = False       ; Don't advance frame

;
; Don't use gsnMaximize. Instead, we will control 
; exactly where plot goes.
;
  res@vpXF              = 0.14
  res@vpYF              = 0.9
  res@vpWidthF          = 0.5
  res@vpHeightF         = 0.5

;---Some calculations to get min/max for Y axis.

;---Resources for X/Y axes.
  res@trYMinF            =  0.               ; min value on y-axis
  res@trYMaxF    = 3.5
 ; res@trYMaxF            =  3.5               ; max value on y-axis
  res@trXMinF            =  0.                  ; min value on x-axis
  res@trXMaxF            =  0.6                  ; max value on x-axis

;---XY marker resources
  res@xyMarkLineMode     = "Markers"
  res@xyMarkerSizeF      = 0.013
  res@xyMarkerThicknessF = 4.
  res@xyMonoMarkerColor  = False          ; Allow different colors for markers
;res@xyMarkerColors     = (/22,22,20,19,19,18,18,18,17,17,17,16,15,15,14,14, \
;                         13,13,13,12,11,10,9,8,7,7,6,6,6,5,3,3,2,2,1,1,1,1,1/)
;res@xyMarkers          = (/6,4,4,6,4,7,6,4,7,6,4,4,6,4,6,4,7,6,4,4,4,4,4,4, \
;                         6,4,7,6,4,4,6,4,6,4,12,5,7,6,4/)
  res@xyMarkerColors    = colors
  res@xyMarkers         = markers

;---Tickmark resources
  res@tmXBMode         = "Manual"
  res@tmYLMode         = "Manual"
  res@tmXBTickSpacingF = .1
  res@tmYLTickSpacingF = .5

;---Title resources
  res@tiMainOffsetXF     = -0.058
  res@tiMainFontHeightF  = 0.015
  res@tiMainString       = ""
  res@tiXAxisString      = "Trend in TLT (~S~o~N~C/decade)"
  res@tiYAxisString      = "Trend in Precipitable Water (%/decade)"
  res@tmYLLabelFontHeightF  = 0.015
  res@tmXBLabelFontHeightF  = 0.015
  res@tiXAxisFontHeightF    = 0.015
  res@tiYAxisFontHeightF    = 0.015
  res@tmXBMajorOutwardLengthF = 0.
  res@tmYLMajorOutwardLengthF = 0.
  res@tmXBMinorOutwardLengthF = 0.
  res@tmYLMinorOutwardLengthF = 0.
  res@tmXBMajorLengthF = 0.006
  res@tmYLMajorLengthF = 0.006
  res@tmXBMinorLengthF = 0.003
  res@tmYLMinorLengthF = 0.003

;----------------------------------------------------------------------
; Create the scatter plot, but don't draw it yet.
;----------------------------------------------------------------------
  plot0 = gsn_csm_xy (wks,transpose((/tr_tlt, tr_tlt/)), \
                         transpose((/tr_pr, tr_pr/)),res)

;----------------------------------------------------------------------
; Draw time stamp first, then the plot.
;----------------------------------------------------------------------
 ;; Regression line needs at least 3 models
 rc  = regline(tr_tlt,tr_pr)
 ;;print(rc)
 x = (-0.5-rc@yintercept)/rc
 y = rc*0.6 + rc@yintercept

 xx = (/x,0.6/)
 yy = (/-0.5,y/)
 res_lines                   = True        ; polyline mods desired
 res_lines@gsLineColor       = "black"   ; line color
 res_lines@gsLineThicknessF  = 2.           ; line thicker
 res_lines@gsLineDashPattern = 0.   ; dash pattern
 res_lines@tfPolyDrawOrder = "PreDraw"
 dum0 = gsn_add_polyline(wks,plot0,xx,yy,res_lines)
 
 draw(plot0)

;----------------------------------------------------------------------
; Draw some markers and text.
;----------------------------------------------------------------------
  ;markercolors     = (/22,22,20,19,19,18,18,18,17,17,17,16,15,15,14,14, \
  ;                  13,13,13,12,11,10,9,8,7,7,6,6,6,5,3,3,2,2,1,1,1,1,1/)
  ;markers1          = (/6,4,4,6,4,7,6,4,7,6,4,4,6,4,6,4,7,6,4,4,4,4,4,4,\
  ;                      6,4,7,6,4,4,6,4,6,4,12,5,7,6,4/)
  markercolors = colors;(::-1)
  markers1 = markers;(::-1)
  m_names = model_names;(::-1)
  ;ypos2    = fspan(0.39,0.89,20)
  ypos2 = fspan(0.89, 0.39, 20)

  mkres               = True         ; Marker resources
  txres               = True         ; Text resources
  txres@txFontHeightF = 0.01
  txres@txJust        = "CenterLeft"

  half = toint(nmodel/2)

  do i = 0,half
     mkres@gsMarkerThicknessF = 4.
     mkres@gsMarkerSizeF      = 0.013
     mkres@gsMarkerIndex      = markers1(i)
     mkres@gsMarkerColor      = markercolors(i)  
     gsn_polymarker_ndc(wks,           0.69,ypos2(i),mkres)
     gsn_text_ndc      (wks,m_names(i),0.71,ypos2(i),txres)
  end do
  do i = half+1,nmodel-1
     mkres@gsMarkerThicknessF = 3.5
     mkres@gsMarkerSizeF      = 0.012
     mkres@gsMarkerIndex      = markers1(i)
     mkres@gsMarkerColor      = markercolors(i)  
     gsn_polymarker_ndc(wks,           0.85,ypos2(i-half-1),mkres)
     gsn_text_ndc      (wks,m_names(i),0.87,ypos2(i-half-1),txres)
  end do

  txres@txFontHeightF = 0.016
  ;gsn_text_ndc (wks,"Slope = 5.7% / ~S~o~N~C",0.45,0.45,txres)
  gsn_text_ndc (wks,"Slope = " + sprintf("%2.1f", rc) + "% / ~S~o~N~C", \
                0.45,0.45,txres)

  frame(wks)    ; Now advance the frame.

    ; =========================================================================
    ; add meta data to plot (for reporting)

    caption = "Scatter plot of decadal trands in tropical " \
              + variables(0) + " as a function of " + variables(1) \
              + " over the world's oceans (similar to IPCC ch. 9 fig. 9.9)."
    id = diag_script

    ESMValMD(outfile, alltags, caption, id, variables, models@name, \
             climofiles, diag_script, authors)
    ; =========================================================================

end
