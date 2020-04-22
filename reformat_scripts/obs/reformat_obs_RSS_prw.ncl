;;#############################################################################
;; REFORMAT SCRIPT FOR RSS prw REANALYISIS DATA
;;#############################################################################
;;
;; Tier
;;    Tier 2: other freely available dataset.
;;
;; Source
;;    ftp://ftp.remss.com/vapor/monthly_1deg/
;;
;; Last access
;;    20170406
;;
;; Download and processing instructions
;;    A registration is required for downloading the data, but no licence
;;    agreement necessary. Download vapor, here ncdf3 file is used.     
;;
;; Caveats
;;
;; Modification history
;;    20170419-A_gier_be: written.
;;
;;#############################################################################

begin

    ;; Source name
    OBSNAME = "RSS"

    ;; Tier
    TIER = 2

    ;; Input dir (raw data)
    INDIR = getenv("ESMValTool_RAWOBSPATH") + "/Tier" + \
        TIER + "/" + OBSNAME + "/"
    
    ;; Output dir (CMOR-ized data)
    OUTDIR = getenv("ESMValTool_OBSPATH") + "/Tier" + \
        TIER + "/" + OBSNAME + "/"

    ;; Period
    YEAR1 = 1988
    YEAR2 = 2016

    ;; Selected variable (standard name)
    VAR = "prw"
     
    ;; Name in the raw data
    NAME = "precipitable_water"

    ;; Initialize global variable
    FIELD = "T2Ms"

end

load "./interface_scripts/messaging.ncl"
load "./reformat_scripts/constants.ncl"
load "./reformat_scripts/obs/reformat_obs_func.ncl"

begin
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    diag_script = "reformat_obs_RSS.ncl"
    enter_msg(diag_script, "", 1)

    ;; Create output dir
    if (.not.isfilepresent(OUTDIR)) then
        system("mkdir -p " + OUTDIR)
    end if

    f = addfile(INDIR + "tpw_v07r01_198801_201702.nc3.nc", "r")
    output = f->$NAME$

    ;; Cut off the first 2 months of 2017
    output_new = output(:dimsizes(output&time)-3, :, :)
    delete(output)
    output = output_new

    ;; Calculate days per month
    date = cd_calendar(output&time, 0)
    dpm = days_in_month(toint(date(:, 0)), toint(date(:, 1)))
    dpmc = conform(output, dpm, 0)
    
    info_output(dimsizes(date(:, 0)) + " and " + 12 * (YEAR2 - YEAR1 + 1), verbosity, 1)
    ;; Check time range
    if (dimsizes(date(:, 0)).ne.12 * (YEAR2 - YEAR1 + 1)) then
        error_msg("f", diag_script, "", "incorrect number of timesteps")
    end if

    ;; Values of -500 are ice -> set these to missing values
    output = where(output.eq.(-500), output@_FillValue, output)

    ;; Set variable attributes
    tmp = var_attrib(output, VAR)
    delete(output)
    output = tmp
    delete(tmp)

    ;; Format time coordinate
    ctime = time_attrib(output&time)
    delete(output&time)
    output&time = ctime
    delete(ctime)

    ;; Format latitude coordinate
    output!1 = "lat"
    if (isMonotonic(output&lat) .eq. 0) then
        error_msg("f", diag_script, "", \
                  "non-monotonic latitude coordinate")
    end if
    if (isMonotonic(output&lat) .eq. -1) then
        output = output(:, ::-1, :)
    end if
    clat = lat_attrib(output&lat)
    delete(output&lat)
    output&lat = clat
    delete(clat)

    ;; Format longitude coordinate
    output!2 = "lon"
    if (isMonotonic(output&lon) .eq. 0) then
        error_msg("f", diag_script, "", \
                  "non-monotonic longitude coordinate")
    end if
    if (any(output&lon.lt.0.)) then
        output = lonFlip(output)
    end if
    clon = lon_attrib(output&lon)
    delete(output&lon)
    output&lon = clon
    delete(clon)

    info_output("  Range: " + min(output) + "; " + max(output), verbosity, 1)

    ;; Set global attributes
    gAtt = True
    gAtt@history       = "Created on " + systemfunc("date")
    gAtt@host          = systemfunc("echo $HOST")
    gAtt@user          = systemfunc("echo $USER")
    gAtt@period        = YEAR1 + "-" + YEAR2
    gAtt@field         = FIELD
    gAtt@tier          = TIER
    gAtt@source        = "ftp://ftp.remss.com/vapor/monthly_1deg/"
    gAtt@reference = "Remote Sensing Systems, 2012, updated 2017-02. " \
                     + "Monthly Mean Total Precipitable Water Data Set " \
                     + "on a 1 degree grid made from Remote Sensing " \
                     + "Systems Version-7 Microwave Radiometer Data, " \
                     + "[accessed on 2017-04-06]. Santa Rosa, CA, USA. " \
                     + "Available at www.remss.com"

    gAtt@title = OBSNAME + \
        " reanalysis data reformatted for the ESMValTool"
    gAtt@conventions = "CF/CMOR"
    gAtt@version = "v7.1"

    ;; Outfile
    fout = OUTDIR + "OBS_" + OBSNAME + "_reanaly_1_" + FIELD + \
        "_" + VAR + "_" + YEAR1 + "01-" + YEAR2 + "12.nc"

    ;; Write variable
    write_nc(fout, VAR, output, gAtt)
    delete(gAtt)
    delete(output)

    leave_msg(diag_script, "", 1)

end     
