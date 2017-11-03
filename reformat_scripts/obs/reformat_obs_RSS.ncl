;;#############################################################################
;; REFORMAT SCRIPT FOR RSS TLT  REANALYISIS DATA
;;#############################################################################
;;
;; Tier
;;    Tier 2: other freely available dataset.
;;
;; Source
;;    http://www.remss.com/measurements/upper-air-temperature
;;    ftp://ftp.remss.com/msu/data/netcdf/
;;
;; Last access
;;    20170413
;;
;; Download and processing instructions
;;    A registration is required for downloading the data, but no licence
;;    agreement necessary. Download netcdf FTP files.     
;;
;; Caveats
;;
;; Modification history
;;    20170413-A_gier_be: written.
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
    YEAR1 = 1978
    YEAR2 = 2017

    ;; Selected variable (standard name)
    VAR = "TLT"
     
    ;; Name in the raw data
    NAME = "brightness_temperature"

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

    f = addfile(INDIR + "RSS_Tb_Maps_ch_TLT_V3_3.nc", "r")
    output = f->$NAME$

    ;; Calculate days per month
    date = cd_calendar(output&months, 0)
    dpm = days_in_month(toint(date(:, 0)), toint(date(:, 1)))
    dpmc = conform(output, dpm, 0)
    
    info_output(dimsizes(date(:, 0)) + " and " + 12 * (YEAR2 - YEAR1 + 1), verbosity, 1)
    ;; Check time range
    if (dimsizes(date(:, 0)).ne.12 * (YEAR2 - YEAR1 + 1)) then
        error_msg("f", diag_script, "", "incorrect number of timesteps")
    end if

    ;; Set variable attributes
    tmp = var_attrib(output, VAR)
    delete(output)
    output = tmp
    delete(tmp)

    ;; Format time coordinate
    ctime = time_attrib(output&months)
    delete(output&months)
    output!0 = "time"
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
    gAtt@source        = "ftp://ftp.remss.com/msu/data/netcdf/"
    gAtt@reference = "Mears, C. A. and F. J. Wentz, (2009) Construction of" + \
                     " the RSS V3.2 Lower Tropospheric Dataset From the" + \
                     " MSU and AMSU Microwave Sounders, Journal of" + \
                     " Atmospheric and Oceanic Technology, 26, 1493-1509."
    gAtt@title = OBSNAME + \
        " reanalysis data reformatted for the ESMValTool"
    gAtt@conventions = "CF/CMOR"
    gAtt@version = "v3.3"

    ;; Outfile
    fout = OUTDIR + "OBS_" + OBSNAME + "_reanaly_1_" + FIELD + \
        "_" + VAR + "_" + YEAR1 + "01-" + YEAR2 + "12.nc"

    ;; Write variable
    write_nc(fout, VAR, output, gAtt)
    delete(gAtt)
    delete(output)

    leave_msg(diag_script, "", 1)

end     
