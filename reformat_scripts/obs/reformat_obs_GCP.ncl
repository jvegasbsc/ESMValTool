;;#############################################################################
;; REFORMAT SCRIPT FOR ESRL OBSERVATIONAL DATA
;;#############################################################################
;;
;; Tier
;;    Tier 2: other freely-available dataset.
;;
;; Source
;;    http://www.globalcarbonproject.org/carbonbudget/
;;    http://cdiac.ornl.gov/GCP/
;;
;; Last access
;;    20170328
;;
;; Download and processing instructions
;;    Select the following filters:
;;        Category: Carbon Budget.
;;        Parameter name: fossil fuel and cement emissions
;;       	          land-use change emissions	
;;                        atmospheric growth	
;;                        ocean sink	
;;                        land sink	
;;        Type: Reanalysis.
;;        Frequency: Annual.
;;    Check the period covered by the data and give it below as YEAR1 and
;;    YEAR2.
;;
;; Modification history
;;    20170329-A_gier_be: updated to new data, now gives both nbp and fgco2
;;    20151001-A_wenz_sa: written.
;;
;;#############################################################################

begin

    ;; Source name
    OBSNAME = "GCP"

    ;; Tier
    TIER = 2

    ;; Input dir (raw data)
    INDIR = getenv("ESMValTool_RAWOBSPATH") + "/Tier" + \
        TIER + "/" + OBSNAME + "/"
    ;INDIR = "obs/GCP/"

    ;; Output dir (CMOR-ized data)
    OUTDIR = getenv("ESMValTool_OBSPATH") + "/Tier" + \
        TIER + "/" + OBSNAME + "/"
    ;OUTDIR = "obs/GCP/"

    ;; Period
    YEAR1 = 1959
    YEAR2 = 2015

    ;; Selected variable (standard name)
    VARS = (/"fgco2", "nbp"/)

    ;; Field
    FLDARR = (/"T0M", "T0M"/)

    ;; Initialize global variable
    FIELD = ""

end

;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "./interface_scripts/messaging.ncl"
load "./interface_scripts/constants.ncl"
load "./reformat_scripts/constants.ncl"
load "./reformat_scripts/obs/reformat_obs_func.ncl"

begin

    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    diag_script = "reformat_obs_GCP.ncl"
    enter_msg(diag_script, "", 1)

    ;; Create output dir
    if (.not.isfilepresent_esmval(OUTDIR)) then
        system("mkdir -p " + OUTDIR)
    end if

    ;; Create time coordinate
    timec = create_timec(YEAR1, YEAR2)
    datec = ut_calendar(timec, 1)

    ;; Reading data info
    info_output("Reading data information", verbosity, 1)

    infile = INDIR + OBSNAME + "_2016_v1.0.csv"
    ;;head = readAsciiHead(infile, "#Year")
    data = readAsciiTable(infile, 6, "float", "#Year")
    ;;delete(head)

    ;; Loop over variables
    do vID = 0, dimsizes(VARS) - 1

        ;; Set field
        FIELD = FLDARR(vID)

        info_output("", verbosity, 1)
        info_output("Processing " + VARS(vID) + " (" + FIELD + ") ", \
                    verbosity, 1)

        ;; Read time and extract full years only
        year = data(:, 0)

        if VARS(vID).eq."fgco2" then
            var = data(:,4)
        end if
        if VARS(vID).eq."nbp" then
            var = data(:,5) - data(:,2)
        end if

        ;; Check time range
        if (min(toint(year)).lt.YEAR1 .or. max(toint(year)).gt.YEAR2) then
            error_msg("f", diag_script, "", "given time range (" + \
                                       YEAR1 + "-" + YEAR2 + ")" + \
                                   " not consistent with data (" + \
                       min(toint(year)) + "-" + max(toint(year)) + ")")
        end if
        delete(year)

        ;; Convert units
        ;; [PgC y-1] --> [kg m-2 s-1]
        var = var / (3600. * 24. * 365 / 1e12)
        conversion = True
        
        ;; Create output variable
        output = new(dimsizes(timec), float)
        output!0 = "time"
        output&time = timec
        output@_FillValue = FILL

        do ii = 0, dimsizes(var) - 1
            output(6+ii*12) = (/var(ii)/)
        end do
        delete(var)

        ;; Set variable attributes (CMOR)
        tmp = var_attrib(output, VARS(vID))
        delete(output)
        output = tmp
        delete(tmp)

        ;; Reformat time coordinate
        ctime = time_attrib(output&time)
        delete(output&time)
        output&time = ctime
        delete(ctime)

        ;; Set global attributes
        gAtt = True
        gAtt@history = "Created on " + systemfunc("date")
        gAtt@host = systemfunc("echo $HOST")
        gAtt@user = systemfunc("echo $USER")
        gAtt@period = YEAR1 + "-" + YEAR2
        gAtt@field = FIELD
        gAtt@tier = TIER
        gAtt@source = "http://www.globalcarbonproject.org/carbonbudget/"
        gAtt@title = OBSNAME + \
        " data reformatted for the ESMValTool"
        gAtt@conventions = "CF/CMOR"
        
        ;; Outfile
        fout = OUTDIR + "OBS_" + OBSNAME + "_reanaly_1_" + FIELD + \
            "_" + VARS(vID) + "_" + YEAR1 + "01-" + YEAR2 + "12.nc"
        
        ;; Write variable
        write_nc(fout, VARS(vID), output, gAtt)
        delete(gAtt)
        
    end do

    leave_msg(diag_script, "", 1)
    
end
