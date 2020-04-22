;#######################################################################
;# Model-specific fixes for reformatDEFAULT_main                       #
;# Author: Ruth Lorenz (ETH-Switzerland)                               #
;#######################################################################
;# PROJECT:     CMIP5                                                  #
;# ENSEMBLE:    [all]                                                  #
;# MODEL:       CESM1-WACCM                                            #
;# VARIABLES:   sos                                                    #
;# FIELDS:      [all]                                                  #
;# DESCRIPTION: inconsistent units                                     #
;#######################################################################

undef("fixing_func")
function fixing_func(var:numeric,name:string)
local f, files, tmp
begin

    ret = 1

    if (name.eq."sos") then
        var = where(var .lt. 1., var * 1000., var)
        var@units = "psu"
        ret = 0
    end if

    return(ret)

end
