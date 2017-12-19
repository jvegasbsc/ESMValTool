;#######################################################################
;# Model-specific fixes for reformatDEFAULT_main                       #
;# Author: Ruth Lorenz (ETH-Switzerland)                               #
;#######################################################################
;# PROJECT:     CMIP5                                                  #
;# ENSEMBLE:    [all]                                                  #
;# MODEL:       FGOALS-s2                                              #
;# VARIABLES:   evspsbl                                                #
;# FIELDS:      all                                                    #
;# DESCRIPTION: flip sign                                              #
;#######################################################################

undef("fixing_func")
function fixing_func(var:numeric,name:string)

begin

    ret = 1
    if (name.eq."evspsbl") then
        var = -1. * var
        ret = 0
    end if

    return(ret)

end

