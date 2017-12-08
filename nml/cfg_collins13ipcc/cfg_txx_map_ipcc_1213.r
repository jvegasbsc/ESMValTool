diag_script_info<-new()

## Used by: map_diff_mmm_stippWilcox_ipcc12.ncl
## map projection, e.g., Mollweide, Mercator
diag_script_info[["projection"]]<-"Robinson"
diag_script_infos[["styleset"]]<-"CMIP5"        #"CMIP5", "DEFAULT"
diag_script_info[["scenarios"]]<-c("rcp85") #  scenarios
diag_script_info[["periods"]]<-2081  # start year in time periods
## max allowed number of rows on a panel page (vertical), number of scenarios
diag_script_info[["max_vert"]]<-1
## max allowed number of columns on a panel page (horizontal), number of periods
diag_script_info[["max_hori"]]<-1
diag_script_info[["title"]]<-"Annual maximum surface air temperature change"
diag_script_info[["label"]]<-c("RCP8.5: 2081-2100")
diag_script_info[["units"]]<-c("(~F35~J~F~C)")
diag_script_info[["plotmask"]]<-c("ocean")
diag_script_info[["sig"]]<-False
diag_script_info[["figure_nr"]]<-"12.13"

## enable to output to netCDF; either use "default" or give a full file name
diag_script_info[["ncdf"]]<-"txx_mmm_sig_ipcc_1213.nc"
