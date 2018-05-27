EMBRACE diagnostics for monsoon, coupled equatorial climate and Southern Ocean
==============================================================================

Overview
--------

The namelist namelist_lauer18esd.xml combines different diagnostics from the namelists "clouds" (lauer13jclim; section :numref:`nml_clouds`), "South Asian monsoon" (SAMonsoon; section :numref:`nml_sam`), "Southern Hemisphere" (SouthernHemisphere; section :numref:`nml_sh`), "tropical variability" (TropicalVariability; section :numref:`nml_tropicalvariability`), and "West African monsoon" (WAMonsoon, WAMonsoon_daily; section :numref:`nml_wam`) to calculate the diagnostics used in Lauer et al. (2018). The diagnostics are applied to assess the performance of the four earth system models CNRM, EC-Earth, HadGEM, and MPI-ESM updated during the European Commission's 7th Framework Programme (FP7) project EMBRACE (Earth system Model Bias Reduction and assessing Abrupt Climate change) in comparison to their CMIP5 predecessor versions and against observations. Topics included are South Asian (SAM) and West African Monsoons (WAM), coupled equatorial climate, and Southern Ocean clouds and radiation.


Available namelists and diagnostics
-----------------------------------

Namelists are stored in nml/

    * namelist_lauer18esd.xml

Diagnostics are stored in diag_scripts/

    * clouds.ncl: global maps of (multi-year) annual and seasonal means + bias
    * SAMonsoon_precip_basic.ncl: maps of SAM rainfall, temperature
    * SAMonsoon_precip_seasonal.ncl: SAM seasonal cycle precipitation
    * SAMonsoon_wind_basic.ncl: maps of SAM wind and divergence
    * SAMonsoon_wind_seasonal.ncl: SAM seasonal cycle wind
    * SouthernHemisphere.py: zonal means
    * SouthernHemisphere_scatter.py: scatter plots cloud fraction vs. radiation
    * TropicalVariability.py: latitude cross sections
    * TropicalVariability_EQ.py: equatorial mean longitude cross sections
    * TropicalVariability_wind.py: equatorial mean longitude cross sections (divergence)
    * WAMonsoon_10W10E_1D_basic.ncl: WAM zonal means
    * WAMonsoon_contour_basic.ncl: maps of WAM rainfall, temperature
    * WAMonsoon_isv_filtered.ncl: maps of 3-10-day bandpass filtered precipitation variance
    * WAMonsoon_wind_basic.ncl: maps of WAM wind


User settings
-------------

User setting files (cfg files) are stored in nml/cfg_embrace/

* clouds.ncl: see section :numref:`nml_clouds`
* SAMonsoon_precip_basic.ncl, SAMonsoon_precip_seasonal.ncl, SAMonsoon_wind_basic.ncl, SAMonsoon_wind_seasonal.ncl: see section :numref:`nml_sam`
* SouthernHemisphere.py, SouthernHemisphere_scatter.py: see section :numref:`nml_sh`
* TropicalVariability.py, TropicalVariability_EQ.py, TropicalVariability_wind.py: see section :numref:`nml_tropicalvariability`
* WAMonsoon_10W10E_1D_basic.ncl, WAMonsoon_contour_basic.ncl, WAMonsoon_isv_filtered.ncl, WAMonsoon_wind_basic.ncl: see section :numref:`nml_wam`


Variables
---------

* clivi (atmos, monthly mean, longitude latitude time)
* clt (atmos, monthly mean, longitude latitude time)
* clwvi (atmos, monthly mean, longitude latitude time)
* pr (atmos, monthly mean / day, longitude latitude time)
* psl (atmos, monthly mean, longitude latitude time)
* rlds (atmos, monthly mean, longitude latitude time)
* rlut, rlutcs (atmos, monthly mean, longitude latitude time)
* rsds (atmos, monthly mean, longitude latitude time)
* rsut, rsutcs (atmos, monthly mean, longitude latitude time)
* tas (atmos, monthly mean, longitude latitude time)
* ts (atmos, monthly mean, longitude latitude time)
* ua (atmos, monthly mean, level longitude latitude time)
* va (atmos, monthly mean, level longitude latitude time)


Observations and reformat scripts
---------------------------------

*Note: (1) obs4mips data can be used directly without any preprocessing; (2) see headers of reformat scripts for non-obs4mips data for download instructions.*

* CERES-EBAF (rlds, rlut, rlutcs, rsds, rsut, rsutcs -- obs4mips)
* CMAP (pr -- reformat_scripts/obs/reformat_obs_CMAP.ncl)
* CRU (tas -- reformat_scripts/obs/reformat_obs_CRU.ncl)
* ERA-Interim (tas, ua, va, psl, clt, lwp, clivi -- reformat_scripts/obs/reformat_obs_ERA-Interim.ncl)
* GPCP-SG (pr -- obs4mips)
* GPCP-1DD (pr -- obs4mips)
* HadISST (ts -- reformat_scripts/obs/reformat_obs_HadISST.ncl)
* MODIS-L3-C6 (clt, clivi -- reformat_scripts/obs/reformat_obs_MODIS-L3-C6.ncl)
* NCEP (tas, ua, va -- reformat_scripts/obs/reformat_obs_NCEP.ncl)
* TRMM-L3 (pr -- obs4mips)
* UWisc (lwp -- reformat_scripts/obs/reformat_obs_UWisc.ncl)


References
----------

* Lauer, A., C. Jones, V. Eyring, M. Evaldsson, S. Hagemann, G. Martin, R. Roehrig, and Shiyu Wang, Process-level improvements in CMIP5 models and their impact on tropical variability, Southern Ocean and monsoons, Earth System Dynamics (accepted).


Example plots
-------------

.. figure:: /namelists/figures/embrace/figure01.png
   :width: 70 %
   
   Bias in annual mean near-surface air temperature (Lauer et al. (2018), figure 1).

.. figure:: /namelists/figures/embrace/figure03.png
   :width: 90 %
   
   Seasonal mean precipitation for JJAS (left) and differences relative to TRMM (right) (Lauer et al. (2018), figure 3).

.. figure:: /namelists/figures/embrace/figure05.png
   :width: 90 %
   
   Seasonal mean zonal wind speed at 850 hPa for JJAS (leftmost two columns) and differences relative to ERA-Interim (rightmost two columns) (Lauer et al. (2018), figure 5).

.. figure:: /namelists/figures/embrace/figure06.png
   :width: 50 %
   
   Mean annual cycle plots averaged over 5°-30°N, 65°-95°W (a), Webster and Yang Monsoon Index (b), and Goswami Monsoon Index (c) (Lauer et al. (2018), figure 6).

.. figure:: /namelists/figures/embrace/figure07.png
   :width: 90 %
   
   Seasonal mean precipitation for JJAS (leftmost two columns) and differences relative to TRMM (rightmost two columns) (Lauer et al. (2018), figure 7).

.. figure:: /namelists/figures/embrace/figure09.png
   :width: 90 %
   
   Seasonal mean wind speed at 925 hPa for JJAS (leftmost two columns) and differences relative to ERA-Interim (rightmost two columns) (Lauer et al. (2018), figure 9).

.. figure:: /namelists/figures/embrace/figure10.png
   :width: 90 %
   
   10°W-10°E zonal average JJAS mean values as a function of latitude (Lauer et al. (2018), figure 10).

.. figure:: /namelists/figures/embrace/figure11.png
   :width: 50 %
   
   JJAS average 3-10 day band-pass filtered precipitation variance calculated from daily precipitation fields (Lauer et al. (2018), figure 11).

.. figure:: /namelists/figures/embrace/figure14.png
   :width: 70 %
   
   Latitude cross-section of DJF SST (left) and precipitation (right) (Lauer et al. (2018), figure 14).

.. figure:: /namelists/figures/embrace/figure15.png
   :width: 70 %
   
   Equatorial mean (2.5°N-2.5°S) values plotted for the Pacific between 120°E and 80°W (Lauer et al. (2018), figure 15).

.. figure:: /namelists/figures/embrace/figure16.png
   :width: 70 %
   
   Latitude cross-section of DJF zonal means (Lauer et al. (2018), figure 16).

.. figure:: /namelists/figures/embrace/figure17.png
   :width: 90 %
   
   Scatterplot of monthly mean TOA SWUP versus total cloud cover for the Southern Ocean region 30°S-65°S and season DJF (top row), surface SWDN versus total cloud cover (bottom row), and fractional occurence of monthly mean cloud cover over this region (middle row) (Lauer et al. (2018), figure 17).

