#!/usr/bin/env bash
set -e

#for ip2f in ../../test_suites/easytest/namelist_collins*xml
for ip2f in ../../namelist_GlobalOcean*xml
do
    echo $ip2f
    ifile=$(echo ${ip2f}|cut -d'/' -f3)
    cp ${ip2f} ./${ifile}
    python /pf/b/$USER/util/translate2ESGF/translate2ESGF.py -n ${ifile}
    touch test.xml
    mv test.xml ${ifile}
done
