"""
call like:
    python extract.py FOLDER
to check if imagefiles in FOLDER has metadata tag.
use python extract.py FOLDER -c 
to check for tags

"""
from PIL import Image
from PIL import PngImagePlugin
from PIL.ExifTags import TAGS, GPSTAGS
import argparse
import os
from os.path import join, getsize
import xml.etree.ElementTree as ET



parser = argparse.ArgumentParser(description='Extract Metadata')
parser.add_argument('indir', metavar='INDIR', type=str, nargs=1,
                                    help='indir containing files with metadata')
parser.add_argument('-c','--check', action='store_true',
                            help='activate checking for tags; default is deactive')


args = parser.parse_args()
indir= args.indir[0]
checkTags = args.check

tags ={ 
'namelist_': 'Namelist',
'P_': 'Projects',
'R_': 'CMIP6 Realms',
'T_': 'Themes',
'DM_': 'Domain',
'PT_': 'Plot Type',
'ST_': 'Statistics',
'D_': 'References',
'V_': 'Variables',
'A_': 'Author/Contributor List',
'M_': 'Models'}

if checkTags:
    print("{0}".format(",".join([v for k, v in tags.iteritems()])))
    print("{0}".format(",".join([item for item in tags])))
    print("{0}".format("|".join([item[0] for item in tags])))

list_info = []
for root, dirs, files in os.walk(indir):
    for f in files:
        image = Image.open("{0}/{1}".format(root,f))
        info = image.info
        for tag, value in info.items():
            if 'replace' not in dir(value):
                continue
            try:
                key = TAGS.get(tag, tag)
                if checkTags:
                    tmp = []
                    for k, v in tags.iteritems():
                        if "|{0}".format(k) in value.replace(' ', '') or ">{0}".format(k) in value.replace(' ', ''):
                            #print("Contains tag for {0}".format(v))
                            tmp.append(str(1))
                        else:
                            #print("Contains no tag for {0}".format(v))
                            tmp.append(str(0))
                    print("{1} for file: {0}".format(f, "|".join(tmp)))
                else:
                    print(f + ":" + key + " " + str(value).replace('\n', ' ').replace('\r', ' '))
            except:
                print(f + " NO KEY")

        list_info.append([f, info]) 


file = open("../../../doc/sphinx/source/print_report2.rst","w") 
file.write("DIAGNOSTIC\n")
file.write("==========\n\n")
file.write(".. centered:: |fig1|\n\n")
file.write(".. |fig1| image:: " + indir + "/" + list_info[0][0]+"\n   :width: 45%")
file.close() 
#../../../doc/sphinx/source/skel_report.rst
