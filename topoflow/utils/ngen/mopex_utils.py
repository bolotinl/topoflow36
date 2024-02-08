
# Copyright (c) 2023-2024, Scott D. Peckham
#
# Jan 2024. Renamed all "ngen/utils" files to end in "_utils.py"
#           instead of "_tools.py"
#           Modified to use new data_utils.py.
# Jun 2023. Started from gages2_utils.py.
#
#---------------------------------------------------------------------
#
#  % conda activate tf36  (has gdal package)
#  % python
#  >>> from topoflow.utils.ngen import mopex_utils as mt
#  >>> mt.create_tsv( nf_max=500 )
#
#---------------------------------------------------------------------
#
#  create_tsv()
#  sort_by_site_code()
#
#---------------------------------------------------------------------

import numpy as np
from osgeo import ogr, osr
import json, sys, time

from topoflow.utils.ngen import data_utils as dtu
from topoflow.utils.ngen import shape_utils as su

#---------------------------------------------------------------------
def create_tsv( data_dir=None, new_dir=None, nf_max=50,
                shape_file='MOPEX431_basins.shp',                    
                prj_file  ='MOPEX431_basins.prj',
                tsv_file='new_mopex431.tsv',
                SWAP_XY=False,  # False works for MOPEX.
                REPORT=True):

    if (data_dir is None):
        data_dir = dtu.get_data_dir( 'MOPEX' )
    if (new_dir is None):
        new_dir  = dtu.get_new_data_dir( 'MOPEX' )
    insert_key = 'Longitude'

    su.create_tsv_from_shapefile(data_dir=data_dir, new_dir=new_dir,
                      shape_file=shape_file, prj_file=prj_file,
                      tsv_file=tsv_file, nf_max=nf_max,
                      SWAP_XY=SWAP_XY, REPORT=REPORT,
                      ADD_BOUNDING_BOX=True, insert_key=insert_key,
                      filter_key=None, filter_value=None)

#   create_tsv()
#---------------------------------------------------------------------
def sort_by_site_code(tsv_file1='new_mopex431.tsv',
                      tsv_file2='new_mopex431_sorted.tsv'):

    new_dir   = dtu.get_new_data_dir( 'MOPEX' )
    tsv_path1 = new_dir + tsv_file1
    tsv_path2 = new_dir + tsv_file2
    
    tsv_unit1 = open(tsv_path1, 'r')
    lines = tsv_unit1.readlines()
    tsv_unit1.close()

    #-------------------------------
    # Swap "SiteID" and "SiteCode"
    #-------------------------------
    tab = '\t'
    lines2 = list()
    for line in lines:
        pos_tab1  = line.find( tab )
        site_id   = line[:pos_tab1]
        line      = line[pos_tab1 + 1:]
        #---------------------------------
        pos_tab2  = line.find(tab)
        site_code = line[:pos_tab2]
        line      = line[pos_tab2 + 1:]
        #---------------------------------
        line2 = site_code + tab + site_id + tab + line
        lines2.append(line2)
        
    #---------------------------------
    # Sort all lines (by "SiteCode")
    #---------------------------------
    header = lines2[0]   # header line
    lines2 = lines2[1:]    
    lines2.sort()
    tsv_unit2 = open(tsv_path2, 'w')
    tsv_unit2.writelines( [header] + lines2 )
    tsv_unit2.close()

#   sort_by_site_code()
#---------------------------------------------------------------------

