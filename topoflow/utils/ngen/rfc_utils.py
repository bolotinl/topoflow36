
# Copyright (c) 2023-2024, Scott D. Peckham
#
# Jan 2024. Wrote create_combo_tsv().
#           Wrote get_rfc_dict1(), get_rfc_dict2().
#           Wrote create_rfc_tsv(); renamed create_tsv() to
#              create_hads_tsv().
#           Wrote get_usgs_rfc_dict_v1().
#           Renamed all "ngen/utils" files to end in "_utils.py"
#           instead of "_tools.py"
#           Modified to use new data_utils.py.
# Dec 2023. Wrote get_site_data_from_api() & create_tsv_via_api().
#           Wrote get_usgs_rfc_dict_v0().
# Nov 2023. Incorporated info in new RFC data sets.  Wrote
#           merge_rfc_basin_info().  Wrote convert_json_to_csv()
#           to work around issues with WGRFC shapefile.
#           Can now include RFC_name and basin bounding box
#           for 7759 basins in 13 NOAA RFCs.
#           Wrote:  cull_usa_dcp_data_rows().
#           Wrote:  get_site_data_from_api().
#           Still need to update create_tsv().
# Oct 2023. create_tsv() and all supporting functions working. 
# Sep 2023. Started.
#
#---------------------------------------------------------------------
#
#  % conda activate tf36
#  % python
#  >>> from topoflow.utils.ngen import rfc_utils as rfc
#  >>> rfc.create_rfc_tsv( nf_max=10000 )
#  >>> rfc.create_hads_tsv()
#  >>> rfc.sort_by_site_code()
#
#---------------------------------------------------------------------
#
#  get_nws_site_url()
#
#  get_rfc_dict1()   # Info from All_USGS-HADS_SITES.tsv
#  get_rfc_dict2()   # Info from new_rfc_info.tsv
#  get_equivalent_nws_ids2()
#  get_equivalent_nws_ids()
#  create_combo_tsv()
#
#  get_usgs_rfc_dict_v0()
#  get_usgs_rfc_dict_v1()
#  get_hydrograph_type_dict()
#  cull_usa_dcp_data_rows()
#  convert_json_to_csv()
#  merge_rfc_basin_info()
#
#  create_rfc_tsv()
#  create_hads_tsv()
#  sort_by_site_code()
#
#  get_site_data_from_api()
#  create_tsv_via_api()
#
#---------------------------------------------------------------------

import numpy as np
import os, os.path
import json, requests, time

import pickle    # to save usgs site info map
from topoflow.utils.ngen import data_utils as dtu
from topoflow.utils.ngen import shape_utils as su
from topoflow.utils.ngen import usgs_utils as usgs

#---------------------------------------------------------------------
def get_nws_site_url( wfo='SGF', nws_id='RBUM7' ):

    site_url = 'https://water.weather.gov/ahps2/hydrograph.php'
    site_url += '?wfo='  + wfo.lower()
    site_url += '&gage=' + nws_id.lower()
    return site_url

#   get_nws_site_url()
#---------------------------------------------------------------------
def get_rfc_dict1( SAVE_TO_FILE=True ):

    #-------------------------------------------------------------    
    # Construct a site_info dictionary, where NWS loc ID
    # is the key used to get dictionary info for a site.
    # This version uses the file: ALL_USGS-HADS_SITES.txt,
    # which was obtained from:
    #   https://hads.ncep.noaa.gov/USGS/ALL_USGS-HADS_SITES.txt
    # and has info for 10510 sites/basins.
    # Note that this file is sorted by 5-char NWS location ID.
    #-------------------------------------------------------------
    # Note: None of the NWS IDs in the set of all USGS-HADS
    # sites is numeric.
    #-------------------------------------------------------------
    new_dir   = dtu.get_new_data_dir( 'NOAA_RFC' )
    data_dir  = dtu.get_data_dir( 'NOAA_HADS' )
    #------------------------------------------------
    dict_file = 'RFC_site_info_dict1.pkl'
    file_path = new_dir + dict_file
    if (os.path.exists( file_path )):
        print('Reading saved RFC site info dictionary...')
        file_unit = open( file_path, 'rb')
        site_info = pickle.load( file_unit )
        file_unit.close()
        print('Finished.')
        print()
        return site_info

    print('Getting RFC site info as dictionary...')
    info_file = 'All_USGS-HADS_SITES.tsv'
    info_path = data_dir + info_file
    info_unit = open( info_path, 'r' )
    delim = '\t'
    
    #--------------------------------
    # Skip over three header lines
    #--------------------------------
    for j in range(3):
        line = info_unit.readline()
    
    site_info = dict()
    k = 0

    while (True):
        info_line = info_unit.readline()
        if (info_line == ''):
            break  # (reached end of file)
        vals = info_line.split( delim )
        #-------------------------------------
        # Note that RFC abbrev. is not given
        #-------------------------------------
        nws_id    = vals[0].strip()
        usgs_id   = vals[1].strip()
        goes_id   = vals[2].strip()
        hsa_id    = vals[3].strip()
        lat_dms   = vals[4].strip()
        lon_dms   = vals[5].strip()
        usgs_name = vals[6].strip()
        #------------------------------------------
        # Note: Only one of the longitudes starts
        # with a minus sign, so need to add this!
        #------------------------------------------
        if (lon_dms[0] != '-'):
            lon_dms  = '-' + lon_dms
        #--------------------------------------------
        lon = dtu.convert_dms_to_dec_deg( lon_dms )
        lat = dtu.convert_dms_to_dec_deg( lat_dms )

        lon_str = '{x:.5f}'.format(x=lon)
        lat_str = '{x:.5f}'.format(x=lat)
                     
        ## url = usgs.get_usgs_site_url( usgs_id )
        site_info[ nws_id ] = \
            {'usgs_id':usgs_id, 'usgs_name':usgs_name,
             'goes_id':goes_id, 'hsa_id':hsa_id,
             'lat':lat_str, 'lon':lon_str } 
        k += 1

    print('Processed', k, 'RFC sites.')

    if (SAVE_TO_FILE):
        file_unit = open( file_path, 'wb' )
        pickle.dump( site_info, file_unit, protocol=pickle.HIGHEST_PROTOCOL)
        file_unit.close()
        print('Saved RFC site info dictionary to file:')
        print('  ' + file_path)
        print()

    return site_info

#   get_rfc_dict1()
#---------------------------------------------------------------------
def get_rfc_dict2( SAVE_TO_FILE=True ):

    #-------------------------------------------------------------    
    # Construct a site_info dictionary, where NWS loc ID
    # is the key used to get dictionary info for a site.
    # This version uses the file: new_rfc_info.tsv
    # which was extracted from the shapefile:
    #   NOAA_RFC/Data/ba12my15.shp
    # and has info for 9370 sites/basins.
    # The attributes differs from what is in rfc_dict1.
    #-------------------------------------------------------------
    new_dir = dtu.get_new_data_dir( 'NOAA_RFC' )
    #------------------------------------------------
    dict_file = 'RFC_site_info_dict2.pkl'
    file_path = new_dir + dict_file
    if (os.path.exists( file_path )):
        print('Reading saved RFC site info dictionary...')
        file_unit = open( file_path, 'rb')
        site_info = pickle.load( file_unit )
        file_unit.close()
        print('Finished.')
        print()
        return site_info

    print('Getting RFC site info as dictionary...')
    info_file = 'new_rfc_info.tsv'
    info_path = new_dir + info_file
    info_unit = open( info_path, 'r' )
    delim = '\t'
    
    #----------------------------
    # Skip over one header line
    #----------------------------
    line = info_unit.readline()
    
    site_info = dict()
    k = 0

    while (True):
        info_line = info_unit.readline()
        if (info_line == ''):
            break  # (reached end of file)
        vals = info_line.split( delim )
        #-------------------------------------
        # Note that RFC abbrev. is not given
        #-------------------------------------
        object_id = vals[0].strip()
        nws_id    = vals[1].strip()
        name      = vals[2].strip()
        cwa_id    = vals[3].strip()
        rfc_id    = vals[4].strip()
        lon_str   = vals[5].strip()
        lat_str   = vals[6].strip()
        minlon    = vals[7].strip()
        maxlon    = vals[8].strip()
        minlat    = vals[9].strip()
        maxlat    = vals[10].strip()

        lon = np.float64( lon_str )
        lat = np.float64( lat_str )
        lon_str = '{x:.5f}'.format(x=lon)
        lat_str = '{x:.5f}'.format(x=lat)
        
        site_info[ nws_id ] = \
            {'rfc_id':rfc_id, 'cwa_id':cwa_id,
             'lon':lon_str, 'lat':lat_str, 'minlon':minlon,
             'maxlon':maxlon, 'minlat':minlat, 'maxlat':maxlat }
        k += 1

    print('Processed', k, 'RFC sites.')

    if (SAVE_TO_FILE):
        file_unit = open( file_path, 'wb' )
        pickle.dump( site_info, file_unit, protocol=pickle.HIGHEST_PROTOCOL)
        file_unit.close()
        print('Saved RFC site info dictionary to file:')
        print('  ' + file_path)
        print()

    return site_info

#   get_rfc_dict2()
#---------------------------------------------------------------------
def get_rfc_dict3( SAVE_TO_FILE=True ):

    #-------------------------------------------------------------    
    # Construct a site_info dictionary, where NWS loc ID
    # is the key used to get dictionary info for a site.
    # This version uses the file:
    #   USGS_NWIS_Web_Site_Info_via_API.tsv
    # which was obtained from a new NOAA API.
    # and has info for 7136 sites/basins that have USGS
    # site IDs in the 'USGS_NWIS_Web' dataset.
    # The attributes differ from what is in rfc_dict1
    # and rfc_dict2.
    #-------------------------------------------------------------
    new_dir1 = dtu.get_new_data_dir( 'NOAA_RFC' )
    new_dir2 = dtu.get_new_data_dir( 'NOAA_via_API' )
    #------------------------------------------------
    dict_file = 'RFC_site_info_dict3.pkl'
    file_path = new_dir1 + dict_file
    if (os.path.exists( file_path )):
        print('Reading saved RFC site info dictionary...')
        file_unit = open( file_path, 'rb')
        site_info = pickle.load( file_unit )
        file_unit.close()
        print('Finished.')
        print()
        return site_info

    print('Getting RFC site info as dictionary...')
    info_file = 'USGS_NWIS_Web_Site_Info_via_API.tsv'
    info_path = new_dir2 + info_file
    info_unit = open( info_path, 'r' )
    delim = '\t'
    
    #----------------------------
    # Skip over one header line
    #----------------------------
    line = info_unit.readline()
    
    site_info = dict()
    k = 0

    while (True):
        info_line = info_unit.readline()
        if (info_line == ''):
            break  # (reached end of file)
        vals = info_line.split( delim )
        #-------------------------------------
        usgs_id    = vals[0].strip()
        nws_id     = vals[1].strip()
        nws_up_id  = vals[2].strip()
        nws_dn_id  = vals[3].strip()
        reach_id   = vals[4].strip()  # NWM reach ID
        rfc_id     = vals[5].strip()
        wfo_id     = vals[6].strip()
        usgs_name  = vals[7].strip()
        descript   = vals[8].strip()
        state_code = vals[9].strip()
        county     = vals[10].strip()
        lon_str    = vals[11].strip()
        lat_str    = vals[12].strip()
        elev       = vals[13].strip()
        upd_time   = vals[14].strip()
        in_service = vals[15].strip()
        pedts_obs  = vals[16].strip()
        pedts_pred = vals[17].strip()

        lon = np.float64( lon_str )
        lat = np.float64( lat_str )
        lon_str = '{x:.5f}'.format(x=lon)
        lat_str = '{x:.5f}'.format(x=lat)
        
        site_info[ nws_id ] = {
            'usgs_id':usgs_id, 'nws_up_id':nws_up_id, 'nws_dn_id':nws_dn_id,
            'reach_id':reach_id, 'rfc_id':rfc_id, 
            ### 'wfo_id',  # always blank
            'usgs_name':usgs_name, 'description':descript,
            'state_code':state_code, 'county':county,
            'lon':lon_str, 'lat':lat_str, 'elev':elev,
            'update_time':upd_time, 'in_service':in_service,
            'pedts_obs':pedts_obs, 'pedts_pred':pedts_pred }
        k += 1

    print('Processed', k, 'RFC sites.')

    if (SAVE_TO_FILE):
        file_unit = open( file_path, 'wb' )
        pickle.dump( site_info, file_unit, protocol=pickle.HIGHEST_PROTOCOL)
        file_unit.close()
        print('Saved RFC site info dictionary to file:')
        print('  ' + file_path)
        print()

    return site_info

#   get_rfc_dict3()
#---------------------------------------------------------------------
def get_equivalent_nws_ids2():

    #-------------------------------------------------------
    # The MBRFC uses two different NWS loc IDs for most
    # or all basins.  One id is a 3 or 4-digit number,
    # while the other is a standard 5+ character NWS
    # location ID.  Both IDs occur in the basin shapefile
    # ba12my15.shp.  This function uses lon and lat to
    # identify which IDs refer to the same MBRFC basin.
    #-------------------------------------------------------
    # The numeric NWS IDs that don't match to a standard
    # 5-char NWS ID seem to be parts of a larger, narrow
    # basin in many cases.
    #-------------------------------------------------------

    #------------------------
    # Get Python dictionary
    #------------------------
    rfc_dict2 = get_rfc_dict2()  # info from RFC shapefile

    #-------------------------------------------------
    # Both of these ID pairs are in rfc_dict2 and 
    # they have the same lon and lat and other info.
    #-------------------------------------------------
    # "MISSOURI RIVER AT ATCHISON KS" (MBRFC, EAX)
    # Both of these have the same name, but 1814 seems
    # to plot beneath ATCK1.
    #----------------------------------------------------
    print("rfc_dict2[ 'ATCK1' ] =", rfc_dict2['ATCK1'])
    print()
    print("rfc_dict2[ '1814' ]  =", rfc_dict2['1814'])
    print()
    #-----------------------------------------------------
    # "ENNIS LAKE MT"  (MBRFC, ???)
    # Both of these have the same name, but 2843 seems
    # to plot beneath ATCK1. Checked w/ QGIS "info" tool
    # that shapefile has them as 2 overlapping polygons.
    #-----------------------------------------------------
    print("rfc_dict2[ 'ELMM8' ] =", rfc_dict2['ELMM8'])
    print()
    print("rfc_dict2[ '2843' ]  =", rfc_dict2['2843'])
    print()
    #----------------------------------------------------------
    # "EAST NISHNABOTNA RIVER NEAR ATLANTIC IA"  (MBRFC, DMX)
    # Both of these have the same name, but 1604 seems
    # to plot beneath ATCI4. Checked w/ QGIS "info" tool
    # that shapefile has them as 2 overlapping polygons.
    #----------------------------------------------------------
    # NOTE! Adjacent basin with ID=1630 also appears twice
    # but with that same ID number, according to info tool.
    #----------------------------------------------------------    
    print("rfc_dict2[ 'ATCI4' ] =", rfc_dict2['ATCI4'])
    print()
    print("rfc_dict2[ '1604' ]  =", rfc_dict2['1604'])
    print()
        
    #---------------------------------------------    
    # Get list of unique NWS IDs from both lists
    #---------------------------------------------
    nws_ids2 = list( rfc_dict2.keys() )

    #---------------------------------------------    
    # Get lists for all longitudes and latitudes
    #---------------------------------------------
    lon_list = list()
    lat_list = list()
    n_numeric_ids = 0
    for nws_id in nws_ids2:
        rec = rfc_dict2[ nws_id ]
        lon_list.append( rec['lon'] )
        lat_list.append( rec['lat'] )

    #---------------------------------------
    # Convert lists to numpy string arrays
    #---------------------------------------
    lon_list2 = np.array( lon_list )
    lat_list2 = np.array( lat_list )
    id_list2  = np.array( nws_ids2 )
    #---------------------------------------
    s = np.argsort( lon_list2 )
    lon_list2 = lon_list2[ s ]
    lat_list2 = lat_list2[ s ]
    id_list2  = id_list2[ s ]
   
    last_lon = '-999'
    last_lat = '-999'
    last_id  = 'XXXX'
    list1    = list()
    list2    = list()
    for k in range( len(lon_list2) ):
        lon = lon_list2[k]
        lat = lat_list2[k]
        id  = id_list2[k]
        if (id.isnumeric()):
            n_numeric_ids += 1
        if (lon == last_lon) and (lat == last_lat):
            print('Found equivalent IDs:')
            print('   ID1 =', last_id, ', ID2 =', id)
            if (last_id.isnumeric()):
                list1.append( last_id )
                list2.append( id)
            else:
                list1.append( id )
                list2.append( last_id )
        last_lon = lon
        last_lat = lat
        last_id  = id
 
    a  = np.array( list1 )
    b  = np.array( list2 )
    s1 = np.argsort( a )
    a  = a[ s1 ]
    b  = b[ s1 ]
    a_list = list(a)
    b_list = list(b)
    print('Total number of numeric IDs =', n_numeric_ids)
    print('Number of equivalent IDs =', len(a_list))
    print()

    return a_list, b_list

#   get_equivalent_nws_ids2()
#---------------------------------------------------------------------
def get_equivalent_nws_ids():

    #-------------------------------------------------------
    # The MBRFC uses two different NWS loc IDs for most
    # or all basins.  This function uses lon and lat to
    # check whether the NWS IDs from the basin shapefile
    # ba12my15.shp are equivalent to any of the ones from
    # the USGS HADS crosswalk.
    #-------------------------------------------------------

    #-------------------------------
    # Get some Python dictionaries
    #-------------------------------
    rfc_dict1 = get_rfc_dict1()  # info from HADS crosswalk
    rfc_dict2 = get_rfc_dict2()  # info from RFC shapefile
   
    #---------------------------------------------    
    # Get list of unique NWS IDs from both lists
    #---------------------------------------------
    nws_ids1 = list( rfc_dict1.keys() )
    nws_ids2 = list( rfc_dict2.keys() )

    #-----------------------------------------    
    # Build a dictionary with lon as the key
    #-----------------------------------------
    lon_dict1 = dict()
    lon_list1 = list()
    for nws_id in nws_ids1:
        rec = rfc_dict1[ nws_id ]
        lon = rec[ 'lon' ]  # string
        lat = rec[ 'lat' ]
        lon_dict1[ lon ] = {'lat':lat, 'nws_id':nws_id }
        lon_list1.append( lon )
        
    lon_dict2 = dict()
    lon_list2 = list()
    for nws_id in nws_ids2:
        rec = rfc_dict2[ nws_id ]
        lon = rec[ 'lon' ]  # string
        lat = rec[ 'lat' ]
        lon_dict2[ lon ] = {'lat':lat, 'nws_id':nws_id }
        lon_list2.append( lon )

    #------------------------------------------         
    # Create a sorted list of all unique lons
    #------------------------------------------
    lon_list = list( set( lon_list1 + lon_list2 ) )
    ## print('lon_list =', lon_list) 
    lon_list.sort()
     
    id_list1 = list()
    id_list2 = list()
    for lon in lon_list:
        if (lon in lon_dict1) and (lon in lon_dict2):
            rec1 = lon_dict1[ lon ]
            rec2 = lon_dict2[ lon ]
            lat1 = rec1[ 'lat' ]
            lat2 = rec2[ 'lat' ]
            if (lat1 == lat2):
                id1 = rec1[ 'nws_id' ]
                id2 = rec2[ 'nws_id' ]
                id_list1.append( id1 )
                id_list2.append( id2 )
  
    n_equivalent = len( id_list1 )
    print('n_equivalent =', n_equivalent)
    
    return id_list1, id_list2
    
#   get_equivalent_nws_ids()
#---------------------------------------------------------------------
def create_combo_tsv( tsv_file='new_combo_rfc_info.tsv',
                      REPORT=True):

    #-----------------------------
    # Open new TSV file to write
    #-----------------------------
    new_dir  = dtu.get_new_data_dir( 'NOAA_RFC' )
    tsv_path = new_dir + tsv_file
    tsv_unit = open( tsv_path, 'w')

    #--------------------------------
    # Write header for new TSV file
    #--------------------------------
    delim  = '\t'  # tab character  
    header = ''
    header += 'USGS_ID'      + delim
    header += 'USGS_name'    + delim
    header += 'Site_type'    + delim
    header += 'Longitude'    + delim
    header += 'Latitude'     + delim
    header += 'Alt_Lon'      + delim
    header += 'Alt_Lat'      + delim
    header += 'Elevation'    + delim  # What are units?
    header += 'Minlon'       + delim
    header += 'Maxlon'       + delim
    header += 'Minlat'       + delim
    header += 'Maxlat'       + delim
    header += 'RFC_ID'       + delim
    header += 'NWS_ID'       + delim
    header += 'NWS_alt_ID'   + delim  ###
    header += 'HSA_ID'       + delim
    header += 'CWA_ID'       + delim
    header += 'GOES_ID'      + delim
    header += 'HUC8'         + delim
    header += 'State_Code'   + delim
    header += 'NWM Reach ID' + delim
    header += 'Active'       + delim
    header += 'PEDTS obs'    + delim
    header += 'PEDTS_pred'   + delim
    header += 'Hgraph_type'  + '\n'  # newline at end  
    tsv_unit.write( header )   
        
    #-------------------------------
    # Get some Python dictionaries
    #-------------------------------
    rfc_dict1 = get_rfc_dict1()
    rfc_dict2 = get_rfc_dict2()
    rfc_dict3 = get_rfc_dict3()

    #--------------------------------------------    
    # Get list of unique NWS IDs from all lists
    #--------------------------------------------
    nws_ids1 = list( rfc_dict1.keys() )
    nws_ids2 = list( rfc_dict2.keys() )
    nws_ids3 = list( rfc_dict3.keys() )
    ## nws_ids  = list( set(nws_ids1 + nws_ids2) )
    nws_ids  = list( set(nws_ids1 + nws_ids2 + nws_ids3) )
    nws_ids.sort()
    print('Number of unique NWS IDs =', len(nws_ids))

    #-------------------------------    
    # Get more Python dictionaries
    #-------------------------------
    rfc_map   = get_usgs_rfc_dict_v0()
    htype_map = get_hydrograph_type_dict()
    site_type_dict = usgs.get_site_type_dict()   ####

    alt_id_list1, alt_id_list2 = get_equivalent_nws_ids2()
    
    k = 0
    n_lons_differ = 0
    n_lats_differ = 0
    n_good_lons_differ = 0
    n_good_lats_differ = 0
    for nws_id in nws_ids:
        usgs_id   = '-'
        usgs_name = '-'
        goes_id   = '-'
        hsa_id    = '-'
        lon1      = '-'
        lat1      = '-'
        #-----------------
        rfc_id    = '-'
        cwa_id    = '-'
        lon2      = '-'
        lat2      = '-'
        minlon    = '-'
        maxlon    = '-'
        minlat    = '-'
        maxlat    = '-'
        #------------------
        rfc_id2  = '-'
        huc8     = '-'
        lon3     = '-'
        lat3     = '-'
        #------------------
        alt_nws_id = '-' 
        alt_lon_id = '-'
        alt_lat_id = '-'
        #------------------
        reach_id   = '-'
        state_code = '-'
        elev       = '-'
        active     = '-'
        pedts_obs  = '-'
        pedts_pred = '-'
        #---------------------------------------      
        if (nws_id in rfc_dict1):
            rec1 = rfc_dict1[ nws_id ]
            usgs_id   = rec1['usgs_id']
            usgs_name = rec1['usgs_name']
            goes_id   = rec1['goes_id']
            hsa_id    = rec1['hsa_id']
            lon1      = rec1['lon']
            lat1      = rec1['lat']
        #---------------------------------------                
        if (nws_id in rfc_dict2):
            rec2 = rfc_dict2[ nws_id ]
            rfc_id = rec2['rfc_id']
            cwa_id = rec2['cwa_id']
            lon2   = rec2['lon']
            lat2   = rec2['lat']
            minlon = rec2['minlon'] 
            maxlon = rec2['maxlon']           
            minlat = rec2['minlat']             
            maxlat = rec2['maxlat']
        #---------------------------------------                
        if (nws_id in rfc_dict3):
            rec3 = rfc_dict3[ nws_id ]
            if (usgs_id == '-'):
                usgs_id = rec3['usgs_id']
            reach_id   = rec3['reach_id']
            if (rfc_id == '-'):
               rfc_id = rec3['rfc_id']
            # usgs_name = rec3['usgs_name']
            state_code = rec3['state_code']
            # county     = rec3['county']
            lon3       = rec3['lon']
            lat3       = rec3['lat']
            elev       = rec3['elev']
            # updated_time = rec3['updated_time']
            active     = rec3['in_service'] 
            pedts_obs  = rec3['pedts_obs']
            pedts_pred = rec3['pedts_pred']          
        #----------------------------------------
        if (usgs_id in rfc_map):
            rec4 = rfc_map[ usgs_id ]
            rfc_id2  = rec4['rfc_name']
            huc8     = rec4['huc8']
            lon4     = rec4['lon']
            lat4     = rec4['lat']

        #----------------------------------
        # Change 'None' and '-9999', etc.
        #----------------------------------
        if (cwa_id == 'None'):
            cwa_id = '-'
        if (elev == '-9999.0') or (elev == '-9999.00000'):
            elev = '-'

        #-------------------
        # Process the lons
        #-------------------
        if (lon1 != '-'):
            lon_str     = lon1
            alt_lon_str = lon2   # valid or '-'
        else:
            # lon1 == '-'
            if (lon2 != '-'):
                lon_str     = lon2
                alt_lon_str = lon1  # '-'
        #---------------------------------------------------------      
        if (lon1 != lon2) and (lon1 != '-') and (lon2 != '-'):
            print('WARNING: Different lons for NWS ID =', nws_id)
            print('   lon1 =', lon1, ', lon2 =', lon2)
            n_good_lons_differ += 1
            
        #-------------------
        # Process the lats
        #-------------------
        if (lat1 != '-'):
            lat_str     = lat1
            alt_lat_str = lat2   # valid or '-'
        else:
            # lat1 == '-'
            if (lat2 != '-'):
                lat_str     = lat2
                alt_lat_str = lat1  # '-'
        #---------------------------------------------------------      
        if (lat1 != lat2) and (lat1 != '-') and (lat2 != '-'):
            print('WARNING: Different lats for NWS ID =', nws_id)
            print('   lat1 =', lat1, ', lat2 = ', lat2)
            n_good_lats_differ += 1
               
        #-------------------------------------
        # Get hydrograph type classification
        # (not available for all sites)
        #-------------------------------------
        if (usgs_id in htype_map):
            rec4 = htype_map[ usgs_id ]
            hgraph_type = rec4[ 'htype' ]
        else:
            hgraph_type = 'unknown'

        #---------------------------          
        # Get the site type string
        #---------------------------
        if (usgs_id in site_type_dict):
            site_type = site_type_dict[ usgs_id ] 
        else:
            site_type = 'unknown'
          
        #-------------------------------
        # Write line to new TSV file ?
        #-------------------------------
        IN_LIST1  = (nws_id in alt_id_list1)
        IN_LIST2  = (nws_id in alt_id_list2)
        NO_ALT_ID = not(IN_LIST1) and not(IN_LIST2)
        if (IN_LIST2):
            index2 = alt_id_list2.index( nws_id )

        if (NO_ALT_ID or IN_LIST2):
            if (IN_LIST2):
                alt_nws_id = alt_id_list1[ index2 ]
            else:
                alt_nws_id = '-'
            line = ''
            line += usgs_id     + delim  # should be first
            line += usgs_name   + delim
            line += site_type   + delim
            line += lon_str     + delim
            line += lat_str     + delim
            line += alt_lon_str + delim   ###
            line += alt_lat_str + delim   ###
            line += elev        + delim   # from rfc_dict3
            line += minlon      + delim
            line += maxlon      + delim
            line += minlat      + delim
            line += maxlat      + delim
            line += rfc_id      + delim
            line += nws_id      + delim
            line += alt_nws_id  + delim   ###
            line += hsa_id      + delim
            line += cwa_id      + delim
            line += goes_id     + delim
            line += huc8        + delim
            #-----------------------------
            line += state_code  + delim
            line += reach_id    + delim
            line += active      + delim
            line += pedts_obs   + delim
            line += pedts_pred  + delim
            line += hgraph_type + '\n'   # add newline at end
            tsv_unit.write( line )       
            k += 1

    #------------------ 
    # Close the files
    #------------------
    print('n_good_lons_differ =', n_good_lons_differ)
    print('n_good_lats_differ =', n_good_lats_differ)
    tsv_unit.close()
    print('Finished creating combo TSV file for RFC basins.')
#     print('  lon_diff_count =', lon_diff_count)
#     print('  lat_diff_count =', lat_diff_count)
    print()

#   create_combo_tsv()
#---------------------------------------------------------------------
def get_usgs_rfc_dict_v0( SAVE_TO_FILE=True ):

    #------------------------------------------------------------    
    # Construct a site_info dictionary, where USGS site ID
    # is the key used to get dictionary info for a site.
    # This version uses the file: usgsMasterMeta.txt, which was
    # obtained from Wanru Wu and has info for 8170 basins.
    # Don't know the original source of this file.
    # Note that this file is sorted by USGS site ID.
    #------------------------------------------------------------
    # See get_usgs_rfc_dict_v1(), which may be better.
    #------------------------------------------------------------    
    new_dir   = dtu.get_new_data_dir( 'NOAA_HADS' )
    data_dir  = dtu.get_data_dir( 'NOAA_HADS' )
    #------------------------------------------------
    dict_file = 'USGS_RFC_site_info_dict0.pkl'
    file_path = new_dir + dict_file
    if (os.path.exists( file_path )):
        print('Reading saved RFC site info dictionary...')
        file_unit = open( file_path, 'rb')
        site_info = pickle.load( file_unit )
        file_unit.close()
        print('Finished.')
        print()
        return site_info

    print('Getting USGS/RFC site info as dictionary...')
    info_file = 'usgsMasterMeta.txt'
    info_path = new_dir + info_file
    info_unit = open( info_path, 'r' )
    delim     = ','

    #--------------------------------
    # Skip over the one header line
    #--------------------------------
    line = info_unit.readline()
    
    site_info = dict()
    k = 0

    while (True):
        info_line = info_unit.readline()
        if (info_line == ''):
            break  # (reached end of file)
        vals = info_line.split( delim )
        sid  = vals[0].strip()
        lat  = vals[1].strip()
        lon  = vals[2].strip()
        huc8 = vals[6].strip()
        rfc  = vals[7].strip()
        if (rfc.strip() == ''):
            rfc = 'unknown'
        ## url  = usgs.get_usgs_site_url( sid )
        site_info[ sid ] = \
            {'rfc_name':rfc, 'huc8':huc8, 'lat':lat, 'lon':lon }
        k += 1

    print('Processed', k, 'USGS/RFC sites.')

    if (SAVE_TO_FILE):
        file_unit = open( file_path, 'wb' )
        pickle.dump( site_info, file_unit, protocol=pickle.HIGHEST_PROTOCOL)
        file_unit.close()
        print('Saved USGS/RFC site info dictionary to file:')
        print('  ' + file_path)
        print()

    return site_info
                    
#   get_usgs_rfc_dict_v0()
#---------------------------------------------------------------------
def get_usgs_rfc_dict_v1( SAVE_TO_FILE=True ):

    #-------------------------------------------------------------    
    # Construct a site_info dictionary, where USGS site ID
    # is the key used to get dictionary info for a site.
    # This version uses the file: ALL_USGS-HADS_SITES.txt,
    # which was obtained from:
    #   https://hads.ncep.noaa.gov/USGS/ALL_USGS-HADS_SITES.txt
    # and has info for 10510 sites/basins.
    # Note that this file is sorted by 5-char NWS location ID.
    #-------------------------------------------------------------
    new_dir   = dtu.get_new_data_dir( 'NOAA_HADS' )
    data_dir  = dtu.get_data_dir( 'NOAA_HADS' )
    #------------------------------------------------
    dict_file = 'USGS_RFC_site_info_dict.pkl'
    file_path = new_dir + dict_file
    if (os.path.exists( file_path )):
        print('Reading saved RFC site info dictionary...')
        file_unit = open( file_path, 'rb')
        site_info = pickle.load( file_unit )
        file_unit.close()
        print('Finished.')
        print()
        return site_info

    print('Getting USGS/RFC site info as dictionary...')
    info_file = 'All_USGS-HADS_SITES.tsv'
    info_path = data_dir + info_file
    info_unit = open( info_path, 'r' )
    delim = '\t'
    
    #--------------------------------
    # Skip over three header lines
    #--------------------------------
    for j in range(3):
        line = info_unit.readline()
    
    site_info = dict()
    k = 0

    while (True):
        info_line = info_unit.readline()
        if (info_line == ''):
            break  # (reached end of file)
        vals = info_line.split( delim )
        #-------------------------------------
        # Note that RFC abbrev. is not given
        #-------------------------------------
        nws_id   = vals[0].strip()
        usgs_id  = vals[1].strip()
        goes_id  = vals[2].strip()
        hsa_id   = vals[3].strip()
        lat_dms  = vals[4].strip()
        lon_dms  = vals[5].strip()
        loc_name = vals[6].strip()
        #--------------------------------------------
        lat = dtu.convert_dms_to_dec_deg( lat_dms )
        lon = dtu.convert_dms_to_dec_deg( lon_dms )
                
        ## url = usgs.get_usgs_site_url( usgs_id )
        site_info[ usgs_id ] = \
            {'nws_id':nws_id, 'goes_id':goes_id, 'hsa_id':hsa_id,
             'lat':lat, 'lon':lon }  ### 'loc_name':loc_name }
        k += 1

    print('Processed', k, 'USGS/RFC sites.')

    if (SAVE_TO_FILE):
        file_unit = open( file_path, 'wb' )
        pickle.dump( site_info, file_unit, protocol=pickle.HIGHEST_PROTOCOL)
        file_unit.close()
        print('Saved USGS/RFC site info dictionary to file:')
        print('  ' + file_path)
        print()

    return site_info

#   get_usgs_rfc_dict_v1()
#---------------------------------------------------------------------
def get_hydrograph_type_dict( SAVE_TO_FILE=True ):

    new_dir   = dtu.get_new_data_dir( 'NOAA_HADS' )
    dict_file = 'RFC_hydrograph_type_dict.pkl'
    pkl_path  = new_dir + dict_file
    if (os.path.exists( pkl_path )):
        print('Reading saved RFC hydrograph type dictionary...')
        pkl_unit = open( pkl_path, 'rb')
        htype_dict = pickle.load( pkl_unit )
        pkl_unit.close()
        print('Finished.')
        print()
        return htype_dict
    
    print('Getting RFC hydrograph type as dictionary...')

    #-------------------------------------------------    
    # Open each of the hydrograph type files to read
    # Each file is sorted by USGS site ID
    #-------------------------------------------------
    delim = ','
    htype_dict = dict()

    #--------------------------------    
    # Make a list of CSV file names
    #--------------------------------
    file_names = list()
    file_names.append( 'NWMv3_basin_type_flashy.txt' ) 
    file_names.append( 'NWMv3_basin_type_slow.txt' )
    file_names.append( 'NWMv3_basin_type_regulation.txt' )
    file_names.append( 'NWMv3_basin_type_snow.txt' )

    #----------------------------------    
    # Make a list of hydrograph types
    #----------------------------------
    htypes = ['flashy', 'slow', 'regulated', 'snow-dom']
    
    #---------------------------------------------   
    # Create a dictionary with hydrograph types
    # and other info for each USGS site ID found
    #---------------------------------------------
    for k in range(len(file_names)):
        htype = htypes[k]
        file_path = new_dir + file_names[k]
        file_unit = open( file_path, 'r' )
        line = file_unit.readline()  # skip 1 header line
        while (True):
            line = file_unit.readline()
            if (line == ''):
                break  # (reached end of file)
            vals = line.split( delim )
            sid  = vals[0].strip()
            lat  = vals[1].strip()
            lon  = vals[2].strip()
            huc8 = vals[6].strip()
            rfc  = vals[7].strip()
            ## url  = usgs.get_usgs_site_url( sid )
            if (lat != 'NA'):  ############## ADD MISSING INFO LATER #######
                htype_dict[ sid ] = \
                    {'htype':htype, 'rfc_name':rfc, 'huc8':huc8,
                     'lat':lat, 'lon':lon }       
        file_unit.close() 

    if (SAVE_TO_FILE):
        pkl_unit = open( pkl_path, 'wb' )
        pickle.dump( htype_dict, pkl_unit, protocol=pickle.HIGHEST_PROTOCOL)
        pkl_unit.close()
        print('Saved RFC hydrograph info dictionary to file:')
        print('  ' + pkl_path)
        print()
        
    return htype_dict 
    
#   get_hydrograph_type_dict()
#---------------------------------------------------------------------
def cull_usa_dcp_data_rows( in_csv_file=None, out_csv_file=None):

    #--------------------------------------------------------------
    #  Note: This website at Iowa State Univ. provides DCP data
    #        for many networks that can be downloaded as CSV
    #        or as a shapefile:
    #  https://mesonet.agron.iastate.edu/sites/networks.php
    #--------------------------------------------------------------
    # Downloaded a CSV with DCP info for all All Networks.
    # This function culls just the ones for US states.
    #--------------------------------------------------------------
    # If viewing HTML Table for a state, it has an "Archive Ends"
    # column for inactive DCPs.  However, if saving as CSV, only
    # ACTIVE DCPs are included and this column is not provided.
    #--------------------------------------------------------------
    # It is likely that many of these DCPs are not of the type
    # "Stream", vs. "Atmosphere", "Well", etc.
    #--------------------------------------------------------------  
    # Wrote an email to Daryl Herzmann (akrherz@iastate.edu) 
    # on 2023-11-09 to request more info.
    # See: https://mesonet.agron.iastate.edu/info/contacts.php
    #--------------------------------------------------------------
    data_dir  = dtu.get_data_dir( 'NOAA_DCP_via_IEM' )          
    if (in_csv_file is None):
        in_csv_file = 'All_Networks_Active_DCP_Data.csv'
    in_csv_path = data_dir + in_csv_file
    #-----------------------------------------------------
    new_dir = dtu.get_new_data_dir( 'NOAA_DCP_via_IEM' )
    if (out_csv_file is None):
        out_csv_file = 'US_State_Active_DCP_Data.csv'
    out_csv_path = new_dir + out_csv_file
    #-----------------------------------------------------    
    in_csv_unit  = open( in_csv_path,  'r' )
    out_csv_unit = open( out_csv_path, 'w' )  
    comma  = ','
    n_rows_in   = 0
    n_rows_out  = 0
    n_rows_skip = 0
    #------------------------------------------------
    # Skip these territories, etc.
    # Keep PR (Puerto Rico) and VI (Virgin Islands)
    # GU_DCP is Guam.
    #------------------------------------------------
    skip_list  = ['GU_DCP','P1_DCP','P2_DCP','P3_DCP','P4_DCP']
    state_dict = dict()
    #---------------------------------------
    # Copy header from in_file to out_file
    #---------------------------------------
    header = in_csv_unit.readline()
    out_csv_unit.write( header ) 
            
    #--------------------------------------
    # Read each data line in the CSV File
    # All US DCP networks end with "_DCP"
    #--------------------------------------
    # Canadian networks have names like:
    #    CA_AB_DCP (Canada, Alberta)
    # Mexican networks have names like:
    #    MX_BJ_DCP (Mexico, Baja)
    # Puerto Rico is included:  PR_DCP
    # Guam is included: GU_DCP
    # Barbuda is included: AG__DCP
    # Others are BD_DCP, BM_DCP (Bermuda)
    #--------------------------------------                                       
    while (True):
        line = in_csv_unit.readline()
        if (line == ''):
            break  # (reached end of file)
        n_rows_in += 1
        #-------------------------------------
#         if (line.endswith('_DCP\n')):
#             out_csv_unit.write( line )
#             n_rows_out += 1
        #-------------------------------------        
        vals = line.split( comma )
        name = vals[-1].strip()  # last entry
        KEEP = not(name in skip_list)
        if (KEEP and name.endswith('_DCP') and (len(name) == 6)):
            #--------------------------------------
            # e.g. AK_DCP, CO_DCP, GU_PCP, PR_DCP
            # Could exclude Guam, others.
            #--------------------------------------
            out_csv_unit.write( line )
            n_rows_out += 1
            #--------------------------------------
            state_code = name[0:2]
            try:
                state_dict[ state_code ] += 1
            except:
                state_dict[ state_code ] = 0
        else:
            #--------------------------------------
            # Skip Canada, Mexico, & Barbuda DCPs
            #--------------------------------------
            n_rows_skip += 1  
    
    in_csv_unit.close()
    out_csv_unit.close()
 
    #-------------------------------------   
    # Print the count of each state code
    #-------------------------------------
    key_list = list( state_dict.keys())
    key_list.sort()
    for key in key_list:
        print('# ' + key, state_dict[key])
    print()

    print('Read   ', n_rows_in,  'rows.')
    print('Skipped', n_rows_skip,'rows.')
    print('Wrote  ', n_rows_out, 'rows.')
    print('Finished.')
    print()

#   cull_usa_dcp_data_rows()
#---------------------------------------------------------------------
def convert_json_to_csv( json_file=None, csv_file=None, 
                         rfc_name='WGRFC', REPORT=True ):

    #-------------------------------------------------------------
    # Note: The shapefile for the West Gulf RFC (WGRFC) had an
    #       empty attribute table.  However, was able to find
    #       a JSON file on an ESRI website with the attributes.
    #       The URL was:
    #       https://www.arcgis.com/home/item.html?
    #           id=fcbec367432d4b5a829a250906189fd6
    #       Then clicked on "Source: Feature Collection",
    #       downloaded to text and saved with JSON extension.
    #       This JSON file has info for 624 basins.
    #-------------------------------------------------------------
    #       It may be possible to access info for other RFCs in
    #       this same way, which gives additional attributes.
    #-------------------------------------------------------------
    new_dir  = dtu.get_new_data_dir( 'NOAA_RFC' )
    new_dir += 'CSVs_for_each_RFC/'
    #-------------------------------------               
    if (json_file is None):
        json_file = rfc_name + '_data_from_ESRI_online.json'
    json_path = new_dir + json_file
    #-------------------------------------
    if (csv_file is None):
        csv_file = json_file.replace('.json', '.csv')
    csv_path = new_dir + csv_file
    #-------------------------------------    
    json_unit = open( json_path, 'r' )
    csv_unit  = open( csv_path,  'w' )
    #-------------------------------------     
    data_dict    = json.load( json_unit )
    layer_list   = data_dict['layers']
    layer0_dict  = layer_list[0]
    featSet_dict = layer0_dict['featureSet']
    feature_list = featSet_dict['features']

    if (REPORT):
        print('data_dict keys    =', list(data_dict.keys()))
        print('len(layer_list)   =', len(layer_list))  # Is "1" here.
        print('layer0_dict keys  =', list(layer0_dict.keys()))
        print('featSet_dict keys =', list(featSet_dict.keys()))
        print('len(feature_list) =', len(feature_list))  # Is "624" here.        
        print('')

    #--------------------------------
    # Write header for new CSV file
    #--------------------------------
    delim  = ',' 
    header = ''
    header += 'RFC_ID'     + delim
    header += 'NWS_loc_ID' + delim
    header += 'MINLON'     + delim    # MIN_X_AXIS
    header += 'MAXLON'     + delim    # MAX_X_AXIS
    header += 'MINLAT'     + delim    # MIN_Y_AXIS
    header += 'MAXLAT'     + delim    # MAX_Y_AXIS
    header += 'NAME'       + '\n'     # newline at end
    csv_unit.write( header )   

    n_basins = 0
    for feature in feature_list:
        #--------------------------------------
        # Each feature is a Python dictionary
        # Keys are: 'geometry', 'attributes'.
        #--------------------------------------
        att_dict = feature['attributes']
        if (REPORT and (n_basins == 0)):
            print('feature_dict keys =', list(feature.keys()) )
            print('att_dict keys     =', list(att_dict.keys()) )
            print()

        name   = att_dict['NAME']        # string
        nws_id = att_dict['CH5_ID']      # string
        minlon = att_dict['MIN_X_AXIS']  # float
        maxlon = att_dict['MAX_X_AXIS']
        minlat = att_dict['MIN_Y_AXIS']
        maxlat = att_dict['MAX_Y_AXIS']
        n_basins += 1      

        #---------------------------------
        # Format the lat and lon strings
        #---------------------------------
        minlon_str = '{x:.5f}'.format(x=minlon)
        maxlon_str = '{x:.5f}'.format(x=maxlon)
        minlat_str = '{x:.5f}'.format(x=minlat)
        maxlat_str = '{x:.5f}'.format(x=maxlat)
            
        #-----------------------------       
        # Write info to new CSV file
        #-----------------------------
        csv_line = ''
        csv_line += rfc_name   + delim
        csv_line += nws_id     + delim
        csv_line += minlon_str + delim
        csv_line += maxlon_str + delim
        csv_line += minlat_str + delim
        csv_line += maxlat_str + delim
        csv_line += name + '\n'      # newline at end
        csv_unit.write( csv_line )
        
    json_unit.close()     
    csv_unit.close()
    print('Wrote info for', n_basins, 'basins.')
    print('Finished.')
    print()
   
#   convert_json_to_csv()
#---------------------------------------------------------------------
def merge_rfc_basin_info( tsv_file='new_rfc_info.tsv', REPORT=True):

    #----------------------------------------------------------------
    # Note: Found a site that provides a shapefile of basins for
    #       each of the 13 NOAA RFCs at:
    #          https://www.nohrsc.noaa.gov/gisdatasets/
    #       Exported shapefile attributes to CSV from QGIS.
    #       The 13 shapefiles do not provide the same set of
    #       attributes but they have several in common, such as
    #       the bounding box info for each basin.
    #       This also allows us to map 5-char NWS IDs to RFC names.
    #----------------------------------------------------------------
    #       The West Gulf RFC shapefile has no attributes and also
    #       seems to not have all of the basin boundaries.
    #----------------------------------------------------------------
    #       The Missouri Basin RFC has numbers instead of 5-char
    #       NWS location IDs for most basins.  On this page:
    #       https://water.weather.gov/ahps/region.php?rfc=mbrfc
    #       after clicking the radio button "River Observations",
    #       it says there are 1185 total gauges.  But the shapefile
    #       has info for 1412 basins.
    #----------------------------------------------------------------
    new_dir  = dtu.get_new_data_dir( 'NOAA_RFC' )
    csv_dir  = new_dir + 'CSVs_for_each_RFC/'

    rfc_list = [
    'abrfc', 'aprfc', 'cbrfc', 'cnrfc', 'lmrfc',
    'marfc', 'mbrfc', 'ncrfc', 'nerfc', 'nwrfc',
    'ohrfc', 'serfc']
    #### 'ohrfc', 'serfc', 'wgrfc']   # wgrfc shapefile is corrupt
    prefix_list = ['b_' + s for s in rfc_list ]
    csv_paths = list()
    
    for prefix in prefix_list:
        csv_path = csv_dir + prefix + '/' + prefix + '.csv'
        csv_paths.append( csv_path )

    #-----------------------------
    # Open new TSV file to write
    #-----------------------------
    tsv_path = new_dir + tsv_file
    tsv_unit = open( tsv_path, 'w')

    #--------------------------------
    # Write header for new TSV file
    #--------------------------------
    tab    = '\t'  # tab character  
    header = ''
    header += 'RFC_ID'     + tab
    header += 'NWS_loc_ID' + tab
    header += 'MINLON'     + tab    # MIN_X_AXIS_
    header += 'MAXLON'     + tab    # MAX_X_AXIS_
    header += 'MINLAT'     + tab    # MIN_Y_AXIS_
    header += 'MAXLAT'     + tab    # MAX_Y_AXIS_
    header += 'NAME'       + '\n'   # newline at end
    tsv_unit.write( header )   
       
    k = 0
    comma = ','
    n_total = 0
 
    for prefix in prefix_list:
        rfc_name = rfc_list[k].upper()
        csv_unit = open( csv_paths[k], 'r' )
        if (REPORT):
            print('For RFC = ' + rfc_name + ':')
            
        #--------------------------------------
        # Get info from the CSV header line
        #--------------------------------------
        # Find col nums of certain attributes
        # NCRFC has many other attributes.
        # Several have centroid lon, lat:
        #   e.g. NCRFC, NERFC, 
        #--------------------------------------
        n_basins = 0
        header   = csv_unit.readline()
        headings = header.split( comma )
        if (REPORT):
            print('   number of columns =', len(headings) )
        if (len(headings) <= 1):
            print('   ### ERROR: No data found for this RFC.')
            break
        name_col   = headings.index('NAME')
        minlon_col = headings.index('MIN_X_AXIS_')
        maxlon_col = headings.index('MAX_X_AXIS_')
        minlat_col = headings.index('MIN_Y_AXIS_')
        maxlat_col = headings.index('MAX_Y_AXIS_')
        if ('NWS_ID' in headings):
            nwsid_col = headings.index('NWS_ID')
        elif ('NWSID' in headings):
            nwsid_col = headings.index('NWSID')  # LMRFC
        elif ('CH5_ID' in headings):
            nwsid_col = headings.index('CH5_ID')
        elif ('BASIN' in headings):
            nwsid_col = headings.index('BASIN')  # CNRFC
        elif ('NAME' in headings):
            # For MBRFC, NAME is NWS_ID + "UPR" or "LWR"
            # It also has SHAPE_LENG, SHAPE_AREA,
            #  CENTROID_L (lat), CENTROID_1 (lon)
            nwsid_col = headings.index('NAME')  # MBRFC
        else:
            break  # WGRFC (West Gulf) has no attributes.

        #--------------------------------------
        # Read each data line in the CSV File
        #--------------------------------------                                   
        while (True):
            csv_line = csv_unit.readline()
            if (csv_line == ''):
                break  # (reached end of file)
            vals = csv_line.split( comma )
                                                              
            nws_loc_id = vals[ nwsid_col ].strip()
            name       = vals[ name_col ].strip()
            if (len(nws_loc_id) == 0):
                nws_loc_id = name  # this occurs for OHRFC
            #----------------------------------------
            # Get geographic bounding box for basin
            #----------------------------------------
            minlon = np.float64( vals[ minlon_col ].strip() )
            maxlon = np.float64( vals[ maxlon_col ].strip() )
            minlat = np.float64( vals[ minlat_col ].strip() )
            maxlat = np.float64( vals[ maxlat_col ].strip() )
            #---------------------------------
            # Format the lat and lon strings
            #---------------------------------
            minlon_str = '{x:.5f}'.format(x=minlon)
            maxlon_str = '{x:.5f}'.format(x=maxlon)
            minlat_str = '{x:.5f}'.format(x=minlat)
            maxlat_str = '{x:.5f}'.format(x=maxlat)
                
            #-----------------------------       
            # Write info to new TSV file
            #-----------------------------
            tsv_line = ''
            tsv_line += rfc_name   + tab
            tsv_line += nws_loc_id + tab
            tsv_line += minlon_str + tab
            tsv_line += maxlon_str + tab
            tsv_line += minlat_str + tab
            tsv_line += maxlat_str + tab
            tsv_line += name + '\n'      # newline at end
            tsv_unit.write( tsv_line )
            n_basins += 1
            n_total  += 1

        if (REPORT):
            print('   number of basins  =', n_basins)
        csv_unit.close()
        k += 1

    #-----------------------------------------------
    # Add info for the WGRFC by a different method
    #-----------------------------------------------
    n_basins = 0
    rfc_name       = 'WGRFC'
    wgrfc_csv_file = 'WGRFC_data_from_ESRI_online.csv'
    wgrfc_csv_path = csv_dir + wgrfc_csv_file
    csv_unit = open( wgrfc_csv_path, 'r' )
    header   = csv_unit.readline()
    if (REPORT):
        print('For RFC = ' + rfc_name + ':')
        print('   by reading info from JSON file')
    #--------------------------------------
    # Read each data line in the CSV File
    #--------------------------------------                                   
    while (True):
        csv_line = csv_unit.readline()
        if (csv_line == ''):
            break  # (reached end of file)
        vals = csv_line.split( comma )
                                                          
        nws_loc_id = vals[1].strip()
        minlon_str = vals[2].strip()  # already formatted
        maxlon_str = vals[3].strip()
        minlat_str = vals[4].strip()
        maxlat_str = vals[5].strip()        
        name       = vals[6].strip()
        #-----------------------------       
        # Write info to new TSV file
        #-----------------------------
        tsv_line = ''
        tsv_line += rfc_name   + tab
        tsv_line += nws_loc_id + tab
        tsv_line += minlon_str + tab
        tsv_line += maxlon_str + tab
        tsv_line += minlat_str + tab
        tsv_line += maxlat_str + tab
        tsv_line += name + '\n'      # newline at end
        tsv_unit.write( tsv_line )
        n_basins += 1
        n_total  += 1

    csv_unit.close()
    tsv_unit.close()
    
    if (REPORT):
        print('   number of basins  =', n_basins)
    print('Total number of RFC basins =', n_total)
    print('Finished merging RFC basin info to create:')
    print('   ' + tsv_path)
    print()
         
#   merge_rfc_basin_info()
#---------------------------------------------------------------------
def create_rfc_tsv( data_dir=None, new_dir=None, nf_max=50,
                shape_file='ba12my15.shp',                    
                prj_file  ='ba12my15.prj',
                tsv_file='new_rfc_info.tsv',
                SWAP_XY=True,   # True works for 'NOAA_RFC'
                REPORT=True):

    if (data_dir is None):
        data_dir = dtu.get_data_dir( 'NOAA_RFC' )
    if (new_dir is None):
        new_dir  = dtu.get_new_data_dir( 'NOAA_RFC' )
    insert_key = 'LAT'

    su.create_tsv_from_shapefile(data_dir=data_dir, new_dir=new_dir,
                      shape_file=shape_file, prj_file=prj_file,
                      tsv_file=tsv_file, nf_max=nf_max,
                      SWAP_XY=SWAP_XY, REPORT=REPORT,
                      ADD_BOUNDING_BOX=True, insert_key=insert_key,
                      filter_key=None, filter_value=None,
                      out_lon_heading='LON', out_lat_heading='LAT')

#   create_rfc_tsv()
#---------------------------------------------------------------------
def create_hads_tsv( data_dir=None,
           tsv_file='new_hads_info.tsv', REPORT=True):

    if (data_dir is None):
        data_dir  = dtu.get_data_dir( 'NOAA_HADS' )
    hads_file = 'ALL_USGS-HADS_SITES.txt'
    hads_path  = data_dir + hads_file
    hads_unit  = open( hads_path, 'r' )
    hads_delim = '|'

    #---------------------------------------
    # Skip over HADS file header (4 lines)
    #---------------------------------------
    for i in range(4):
       line = hads_unit.readline()

    #-----------------------------
    # Open new TSV file to write
    #-----------------------------
    new_dir  = dtu.get_new_data_dir( 'NOAA_HADS' )
    tsv_path = new_dir + tsv_file
    tsv_unit = open( tsv_path, 'w')

    #--------------------------------
    # Write header for new TSV file
    #--------------------------------
    delim  = '\t'  # tab character  
    header = ''
    header += 'USGS_ID'      + delim
    header += 'USGS_name'    + delim
    header += 'Site_type'    + delim
    header += 'Longitude'    + delim
    header += 'Latitude'     + delim
    header += 'RFC_ID'       + delim
    header += 'NWS_loc_ID'   + delim
    header += 'NWS_HSA_ID'   + delim
    header += 'GOES_ID'      + delim
    header += 'HUC8'         + delim
    header += 'Hgraph_type'  + '\n'  # newline at end  
    tsv_unit.write( header )   

    rfc_map   = get_usgs_rfc_dict()      ############### SEE NEW METHOD
    htype_map = get_hydrograph_type_dict()
    null_info = {'rfc_name': 'unknown', 'huc8': 'unknown',
                 'lat':'-999', 'lon':'-999'}
    site_type_dict = usgs.get_site_type_dict()

    k = 0
    lat_diff_count = 0
    lon_diff_count = 0

    while (True):
        hads_line = hads_unit.readline()
        if (hads_line == ''):
            break  # (reached end of file)
        vals = hads_line.split( hads_delim )
        nws_loc_id = vals[0].strip()
        usgs_id    = vals[1].strip()
        goes_id    = vals[2].strip()
        nws_hsa_id = vals[3].strip()
        lat_dms    = vals[4].strip()
        lon_dms    = vals[5].strip()
        if (lon_dms[0] != '-'):
            #--------------------------
            # Add missing minus sign.
            # One entry does have one.
            #--------------------------
            lon_dms = '-' + lon_dms
        loc_name   = vals[6].strip()
        #-----------------------------------
        try:
            rfc_info = rfc_map[ usgs_id ]
        except:
            rfc_info = null_info
        rfc_name   = rfc_info['rfc_name']
        huc8_code  = rfc_info['huc8']
        rfc_lat    = rfc_info['lat']
        rfc_lon    = rfc_info['lon']
        if (rfc_lat == 'NA'):
            rfc_lat = '-999'
        if (rfc_lon == 'NA'):
            rfc_lon = '-999'
        #-----------------------------------
        lat = dtu.convert_dms_to_dec_deg( lat_dms )
#         dms_list = lat_dms.split()  # split on white space
#         D = np.float64( dms_list[0] )
#         M = np.float64( dms_list[1] )
#         S = np.float64( dms_list[2] )
#         if (D < 0):
#             lat = D - (M/60.0) - (S/3600.0)
#         else:
#             lat = D + (M/60.0) + (S/3600.0)
        #-----------------------------------
        lon = dtu.convert_dms_to_dec_deg( lon_dms )
#         dms_list = lon_dms.split()  # split on white space
#         D = np.float64( dms_list[0] )
#         M = np.float64( dms_list[1] )
#         S = np.float64( dms_list[2] )
#         if (D < 0):
#             lon = D - (M/60.0) - (S/3600.0)
#         else:
#             lon = D + (M/60.0) + (S/3600.0)

        #---------------------------------
        # Format the lat and lon strings
        #---------------------------------
        lon_str = '{x:.5f}'.format(x=lon)
        lat_str = '{x:.5f}'.format(x=lat)
        
        #-------------------------------------
        # Compare lon/lat to rfc_lon/rfc_lat
        #-----------------------------------------------------------
        # Note:  rfc_lon & rfc_lat are only available for the 8170
        # basins listed in usgsMasterMeta.txt, but lon & lat from
        # DMS are available for the 10512 basins in the HADS file.
        #-----------------------------------------------------------
        # 1601 USGS IDs in the HADS file do not occur in the
        # USGS_NWIS_all (streams only), and appear to be other
        # site types like atmosphere, canal, lake, and well.
        # See: get_site_type_dict() in usgs_utils.py. 
        #-----------------------------------------------------------
        # See:  https://en.wikipedia.org/wiki/Decimal_degrees
        #-----------------------------------------------------------
        ## tol = 0.0001  # difference over 10 meters
        tol = 0.001   # difference over 100 meters
        if (rfc_lat != '-999'):
            lat_diff = np.abs(lat - np.float64(rfc_lat))
            if (lat_diff > tol):
                lat_diff_count += 1
        if (rfc_lon != '-999'):
            lon_diff = np.abs(lon - np.float64(rfc_lon))
            if (lon_diff > tol):
                lon_diff_count += 1        

        #-------------------------------------
        # Get hydrograph type classification
        # (not available for all sites)
        #-------------------------------------
        try:
            htype_info = htype_map[ usgs_id ]
            htype = htype_info[ 'htype' ]
        except:
            htype = 'unknown'

        #---------------------------          
        # Get the site type string
        #---------------------------
        try:
            site_type = site_type_dict[ usgs_id ] 
        except:
            site_type = 'unknown'
          
        #-----------------------------
        # Write line to new TSV file
        #-----------------------------
        line = ''
        line += usgs_id      + delim
        line += loc_name     + delim
        line += site_type    + delim
        line += lon_str      + delim
        line += lat_str      + delim
        line += rfc_name     + delim
        line += nws_loc_id   + delim
        line += nws_hsa_id   + delim
        line += goes_id      + delim
        line += huc8_code    + delim
        line += htype        + '\n'   # add newline at end
        tsv_unit.write( line )       
        k += 1

    #------------------ 
    # Close the files
    #------------------
    hads_unit.close()
    tsv_unit.close()
    print('Finished creating new TSV file for RFC basins.')
    print('  lon_diff_count =', lon_diff_count)
    print('  lat_diff_count =', lat_diff_count)
    print()
                
#   create_hads_tsv()
#---------------------------------------------------------------------
def sort_by_site_code(tsv_file1='new_hads_info.tsv',
                      tsv_file2='new_hads_info_sorted.tsv'):

    new_dir   = dtu.get_new_data_dir( 'NOAA_HADS' )
    tsv_path1 = new_dir + tsv_file1
    tsv_path2 = new_dir + tsv_file2
    
    tsv_unit1 = open(tsv_path1, 'r')
    lines = tsv_unit1.readlines()
    tsv_unit1.close()
        
    #--------------------------------
    # Sort all lines (by "USGS_ID")
    #--------------------------------
    header = lines[0]   # header line  
    lines2 = lines[1:]    
    lines2.sort()
    tsv_unit2 = open(tsv_path2, 'w')
    tsv_unit2.writelines( [header] + lines2 )
    tsv_unit2.close()

#   sort_by_site_code()
#---------------------------------------------------------------------
def get_site_data_from_api( gauge_id='13334300', REPORT=True,
                            RAW_JSON=False):

    #-----------------------------------------------------------
    # Note: This function uses a new (beta) API from NOAA
    #       to retrieve information for a given gauge.
    #       In the API, gauge_ID can be USGS ID or NWS Loc ID.
    #-----------------------------------------------------------
    # Q:  What is best way to check if this is stream gauge?
    # A:  The PEDTS Code seems to contain this information.
    #-----------------------------------------------------------
    # The "PEDTS Code" is explained in the NWS SHEF Manual.
    # The first 2 letters are "Physical Element" (PE) code.
    #    HG stands for "gauge height".
    #    QR stands for "discharge".
    #    PP stands for "precipitation".
    # The 3rd letter is the Duration Code and can be:
    #    I (instantaneous), H (hourly), D (daily), etc.
    # The 4th letter is the Type Code and can be:
    #    R (reading/observed), F (forecast), etc.
    # The 5th letter is the Source Code and is used to
    #    specify how the data was created or transmitted.
    # Note that "RG" = GOES observation.
    # Note that "FZ" = nonspecific forecast data (default)
    # See Tables 3 and 4 in the NWS SHEF Manual for details.
    #-----------------------------------------------------------
    # 'https://preview-api.water.noaa.gov/v1/gauges/13334300'
    #-----------------------------------------------------------     
    url = 'https://preview-api.water.noaa.gov/v1/gauges/' + gauge_id
    response  = requests.get( url )  # returns JSON as dict
    data_dict = response.json()
    if (RAW_JSON):
        return data_dict

    #------------------------------------------------------------------------------
    # Note: If the API can't find a gauge ID, it returns this:
    # {'code': 5, 'message': '[01016500] could not find gauge ID', 'details': []}
    #------------------------------------------------------------------------------
    if ('code' in data_dict):   # Python 3 syntax
        if (data_dict['code'] == 5):
            if (REPORT):
                print('### SORRY, API could not find gauge_ID:', gauge_id)
            return None

    #------------------------------------------------
    # Extract some info and store in new dictionary
    #------------------------------------------------
    new_dict = dict()
    try:
        new_dict['usgs_id'] = data_dict['usgsId']
    except:
        #------------------------------------
        # Still need this for USGS_NWIS2
        #------------------------------------
        new_dict['usgs_id'] = gauge_id   
    new_dict['nws_loc_id']   = data_dict['lid']
    new_dict['upstrm_lid']   = data_dict['upstreamLid']   
    new_dict['dnstrm_lid']   = data_dict['downstreamLid']
    new_dict['nwm_reach_id'] = data_dict['reachId']
    new_dict['rfc_abbrev']   = data_dict['rfc']['abbreviation']
    new_dict['wfo_abbrev']   = data_dict['wfo']['abbreviation']
    new_dict['usgs_name']    = data_dict['name']
    new_dict['description']  = data_dict['description']
    new_dict['state_code']   = data_dict['state']['abbreviation']
    new_dict['county_name']  = data_dict['county']
    new_dict['longitude']    = data_dict['longitude']
    new_dict['latitude']     = data_dict['latitude']
    new_dict['elevation']    = data_dict['datum']['elevation']
    new_dict['vert_datum']   = data_dict['datums']['vertical']['value']
    new_dict['horiz_datum']  = data_dict['datums']['horizontal']['value']
    new_dict['updated_time'] = data_dict['status']['updatedTime']
    if (new_dict['updated_time'] is None):
        new_dict['updated_time'] = ''
    new_dict['in_service']   = data_dict['inService']['enabled']
    new_dict['pedts_obs']    = data_dict['pedts']['observed'] 
    new_dict['pedts_pred']   = data_dict['pedts']['forecast']
                 
    if (REPORT):
        print('USGS ID        =', new_dict['usgs_id'])
        print('NWS LID        =', new_dict['nws_loc_id'])
        print('Upstream LID   =', new_dict['upstrm_lid'])
        print('Downstream LID =', new_dict['dnstrm_lid'])
        print('NWM Reach ID   =', new_dict['nwm_reach_id'])
        print('RFC name       =', new_dict['rfc_abbrev'])
        print('WFO name       =', new_dict['wfo_abbrev'])
        print('USGS name      =', new_dict['usgs_name'])
        print('Description    =', new_dict['description'])            
        print('US State Code  =', new_dict['state_code'])
        print('US County      =', new_dict['county_name'])
        print('Longitude      =', new_dict['longitude'])
        print('Latitude       =', new_dict['latitude'])
        print('Elevation      =', new_dict['elevation'])
        print('Vert Datum     =', new_dict['vert_datum'])
        print('Horiz Datum    =', new_dict['horiz_datum'])
        print('Last updated   =', new_dict['updated_time'])
        print('InService      =', new_dict['in_service'])
        print('PEDTS (obs)    =', new_dict['pedts_obs'])
        print('PEDTS (pred)   =', new_dict['pedts_pred'])
        print()
    
    return new_dict

#   get_site_data_from_api()
#---------------------------------------------------------------------
def create_tsv_via_api( base_key=None, ADD_INFO=False):
                        ###### max_lines=100 ):

    #---------------------------------------------------------------
    # Note:  This function calls get_site_data_from_api(), above,
    #        to get info for every USGS site ID in one of the
    #        USGS basin collections.
    #        "base_key" can be: USGS_NWIS_Web, USGS_NWIS_WQP1,
    #        USGS_NWIS_WQP2, USGS_NWIS_WQP3, USGS_GAGES2_all, etc.
    #---------------------------------------------------------------
    start_time = time.time()
    if (base_key is None):
        base_key = 'USGS_NWIS_Web'   # 27915 data rows; 
    # base_key = 'USGS_NWIS_WQP1'    # 14707 data rows;
    # base_key = 'USGS_NWIS_WQP2'    # 46760 data rows;
    # base_key = 'USGS_NWIS_WQP3'    # 145375 data rows; 9.56 hrs.
    # base_key = 'USGS_GAGES2_all'   # 9322 data rows

    print('Base key =', base_key)
    print('Working...')
    usgs_id_col = dtu.get_id_column( base_key )
    n_sites    = 0
    n_skipped  = 0
    n_lines    = 0
    out_tsv_file = base_key + '_Site_Info_via_API.tsv'
    #----------------------------------------------------
    new_dir = dtu.get_new_data_dir( 'NOAA_via_API' )
    out_tsv_path = new_dir + out_tsv_file
    out_tsv_unit = open( out_tsv_path, 'w') 

    #----------------------------------------------
    # Get hsa_id and goes_id from this dictionary
    #----------------------------------------------
    if (ADD_INFO):
        usgs_rfc_dict = get_usgs_rfc_dict_v1()
    
    #-------------------------------------------
    # Open the TSV file with USGS IDs and info
    #-------------------------------------------
    delim = '\t'  # tab
    usgs_path = dtu.get_new_tsv_filepath( base_key )   
    usgs_unit = open(usgs_path, 'r')
    #-------------------        
    # Skip header line
    #-------------------
    dtu.skip_header_lines( usgs_unit, key=base_key)

    #-----------------------------------------
    # Write column headings for new TSV file
    #-----------------------------------------
    headings = [
    'usgs_id', 'nws_loc_id', 'upstrm_lid', 'dnstrm_lid',
    'nwm_reach_id', 'rfc_abbrev', 'wfo_abbrev']
    if (ADD_INFO):
        headings = headings + ['hsa_id', 'goes_id']   # from USGS-HADS crosswalk
    headings = headings + [
    'usgs_name', 'description', 'state_code',
    'county_name', 'longitude', 'latitude', 'elevation', 'updated_time',
    'in_service', 'pedts_obs', 'pedts_pred' ]
    header = ''
    for h in headings:
        header += h + delim
    header = header[:-1] + '\n'
    out_tsv_unit.write( header ) 

    while (True):
        usgs_line = usgs_unit.readline()
        if (usgs_line == ''):
            break  # (reached end of file)
        n_lines += 1
#         if (n_lines > max_lines):
#             break

        values = usgs_line.split( delim )   
        usgs_id = values[ usgs_id_col ].strip()
        if (base_key.startswith('USGS_NWIS_WQP')):
            usgs_id = usgs_id.replace('USGS-','')    ###### NEED THIS!
        new_dict = get_site_data_from_api( gauge_id=usgs_id, REPORT=False)
        if (new_dict is None):
            # print('### Skipping USGS ID:', usgs_id)
            n_skipped += 1
            if ((n_skipped % 50) == 0):
                print('status: n_skipped  =', n_skipped)
            continue

        #------------------------------
        # Format some numeric strings
        #------------------------------
        lon      = new_dict['longitude']
        lat      = new_dict['latitude']
        elev      = new_dict['elevation']
        lon_str  = '{x:.5f}'.format(x=lon)
        lat_str  = '{x:.5f}'.format(x=lat)
        elev_str = '{x:.5f}'.format(x=elev)

        #-------------------------------------------------      
        # Get more attributes from a USGS-HADS crosswalk
        #---------------------------------------------------
        # Could also get CWA_ID from:
        # NOAA_RFCs/_New/ba12my15.csv
        # for 9370 basins, but CWA_ID should be WFO_ID.
        # Can use this other file to double-check info.
        #---------------------------------------------------
        if (ADD_INFO):        
            rfc_info = usgs_rfc_dict[ usgs_id ]
            hsa_id   = rfc_info[ 'hsa_id' ]
            goes_id  = rfc_info[ 'goes_id' ]
                 
        #--------------------------------------
        # Write all site data to new TSV file
        #--------------------------------------
        line = ''
        line += new_dict['usgs_id']         + delim
        line += new_dict['nws_loc_id']      + delim
        line += new_dict['upstrm_lid']      + delim 
        line += new_dict['dnstrm_lid']      + delim
        line += new_dict['nwm_reach_id']    + delim
        line += new_dict['rfc_abbrev']      + delim
        line += new_dict['wfo_abbrev']      + delim
        if (ADD_INFO):
            line += hsa_id                      + delim  ####
            line += goes_id                     + delim  ####       
        line += new_dict['usgs_name']       + delim
        line += new_dict['description']     + delim
        line += new_dict['state_code']      + delim
        line += new_dict['county_name']     + delim
        line += lon_str                     + delim
        line += lat_str                     + delim 
        line += elev_str                    + delim
        line += new_dict['updated_time']    + delim
        line += str(new_dict['in_service']) + delim
        line += new_dict['pedts_obs']       + delim 
        line += new_dict['pedts_pred']      + '\n'
        out_tsv_unit.write( line )
        n_sites += 1
        if ((n_sites % 10) == 0):
            print('status: n_sites =', n_sites)
  
    usgs_unit.close()
    out_tsv_unit.close()
    run_time = (time.time() - start_time)
    print('Run time =', run_time, '[secs].')
    print('Wrote info for', n_sites, 'sites.')
    print('Finished.')
    print()
                
#   create_tsv_via_api()
#---------------------------------------------------------------------




