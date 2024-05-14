import os
from topoflow.utils import rti_files
from topoflow.utils import regrid
from topoflow.utils import rts_files
from topoflow.utils import ncgs_files
from topoflow.utils import model_input

#-------------------------------------------------------
# Create individual .nc files for each forcing variable
#-------------------------------------------------------
# Input file
input_file = '/Users/laurenbolotin/basins/Wolverine/__met/TopoFlow_AORC_Alaska_Forcings_full.nc'

# Get DEM grid info
grid_info = rti_files.read_info( '/Users/laurenbolotin/basins/Wolverine/__topo/Wolverine.rti', SILENT=True )

# Forcings var names (AORC)
var_names = ['Lwnet_tavg','Psurf_f_inst','Qair_f_inst','Rainf_tavg','Swnet_tavg','Tair_f_inst','Wind_f_inst']

# Loop over vars and create one .nc file for each
if not(os.path.exists( os.path.join(os.path.dirname(input_file), 'individual_var_nc_files/') )): 
        os.mkdir( os.path.join(os.path.dirname(input_file), 'individual_var_nc_files/') )

# TODO: Need to make sure the start datetime is not hard coded in regrid.py:create_ncgs_forcings_file(), as it currently is
for var_name in var_names:

    # name the output file
    output_file = os.path.join(os.path.dirname(input_file), 'individual_var_nc_files/TopoFlow_AORC_Alaska_Forcings_'+var_name+'.nc')
    
    # create regridded forcings file
    regrid.create_ncgs_forcings_file(var_name=var_name, nc_file_in=input_file, ncgs_file_out=output_file, grid_info=grid_info)

#-------------------------------------------------------
# Convert .nc files to .rts files
#-------------------------------------------------------
var_type = 'Grid_Sequence'
# Doesn't give full paths:
# input_files = os.listdir('/Users/laurenbolotin/basins/Gulkana/__met/individual_var_nc_files')
# Give full paths:
input_files = []
for root, dirs, files in os.walk(os.path.abspath("/Users/laurenbolotin/basins/Wolverine/__met/individual_var_nc_files/")):
    for file in files:
        nc_forcing_file = os.path.join(root, file)
        input_files.append(nc_forcing_file)

file = input_files[0]

for file in input_files:
    print('Now converting "'+file+'" to .RTS format...')
    file_unit = model_input.open_file(var_type, file)
    nc_obj = ncgs_files.ncgs_file()
    nc_obj.open_file(file)
    var_names_nc = nc_obj.get_var_names()

    # See which of the vars is in the file the loop is currently on, and where that var is so it can be extracted
    for x in var_names:
        if x in var_names_nc:
            var_index = var_names_nc.index(x)

    rts = rts_files.rts_file()
    rts_fn = '/Users/laurenbolotin/basins/Wolverine/__met/'+var_names_nc[var_index]+'_nc_to_rts.rts'
    # var_names = ['Lwnet_tavg','Psurf_f_inst','Qair_f_inst','Rainf_tavg','Swnet_tavg','Tair_f_inst','Wind_f_inst']
    if var_names_nc[var_index] == 'Tair_f_inst':
         rts_var_name = 'T_air'
    if var_names_nc[var_index] == 'Lwnet_tavg':
         rts_var_name = 'Qn_SW'
    if var_names_nc[var_index] == 'Psurf_f_inst':
         rts_var_name = 'p0'
    if var_names_nc[var_index] == 'Qair_f_inst':
         rts_var_name = 'Qair'
    # TODO: add conversion of SH to RH
    if var_names_nc[var_index] == 'Rainf_tavg':
         rts_var_name = 'P'
    if var_names_nc[var_index] == 'Swnet_tavg':
         rts_var_name = 'Qn_SW'
    if var_names_nc[var_index] == 'Wind_f_inst':
         rts_var_name = 'uz'
    rts.open_new_file(file_name = rts_fn, info=grid_info, var_name = 'T_air', dtype='float32')
    # TODO: don't have the var_name hard coded above ^^^^

    # Set how many grids (hours) of data there should be:
    time_info = nc_obj.get_time_info()
    mins = time_info.duration
    hours = mins/60
    n_grids = int(hours+1)

    for time_index in range(n_grids):
        var_data = nc_obj.get_grid(var_name=var_names_nc[var_index], time_index=time_index)
        rts.add_grid( var_data)       
    rts.close_file()

#-------------------------------------------------------
# Convert specific humidity to relative humidity
#-------------------------------------------------------
import numpy as np
from topoflow.utils import met_utils as mu

# Do the needed files exist?
q_air_file  = '/Users/laurenbolotin/basins/Wolverine/__met/Qair_f_inst_nc_to_rts.rts'
T_air_file  = '/Users/laurenbolotin/basins/Wolverine/__met/Tair_f_inst_nc_to_rts.rts' 
P_surf_file = '/Users/laurenbolotin/basins/Wolverine/__met/Psurf_f_inst_nc_to_rts.rts' 
HAVE_Q_AIR  = os.path.exists( q_air_file )
HAVE_T_AIR  = os.path.exists( T_air_file )
HAVE_P_SURF = os.path.exists( P_surf_file ) 

if (HAVE_Q_AIR and HAVE_T_AIR and HAVE_P_SURF):
    rts_q_air  = rts_files.rts_file()
    rts_T_air  = rts_files.rts_file()
    rts_P_surf = rts_files.rts_file()

    rts_q_air.open_file( q_air_file )
    rts_T_air.open_file( T_air_file )
    rts_P_surf.open_file( P_surf_file )

    RH_file = '/Users/laurenbolotin/basins/Wolverine/__met/aorc_rh.rts'
    rts_RH  = rts_files.rts_file()   # For new RTS file
    OK = rts_RH.open_new_file( RH_file, info=grid_info, var_name='RH' )
    method = 'BRUTSAERT'  # (or 'SATTERLUND', or 'BOLTON')
    n_grids = rts_q_air.number_of_grids()

    print('Creating RTS Grid Stack for RH...')
    for time_index in range(n_grids):
        q_air  = rts_q_air.read_grid(  time_index, dtype='float32' )
        T_air  = rts_T_air.read_grid(  time_index, dtype='float32' )
        P_surf = rts_P_surf.read_grid( time_index, dtype='float32' )

        # Compute RH via met_utils.py. RH in [0,1]
        RH = mu.relative_humidity( q_air, T_air, P_surf, method=method)
        RH = np.float32( RH )
        rts_RH.add_grid( RH )

    rts_q_air.close_file()
    rts_T_air.close_file()
    rts_P_surf.close_file()
    rts_RH.close_file()
    print('Finished creating RTS Grid Stack for RH.')
    print()
else:
    print('ERROR:  Could not find the required grids for:')
    print('        Q_air, T_air and P_surf.')
    print()