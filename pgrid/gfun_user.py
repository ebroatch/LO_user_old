"""
User-specific code for pgrid.

You would edit the information to reflect whatever grid you are working on.
"""
import numpy as np
from lo_tools import zfun

import sys
from pathlib import Path
pth = Path(__file__).absolute().parent.parent.parent / 'LO' / 'pgrid'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import gfun_utility as gfu
import gfun

# This is the name of the grid that you are working on.
gridname = 'so0'

# default s-coordinate info (could override below)
s_dict = {'THETA_S': 4, 'THETA_B': 2, 'TCLINE': 10, 'N': 30,
        'VTRANSFORM': 2, 'VSTRETCHING': 4}

if gridname in ['ai0','hc0', 'sal0', 'so0']:
    # These are the gridname and tag to feed to use when creating the Ldir paths.
    # They are used for accessing the river tracks, which may be developed for one
    # grid but reused in others.
    base_gridname = 'cas6'
    base_tag = 'v3'

def make_initial_info(gridname=gridname):
    # Add an elif section for your grid.

    if gridname == 'sal0':
        # A Salish Sea grid, used as an example.
        dch = gfun.default_choices()
        aa = [-124, -122, 47, 49]
        res = 600 # target resolution (m)
        Lon_vec, Lat_vec = gfu.simple_grid(aa, res)
        dch['nudging_edges'] = ['north', 'west']
        # Make the rho grid.
        lon, lat = np.meshgrid(Lon_vec, Lat_vec)
        # initialize the bathymetry array
        z = np.nan * lon
        # add bathymetry automatically from files
        for t_fn in dch['t_list']:
            print('\nOPENING BATHY FILE: ' + t_fn.name)
            tlon_vec, tlat_vec, tz = gfu.load_bathy_nc(t_fn)
            tlon, tlat = np.meshgrid(tlon_vec, tlat_vec)
            z_part = zfun.interp2(lon, lat, tlon, tlat, tz)
            # put good values of z_part in z
            z[~np.isnan(z_part)] = z_part[~np.isnan(z_part)]
        if dch['use_z_offset']:
            z = z + dch['z_offset']
            
    elif gridname == 'hc0':
        dch = gfun.default_choices()
        aa = [-123.2, -122.537, 47.3, 47.9]
        res = 100 # target resolution (m)
        Lon_vec, Lat_vec = gfu.simple_grid(aa, res)
        dch['t_list'] = [dch['t_dir'] / 'psdem' / 'PS_27m.nc']
        dch['nudging_edges'] = ['north']
        dch['nudging_days'] = (0.1, 1.0)
        
        # Make the rho grid.
        lon, lat = np.meshgrid(Lon_vec, Lat_vec)
        # initialize the bathymetry array
        z = np.nan * lon
        # add bathymetry automatically from files
        for t_fn in dch['t_list']:
            print('\nOPENING BATHY FILE: ' + t_fn.name)
            tlon_vec, tlat_vec, tz = gfu.load_bathy_nc(t_fn)
            tlon, tlat = np.meshgrid(tlon_vec, tlat_vec)
            z_part = zfun.interp2(lon, lat, tlon, tlat, tz)
            # put good values of z_part in z
            z[~np.isnan(z_part)] = z_part[~np.isnan(z_part)]
        if dch['use_z_offset']:
            z = z + dch['z_offset']
            
    elif gridname == 'ai0':
        dch = gfun.default_choices()
        aa = [-122.82, -122.36, 47.758, 48.18]
        res = 100 # target resolution (m)
        Lon_vec, Lat_vec = gfu.simple_grid(aa, res)
        dch['t_list'] = [dch['t_dir'] / 'psdem_10m' / 'PS_30m.nc']
        dch['nudging_edges'] = ['north', 'south', 'east', 'west']
        dch['nudging_days'] = (0.1, 1.0)
        
        # Make the rho grid.
        lon, lat = np.meshgrid(Lon_vec, Lat_vec)
        # initialize the bathymetry array
        z = np.nan * lon
        # add bathymetry automatically from files
        for t_fn in dch['t_list']:
            print('\nOPENING BATHY FILE: ' + t_fn.name)
            tlon_vec, tlat_vec, tz = gfu.load_bathy_nc(t_fn)
            tlon, tlat = np.meshgrid(tlon_vec, tlat_vec)
            z_part = zfun.interp2(lon, lat, tlon, tlat, tz)
            # put good values of z_part in z
            z[~np.isnan(z_part)] = z_part[~np.isnan(z_part)]
        if dch['use_z_offset']:
            z = z + dch['z_offset']
            
    elif gridname == 'so0':
        # South Sound
        dch = gfun.default_choices()
        
        dch['z_offset'] = -1.3 # NAVD88 is 1.3 m below MSL at Seattle
        
        #dch['min_depth'] = 0.4 # requires WET_DRY in ROMS
        
        dch['excluded_rivers'] = ['skokomish']
        
        #dch['rx0max'] = 0.12 # default is 0.15 (smaller = smoother)
        
        aa = [-123.13, -122.76, 47, 47.42]
        res = 50 # target resolution (m)
        Lon_vec, Lat_vec = gfu.simple_grid(aa, res)
        dch['t_list'] = [dch['t_dir'] / 'srtm15' / 'topo15.nc',
                    dch['t_dir'] / 'psdem_10m' / 'PS_30m.nc']
        dch['nudging_edges'] = ['east']
        dch['nudging_days'] = (0.1, 1.0)
        
        # Make the rho grid.
        lon, lat = np.meshgrid(Lon_vec, Lat_vec)
        # initialize the bathymetry array
        z = np.nan * lon
        # add bathymetry automatically from files
        for t_fn in dch['t_list']:
            print('\nOPENING BATHY FILE: ' + t_fn.name)
            tlon_vec, tlat_vec, tz = gfu.load_bathy_nc(t_fn)
            tlon, tlat = np.meshgrid(tlon_vec, tlat_vec)
            z_part = zfun.interp2(lon, lat, tlon, tlat, tz)
            # put good values of z_part in z
            z[~np.isnan(z_part)] = z_part[~np.isnan(z_part)]
        if dch['use_z_offset']:
            z = z + dch['z_offset']
    
        
        
    else:
        print('Error from make_initial_info: unsupported gridname')
        return
        
    return lon, lat, z, dch
    

