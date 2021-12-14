"""
This is the one place where you set the path structure of the LO code.
The info is stored in the dict Ldir.

All paths are pathlib.Path objects

This program is meant to be loaded as a module by Lfun which then adds more
entries to the Ldir dict based on which model run you are working on.

Users should create LO_user/get_lo_info.py, edit as needed, and make it into
their own GitHub repo.

"""
import os
from pathlib import Path

# defaults that should work on all machines
parent = Path(__file__).absolute().parent.parent
LO = parent / 'LO'
LOo = parent / 'LO_output'
LOu = parent / 'LO_user'
data = parent / 'LO_data'

# This is where the ROMS source code, makefiles, and executables are
roms_code = parent / 'LiveOcean_roms'

# These are places where the ROMS history files are kept
roms_out = parent / 'LO_roms'
roms_out1 = parent / 'BLANK' # placeholder
roms_out2 = parent / 'BLANK' # placeholder
roms_out3 = parent / 'BLANK' # placeholder
roms_out4 = parent / 'BLANK' # placeholder

# default for linux machines
which_matlab = '/usr/local/bin/matlab'

HOME = Path.home()
try:
    HOSTNAME = os.environ['HOSTNAME']
except KeyError:
    HOSTNAME = 'BLANK'

if str(HOME) == '/Users/pm8':
    lo_env = 'pm_mac'
    which_matlab = '/Applications/MATLAB_R2020a.app/bin/matlab'
    roms_out2 = parent / 'LiveOcean_roms' / 'output'    # LiveOcean past I
    roms_out3 = parent / 'LiveOcean_roms' / 'output'    # LiveOcean past II (same as I)

elif (str(HOME) == '/home/parker') & ('boiler' in HOSTNAME):
    lo_env = 'pm_boiler'
    roms_out1 = Path('/data1/parker/LiveOcean_roms/output')     # LiveOcean current
    roms_out2 = Path('/pgdat1/parker/LiveOcean_roms/output')    # LiveOcean past I
    roms_out3 = Path('/pgdat2/parker/LiveOcean_roms/output')    # LiveOcean past II

elif (str(HOME) == '/home/parker') & ('perigee' in HOSTNAME):
    lo_env = 'pm_perigee'
    roms_out1 = Path('/boildat1/parker/LiveOcean_roms/output')  # LiveOcean current
    roms_out2 = Path('/data1/parker/LiveOcean_roms/output')     # LiveOcean past I
    roms_out3 = Path('/data2/parker/LiveOcean_roms/output')     # LiveOcean past II
    roms_out4 = Path('/boildat1/parker/LO_roms')                # LO boiler
    
elif (str(HOME) == '/home/parker') & ('apogee' in HOSTNAME):
    lo_env = 'pm_apogee'
    roms_out1 = Path('/boildat/parker/LiveOcean_roms/output')  # LiveOcean current
    roms_out2 = Path('/pgdat1/parker/LiveOcean_roms/output')     # LiveOcean past I
    roms_out3 = Path('/pgdat2/parker/LiveOcean_roms/output')     # LiveOcean past II
    roms_out4 = Path('/boildat/parker/LO_roms')                # LO boiler
    
elif (str(HOME) == '/usr/lusers/pmacc'):
    lo_env = 'pm_mox'

elif (str(HOME) == '/mmfs1/home/pmacc'):
    lo_env = 'pm_klone'

if str(HOME) == '/Users/erinbroatch':
    lo_env = 'eb_mac'
    which_matlab = '/Applications/MATLAB_R2020b.app/toolbox/matlab/general/matlab.m' #not sure about this
    #roms_out2 = parent / 'LiveOcean_roms' / 'output'    # LiveOcean past I
    #roms_out3 = parent / 'LiveOcean_roms' / 'output'    # LiveOcean past II (same as I)

elif (str(HOME) == '/home/ebroatch') & ('boiler' in HOSTNAME):
    lo_env = 'eb_boiler'
    roms_out1 = Path('/data1/parker/LiveOcean_roms/output')     # LiveOcean current
    roms_out2 = Path('/pgdat1/parker/LiveOcean_roms/output')    # LiveOcean past I
    roms_out3 = Path('/pgdat2/parker/LiveOcean_roms/output')    # LiveOcean past II

elif (str(HOME) == '/home/ebroatch') & ('perigee' in HOSTNAME):
    lo_env = 'eb_perigee'
    roms_out1 = Path('/boildat1/parker/LiveOcean_roms/output')  # LiveOcean current
    roms_out2 = Path('/data1/parker/LiveOcean_roms/output')     # LiveOcean past I
    roms_out3 = Path('/data2/parker/LiveOcean_roms/output')     # LiveOcean past II
    roms_out4 = Path('/boildat1/parker/LO_roms')                # LO boiler
  
Ldir0 = dict()
Ldir0['lo_env'] = lo_env
Ldir0['parent'] = parent
Ldir0['LO'] = LO
Ldir0['LOo'] = LOo
Ldir0['LOu'] = LOu
Ldir0['data'] = data
Ldir0['roms_code'] = roms_code
Ldir0['roms_out'] = roms_out
Ldir0['roms_out1'] = roms_out1
Ldir0['roms_out2'] = roms_out2
Ldir0['roms_out3'] = roms_out3
Ldir0['roms_out4'] = roms_out4
Ldir0['which_matlab'] = which_matlab

