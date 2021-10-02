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
roms_out2 = parent / 'BLANK'
roms_out3 = parent / 'BLANK'
roms_out4 = parent / 'BLANK'

# these are for mox and klone, other hyak mackines
remote_user = 'BLANK'
remote_machine = 'BLANK'
remote_dir0 = 'BLANK'
local_user = 'BLANK'

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
    remote_user = 'parker'
    remote_machine = 'perigee.ocean.washington.edu'
    remote_dir0 = '/data1/parker'
    local_user = 'pmacc'

elif (str(HOME) == '/mmfs1/home/pmacc'):
    lo_env = 'pm_klone'
    remote_user = 'parker'
    remote_machine = 'apogee.ocean.washington.edu'
    remote_dir0 = '/dat1/parker'
    local_user = 'pmacc'
  
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
#
Ldir0['remote_user'] = remote_user
Ldir0['remote_machine'] = remote_machine
Ldir0['remote_dir0'] = remote_dir0
Ldir0['local_user'] = local_user


