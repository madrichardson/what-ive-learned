# regridding a batch of files

import os
import subprocess
import shutil
from datetime import datetime

# set software locations
cwutl_dir = '"/Applications/CoastWatch Utilities/bin"'

# set up local directories
root_dir = rootdir = os.path.expanduser("~/netpp")
data_dir = os.path.join(root_dir, 'data', 'noaa20')
sst_in_dir = os.path.join(data_dir, 'sst_2k')
sst_out_dir = os.path.join(data_dir, 'sst')
work_dir = os.path.join(root_dir, 'bin')
res_dir = os.path.join(root_dir, 'resources')

# ensure directories exist
os.makedirs(data_dir, exist_ok=True)
os.makedirs(sst_in_dir, exist_ok=True)
os.makedirs(sst_out_dir, exist_ok=True)
os.makedirs(work_dir, exist_ok=True)
os.makedirs(res_dir, exist_ok=True)

# name the output file for cwregister2
cwutil_ofile_hdf = 'nasa_4k.hdf'
# name the output file for cwexport
cwutil_ofile_nc = 'nasa_4k.nc'

# template for output files
master_file = 'nasa_4k_template.nc'

# list of 2km files to be regridded
ifiles = [
    '20230601120000-STAR-L3S_GHRSST-SSTsubskin-LEO_Daily-ACSPO_V2.81-v02.0-fv01.0_sstcelsius.nc',
    '20230602120000-STAR-L3S_GHRSST-SSTsubskin-LEO_Daily-ACSPO_V2.81-v02.0-fv01.0_sstcelsius.nc',
    '20230603120000-STAR-L3S_GHRSST-SSTsubskin-LEO_Daily-ACSPO_V2.81-v02.0-fv01.0_sstcelsius.nc',
    '20230604120000-STAR-L3S_GHRSST-SSTsubskin-LEO_Daily-ACSPO_V2.81-v02.0-fv01.0_sstcelsius.nc',
    '20230605120000-STAR-L3S_GHRSST-SSTsubskin-LEO_Daily-ACSPO_V2.81-v02.0-fv01.0_sstcelsius.nc',
    '20230606120000-STAR-L3S_GHRSST-SSTsubskin-LEO_Daily-ACSPO_V2.81-v02.0-fv01.0_sstcelsius.nc',
    '20230607120000-STAR-L3S_GHRSST-SSTsubskin-LEO_Daily-ACSPO_V2.81-v02.0-fv01.0_sstcelsius.nc',
    '20230608120000-STAR-L3S_GHRSST-SSTsubskin-LEO_Daily-ACSPO_V2.81-v02.0-fv01.0_sstcelsius.nc',
    '20230609120000-STAR-L3S_GHRSST-SSTsubskin-LEO_Daily-ACSPO_V2.81-v02.0-fv01.0_sstcelsius.nc',
    '20230610120000-STAR-L3S_GHRSST-SSTsubskin-LEO_Daily-ACSPO_V2.81-v02.0-fv01.0_sstcelsius.nc',
]

for ifile in ifiles:
    # regrid with cwregister2
    myCmd = ' '.join([os.path.join(cwutl_dir, 'cwregister2'),
                  '--clobber',
                  '--match=' + 'sea_surface_temperature',
                  '--master=' + os.path.join(res_dir,
                                             master_file),
                  os.path.join(sst_in_dir, ifile),
                  os.path.join(work_dir, cwutil_ofile_hdf)
                  ])
    print(myCmd)
    print('Regrid NASA file',
        subprocess.call(myCmd, shell=True))
    
    # add lat/lon with cwangles
    cw_cmd = ' '.join([os.path.join(cwutl_dir, 'cwangles'),
                   '--location --float --units=deg',
                   os.path.join(work_dir, cwutil_ofile_hdf)
                   ])
    print('Add lat/lon', subprocess.call(cw_cmd, shell=True))

    # convert hdf to netcdf with cw utilities
    myCmd = ' '.join([os.path.join(cwutl_dir, 'cwexport'), '-v ',
                  os.path.join(work_dir, cwutil_ofile_hdf),
                  os.path.join(work_dir, cwutil_ofile_nc)
                  ])
    print('convert to netCDF',
        subprocess.call(myCmd, shell=True))
    
    # pull the date from input file
    date_str = ifile.split('.')[0][:8]
    date_obj = datetime.strptime(date_str, '%Y%m%d')

    # create output file name using date object
    ofile_template = "NOAA_LEO.{date}.L3.DAY.SST.sst.4km.nc"
    ofile_name = ofile_template.format(date=date_obj.strftime('%Y%m%d'))

    # rename and move regridded file to sst_out_dir
    source_fpath = os.path.join(work_dir, cwutil_ofile_nc)
    new_fpath = os.path.join(sst_out_dir, ofile_name)

    shutil.move(source_fpath, new_fpath)

    print(f"File has been renamed to {ofile_name} and moved to {sst_out_dir}")