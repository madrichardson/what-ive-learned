# Import necessary packages
import argparse
from dateutil import rrule
from dateutil.parser import parse
from datetime import timezone
from netCDF4 import Dataset
import numpy.ma as ma
import os
import pandas as pd
import subprocess
import xarray as xr


# Create global variables
ROOT_DIR = ("/Users/madisonrichardson/netpp")
VIIRS_DIR_T = os.path.join(ROOT_DIR, "data/{}/monthly_noaa20")
MODIS_DIR = os.path.join(ROOT_DIR, "data/monthly_modis")
WORK_DIR = os.path.join(ROOT_DIR, 'work')
BIN_DIR = os.path.join(ROOT_DIR, 'bin')
RES_DIR = os.path.join(ROOT_DIR, 'resources')
NC_OUT_DIR_T = os.path.join(ROOT_DIR, "data/{}/matrix")
CDL_IN_FILE = 'psi_nasa_9k_mr_dr.cdl'
TEMP_OUT_FILE = 'tempoutfile.nc'
NCO_DIR = '/Users/madisonrichardson/miniforge3/bin/'


# Create function to make NetCDF file from CDL file
def make_ncfile_from_cdl_dr(path_cdl, path_ncfile, nco_path):
    """
    Make a NetCDF file from a CDL file using ncgen. The function
    also prints a message to indicate if the subprocess call was
    successful (prints 0 for success).
    """
    myCmd = ' '.join([os.path.join(nco_path, 'ncgen'),
                      '-o', path_ncfile,
                      path_cdl])

    print("Generated NetCDF template",
          subprocess.call(myCmd, shell=True))


# Function to calculate unbiased relative difference (psi) using NetPP
def calculate_psi_dr(minuend_data, subtrahend_data):
    """
    Calculate the unbiased relative difference (psi) between the two datasets.
    """
    # Calculate the average of input datasets
    avg = ma.divide(ma.add(minuend_data, subtrahend_data), 2)

    # Calculate unbiased relative difference (psi)
    psi_values = ma.divide(ma.subtract(minuend_data, subtrahend_data), avg)

    return psi_values


# Function to extract variable data from NetCDF file
def get_nc_var_data_dr(file1_path, my_var):
    """Extract the data from the specified NetCDF variable"""
    file1 = Dataset(file1_path, 'r')
    netpp_file1 = file1.variables[my_var][0, :, :]
    file1.close()
    return netpp_file1


# Main function
def main():
    """Main function to process NetCDF files and calculate psi values"""

    # Set up argparse
    doc_formatter = argparse.RawDescriptionHelpFormatter
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=doc_formatter
        )

    parser.add_argument("-a", "--start",
                        type=str,
                        required=True,
                        help="Start date in YYYY-MM format")
    parser.add_argument("-z",
                        "--end",
                        type=str,
                        required=True,
                        help="End date in YYYY-MM format")
    parser.add_argument("-s",
                        "--sensor",
                        type=str,
                        required=True,
                        choices=['noaa20', 'noaa21'],
                        help="Sensor to use (noaa20 or noaa21)")
    parser.add_argument("-o", "--overwrite", required=False,
                        action='store_true',
                        help="set to overwrite a netpp file that exists")

    args = parser.parse_args()

    # Parse the start and end dates and center them on the 16th day at 00:00:00
    start_date = parse(args.start)
    end_date = parse(args.end)
    start_date = start_date.replace(day=16, hour=0, minute=0, second=0, microsecond=0)
    end_date = end_date.replace(day=16, hour=0, minute=0, second=0, microsecond=0)
    print(start_date, end_date)

    yr_mo_to_bin = rrule.rrule(
        rrule.MONTHLY,
        dtstart=start_date,
        until=end_date
    )
    yr_mo_to_bin = list(yr_mo_to_bin)
    # yr_mo_to_bin contains list of monthly date objects
    # e.g. [datetime.datetime(2018, 1, 16, 0, 0),
    # datetime.datetime(2018, 2, 16, 0, 0),...]

    # Create dynamic directories and verify
    VIIRS_DIR = VIIRS_DIR_T.format(args.sensor)
    NC_OUT_DIR = NC_OUT_DIR_T.format(args.sensor)

    DIR_LIST = [ROOT_DIR,
                WORK_DIR,
                RES_DIR,
                VIIRS_DIR,
                MODIS_DIR,
                NC_OUT_DIR,
                NCO_DIR
                ]

    for dr in DIR_LIST:
        os.makedirs(dr, exist_ok=True)
    print(len(DIR_LIST), 'directories validated')

    # Define the final output file path
    nc_filename = 'psi_3D_matrix_2018to2022.nc'
    nc_file_path = os.path.join(NC_OUT_DIR, nc_filename)

    # Generate the initial NetCDF template from CDL for all years
    make_ncfile_from_cdl_dr(os.path.join(RES_DIR, CDL_IN_FILE),
                            os.path.join(WORK_DIR, nc_file_path),
                            NCO_DIR
                            )

    # Name of variable in NetPP file
    netpp_var = 'productivity'

    nc_file = Dataset(nc_file_path, 'a', format='NETCDF4')

    time_index = 0
    for dt in yr_mo_to_bin:
        file_date = '{0:%Y%m}'.format(dt)
        print("Processing", file_date)

        viirs_file_path = os.path.join(
            VIIRS_DIR, f"netpp_viirs_noaa20_monthly_{file_date}.nc"
        )
        modis_file_path = os.path.join(
            MODIS_DIR, f"productivity_month_modis_{file_date}_9km.nc"
        )

        # Check if both VIIRS and MODIS files exist
        if not os.path.exists(viirs_file_path) or not os.path.exists(modis_file_path):
            print(f"Files for {file_date} are missing! Skipping...")
            continue

        # Extract data from VIIRS and MODIS NetCDF files
        viirs_npp = get_nc_var_data_dr(viirs_file_path, netpp_var)
        modis_npp = get_nc_var_data_dr(modis_file_path, netpp_var)

        # Check if the dimensions match
        if viirs_npp.shape != modis_npp.shape:
            err_msg = ("Dimension mismatch. VIIRS shape: "
                       f"{viirs_npp.shape}, MODIS shape: "
                       f"{modis_npp.shape}")
            raise ValueError(err_msg)
            print("Skipping...", file_date)
            continue

        # Calculate psi
        psi = calculate_psi_dr(viirs_npp, modis_npp)
        print('psi calculated for', file_date)

        # Ensure the date is set to the 16th day at 00:00:00 UTC
        dt = dt.replace(day=16, hour=0, minute=0, second=0, microsecond=0, tzinfo=timezone.utc)

        # Get the timestamp (seconds since 1970-01-01 00:00:00 UTC)
        timestamp = int(dt.timestamp())

        # Save the timestamp and psi values in the NetCDF file
        nc_file['time'][time_index] = timestamp
        nc_file['psi'][time_index, :, :] = psi[:, :]
        nc_file.sync()
        print(file_date, 'data saved to nc file with timestamp:', timestamp)

        time_index += 1

    nc_file.close()


# Call the main function
main()

# Find monthly means, minimums, and maximums of psi using xarray
matrix_file = '/Users/madisonrichardson/netpp/data/noaa20/matrix/psi_3D_matrix_2018to2022.nc'

ds = xr.open_dataset(matrix_file)
da = ds.psi

psi_mean3 = da.mean(dim=["latitude", "longitude"], skipna=True)
psi_max = da.max(dim=["latitude", "longitude"], skipna=True)
psi_min = da.min(dim=["latitude", "longitude"], skipna=True)

my_date = pd.to_datetime(da.time.values)

for i in range(0, my_date.size):
    print('For', '{0:%B %Y}'.format(my_date[i]))
    print(
        f"Mean psi {psi_mean3[0]:.2f}",
        f"Min psi {psi_min[0]:.2f}",
        f"Max psi {psi_max[0]:.2f}"
    )