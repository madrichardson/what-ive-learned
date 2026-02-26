# Import necessary packages
import argparse
from dateutil import rrule
from dateutil.parser import parse
from datetime import timezone, datetime
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
CDL_IN_FILE = 'delta_nasa_9k.cdl'
TEMP_OUT_FILE = 'tempoutfiledelta.nc'
NCO_DIR = '/Users/madisonrichardson/miniforge3/bin/'

# Create functions
# Create function to make NetCDF file from CDL file


def make_ncfile_from_cdl_dr(path_cdl, path_ncfile, nco_path):
    """
    Make a NetCDF file from a CDL (Common Data Language) file
    template using the ncgen tool.

    Args:
        path_cdl (str): The path to the input CDL file that defines
                        the structure of the NetCDF file.
        path_ncfile (str): The path to the output NetCDF file that
                            will be created by ncgen.
        nco_path (str): The path to the directory containing the
                        NCO (NetCDF Operators) tools such as ncgen
                        to make the NetCDF file from the CDL file.
    """
    myCmd = ' '.join([os.path.join(nco_path, 'ncgen'),
                      '-o', path_ncfile,
                      path_cdl])

    print("Generated NetCDF template",
          subprocess.call(myCmd, shell=True))
    # prints zero if myCmd is successful


# Create function to find the relative difference (delta) between
# NOAA20 and MODIS NetPP datasets


def calculate_relative_diff(minuend_data, subtrahend_data):
    """
    Calculate the relative difference (delta) between the two datasets.

    Args:
        minuend_data (ma.MaskedArray): The dataset where values will be subtracted.
        subtrahend_data (ma.MaskedArray): The dataset that will be subtracted
                                          from 'minuend_data'.

    Returns:
        ma.MaskedArray: The array of relative difference values
                        where missing data remains masked.
    """

    # Calculate relative difference (delta)
    delta = ma.divide(
        ma.subtract(minuend_data, subtrahend_data),
        subtrahend_data
    )

    return delta


# Create function to extract variable data from NetCDF file


def get_nc_var_data_dr(file1_path, my_var):
    """
    Extract data from specified variable in a NetCDF file.

    Args:
        file1_path (str): The full path to the NetCDF file from where the data is extracted.
        my_var (str): The name of the variable in the NetCDF file to extract.

    Returns:
        netpp_file1 (numpy.ndarray): A 2D array (latitude, longitude)
                                     containing the data for the first time step
                                     of the specified variable.
    """

    file1 = Dataset(file1_path, 'r')

    # Extract NetPP data
    netpp_file1 = file1.variables[my_var][0, :, :]

    # Close the file
    file1.close()

    return netpp_file1


# Run main function


def main():
    """
    Runs main function.
    Raises:
        ValueError: If the dimensions of the VIIRS and MODIS NetPP datasets do not match.
    """

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
                        help="Set to overwrite a NetPP file that exists")

    args = parser.parse_args()

    # Parse the start and end dates
    start_date = parse(args.start)
    end_date = parse(args.end)
    start_date = start_date.replace(day=16)
    end_date = end_date.replace(day=16)
    print(start_date, end_date)

    yr_mo_to_bin = rrule.rrule(
        rrule.MONTHLY,
        dtstart=start_date,
        until=end_date
    )
    yr_mo_to_bin = list(yr_mo_to_bin)

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
    nc_filename = 'delta_3D_matrix_2018to2022.nc'
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
            VIIRS_DIR,
            f"netpp_viirs_noaa20_monthly_{file_date}.nc"
        )
        modis_file_path = os.path.join(
            MODIS_DIR,
            f"productivity_month_modis_{file_date}_9km.nc"
        )

        # Check if both VIIRS and MODIS files exist
        if not os.path.exists(viirs_file_path) or not os.path.exists(modis_file_path):
            print(f"Files for {file_date} are missing! Skipping...")
            continue

        viirs_npp = get_nc_var_data_dr(viirs_file_path, netpp_var)
        modis_npp = get_nc_var_data_dr(modis_file_path, netpp_var)

        if viirs_npp.shape != modis_npp.shape:
            err_msg = ("Dimension mismatch. VIIRS shape: "
                       f"{viirs_npp.shape}, MODIS shape: "
                       f"{modis_npp.shape}"
                       )
            raise ValueError(err_msg)
            print("Skipping...", file_date)
            continue

        delta = calculate_relative_diff(viirs_npp, modis_npp)
        print('Relative difference calculated for', file_date)

        # Ensure the date is centered on the 16th of the month at 00:00:00
        dt = dt.replace(day=16, hour=0, minute=0, second=0, microsecond=0)
        
        # Convert the datetime to a timestamp (seconds since 1970-01-01 00:00:00)
        timestamp = int(dt.replace(tzinfo=timezone.utc).timestamp())

        # Save the timestamp in your NetCDF file
        nc_file['time'][time_index] = timestamp
        nc_file['delta'][time_index, :, :] = delta[:, :]
        nc_file.sync()
        print(file_date, 'data saved to nc file with timestamp:', timestamp)

        time_index += 1

    nc_file.close()


main()
