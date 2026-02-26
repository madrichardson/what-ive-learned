# Finding Unbiased Unbiased Relative Difference (Delta)
# Between 8-DAY MODIS and VIIRS (NOAA21)

# Import packages
import os
import subprocess
import numpy.ma as ma
from netCDF4 import Dataset
import argparse
from dateutil.parser import parse
from datetime import timedelta, datetime

# Create global variables
ROOT_DIR = ("/Users/madisonrichardson/netpp")
VIIRS_DIR_T = os.path.join(ROOT_DIR, "data/{}/netpp/8_day_netpp")
SNPP_DIR = os.path.join(ROOT_DIR, "data/snpp_viirs/8_day/netpp")
MODIS_DIR = os.path.join(ROOT_DIR, "data/aqua_modis/8_day/netpp")
WORK_DIR = os.path.join(ROOT_DIR, 'work')
BIN_DIR = os.path.join(ROOT_DIR, 'bin')
RES_DIR = os.path.join(ROOT_DIR, 'resources')
NC_OUT_DIR_T = os.path.join(ROOT_DIR, "data/{}/matrix/delta")
CDL_IN_FILE = 'delta_nasa_8day_9k.cdl'
TEMP_OUT_FILE = 'tempoutfiledelta8day.nc'
NCO_DIR = '/Users/madisonrichardson/miniforge3/bin/'

# Create function to generate a NetCDF file from CDL file


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


# Create function to find the unbiased relative difference (psi)
# between NOAA21 and MODIS NetPP datasets


def calculate_relative_diff(minuend_data, subtrahend_data):
    """
    Calculate the relative difference (delta) between the two
    datasets.

    Nerd Notation
    minuend is math-speak for the value substracted from.

    subtrahend is math-speak for the value that is substracted.

    Args:
        minuend_data (ma.MaskedArray): The dataset where values will be
                                    subtracted.
        subtrahend_data (ma.MaskedArray): The dataset that will be subtracted
                                        from 'minuend_data'.

    Returns:
        ma.MaskedArray: The array of relative difference values
                        where missing data remains masked.
    """

    # Calculate relative difference (delta)
    delta = ma.divide(
        ma.subtract(minuend_data, subtrahend_data), subtrahend_data
    )

    return delta


# Create function to populate a NetCDF file
# with a 3D matrix of the Delta values


def get_nc_var_data_dr(file1_path, my_var):
    """
    Extract data from specified variable in a NetCDF file.

    Args:
        file1_path (str): The full path to the NetCDF file from
                            where the data is extracted.
        my_var (str): The name of the variable in the NetCDF file
                        to extract (e.g., 'productivity).

    Returns:
        netpp_file1 (numpy.ndarray): A 2D array (latitude, longitude)
        containing the data for the first time step of the
        specified variable.
    """

    file1 = Dataset(file1_path, 'r')

    # Extract NetPP data
    netpp_file1 = file1.variables[my_var][0, :, :]

    # Close the VIIRS and MODIS files
    file1.close()

    return netpp_file1


# Run main function


def main():
    """
    Runs main function.

    Raises:
        ValueError: If the dimensions of the VIIRS and MODIS
                    NetPP datasets do not match.
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
                        help="End date in YYYY-MM format")
    parser.add_argument("-c",
                        "--comparison_sensor",
                        type=str,
                        required=True,
                        choices=['modis', 'snpp_viirs'],
                        help="End date in YYYY-MM format")
    parser.add_argument("-o", "--overwrite", required=False,
                        action='store_true',
                        help="set to overwrite a netpp file that exists")

    args = parser.parse_args()

    # Parse the start and end dates
    start_date = parse(args.start)
    end_date = parse(args.end)

    print(start_date, end_date)

    # Define 8-day intervals
    day_step = timedelta(days=8)
    dates_to_bin = []
    current_date = start_date

    while current_date <= end_date:
        next_date = current_date + day_step

        # Check if the interval is Dec 27-31 and skip
        if (
            current_date.year == 2023
            and current_date.month == 12
            and current_date.day == 27
        ):
            print(f"Skipping interval from {current_date} to Dec 31")
            current_date = datetime(2024, 1, 1)
            continue

        # Add regular 8-day interval adjusting if it goes past Dec 31
        if next_date.year > current_date.year and current_date.month == 12:
            # Adjust to end on Dec 26 if it goes past the year
            end_dt = datetime(current_date.year, 12, 26)
        else:
            end_dt = next_date - timedelta(days=1)

        # Add the interval (start_date, end_date) to the list
        dates_to_bin.append((current_date, end_dt))

        # Move to the next interval
        current_date = next_date

    # Create dynamic directories and verify
    VIIRS_DIR = VIIRS_DIR_T.format(args.sensor)
    NC_OUT_DIR = NC_OUT_DIR_T.format(args.sensor)

    DIR_LIST = [ROOT_DIR,
                WORK_DIR,
                RES_DIR,
                VIIRS_DIR,
                SNPP_DIR,
                MODIS_DIR,
                NC_OUT_DIR,
                NCO_DIR
                ]

    for dr in DIR_LIST:
        os.makedirs(dr, exist_ok=True)
    print(len(DIR_LIST), 'directories validated')

    # Define the final output file path
    nc_filename = f'delta_8DAY_3D_matrix_{args.comparison_sensor}_2023to2024.nc'
    nc_file_path = os.path.join(NC_OUT_DIR, nc_filename)

    # Add logic to not overwrite existing files
    if os.path.isfile(nc_file_path):
        if not args.overwrite:
            print(f'{nc_filename} already exists for {args.sensor}')
        else:
            print(f'Overwriting {nc_filename} for {args.sensor}')

    # Generate the initial NetCDF template from CDL for all years
    make_ncfile_from_cdl_dr(os.path.join(RES_DIR, CDL_IN_FILE),
                            os.path.join(WORK_DIR, nc_file_path),
                            NCO_DIR
                            )

    # Name of variable in NetPP file
    netpp_var = 'productivity'

    nc_file = Dataset(nc_file_path, 'a', format='NETCDF4')

    time_index = 0
    for start_dt, end_dt in dates_to_bin:
        midpoint_dt = start_dt + timedelta(days=4)

        file_date = '{0:%Y%m%d}_{1:%Y%m%d}'.format(start_dt, end_dt)
        year = start_dt.year
        print("Processing", file_date)

        viirs_file_path = os.path.join(
            VIIRS_DIR,
            str(year),
            f"netpp_viirs_noaa21_8day_{file_date}.nc"
        )

        # Determine comparison file path based on selected comparison_sensor
        if args.comparison_sensor == 'modis':
            comparison_dir = MODIS_DIR
            comparison_file_path = os.path.join(
                comparison_dir,
                str(year),
                f"netpp_aqua_modis_8day_{file_date}.nc"
            )
        else:
            comparison_dir = SNPP_DIR
            comparison_file_path = os.path.join(
                comparison_dir,
                str(year),
                f"netpp_snpp_viirs_8day_{file_date}.nc"
            )

        if (
            not os.path.exists(viirs_file_path)
            or not os.path.exists(comparison_file_path)
        ):
            print(f"Files for {file_date} are missing! Skipping...")
            continue

        viirs_npp = get_nc_var_data_dr(viirs_file_path, netpp_var)
        comparison_npp = get_nc_var_data_dr(comparison_file_path, netpp_var)

        if viirs_npp.shape != comparison_npp.shape:
            err_msg = ("Dimension mismatch. Viirs shape: "
                       f"{viirs_npp.shape}, Legacy shape: "
                       f"{comparison_npp.shape}"
                       )
            raise ValueError(err_msg)
            print("Skipping...", file_date)
            continue

        delta = calculate_relative_diff(viirs_npp, comparison_npp)
        print('relative difference calculated for', file_date)

        time_stamp = midpoint_dt.replace(hour=0, minute=0, second=0).timestamp()
        nc_file['time'][time_index] = time_stamp
        nc_file['delta'][time_index, :, :] = delta[:, :]
        nc_file.title = (
            f"Pixel by pixel Delta, VIIRS {args.sensor.upper()} "
            f"minus {args.comparison_sensor.upper()}"
        )
        nc_file.summary = (
            f"The relative difference (Delta) between primary productivity "
            f"(NetPP) calculated using VIIRS {args.sensor.upper()} data and "
            f"{args.comparison_sensor.upper()} data. For each pixel, Delta is "
            f"calculated as netPP(VIIRS {args.sensor.upper()}) - "
            f"netPP({args.comparison_sensor.upper()}) divided by "
            f"netPP({args.comparison_sensor.upper()}). Primary productivity "
            f"was calculated as described by Behrenfeld and Falkowski 1997. "
            f"The data is a 1-year mean from Apr 2023 to May 2024 and "
            f"at 9km resolution. Input data for primary productivity were "
            f"obtained from NASA and included chlorophyll_a, sea surface "
            f"temperature, and photosynthetically active radiation from "
            f"either {args.comparison_sensor.upper()} or "
            f"{args.sensor.upper()}."
        )
        nc_file.sync()
        print(file_date, 'data saved to nc file')

        time_index += 1

    nc_file.close()

    # Compress file with nccopy
    compressed_nc_file_path = os.path.join(
        NC_OUT_DIR,
        'compressed_' + nc_filename
    )

    compressed_cmd = ' '.join(
        [
            'nccopy',
            '-d4',
            nc_file_path,
            compressed_nc_file_path
        ]
    )
    print('Compress ofile', subprocess.call(compressed_cmd, shell=True))

    archive_cmd = ' '.join(['mv', compressed_nc_file_path, nc_file_path])
    print('Archive ofile', subprocess.call(archive_cmd, shell=True))

    print(f"NetCDF file '{nc_filename}' archived at {nc_file_path}")


if __name__ == '__main__':
    main()
