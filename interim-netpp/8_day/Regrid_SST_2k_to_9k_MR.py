"""
    Regridding NOAA ACSPO Daily Global 0.02° Gridded Super-collated SST data
    from 2km to match NASA 9km grid.
"""

# Import necessary packages
import argparse
import os
import subprocess
from netCDF4 import Dataset
from datetime import timedelta
import sys
from dateutil.parser import parse
import warnings
warnings.filterwarnings('ignore')

# Create global variables
ROOT_DIR = ("/Users/madisonrichardson/netpp")
SST_2k_DIR_T = os.path.join(ROOT_DIR, "data/{}/sst_2k")
SST_DIR_T = os.path.join(ROOT_DIR, "data/{}/sst")
WORK_DIR = os.path.join(ROOT_DIR, 'work')
BIN_DIR = os.path.join(ROOT_DIR, 'bin')
RES_DIR = os.path.join(ROOT_DIR, 'resources')
CWUTL_DIR = '"/Applications/CoastWatch Utilities/bin"'
CWUTIL_OFILE_HDF = 'nasa_9k.hdf'
CWUTIL_OFILE_NC = 'nasa_9k.nc'
MASTER_FILE = 'master_file_sst.nc'
CDL_IN_FILE = 'temp21.cdl'
TEMP_OUT_FILE = 'temp21.nc'

# Create functions to source data


def generate_url(date_obj):
    """Create a URL for downloading a data file.

    Creates a URL for downloading a data file from
    NOAA ACSPO Daily Global 0.02 degree ridded Super-collated SST
    data on a specific date.

    Args_:
        date_obj (datetime): A datetime object representing the specified date
                            for which the SST (sea surface temperature) data
                            URL is generated.

    Returns_:
        tuple: A tuple containing:
            - url (str): The complete URL to download the specified data file.
            - file_name (str): The name of the data file to be downloaded.
    """
    year = '{0:%Y}'.format(date_obj)
    doy = '{0:%j}'.format(date_obj)

    date_str = '{0:%Y%m%d}'.format(date_obj)

    base_url = (
        "https://coastwatch.noaa.gov/erddap/files/"
        "noaacwLEOACSPOSSTL3SKDaily/"
    )

    file_name = (
        f"{date_str}120000-STAR-L3S_GHRSST-SSTsubskin-"
        "LEO_Daily-ACSPO_V2.81-v02.0-fv01.0.nc"
    )

    url = f"{base_url}{year}/{doy}/{file_name}"

    return url, file_name


def download_sst_file(url, SST_2k_DIR):
    """Download a file using generate_url_batch output.

    Downloads a file from the specified URL that is generated from
    generate_url_batch function and saves it to the specified
    directory.

    Args_:
        url (str): The URL of the file to be downloaded.
        SST_2k_DIR (str): The directory where the downloaded
        file will be saved.

    """
    command = f"wget -P {SST_2k_DIR} {url}"

    subprocess.call(command, shell=True)

    print(f"Downloaded: {url} to {SST_2k_DIR}")


# Run main function to source data,
# regrid to 4km, and convert to Celsius


def main():
    """Run main function."""

    # Simulated argparse argument variables
    # Set up argparse
    doc_formatter = argparse.RawDescriptionHelpFormatter
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=doc_formatter
        )

    parser.add_argument("-a", "--start",
                        type=str,
                        required=True,
                        help="Start date in YYYY-MM-DD format")
    parser.add_argument("-z",
                        "--end",
                        type=str,
                        required=True,
                        help="End date in YYYY-MM-DD format")
    parser.add_argument("-s",
                        "--sensor",
                        type=str,
                        required=True,
                        choices=['noaa20', 'noaa21', 'both'],
                        help="End date in YYYY-MM-DD format")
    parser.add_argument("-o", "--overwrite", required=False,
                        action='store_true',
                        help="set to overwrite a netpp file that exists")

    args = parser.parse_args()

    # Parse the start and end dates
    start_date = parse(args.start)
    end_date = parse(args.end)

    sensors = ['noaa20', 'noaa21'] if args.sensor == 'both' else [args.sensor]

    # Determine range of years
    started_year = start_date.year
    ended_year = end_date.year

    for args.sensor in sensors:  # Loop over sensors
        current_date = start_date

        # Create directories for all years between start_year and end_year
        for year in range(started_year, ended_year + 1):
            # Create dynamic directories based on sensor and year range
            SST_2k_DIR = os.path.join(SST_2k_DIR_T.format(args.sensor))
            SST_DIR = os.path.join(SST_DIR_T.format(args.sensor))

            DIR_LIST = [ROOT_DIR,
                        WORK_DIR,
                        SST_2k_DIR,
                        SST_DIR,
                        ]

            for dr in DIR_LIST:
                os.makedirs(dr, exist_ok=True)
            print(len(DIR_LIST), 'directories validated')

        # Check for correct date values
        if start_date > end_date:
            sys.exit("start date must be < end date")

        while current_date <= end_date:
            print('Processing this date', current_date)

            # Define the output NetCDF file for the current date
            formatted_date = current_date.strftime('%Y%m%d')
            nc_filename = f'sst_leo_9km_{formatted_date}_daily.nc'

            odir = os.path.join(SST_DIR, str(current_date.year))
            os.makedirs(odir, exist_ok=True)

            nc_file_path = os.path.join(odir, nc_filename)

            # Add logic to not overwrite existing files unless argparse -o
            if os.path.isfile(nc_file_path):
                if not args.overwrite:
                    print(nc_filename, 'already exists')
                    current_date += timedelta(days=1)
                    continue
                else:
                    print('overwriting', nc_filename)

            date_noon = current_date.replace(hour=12, minute=0, second=0)
            timestamp = date_noon.timestamp()

            odir2 = os.path.join(SST_2k_DIR, str(current_date.year))
            os.makedirs(odir2, exist_ok=True)

            # Use a try to catch when data downloads fail
            try:
                # SST download and processing
                sst_url, sst_file = generate_url(current_date)
                download_sst_file(sst_url, odir2)
            except Exception as e:
                print('The NOAA file did not download', e)
                current_date += timedelta(days=1)
                continue

            # Regrid steps here
            formatted_date = current_date.strftime('%Y%m%d')
            ifile = os.path.join(
                odir2,
                f'{formatted_date}120000-STAR-L3S_GHRSST-SSTsubskin-LEO_Daily'
                '-ACSPO_V2.81-v02.0-fv01.0.nc'
            )

            myCmd = ' '.join([os.path.join(CWUTL_DIR, 'cwregister2'),
                              '--clobber',
                              '--match=' + 'sea_surface_temperature',
                              '--master=' + os.path.join(RES_DIR, MASTER_FILE),
                              os.path.join(odir2, ifile),
                              os.path.join(WORK_DIR, CWUTIL_OFILE_HDF)
                              ])
            print(myCmd)
            print('Regrid NASA file',
                  subprocess.call(myCmd, shell=True))

            cw_cmd = ' '.join([os.path.join(CWUTL_DIR, 'cwangles'),
                               '--location --float --units=deg',
                               os.path.join(WORK_DIR, CWUTIL_OFILE_HDF)
                               ])
            print('Add lat/lon', subprocess.call(cw_cmd, shell=True))

            # Convert hdf to netCDF with CW utilities

            myCmd = ' '.join([os.path.join(CWUTL_DIR, 'cwexport'), '-v',
                              os.path.join(WORK_DIR, CWUTIL_OFILE_HDF),
                              os.path.join(WORK_DIR, CWUTIL_OFILE_NC)
                              ])
            print('Convert  to NetCDF',
                  subprocess.call(myCmd, shell=True))

            # Load datasets
            sst_ds = Dataset(os.path.join(WORK_DIR, CWUTIL_OFILE_NC), 'r')

            sst = sst_ds['sea_surface_temperature'][0, 0, :, :]
            sst = sst - 273.15

            sst_ds.close()

            # Generate output template file from cdl file
            myCmd = ' '.join(['ncgen',
                              '-o',
                              os.path.join(WORK_DIR, TEMP_OUT_FILE),
                              os.path.join(RES_DIR, CDL_IN_FILE)])
            print('Run ncgen',
                  subprocess.call(myCmd, shell=True))

            # Open output template file in append mode
            nc_file = Dataset(os.path.join(WORK_DIR, TEMP_OUT_FILE),
                              'a',
                              format='NETCDF4')

            # Write sst, chlorophyll, par, and PPeu data to the netCDF file
            nc_file['sea_surface_temperature'][0, :, :] = sst[:, :]
            nc_file['time'][0] = timestamp

            # Modify metadata
            formatted_date_start = current_date.strftime('%Y-%m-%dT00:00:00Z')
            formatted_date_end = current_date.strftime('%Y-%m-%dT23:59:59Z')
            nc_file.time_coverage_start = formatted_date_start
            nc_file.time_coverage_end = formatted_date_end
            nc_file.title = ', '.join([
                "NOAA LEO SST"
                ])

            # Close netCDF file for this date
            nc_file.close()

            # Compress and archive
            myCmd = ' '.join(['nccopy',
                              '-d6',
                              os.path.join(WORK_DIR, TEMP_OUT_FILE),
                              os.path.join(WORK_DIR, nc_filename)
                              ])
            print('Compress ofile',
                  subprocess.call(myCmd, shell=True))

            myCmd = ' '.join(['mv',
                              os.path.join(WORK_DIR, nc_filename),
                              nc_file_path
                              ])
            print('Archive ofile',
                  subprocess.call(myCmd, shell=True))

            # Print where netCDF files were saved
            print(f"NetCDF file '{nc_filename}' archived at {nc_file_path}")

            # Add code to send to ERDDAP

            current_date += timedelta(days=1)


if __name__ == '__main__':
    main()
