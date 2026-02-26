"""Create primary productivity satellite-based products.

This script generates primary productivity fields from chlorophyll, SST, PAR
satellite data using the method of Behrenfeld and Falkowski 1997. It accepts
source satellite data that has been gridded to the NASA 9Km SMI.
NOAA CoastWatch standard ocean color grid. As written, the script is tailored
for pairing the following global, monthly input data to produce
primary productivity fields:
    * NASA NOAA20-VIIRS chlorophyll and PAR
    * NASA NOAA21-VIIRS chlorophyll and PAR
    * And SST frome the NOAA Gridded Super-collated product

Users set the satellite input data, start dates and stop date
via command-line arguments.

The output files are created by generating a template netCDF file from a cdl
file that is prepopulate with metadata but has only latitude and longitude
data. The template file is populated with the primary productivity data
and then renamed.

The script can be repurposed for different input data with these modifications:
    * Gridding any input chlorophyll, PAR, and SST data to a common grid
    * Creating an appropriate cdl file for the dataset
    * Adjusting the logic to make directory paths
    * Adjust output file naming and input file search pattern.
"""

# Import necessary packages
import argparse
from calendar import monthrange
from dateutil.relativedelta import relativedelta
import numpy.ma as ma
import numpy as np
import os
import subprocess
import pandas as pd
from netCDF4 import Dataset
from datetime import datetime, timezone
import sys
from dateutil.parser import parse
import warnings

warnings.filterwarnings("ignore")

# Create your global variables
ROOT_DIR = "/home/madison/projects/netpp"
CHL_DIR_T = os.path.join(ROOT_DIR, "data/{}/chl/monthly_chl")
PAR_DIR_T = os.path.join(ROOT_DIR, "data/{}/par/monthly_par")
SST_DIR_T = os.path.join(ROOT_DIR, "data/{}/monthly_sst")
WORK_DIR = os.path.join(ROOT_DIR, "work")
BIN_DIR = os.path.join(ROOT_DIR, "bin")
RES_DIR = os.path.join(ROOT_DIR, "resources")
NC_OUT_DIR_T = os.path.join(ROOT_DIR, "data/{}/monthly_netpp")
# Add the ncgen input cdl and output nc file names
CDL_IN_FILE = "MonthlyNetpp.cdl"
TEMP_OUT_FILE = "temp_out_file.nc"

# Create functions


def download_nasa_data(date, data_type, directory):
    """Download NASA data files.

    Downloads chlorophyll and PAR data from
    L3 NOAA20-VIIRS on NASA data search website.

    Args_:
        date (datetime): The date format %Y-%m-%d for which the data is
                        downloaded.
        data_type (str): The type of data to download, either
                        'chlorophyll' or 'par'.
        directory (str): The directory where the downloaded data will be saved.

    Raises_:
        ValueError: If the 'data_type is not 'chlorophyll' or 'par'.
    """

    dtid = "1197"
    if data_type == "chlorophyll":
        prod_id = "chlor_a"
        product_prefix = "CHL"
    elif data_type == "par":
        prod_id = "par"
        product_prefix = "PAR"
    else:
        raise ValueError("data_type must be chlorophyll or par")

    wget_template = " ".join(
        [
            "wget",
            "-q",
            '--post-data="results_as_file=1&sensor_id=33&dtid={}&sdate={}&edate={}',
            '&subType=1&addurl=1&prod_id={}&resolution_id=9km&period=MO"',
            "-O {}",
            "https://oceandata.sci.gsfc.nasa.gov/api/file_search",
        ]
    )

    sdate = date.strftime("%Y-%m-%d 00:00:00")
    edate = date.strftime("%Y-%m-%d 23:59:59")

    # File name for list of files
    file_list_name = f"{data_type}_nasa_file_list" f'{date.strftime("%Y%m")}.txt'
    file_list_path = os.path.join(directory, file_list_name)

    # Format wget command for retrieving file list
    nasa_wget = wget_template.format(dtid, sdate, edate, prod_id, file_list_path)

    subprocess.call(nasa_wget, shell=True)

    df = pd.read_csv(file_list_path, names=["wget_url"])
    files = list(df.wget_url)
    file_url = next((ln for ln in files if "NRT" not in ln), "default_value")

    # Debugging output: Print full wget command and the dates
    print(f"Start Date: {sdate}")
    print(f"End Date: {edate}")
    print(f"Generated wget command: {nasa_wget}")

    # Command to download file
    myCmd_template = " ".join(
        [
            "wget",
            "--load-cookies ~/.urs_cookies",
            "--save-cookies ~/.urs_cookies",
            "--auth-no-challenge=on",
            "--content-disposition",
            "-P {}",
            "{}",
        ]
    )
    subprocess.call(myCmd_template.format(directory, file_url), shell=True)

    # Verify downloaded file exists
    downloaded_file_name = os.path.basename(file_url)
    downloaded_file_path = os.path.join(directory, downloaded_file_name)

    if os.path.exists(downloaded_file_path):
        print(
            f"Downloaded {data_type} data for"
            f"{date.strftime("%Y%m")} to"
            f"{downloaded_file_path}"
        )

        # Make new filename in YYYYMM format
        new_file_name = (
            f"JPSS1_VIIRS.{date.strftime("%Y%m")}.L3m.MO."
            f"{product_prefix}.{prod_id}.9km.nc"
        )
        new_file_path = os.path.join(directory, new_file_name)

        # Rename downloaded file
        os.rename(downloaded_file_path, new_file_path)

        print(f"Renamed file to: {new_file_path}")

    else:
        print(f"Failed to download {data_type} data for" f"{date.strftime("%Y%m")}")


def daylength_month(month, lat):
    """Determine the length of the daylight period for a given month.

    Computes the length of the day (the time between sunrise and
    sunset) for the middle day of the month (16th day) given
    the latitude of the location.
    Function uses the Brock model for the computations.
    For more information see, for example,
    Forsythe et al., "A model comparison for daylength as a
    function of latitude and day of year", Ecological Modelling,
    1995. Modified (vectorize) from by Dale Robinson:
    https://gist.github.com/anttilipp/ed3ab35258c7636d87de6499475301ce


    Args:
        month (int): The month of the year, where 1 corresponds to January
        and 12 corresponds to December.
        lat (ndarray): A numpy array of latitude of the pixel
        location in decimal degrees. Positive values for north and negative
        for south.


    Returns:
        day_len (ndarray): A numpy array of daylength in decimal hours,
        e.g., 12:30pm is 12.5.
    """

    # Define the approximate day of the year for the 16th day of each month
    dayOfYear_month = {
        1: 16,
        2: 47,
        3: 75,
        4: 106,
        5: 136,
        6: 167,
        7: 197,
        8: 228,
        9: 259,
        10: 289,
        11: 320,
        12: 350,
    }

    # Get the approximate day of the year for the 16th day of the given month
    dayOfYear = dayOfYear_month.get(month, 16)

    # Correct for leap year
    if dayOfYear == 366:
        dayOfYear = 365

    latInRad = np.deg2rad(lat)
    declinationOfEarth = 23.45 * np.sin(np.deg2rad(360.0 * (283.0 + dayOfYear) / 365.0))

    cos_hour_angle = -np.tan(latInRad) * np.tan(np.deg2rad(declinationOfEarth))
    cos_hour_angle = np.clip(cos_hour_angle, -1.0, 1.0)

    hourAngle = np.rad2deg(np.arccos(cos_hour_angle))

    day_len = 2.0 * hourAngle / 15.0
    day_len = np.where(cos_hour_angle <= -1.0, 24, day_len)
    day_len = np.where(cos_hour_angle >= 1.0, 0, day_len)

    return day_len


def calculate_PbOpt(sst_data_mod):
    """Calculate the maximum chlorophyll fixation rate.

    Calculates the maximum chlorophyll fixation rate (PbOpt) within the
    water column based on modified SST data using a seventh-order
    polynomial equation.

    Args_:
        sst_data_mod (ndarray or float): Modified SST where SST < -1 equal
            to the value at -1 and to set all values where SST > 29 equal to
            the value at 29. This can either be a single value or a numpy
            array of SST values.

    Returns_:
        ndarray or float: The calculated maximum chlorophyll fixation rate
                        (PbOpt) (units: mg C (mg chl)^-1 h^-1) corresponding
                        to the input SST data. The return matches the input
                        data meaning if a single value is inputted, a float
                        is returned, or if the input is an array, a numpy
                        array is returned.

    """
    return (
        -3.27e-8 * sst_data_mod**7
        + 3.4132e-6 * sst_data_mod**6
        - 1.348e-4 * sst_data_mod**5
        + 2.462e-3 * sst_data_mod**4
        - 0.0205 * sst_data_mod**3
        + 0.0617 * sst_data_mod**2
        + 0.2749 * sst_data_mod
        + 1.2956
    )


def calculate_Z_eu(chl):
    """Calculate the euphotic depth.

    Calculates the euphotic depth where light is 1% of that at the surface
    (Z_eu) based on chlorophyll-a concentration (CHL_eu) using the Case I
    models of Morel and Berthon (1989).

    Args_:
        chl (xarray.DataArray or ndarray): Chlorophyll-a concentration
                                                (mg m^-3). The input can be
                                                either an xarray.DataArray or
                                                a numpy array.

    Return_:
        xarray.DataArray or ndarray: The calculated euphotic depth (Z_eu) in
                                    meters, where the return type matches the
                                    input type.

    """
    chl_eu = ma.where(chl > 1.0, 40.2 * chl**0.5070, 38.0 * chl**0.4250)

    Z_eu = ma.where(chl_eu > 10.0, 568.2 * chl_eu**-0.746, 200.0 * chl_eu**-0.293)

    return Z_eu


def calculate_PPeu(chl, Pbopt, Z_eu, par, day_len_2d):
    """Calculate the daily depth-integrated primary production.

    Calculates the daily depth-integrated primary production (PP_eu) using
    a Vertically Generalized Production Model (VGPM).

    Args_:
        chl (xarray.DataArray or ndarray): Chlorophyll-a concentration
            (mg m^-3).
        Pbopt (xarray.DataArray or ndarray): Maximum chlorophyll fixation rate
            (PbOpt) (mg C (mg chl)^-1 h^-1).
        Z_eu (xarray.DataArray or ndarray): Euphotic depth (Z_eu) in meters.
        par (xarray.DataArray or ndarray): Photosynthetically ACtive
        Radiation (PAR) (Einstein m^-2 d^-1).
        day_len_2d (xarray.DataArray or ndarray): Length of the daylight
            period expanded into a 2D array to match the dimensions
            of the other inputs.

    Returns_:
        xarray.DataArray or ndarray: The calculated primary production (PP_eu)
            (mg C m^-2 d^-1), with the return type matching the input types.

    """
    par_ratio = par / (par + 4.1)
    PPeu = 0.66125 * Pbopt * par_ratio * Z_eu * chl * day_len_2d
    return PPeu


def main():
    """Run main function."""
    # Verify and create all directories at once

    # Main code for downloading data and calculations for NetPP algorithm
    # Simulated argparse argument variables
    # Set up argparse
    doc_formatter = argparse.RawDescriptionHelpFormatter
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=doc_formatter
    )

    parser.add_argument(
        "-a",
        "--start",
        type=str,
        required=True,
        help="Start date in YYYY-MM-DD format"
    )
    parser.add_argument(
        "-z",
        "--end",
        type=str,
        required=True,
        help="End date in YYYY-MM-DD format"
    )
    parser.add_argument(
        "-s",
        "--sensor",
        type=str,
        required=True,
        choices=["noaa20", "noaa21"],
        help="End date in YYYY-MM-DD format",
    )
    parser.add_argument(
        "-o",
        "--overwrite",
        required=False,
        action="store_true",
        help="set to overwrite a netpp file that exists",
    )

    args = parser.parse_args()

    # Parse the start and end dates
    start_date = parse(args.start)
    end_date = parse(args.end)

    now = datetime.now()

    # Dictionary for added metadata
    end_year = {"noaa20": "2017", "noaa21": "2023"}

    # Create dynamic directories and verify
    CHL_DIR = CHL_DIR_T.format(args.sensor)
    PAR_DIR = PAR_DIR_T.format(args.sensor)
    SST_DIR = SST_DIR_T.format(args.sensor)
    NC_OUT_DIR = NC_OUT_DIR_T.format(args.sensor)

    DIR_LIST = [ROOT_DIR, WORK_DIR, CHL_DIR, PAR_DIR, SST_DIR, NC_OUT_DIR]

    for dr in DIR_LIST:
        os.makedirs(dr, exist_ok=True)
    print(len(DIR_LIST), "directories validated")

    # Check for correct date values
    if start_date > end_date:
        sys.exit("start date must be < end date")

    current_date = start_date
    while current_date <= end_date:
        print("Processing this date", current_date)

        # Define the output NetCDF file for the current date
        formatted_date = current_date.strftime("%Y%m")
        nc_filename = (
            f"productivity_month_{args.sensor}_{formatted_date}_9km.nc"
        )

        odir = os.path.join(NC_OUT_DIR, str(current_date.year))
        os.makedirs(odir, exist_ok=True)

        nc_file_path = os.path.join(odir, nc_filename)

        # Add logic to not overwrite existing files unless argparse -o
        if os.path.isfile(nc_file_path):
            if not args.overwrite:
                print(nc_filename, "already exists")
                current_date += relativedelta(months=1)
                continue
            else:
                print("overwriting", nc_filename)

        middle_of_month = current_date.replace(
            day=16, hour=12, minute=0, second=0, tzinfo=timezone.utc
        )
        timestamp = middle_of_month.timestamp()

        # Use a try to catch when data downloads fail
        # Put the NASA downloads here
        try:
            download_nasa_data(current_date, "chlorophyll", CHL_DIR)
            download_nasa_data(current_date, "par", PAR_DIR)
        except Exception as e:
            print("One or more of the NASA files did not download", e)
            current_date += relativedelta(months=1)
            continue

        # Load datasets
        chl_file = os.path.join(
            CHL_DIR, f"JPSS1_VIIRS.{formatted_date}.L3m.MO.CHL.chlor_a.9km.nc"
        )
        par_file = os.path.join(
            PAR_DIR, f"JPSS1_VIIRS.{formatted_date}.L3m.MO.PAR.par.9km.nc"
        )
        sst_file = os.path.join(
            SST_DIR,
            f"sst_leo_9km_{formatted_date}_monthly.nc"
        )

        sst_ds = Dataset(sst_file, "r")
        chl_ds = Dataset(chl_file, "r")
        par_ds = Dataset(par_file, "r")

        sst = sst_ds["sea_surface_temperature"][0, :, :]
        sst = sst - 273.15
        chl = chl_ds["chlor_a"][:, :]
        par = par_ds["par"][:, :]

        sst_ds.close()
        chl_ds.close()
        par_ds.close()

        # Mask all sst values < -2
        sst = ma.masked_where(sst < -2, sst)

        # Adjust sst values outside of range
        sst_data_mod = ma.where(sst < -1, -1, sst)
        sst_data_mod = ma.where(sst_data_mod > 29, 29, sst_data_mod)
        print("sst_data_mod made", sst_data_mod.min(), sst_data_mod.max())

        # Calculate PbOPt and verify
        PbOpt = calculate_PbOpt(sst_data_mod)
        print("PbOpt made", PbOpt.min(), PbOpt.max())

        # Calculate components of the algorithm
        Z_eu = calculate_Z_eu(chl)
        print("Z_eu made", Z_eu.min(), Z_eu.max())

        # Generate output template file from cdl file
        myCmd = " ".join(
            [
                "ncgen",
                "-o",
                os.path.join(WORK_DIR, TEMP_OUT_FILE),
                os.path.join(RES_DIR, CDL_IN_FILE),
            ]
        )
        print("Run ncgen", subprocess.call(myCmd, shell=True))

        # Open output template file in append mode
        nc_file = Dataset(
            os.path.join(WORK_DIR, TEMP_OUT_FILE),
            "a",
            format="NETCDF4"
        )

        # Get lan and lon vectors from temp ofile
        lat_data = nc_file["latitude"][:]
        lon_data = nc_file["longitude"][:]

        # Calculate daylength
        month = current_date.month
        daylength1D = daylength_month(month, lat_data)
        day_len_2d = np.outer(daylength1D, np.ones(len(lon_data)))

        # Generate netPP
        PPeu = calculate_PPeu(chl, PbOpt, Z_eu, par, day_len_2d)
        print("PPeu made", PPeu.min(), PPeu.max())

        # Mask out values < 0 to be sure there are no negative numbers
        PPeu = ma.masked_where(PPeu <= 0, PPeu)

        # Write sst, chlorophyll, par, and PPeu data to the netCDF file
        nc_file["sea_surface_temperature"][0, :, :] = sst[:, :]
        nc_file["chlor_a"][0, :, :] = chl[:, :]
        nc_file["par"][0, :, :] = par[:, :]
        nc_file["productivity"][0, :, :] = PPeu[:, :]
        nc_file["time"][0] = timestamp

        # Get first and last days of the month
        first_day_of_month = current_date.replace(day=1)
        last_day_of_month = current_date.replace(
            day=monthrange(current_date.year, current_date.month)[1]
        )

        # Format the start and end date for the month
        formatted_date_start = first_day_of_month.strftime("%Y-%m-%dT00:00:00Z")
        formatted_date_end = last_day_of_month.strftime("%Y-%m-%dT23:59:59Z")

        # Modify NetCDF metadata
        nc_file.time_coverage_start = formatted_date_start
        nc_file.time_coverage_end = formatted_date_end
        nc_file.date_created = now.isoformat("T", "seconds")
        nc_file.platform = args.sensor.upper()
        nc_file.id = f"productivity_{args.sensor}_month"
        nc_file.acknowledgement = "The project was supported by funding from the Portfolio Management Branch of NESDIS and NOAA CoastWatch."
        nc_file.contributors = "Dale Robinson, Isaac Shroeder, Ryan Vandermeulen, Jonathan Sherman, Jesse Espinoza, & Madison Richardson"
        nc_file.product_name = "VIIRS {} Primary Productivity".format(
            args.sensor.upper()
        )
        nc_file.source = "satellite observations from VIIRS {}".format(
            args.sensor.upper()
        )
        nc_file.product_name = "VIIRS {} SNPP Primary Productivity".format(
            args.sensor.upper()
        )
        nc_file.title = ", ".join(
            [
                "Primary Productivity",
                "VIIRS {}".format(args.sensor.upper()),
                "Science Quality",
                "Global",
                "9km",
                "{}-present (1 Month Composite)".format(end_year[args.sensor]),
            ]
        )
        nc_file.summary = " ".join(
            [
                "The Visible and Infrared Imager/Radiometer",
                "Suite (VIIRS), {}".format(args.sensor.upper()),
                "Primary Productivity product",
                "estimates net carbon fixation by phytoplankton",
                "in oceanic waters using the algorithm of",
                "Behrenfeld and P. G. Falkowski (1997) with",
                "chlorophyll a and Photosynthetically Available",
                "Radiation (PAR) values from VIIRS {} and".format(args.sensor.upper()),
                "Sea Surface Temperature (SST) values from the",
                "NOAA Gridded Super-collated product as inputs."
                "The science quality chlorophyll a, SST, and PAR",
                "data are",
                "included in the dataset. Data are mapped to a",
                "NASA 9km Standard Mapped Image.",
            ]
        )

        # Close netCDF file for this date
        nc_file.close()

        # Compress and archive
        myCmd = " ".join(
            [
                "nccopy",
                "-d6",
                os.path.join(WORK_DIR, TEMP_OUT_FILE),
                os.path.join(WORK_DIR, nc_filename),
            ]
        )
        print("Compress ofile", subprocess.call(myCmd, shell=True))

        myCmd = " ".join(
            ["mv", os.path.join(WORK_DIR, nc_filename), nc_file_path]
        )
        print("Archive ofile", subprocess.call(myCmd, shell=True))

        # Print where netCDF files were saved
        print(f"NetCDF file '{nc_filename}' archived at {nc_file_path}")

        # Add code to send to ERDDAP

        current_date += relativedelta(months=1)


if __name__ == "__main__":
    main()
