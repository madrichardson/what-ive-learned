"""Create primary productivity satellite-based products.

This script generates primary productivity fields from chlorophyll, SST, PAR
satellite data using the method of Behrenfeld and Falkowski 1997. It accepts
source satellite data that has been gridded to the NASA 9Km SMI.
NOAA CoastWatch standard ocean color grid. As written, the script is tailored
for pairing the following global, 8-day composite input data to produce
primary productivity fields:
    * NASA NOAA20-VIIRS chlorophyll and PAR
    * NASA NOAA21-VIIRS chlorophyll and PAR
    * And SST from the NOAA Gridded Super-collated product

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

# Import packages
import os
import ntpath
import argparse
import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset
import pandas as pd
import subprocess
from datetime import timedelta, datetime, timezone
import sys
from dateutil.parser import parse
import warnings

warnings.filterwarnings("ignore")

# Create global variables
ROOT_DIR = "/Users/madisonrichardson/netpp"
CHL_DIR_T = os.path.join(ROOT_DIR, "data/{}/chl/8_day_chl")
PAR_DIR_T = os.path.join(ROOT_DIR, "data/{}/par/8_day_par")
SST_DIR_T = os.path.join(ROOT_DIR, "data/{}/sst/{}")
WORK_DIR = os.path.join(ROOT_DIR, "work")
BIN_DIR = os.path.join(ROOT_DIR, "bin")
RES_DIR = os.path.join(ROOT_DIR, "resources")
NC_OUT_DIR_T = os.path.join(ROOT_DIR, "data/{}/netpp/8_day_netpp")
CDL_IN_FILE = "composite.cdl"
TEMP_OUT_FILE = "composite.nc"

# Create useful functions


def make_n21_filelist(file_date_ranges, prefix, p_id):
    """
    Generate a list of file URLs for NOAA21 VIIRS 8-day composite data.

    This function constructs URLs for NOAA21 VIIRS data files
    using a base URL, date ranges, a specified file prefix,
    and a product ID. The generated URLs are formatted to
    match the file naming convention for 8-day composites.

    Args:
        file_date_ranges (list of str): A list of date ranges formatted as
        strings (e.g.  '20230407-20230414') representing the 8-day periods.
        prefix (str): The prefix for the file type (e.g., 'SST', 'CHL').
        p_id (str): The product ID to specify the type of data (e.g., 'PAR').

    Returns:
        list of str: A list of full URLs to the data files.
    """
    b_url = "http://oceandata.sci.gsfc.nasa.gov/getfile/"
    fl_list = [
        b_url + f"JPSS2_VIIRS.{fdr}.L3m.8D.{prefix}.{p_id}.9km.NRT.nc"
        for fdr in file_date_ranges
    ]

    return fl_list


def download_nasa20_file_list(start_date, end_date, data_type, directory):
    """
    Download list if NASA data files for NOAA20.

    Args_:
        start_date (datetime): First date of time series to process.
        end_date (datetime): Last date of time series to process.
        data_type (str): The type of data to download, either
                        'chlorophyll' or 'par'.
        directory (str): The directory where the downloaded data will be saved.
                         The word dir is recommended

    Returns_:
        Path to the downloaded NASA file list.
    """

    # Validate data type
    if data_type == "chlorophyll":
        prod_id = "chlor_a"
    elif data_type == "par":
        prod_id = "par"
    else:
        raise ValueError("data_type must be 'chlorophyll' or 'par'")

    # Extract the start and end dates for the 8-day composite
    sdate = start_date.strftime("%Y%m%d")
    edate = end_date.strftime("%Y%m%d")

    dtid = "1197"
    wget_template = " ".join(
        [
            "wget",
            "-q",
            '--post-data="results_as_file=1&sensor_id=33&dtid={}&sdate={}&edate={}',
            '&subType=1&addurl=1&prod_id={}&resolution_id=9km&period=8D"',
            "-O {}",
            "https://oceandata.sci.gsfc.nasa.gov/api/file_search",
        ]
    )

    # File name for list of files
    file_list_name = f"{data_type}_nasa_file_list_" f"{sdate}_{edate}.txt"
    file_list_path = os.path.join(directory, file_list_name)

    # Format wget command for retrieving file list
    nasa_wget = wget_template.format(
        dtid, sdate, edate, prod_id, file_list_path
    )
    subprocess.call(nasa_wget, shell=True)

    return file_list_path


def download_nasa_data_dhr(file_url, directory):
    """
    Download NASA data files for specified sensor.

    Downloads chlorophyll and PAR data files from NASA's
    Ocean Biology Processing Group services using a provided
    file URL. It uses wget with cookie-based authentication to
    access protected resources. The downloaded file is saved
    to the specified directory.

    Args_:
        file_url (str): The URL of the file to download.
        directory (str): The target directory where the
        downloaded file will be saved

    Returns_:
        str: The full path to the downloaded file if successful.
        None: If the download fails.
    """

    # File name for the downloaded file
    import ntpath

    target_file_name = ntpath.basename(file_url)

    downloaded_file_path = os.path.join(directory, target_file_name)

    # Command to download the file using wget with cookies
    wget_Cmd = (
        f"wget --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies "
        f"--auth-no-challenge=on --content-disposition "
        f"-O {downloaded_file_path} {file_url}"
    )

    # Print generated wget command and dates
    print(wget_Cmd)
    print("downloaded = zero", subprocess.call(wget_Cmd, shell=True))

    # Verify if the file was downloaded
    if os.path.exists(downloaded_file_path):
        print("Downloaded:", target_file_name)
        return downloaded_file_path
    else:
        print("Failed to download data", target_file_name)


def meanVar(mean, num, obs):
    """
    Adjusts the running mean and number of valid
    observations by incorporating a new observation
    set ('obs'). It handles masked arrays, ensuring
    that missing or invalid values do not affect the
    mean calculation.

    Args:
        mean (numpy.ndarray): The current mean array, to be updated.
        num (numpy.ndarray): An array that keeps track of the number
                            of valid observations per element.
        obs (numpy.ma.MaskedArray): The new set of observations to
                                    update the mean and count.

    Returns:
        tuple:
            - numpy.ndarray: The updated mean array.
            - numpy.ndarray: The updated count of valid
                            observations array.
    """
    numShape = num.shape
    temp = np.subtract(obs, mean, dtype=np.single)
    numAdd = np.ones(numShape, dtype=np.int32)

    my_mask = ma.getmask(obs)
    numAdd = ma.masked_array(numAdd, mask=my_mask)
    numAdd = numAdd.filled(fill_value=0)

    # numAdd[obs.mask] = 0
    num = np.add(num, numAdd, dtype=np.int32)
    tempNum = ma.array(num, mask=(num == 0), dtype=np.int32)
    tNfloat = tempNum.astype("float")

    temp = np.divide(temp, tNfloat, dtype=np.single)
    # temp = np.divide(temp, tempNum.astype('float'), dtype=np.single)
    mean = np.add(mean, temp.filled(0.0), dtype=np.single)

    return (mean, num)


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

    Z_eu = ma.where(
        chl_eu > 10.0, 568.2 * chl_eu**-0.746, 200.0 * chl_eu**-0.293
    )

    return Z_eu


def daylength(dayOfYear, lat):
    """Determine the length of the daylight period.

    Computes the length of the day (the time between sunrise and
    sunset) given the day of the year and latitude of the location.
    Function uses the Brock model for the computations.
    For more information see, for example,
    Forsythe et al., "A model comparison for daylength as a
    function of latitude and day of year", Ecological Modelling,
    1995. Modified (vectorize) from by Dale Robinson:
    https://gist.github.com/anttilipp/ed3ab35258c7636d87de6499475301ce

    Args_:
        dayOfYear (int): The day of the year, where 1 corresponds to
            the 1st of January.

        lat (ndarray): A numpy array of latitude of the pixel location
        in decimal degrees.
            Positive values for north and negative for south.

    Return:_
        day_len (ndarray): A numpy array of daylength in decimal hours,
            e.g 12:30pm is 12.5.
    """
    # Calculate the  center day for the 8-day composite
    # Center day is the 5th day of the 8-day period
    center_dayOfYear = dayOfYear + 4

    # Correct for leap year
    if dayOfYear == 366:
        dayOfYear = 365

    latInRad = np.deg2rad(lat)
    declinationOfEarth = 23.45 * np.sin(
        np.deg2rad(360.0 * (283.0 + center_dayOfYear) / 365.0)
    )

    cos_hour_angle = -np.tan(latInRad) * np.tan(np.deg2rad(declinationOfEarth))
    cos_hour_angle = np.clip(cos_hour_angle, -1.0, 1.0)

    hourAngle = np.rad2deg(np.arccos(cos_hour_angle))

    day_len = 2.0 * hourAngle / 15.0
    day_len = np.where(cos_hour_angle <= -1.0, 24, day_len)
    day_len = np.where(cos_hour_angle >= 1.0, 0, day_len)

    return day_len


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
    """
    Run main function.

    """

    # Main code for downloading data and calculations for NetPP algorithm
    # Simulated argparse argument variables
    # Set up argparse
    doc_formatter = argparse.RawDescriptionHelpFormatter
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=doc_formatter
    )

    parser.add_argument(
        "-a",
        "--start",
        type=str,
        required=True,
        help="Start date in YYYY-MM-DD format",
    )
    parser.add_argument(
        "-z",
        "--end",
        type=str,
        required=True,
        help="End date in YYYY-MM-DD format",
    )
    parser.add_argument(
        "-s",
        "--sensor",
        type=str,
        required=True,
        choices=["noaa20", "noaa21", "both"],
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
    sensors = ["noaa20", "noaa21"] if args.sensor == "both" else [args.sensor]

    if start_date > end_date:
        sys.exit("Error: start_date must be earlier than end_date.")

    now = datetime.now()

    # Dictionary for added metadata
    end_year = {"noaa20": "2017", "noaa21": "2023"}

    # Generates a list of file URLs for the specified start and
    # end dates and saves it in WORK_DIR
    chl_n20_list = download_nasa20_file_list(
        start_date, end_date, "chlorophyll", WORK_DIR
    )

    # File list contains one URL per line so we load it into a DataFrame
    # with a single columm named "wget_chl_noaa20"
    df_chl = pd.read_csv(chl_n20_list, names=["wget_chl_noaa20"])

    # Extracts all URLs from the "wget_chl_noaa20" column into a Python list
    chl_list = df_chl.wget_chl_noaa20.values.tolist()

    # The date range is extracted and printed to verify the data
    file_date_ranges = [ntpath.basename(ln).split(".")[1] for ln in chl_list]
    df_chl["date_ranges"] = file_date_ranges
    df_chl.head()

    # Removes the temporary file list created
    myCmd = "rm " + chl_n20_list
    print("rm chl file list", subprocess.call(myCmd, shell=True))

    # Generates a list of file URLs for the specified start and end dates 
    # and saves it in WORK_DIR
    par_n20_list = download_nasa20_file_list(
        start_date, end_date, "par", WORK_DIR
    )
    print(par_n20_list)

    # File list contains one URL per line so we load it into a DataFrame 
    # with a single columm named "wget_par_noaa20"
    df_par = pd.read_csv(par_n20_list, names=["wget_par_noaa20"])

    # Extracts all URLs from the "wget_par_noaa20" column into a Python list
    par_list = df_par.wget_par_noaa20.values.tolist()

    # The date range is extracted and printed to verify the data
    file_date_ranges = [ntpath.basename(ln).split(".")[1] for ln in par_list]
    df_par["date_ranges"] = file_date_ranges
    df_par.head()

    # Removes the temporary file list created
    myCmd = "rm " + par_n20_list
    print("rm par file list", subprocess.call(myCmd, shell=True))

    # Keep only the rows where the 'date_ranges' column values are present
    # in both DataFrames
    df_common = pd.merge(df_chl, df_par, on="date_ranges", how="inner")

    # Date ranges that will be used to construct URLs for downloading PAR
    # and CHL
    merged_dt_ranges = df_common.date_ranges.values.tolist()

    # Make URL for PAR
    par_wget_url = make_n21_filelist(merged_dt_ranges, "PAR", "par")

    # Make URL for CHL
    chl_wget_url = make_n21_filelist(merged_dt_ranges, "CHL", "chlor_a")

    df_common["wget_chl_noaa21"] = chl_wget_url
    df_common["wget_par_noaa21"] = par_wget_url

    # Iterate over rows in the DataFrame containing common data
    for index, row in df_common.iterrows():
        # Extract the date range for the current composite
        composite_date_range = row["date_ranges"]  # pull date range from df
        composite_start_date = parse(composite_date_range.split("_")[0])
        composite_end_date = parse(composite_date_range.split("_")[1])
        # Center the timestamp on the 5th day of the 8-day composite
        centered_date = composite_start_date + timedelta(days=4)
        # Make the datetime object timezone-aware and set it to UTC
        centered_date = centered_date.replace(tzinfo=timezone.utc)

        # Specify the SST source directory
        SST_DIR = os.path.join(
            SST_DIR_T.format(args.sensor, composite_start_date.strftime("%Y"))
        )
        print("SST_DIR:", SST_DIR)

        # Format both start and end dates as strings
        formatted_start_date = composite_start_date.strftime("%Y%m%d")
        formatted_end_date = composite_end_date.strftime("%Y%m%d")

        # Do a quick check to see if the output files are made
        # If so, stop generation of files and continue to next date
        if not args.overwrite:
            nc_filename = "productivity_viirs_{}_8day_{}_{}.nc"
            opaths = [
                os.path.join(
                    NC_OUT_DIR_T.format(ln),
                    nc_filename.format(
                        ln, formatted_start_date, formatted_end_date
                    ),
                )
                for ln in sensors
            ]

            if sum([not os.path.isfile(fl) for fl in opaths]) == 0:
                print("output files already exist")
                continue
            else:
                print(
                    "make files for:",
                    formatted_start_date,
                    "to",
                    formatted_end_date,
                )

        # Below are steps to calculate the 8-day composite range
        # 1. Make a list of the dates the composite
        formatted_date = (
            pd.date_range(composite_start_date, composite_end_date, freq="d")
            .strftime("%Y%m%d")
            .tolist()
        )
        # Generate the list of files to make
        fileList = [f"sst_leo_9km_{ln}_daily.nc" for ln in formatted_date]

        # 2. Calculate mean if 4 or more of the 8 files are there
        # Less gaps if a few files are missing
        if (
            sum([os.path.isfile(os.path.join(SST_DIR, fl)) for fl in fileList])
            < 4
        ):
            print(
                "less than 4 of 8 files available for",
                formatted_start_date,
                formatted_end_date,
            )
            continue

        # 3. Calculate the 8-day composite
        first_loop = True
        my_var = "sea_surface_temperature"

        for fl in fileList:
            try:
                # nc = Dataset(os.path.join(SST_DIR, fl), 'r')
                with Dataset(os.path.join(SST_DIR, fl)) as nc:
                    variable = nc.variables[my_var][0, :, :]
            except Exception:
                continue

            if first_loop:
                num = np.zeros(
                    (variable.shape[0], variable.shape[1]), dtype=np.int32
                )
                mean = np.zeros(
                    (variable.shape[0], variable.shape[1]), np.single
                )

            mean, num = meanVar(mean, num, variable)
            first_loop = False
            # nc.close()

        # Apply mask to final composite
        sst = ma.masked_where(num == 0, mean).filled(fill_value=-999.0)

        # Mask sst values < -2 C
        sst = ma.masked_where(sst < -2, sst)

        # Adjust sst values outside of range
        sst_data_mod = ma.where(sst < -1, -1, sst)
        sst_data_mod = ma.where(sst_data_mod > 29, 29, sst_data_mod)

        print("sst_data_mod made", sst_data_mod.min(), sst_data_mod.max())

        # Make optimal photosynthetic rate
        PbOpt = calculate_PbOpt(sst_data_mod)
        print("PbOpt made", PbOpt.min(), PbOpt.max())

        for args.sensor in sensors:  # Loop over sensors
            # Set ofile name and directories
            nc_filename = f"productivity_viirs_{args.sensor}_8day_{formatted_start_date}_{formatted_end_date}.nc"
            NC_OUT_DIR = os.path.join(NC_OUT_DIR_T.format(args.sensor))
            odir = os.path.join(NC_OUT_DIR, str(centered_date.year))
            os.makedirs(odir, exist_ok=True)

            nc_file_path = os.path.join(odir, nc_filename)

            # Try to download data
            try:
                chl_file_path = download_nasa_data_dhr(
                    row["wget_chl_" + args.sensor], WORK_DIR
                )
                par_file_path = download_nasa_data_dhr(
                    row["wget_par_" + args.sensor], WORK_DIR
                )
            except Exception as e:
                print(
                    f"One or more of the NASA files did not download for {args.sensor}",
                    e,
                )
                continue

            # Load datasets and extract data
            with Dataset(chl_file_path, "r") as chl_ds:
                chl = chl_ds["chlor_a"][:, :]

            with Dataset(par_file_path, "r") as par_ds:
                par = par_ds["par"][:, :]

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
                os.path.join(WORK_DIR, TEMP_OUT_FILE), "a", format="NETCDF4"
            )

            # Get lan and lon vectors from temp ofile
            lat_data = nc_file["latitude"][:]
            print("Latitude data:", lat_data)

            lon_data = nc_file["longitude"][:]
            print("Longitude data:", lon_data)

            # Calculate daylength
            # dayOfYear = current_date.timetuple().tm_yday
            dayOfYear = int("{0:%j}".format(centered_date))
            day_len = daylength(dayOfYear, lat_data)
            day_len_2d = np.outer(day_len, np.ones(len(lon_data)))

            # Generate NetPP
            PPeu = calculate_PPeu(chl, PbOpt, Z_eu, par, day_len_2d)
            print("PPeu made", PPeu.min(), PPeu.max())

            # There should be no negative numbers, but the source data might
            # have errors. So, mask out values < 0 to be sure
            PPeu = ma.masked_where(PPeu <= 0, PPeu)

            # Write sst, chlorophyll, par, and PPeu data to the netCDF file
            nc_file["sea_surface_temperature"][0, :, :] = sst[:, :]
            nc_file["chlor_a"][0, :, :] = chl[:, :]
            nc_file["par"][0, :, :] = par[:, :]
            nc_file["productivity"][0, :, :] = PPeu[:, :]
            nc_file["time"][0] = int(centered_date.timestamp())

            # Modify metadata
            formatted_date_start = composite_start_date.strftime(
                "%Y-%m-%dT00:00:00Z"
            )
            formatted_date_end = composite_end_date.strftime(
                "%Y-%m-%dT00:00:00Z"
            )
            nc_file.time_coverage_start = formatted_date_start
            nc_file.time_coverage_end = formatted_date_end
            nc_file.date_created = now.isoformat("T", "seconds")
            nc_file.platform = args.sensor.upper()
            nc_file.id = f"productivity_{args.sensor}_8day"
            nc_file.keywords = ", ".join(
                [
                    "chla",
                    "chlor_a",
                    "chlorophyll",
                    "chlorophyll-a",
                    "coastwatch",
                    "Earth Science > Biosphere > Ecological Dynamics > Ecosystem Functions > Primary Production",
                    "Earth Science > Biosphere > Vegetation > Carbon",
                    "Earth Science > Biosphere > Vegetation > Photosynthetically Active Radiation",
                    "Earth Science > Oceans > Ocean Chemistry > Carbon",
                    "Earth Science > Oceans > Ocean Chemistry > Chlorophyll",
                    "Earth Science > Oceans > Ocean Chemistry > Pigments > Chlorophyll",
                    "Earth Science > Oceans > Ocean Optics > Photosynthetically Available Radiation",
                    "Earth Science > Oceans > Ocean Temperature > Sea Surface Temperature",
                    "mass_concentration_of_chlorophyll_in_sea_water",
                    "net_primary_productivity_of_carbon",
                    "noaa, par, production, productivity",
                    "sea_surface_temperature",
                    "sst",
                    "surface_downwelling_photosynthetic_photon_flux_in_air",
                    "viirs",
                    "west coast node",
                ]
            )
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
                    f"VIIRS {args.sensor.upper()}",
                    "NRT",
                    "Global",
                    "9km",
                    "{}-present (8 Day Composite)".format(
                        end_year[args.sensor]
                    ),
                ]
            )
            # Custom summary and history for noaa21
            if args.sensor == "noaa21":
                nc_file.summary = " ".join(
                    [
                        "The Visible and Infrared Imager/Radiometer",
                        "Suite (VIIRS), NOAA21 Primary Productivity product",
                        "provides near real-time (NRT) estimates of net carbon",
                        "fixation by phytoplankton using enhanced",
                        "satellite measurements with the",
                        "algorithm by Behrenfeld and Falkowski (1997).",
                        "This product incorporates NRT chlorophyll a and",
                        "Photosynthetically Available Radiation (PAR) data"
                        "along with high-resolution SST data.",
                        "Mapped to a NASA 9km Standard Mapped Image for accuracy",
                        "in global assessments.",
                    ]
                )
                nc_file.history = " ".join(
                    [
                        "Chlorophyll a, PAR, and SST satellite near",
                        "real-time (NRT) data were applied to the equation",
                        "of Behrenfeld and Falkowski, 1997",
                    ]
                )
            else:
                # Default summary for other args.sensors
                nc_file.summary = " ".join(
                    [
                        "The Visible and Infrared Imager/Radiometer",
                        "Suite (VIIRS), {}".format(args.sensor.upper()),
                        "Primary Productivity product estimates net carbon fixation by phytoplankton",
                        "in oceanic waters using the algorithm of Behrenfeld and P. G. Falkowski (1997) with",
                        "chlorophyll a and Photosynthetically Available Radiation (PAR) values from VIIRS {} and".format(
                            args.sensor.upper()
                        ),
                        "Sea Surface Temperature (SST) values from the",
                        "NOAA Gridded Super-collated product as inputs.",
                        "The science quality chlorophyll a, SST, and PAR data are",
                        "included in the dataset. Data are mapped to a",
                        "NASA 9km Standard Mapped Image.",
                    ]
                )
                nc_file.history = " ".join(
                    [
                        "Chlorophyll a, PAR, and SST satellite",
                        "data were applied to the equation"
                        "of Behrenfeld and Falkowski, 1997",
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

            # uncomment to clean up the work dir
            myCmd = "rm " + chl_file_path
            print("remove chl file", subprocess.call(myCmd, shell=True))
            myCmd = "rm " + par_file_path
            print("remove par file", subprocess.call(myCmd, shell=True))


if __name__ == '__main__':
    main()
