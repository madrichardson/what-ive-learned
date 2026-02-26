"""Create primary productivity satellite-based products.

This script generates primary productivity fields from chlorophyll, SST, PAR
satellite data using the method of Behrenfeld and Falkowski 1997. It accepts
source satellite data that has been gridded to the NASA 9Km SMI.
NOAA CoastWatch standard ocean color grid. As written, the script is tailored
for pairing the following global, 8-day composite input data to produce
primary productivity fields:
    * NASA AQUA-MODIS chlorophyll, PAR, and SST

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
import os
import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset
import subprocess
from datetime import timedelta, datetime
import sys
from dateutil.parser import parse
import warnings
warnings.filterwarnings('ignore')

# Create global variables
# ROOT_DIR = ("/home/madison/projects/netpp")
# CHL_DIR_T = os.path.join(ROOT_DIR, "data/{}/8_day/chl")
# PAR_DIR_T = os.path.join(ROOT_DIR, "data/{}/8_day/par")
# SST_DIR_T = os.path.join(ROOT_DIR, "data/{}/8_day/sst")
ROOT_DIR = ("/Users/madisonrichardson/netpp")
CHL_DIR_T = os.path.join(ROOT_DIR, "data/{}/8_day/chl")
PAR_DIR_T = os.path.join(ROOT_DIR, "data/{}/8_day/par")
SST_DIR_T = os.path.join(ROOT_DIR, "data/{}/8_day/sst")
WORK_DIR = os.path.join(ROOT_DIR, 'work')
BIN_DIR = os.path.join(ROOT_DIR, 'bin')
RES_DIR = os.path.join(ROOT_DIR, 'resources')
# NC_OUT_DIR_T = os.path.join(ROOT_DIR, "data/{}/8_day/netpp")
NC_OUT_DIR_T = os.path.join(ROOT_DIR, "data/{}/8_day/netpp")
CDL_IN_FILE = '8daycomposite.cdl'
TEMP_OUT_FILE = '8daycomposite.nc'

# Create functions


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
    declinationOfEarth = (23.45
                          * np.sin(np.deg2rad(360.0 * (283.0 + center_dayOfYear)
                                              / 365.0))
                          )

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
    chl_eu = ma.where(chl > 1.0,
                      40.2 * chl**0.5070,
                      38.0 * chl**0.4250)

    Z_eu = ma.where(chl_eu > 10.0,
                    568.2 * chl_eu**-0.746,
                    200.0 * chl_eu**-0.293)

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


# Run main function


def main():

    """
    Run main function.

    """

    # Main code for downloading data and calculations for NetPP algorithm
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
                        choices=['aqua_modis', 'snpp_viirs'],
                        help="End date in YYYY-MM-DD format")
    parser.add_argument("-o", "--overwrite", required=False,
                        action='store_true',
                        help="set to overwrite a netpp file that exists")

    args = parser.parse_args()

    # Parse the start and end dates
    start_date = parse(args.start)
    end_date = parse(args.end)

    now = datetime.now()

    # Determine range of years
    started_year = start_date.year
    ended_year = end_date.year

    # Create directories for all years between start_year and end_year
    for year in range(started_year, ended_year + 1):
        # Create dynamic directories based on sensor and year 
        # range (except for NC_OUT)
        CHL_DIR = os.path.join(CHL_DIR_T.format(args.sensor), str(year))
        PAR_DIR = os.path.join(PAR_DIR_T.format(args.sensor), str(year))
        SST_DIR = os.path.join(SST_DIR_T.format(args.sensor), str(year))
        NC_OUT_DIR = os.path.join(NC_OUT_DIR_T.format(args.sensor))

        DIR_LIST = [ROOT_DIR, CHL_DIR, PAR_DIR, NC_OUT_DIR]

        for dr in DIR_LIST:
            os.makedirs(dr, exist_ok=True)
        print(len(DIR_LIST), 'directories validated')

    # Check for correct date values
    if start_date > end_date:
        sys.exit("start date must be < end date")

    current_date = start_date

    while current_date <= end_date:
        print(
            f'Processing 8-day composite starting '
            f'{current_date} for {args.sensor}'
        )

        # Calculate the 8-day composite range
        composite_start_date = current_date
        composite_end_date = current_date + timedelta(days=7)

        # Format both start and end dates
        formatted_start_date = composite_start_date.strftime('%Y%m%d')
        formatted_end_date = composite_end_date.strftime('%Y%m%d')

        # Define the output NetCDF file for the 8-day period
        nc_filename = (
            f'netpp_{args.sensor}_8day_{formatted_start_date}'
            f'_{formatted_end_date}.nc'
        )

        odir = os.path.join(NC_OUT_DIR, str(current_date.year))
        os.makedirs(odir, exist_ok=True)

        nc_file_path = os.path.join(odir, nc_filename)

        # Add logic to not overwrite existing files unless argparse -o
        if os.path.isfile(nc_file_path):
            if not args.overwrite:
                print(f'{nc_filename} already exists for {args.sensor}')
                current_date += timedelta(days=8)
                continue
            else:
                print(f'Overwriting {nc_filename} for {args.sensor}')

        # Calculate the fifth day as the timestamp
        fifth_day = current_date + timedelta(days=4)
        time_stamp = fifth_day.replace(hour=0, minute=0, second=0).timestamp()

        # Load datasets
        chl_file = os.path.join(
            CHL_DIR,
            f'{args.sensor.upper()}.{formatted_start_date}_'
            f'{formatted_end_date}.L3m.8D.CHL.chlor_a.9km.nc'
        )
        par_file = os.path.join(
            PAR_DIR,
            f'{args.sensor.upper()}.{formatted_start_date}_'
            f'{formatted_end_date}.L3m.8D.PAR.par.9km.nc'
        )
        sst_file = os.path.join(
            SST_DIR,
            f'{args.sensor.upper()}.{formatted_start_date}_'
            f'{formatted_end_date}.L3m.8D.SST.sst.9km.nc'
        )

        chl_ds = Dataset(chl_file, 'r')
        par_ds = Dataset(par_file, 'r')
        sst_ds = Dataset(sst_file, 'r')

        chl = chl_ds['chlor_a'][:, :]
        par = par_ds['par'][:, :]
        sst = sst_ds['sst'][:, :]

        # Close the datasets after loading
        chl_ds.close()
        par_ds.close()
        sst_ds.close()

        # Mask sst values < -2 C
        sst = ma.masked_where(sst < -2, sst)

        # Adjust sst values outside of range
        sst_data_mod = ma.where(sst < -1, -1, sst)
        sst_data_mod = ma.where(sst_data_mod > 29, 29, sst_data_mod)

        print('sst_data_mod made', sst_data_mod.min(), sst_data_mod.max())

        # Calculate PbOPt and verify
        PbOpt = calculate_PbOpt(sst_data_mod)
        print('PbOpt made', PbOpt.min(), PbOpt.max())

        # Calculate components of the algorithm
        Z_eu = calculate_Z_eu(chl)
        print('Z_eu made', Z_eu.min(), Z_eu.max())

        # Generate output template file from cdl file
        myCmd = ' '.join(
            [
                'ncgen',
                '-o',
                os.path.join(WORK_DIR, TEMP_OUT_FILE),
                os.path.join(RES_DIR, CDL_IN_FILE)
            ]
        )
        print('Run ncgen', subprocess.call(myCmd, shell=True))

        # Open output template file in append mode
        nc_file = Dataset(os.path.join(WORK_DIR, TEMP_OUT_FILE),
                          'a',
                          format='NETCDF4')

        # Get lan and lon vectors from temp ofile
        lat_data = nc_file['latitude'][:]
        lon_data = nc_file['longitude'][:]

        # Calculate daylength
        dayOfYear = current_date.timetuple().tm_yday
        day_len = daylength(dayOfYear, lat_data)
        day_len_2d = np.outer(day_len, np.ones(len(lon_data)))

        # Generate NetPP
        PPeu = calculate_PPeu(chl, PbOpt, Z_eu, par, day_len_2d)
        print('PPeu made', PPeu.min(), PPeu.max())

        # There should be not negative numbers, but the source data might
        # have errors. So, mask out values < 0 to be sure
        PPeu = ma.masked_where(PPeu <= 0, PPeu)

        # Write sst, chlorophyll, par, and PPeu data to the netCDF file
        nc_file['sea_surface_temperature'][0, :, :] = sst_data_mod[:, :]
        nc_file['chlor_a'][0, :, :] = chl[:, :]
        nc_file['par'][0, :, :] = par[:, :]
        nc_file['productivity'][0, :, :] = PPeu[:, :]
        nc_file['time'][0] = time_stamp

        # Modify metadata
        formatted_date_start = current_date.strftime('%Y-%m-%dT00:00:00Z')
        formatted_date_end = (
            (current_date + timedelta(days=7)).strftime('%Y-%m-%dT23:59:59Z')
        )
        nc_file.time_coverage_start = formatted_date_start
        nc_file.time_coverage_end = formatted_date_end
        nc_file.date_created = now.isoformat("T", "seconds")
        nc_file.platform = args.sensor.upper()
        nc_file.id = f"productivity_{args.sensor}_8day"
        nc_file.keywords = ', '.join([
            "2023-present",
            "chla",
            "chlor_a",
            "chlorophyll",
            "chlorophyll-a",
            "coastwatch",
            "Earth Science > Biosphere > Ecological Dynamics > ",
            "Ecosystem Functions > Primary Production",
            "Earth Science > Biosphere > Vegetation > Carbon",
            "Earth Science > Biosphere > Vegetation > "
            "Photosynthetically Active Radiation",
            "Earth Science > Oceans > Ocean Chemistry > Carbon",
            "Earth Science > Oceans > Ocean Chemistry > Chlorophyll",
            "Earth Science > Oceans > Ocean Chemistry > Pigments > ",
            "Chlorophyll Earth Science > Oceans > Ocean Optics > ",
            "Photosynthetically Available Radiation",
            "Earth Science > Oceans > Ocean Temperature > Sea ",
            "Surface Temperature",
            "mass_concentration_of_chlorophyll_in_sea_water",
            "net_primary_productivity_of_carbon",
            "noaa, par, production, productivity",
            "sea_surface_temperature",
            "sst",
            "surface_downwelling_photosynthetic_photon_flux_in_air",
            "aqua",
            "west coast node"
            ])
        nc_file.product_name = "{} Primary Productivity".format(
            args.sensor.upper())
        nc_file.source = "satellite observations from {}".format(
            args.sensor.upper())
        nc_file.product_name = "{} SNPP Primary Productivity".format(
            args.sensor.upper())
        nc_file.title = ', '.join([
            "Primary Productivity",
            f"{args.sensor.upper()}",
            "Science Quality",
            "Global",
            "9km",
            f"2023-present (8 Day Composite)".format([args.sensor])
            ])
        nc_file.summary = ' '.join([
            f"The {args.sensor.upper()} Primary Productivity product ",
            "provides scence quality estimates of net carbon ",
            "fixation by phytoplankton using enhanced ",
            "satellite measurements with the ",
            "algorithm by Behrenfeld and Falkowski (1997). ",
            "This product incorporates SQ chlorophyll a and ",
            "Photosynthetically Available Radiation (PAR) data "
            "along with high-resolution SST data. ",
            "Mapped to a NASA 9km Standard Mapped Image for accuracy ",
            "in global assessments."
        ])
        nc_file.history = ' '.join([
            "Chlorophyll a, PAR, and SST satellite science ",
            "quality data were applied to the equation ",
            "of Behrenfeld and Falkowski, 1997"
        ])

        # Close netCDF file for this date
        nc_file.close()

        # Compress and archive
        myCmd = ' '.join(
            [
                'nccopy',
                '-d6',
                os.path.join(WORK_DIR, TEMP_OUT_FILE),
                os.path.join(WORK_DIR, nc_filename)
            ]
        )
        print('Compress ofile', subprocess.call(myCmd, shell=True))

        myCmd = ' '.join(
            [
                'mv',
                os.path.join(WORK_DIR, nc_filename),
                nc_file_path
            ]
        )
        print('Archive ofile', subprocess.call(myCmd, shell=True))

        # Print where netCDF files were saved
        print(f"NetCDF file '{nc_filename}' archived at {nc_file_path}")

        # Add code to send to ERDDAP

        # Move to next 8-day period
        current_date += timedelta(days=8)


if __name__ == '__main__':
    main()
