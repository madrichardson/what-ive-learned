"""Calculating Trends and P-Values Between Legacy and Interim NetPP Products

This script generates trends for time series of monthly means of primary
productivity described by Behrenfeld and Falkowski 1997. The trend analysis
follows methods outlined in Melin et al 2017, see the section 2.3 "Trend
estimates and comparison of trends". Trends can also be calculated for
dataset of PAR, chlorophyll, SST, etc.

Users set the satellite input data, resolution, start and stop years,
dataset, data or anomaly time series via command-line arguments.

Users set if trends are calculated for time series that have seasonal
cycles or anomalies.

The output files are created by generating a template netCDF file from a cdl
file that is prepopulate with metadata for pairing the following global,
daily input data:
    * SNPP-VIIRS chlorophyll and PAR, plus ACSPO Gridded Super-collated SST
    * NOAA20-VIIRS chlorophyll and PAR, plus ACSPO Gridded Super-collated SST
    * MODIS Aqua chlorophyll, PAR, and SST

The cdl files have latitude and longitude that have been gridded to the
9km NASA Standard Mapped Image (SMI). The template file is
populated with the primary productivity data and then renamed.

The script can be repurposed for different input data with these modifications:
    * Gridding any input chlorophyll, PAR, and SST data to a common grid
    * Creating an appropriate cdl file for the dataset
    * Adjusting the logic to make directory paths
    * Adjust output file naming and input file search pattern.
"""

# Import necessary packages
import argparse
import xarray as xr
import numpy as np
import numpy.ma as ma
from scipy import stats
import os
import sys
import subprocess
from datetime import datetime
import netCDF4
from dateutil.parser import parse
import warnings

warnings.filterwarnings("ignore")

# Create your global variables
BASE_DIR = "/Users/madisonrichardson/netpp"
ODATA_DIR = os.path.join(BASE_DIR, "data", "trends")
WORK_DIR = os.path.join(BASE_DIR, "work")
RESOURCES_DIR = os.path.join(BASE_DIR, "resources")

# Create functions


def validate_params(start_date, end_date, percent_keep):
    """
    Validates the inputted parameters to ensure
    correct ranges and logical order.

    Args:
        start_date (str): The start year in 'YYYYMM' format.
        end_date (str): The end year in 'YYYYMM' format.
        percent_keep (float): Percentage of grid points to keep (0 to 1)

    Raises:
        ValueError: The start_date is not earlier than end_date.
        ValueError: The percent_keep is outside the range (0 to 1).
    """
    if start_date > end_date:
        raise ValueError("The start year must be earlier than the end year.")
    if not (0.0 <= percent_keep <= 1.0):
        raise ValueError("percent_keep must be between 0 and 1")


def xr_open_ds(
    e_id, e_source="http://localhost:8080/erddap", dap="griddap", var_name=None
):
    """
    Open a remote ERDDAP dataset by constructing
    the full URL for the specified dataset and opens it as
    an xarray Dataset. If a variable name is provided, it attempts
    to extract and return only that variable as an xarray
    DataArray.

    Args:
        e_id (str): The dataset ID for the dataset on the ERDDAP server.
        e_source (str): The base URL of the ERDDAP server. Defaults to
        'http://localhost:8080/erddap'.
        dap (str): The data acess protocol to us (e.g., 'griddap').
        Defaults to 'griddap'.
        var_name (str, optional): The name of the variable to extract from
        the dataset. If None, the entire dataset is returned.

    Returns:
        xarray.Dataset or xarray.DataArray:
            - The entire dataset if 'var_name' is None.
            - The specified variable as a DataArray if 'var_name' is returned.

    Raises:
        KeyError: If the specified variable name does not exist in the dataset.
    """
    # remove any trailing /
    e_source = e_source.rstrip("/")

    erddap_url = "/".join([e_source, dap, e_id])

    ds = xr.open_dataset(erddap_url)

    if var_name:
        if var_name in ds:
            return ds[var_name]
        else:
            print(f"Variable '{var_name}' not found in the dataset.")
            print("Available variables:")
            print(list(ds.variables.keys()))
            raise KeyError(f"Variable  '{var_name}' not found.")

    return ds


def main():
    """Run main function."""

    # Main code for calculating trends for a primary productivity product.
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
        "--source",
        type=str,
        required=True,
        choices=["modis", "snpp", "noaa20"],
        help="sensor to calculate trends for",
    )
    parser.add_argument(
        "-n",
        "--ncvar",
        type=str,
        required=True,
        help="variable to calculate trends from",
    )
    parser.add_argument(
        "-tt",
        "--timeseries_type",
        type=str,
        required=True,
        choices=["data", "anom"],
        help="timeseries type for calculating trends",
    )
    parser.add_argument(
        "-pk",
        "--percent_keep",
        type=float,
        required=True,
        help="percentage of valid points to keep",
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
    start_date = parse(args.start).strftime("%Y-%m-%d")
    end_date = parse(args.end).strftime("%Y-%m-%d")

    # Validate parameters
    try:
        validate_params(start_date, end_date, args.percent_keep)
    except ValueError as e:
        print(f"Parameter validation error: {e}")
        sys.exit(1)

    print(f"Processing data from {start_date} to {end_date} for {args.source}")

    # Put the directories in a list
    DIR_LIST = [
        BASE_DIR,
        ODATA_DIR,
        WORK_DIR,
        RESOURCES_DIR,
    ]

    # Ensure directories exist
    for dr in DIR_LIST:
        os.makedirs(dr, exist_ok=True)
    print(len(DIR_LIST), "directories validated")

    # Set up templates for the output files
    ofile_tmpl = "{}_{}_trend_month_{}_{}_{}_{}_{:03d}percent.nc"

    # List files in output directory
    odir_filelist = os.listdir(ODATA_DIR)

    # Generate output filename
    ofile = ofile_tmpl.format(
        args.ncvar,
        args.timeseries_type,
        args.source,
        "9km",
        int(datetime.strptime(start_date, "%Y-%m-%d").year),
        int(datetime.strptime(end_date, "%Y-%m-%d").year),
        int(args.percent_keep * 100),
    )

    # Skip if output file exists and overwrite is False
    if not args.overwrite and ofile in odir_filelist:
        print(f"{ofile} already exists")
    else:
        print("Creating output file...")

    # Prepare ncgen input and output filenames
    now = datetime.now()
    ncgen_ofile_nc = f"ncgen_trend_ofile{now:%Y%m%d%H%M%S}.nc"

    # Create temporary file to accept output data from a .cdl file
    # cdl is unique for each ncvar
    ncgen_ifile_cdl = f"trends_{args.source}.cdl"

    # Run ncgen command
    myCmd1 = " ".join(
        [
            "ncgen",
            "-o",
            os.path.join(WORK_DIR, ncgen_ofile_nc),
            os.path.join(RESOURCES_DIR, ncgen_ifile_cdl),
        ]
    )
    print(myCmd1)
    print("ncgen", subprocess.call(myCmd1, shell=True))

    # Load ERDDAP dataset
    erddap_url = "http://localhost:8080/erddap"
    e_id = "productivity_viirs_noaa20_monthly"

    print(f"Fetching dataset {e_id} from {erddap_url}")

    ds = xr_open_ds(e_source=erddap_url, e_id=e_id, var_name=args.ncvar)
    print(f"Dataset {e_id} loaded successfully: {ds}")

    # Subset the data to the start and end date
    da = ds.sel(time=slice(start_date, end_date))

    # Use Dask's chunking to quicken compuations
    da = da.chunk({"time": 10, "latitude": 500, "longitude": 500})

    # Filter valid pixels based on percent_keep
    print(
        f"Filtering valid data points using percent_keep = {args.percent_keep}"
    )
    valid_pixels = da.notnull().sum("time") >= (
        args.percent_keep * da.time.size
    )
    da_filtered = da.where(valid_pixels, np.nan)

    # Count valid time points
    n = da_filtered.count(dim="time")
    n = n.where(n > 0)
    print("Unique n values after filtering:", np.unique(n))

    # Convert to hours since 1970-01-01
    print("Convert to hours since 1970-01-01")
    hours = (
        da_filtered["time"].astype("datetime64[h]")
        - np.datetime64("1970-01-01T00:00:00")
    ).astype(int)

    # Compute anomalies if timeseries_type = "anom"
    if args.timeseries_type == "anom":
        print("Calculating with anomalies")
        monthly_climatology = da_filtered.groupby("time.month").mean("time")
        anomalies = da_filtered.groupby("time.month") - monthly_climatology
        anomalies = anomalies.assign_coords(hours=hours)
        anomalies = anomalies.swap_dims({"time": "hours"})
        da_final = anomalies
    else:
        da_final = da_filtered.swap_dims({"time": "hours"})
        print("Calculating with raw data")

    # Fit linear regression
    print("Fitting a linear regression")
    results = da_final.polyfit("hours", deg=1, cov=True)

    # Extract the slope of the regression
    print("Extracting the slope of the regression")
    slope = ma.masked_invalid(
        results.polyfit_coefficients.sel(degree=1).values
    )

    # Compute variance of the time variable
    print("Compute variance of the time variable")
    sigma_X_sq = np.square(da_final.hours.astype(float).std())

    # Calculate standard error
    print("Calculate standard error")
    se = np.sqrt(
        (1 / (n - 2))
        * ((np.square(da_final.std("hours")) / sigma_X_sq) - np.square(slope))
    )

    # Compute the t-statistic
    print("Compute the t-statistic")
    t = ma.masked_invalid(slope / se)

    # Compute p-value
    print("Compute the p-values")
    p = ma.masked_invalid(stats.t.sf(np.abs(t), df=n - 2) * 2)

    # Save results to NetCDF file
    print("Saving the results to a NetCDF file")
    # Define platform name based on source
    if args.source.lower() == "modis":
        platform_name = "MODIS-Aqua"
    elif args.source.lower() == "noaa20":
        platform_name = "VIIRS-NOAA20"
    elif args.source.lower() == "snpp":
        platform_name = "VIIRS-SNPP"
    else:
        platform_name = args.source.upper()

    # Open temporary file and load data into it
    with netCDF4.Dataset(os.path.join(WORK_DIR, ncgen_ofile_nc), "a") as nc:

        # Write in data
        nc["beta"][0, :, :] = slope
        nc["pval"][0, :, :] = p
        nc["n"][0, :, :] = n

        # Set global attributes

        nc.acknowledgement = (
            "The project was supported by funding from the "
            "Portfolio Management Branch of NESDIS and NOAA CoastWatch."
        )
        nc.contributors = (
            "Dale Robinson, Isaac Shroeder, Ryan Vandermeulen, "
            "Jonathan Sherman, Jesse Espinoza, & Madison Richardson"
        )
        nc.date_created = now.isoformat("T", "seconds")
        nc.instrument = f"{args.source.upper()}"
        nc.id = f"trends_{args.source}_monthly_9km"
        nc.platform = f"{args.source.upper}"
        nc.source = f"{args.source.upper()}"
        nc.title = (
            f"Trend coefficients and p values for monthly "
            f"primary productivity from {platform_name} globally at 9km "
            "resolution"
        )
        nc.summary = (
            f"Trends between primary productivity or PAR, "
            f"chlorophyll and SST from {platform_name}. These "
            f"are 9km products generated from time series of monthly "
            f"means. Trends are slopes of linear regressions and the level "
            f"of significance of the trend are computed from a t-test. "
            f"See Melin et al 2017 for more details."
        )

    # Compress and save the temporary file to the final file name
    myCmd = " ".join(
        [
            "nccopy",
            "-d6",
            os.path.join(WORK_DIR, ncgen_ofile_nc),
            os.path.join(ODATA_DIR, ofile),
        ]
    )
    print("nccopy", subprocess.call(myCmd, shell=True))
    print("Processing complete. NetCDF generated:", ofile)

    # Clean up up temporary file after saving the final output
    os.remove(os.path.join(WORK_DIR, ncgen_ofile_nc))


if __name__ == "__main__":
    main()
