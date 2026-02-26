"""Calculating Correlations and P-Values Between Legacy and Interim NetPP
Products

This script generates correlations between two time series of monthly means of
primary productivity described by Behrenfeld and Falkowski 1997. The
correlations are Pearson correlations. Correlations can also be calculated for
dataset of PAR, chlorophyll, SST, etc.

Users set the satellite input data, start and stop years,
dataset, data or anomaly time series via command-line arguments.

Users set if correlations are calculated for time series that have seasonal
cycles or anomalies.

The cdl files have latitude and longitude that have been gridded to the
9km NOAA CoastWatch standard ocean color grid. The template file is
populated with the latitude and longitude data and then renamed.

The script can be repurposed for different input data with these modifications:
    * Gridding any input chlorophyll, PAR, and SST data to a common grid
    * Creating an appropriate cdl file for the dataset
    * Adjusting the logic to make directory paths
    * Adjust output file naming and input file search pattern.
"""

# Import necessary packages
import argparse
import xarray as xr
import xskillscore as xs
import numpy as np
import numpy.ma as ma
import pandas as pd
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
WORK_DIR = os.path.join(BASE_DIR, "work")
RESOURCES_DIR = os.path.join(BASE_DIR, "resources")
CORRELATIONS_DIR = os.path.join(BASE_DIR, "data", "correlations")

# Create some functions


def validate_params(startyear, endyear, percent_keep):
    """
    Validates the inputted parameters to ensure
    correct ranges and logical order.

    Args:
        startyear (str): The start year in 'YYYYMM' format.
        endyear (str): The end year in 'YYYYMM' format.
        percent_keep (float): Percentage of grid points to keep (0 to 1)

    Raises:
        ValueError: The startyear is not earlier than endyear.
        ValueError: The percent_keep is outside the range (0 to 1).
    """
    if startyear > endyear:
        raise ValueError("The start year must be earlier than the end year.")
    if not (0.0 <= percent_keep <= 1.0):
        raise ValueError("percent_keep must be between 0 and 1")


def parse_dates(startyear, endyear):
    """
    Parses the start and end years in 'YYYYMM' format
    and extracts the corresponding years and months as
    separate integer values.

    Args:
        startyear (str): The start date in 'YYYYMM' format.
        endyear (str): The end date in 'YYYYMM' format.

    Returns:
        tuple: A tuple containing four integers:
        (syear, smonth, eyear, emonth).
            - syear (int): The year component of the start date.
            - smonth (int): The month component of the start date.
            - eyear (int): The year component of the end date.
            - emonth (int): The monht component of the end date.
    """
    syear = int(np.floor(int(startyear) / 100))
    smonth = int(startyear) - syear * 100
    eyear = int(np.floor(int(endyear) / 100))
    emonth = int(endyear) - eyear * 100
    return syear, smonth, eyear, emonth


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


def reshape_data_block(data, block_size):
    """
    Reshapes a 2D array's latitude dimension into
    smaller blocks for efficient processing.

    Args:
        data (np.ndarray): A 2D array with shape (time, latitude, longitude)
        where the latitude dimension (axis=1) will be divide into blocks.
        block_size (int): The number of latitude grid points in each block.

    Returns:
        np.ndarray: A 2D array where the latitude indices are reshaped
        into blocks of the specified size, with shape (num_blocks, block_size).
    """
    ny_index = np.arange(data.shape[1])
    num_blocks = int(data.shape[1] / block_size)
    return np.reshape(ny_index, [num_blocks, block_size])


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
        help="Start date in YYYYMM format",
    )
    parser.add_argument(
        "-z",
        "--end",
        type=str,
        required=True,
        help="End date in YYYYMM format",
    )
    parser.add_argument(
        "-s1",
        "--source1",
        type=str,
        required=True,
        choices=["modis", "snpp", "noaa20"],
        help="first sensor to calculate correlations for",
    )
    parser.add_argument(
        "-s2",
        "--source2",
        type=str,
        required=True,
        choices=["modis", "snpp", "noaa20"],
        help="second sensor to calculate correlations for",
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
    parser.add_argument(
        "-ny",
        "--ny_block",
        type=str,
        required=True,
        help="number of latitude grid points per block processing",
    )

    args = parser.parse_args()

    # Parse the start and end dates
    startyear = parse(args.start).strftime("%Y%m")
    endyear = parse(args.end).strftime("%Y%m")

    # Validate parameters
    try:
        validate_params(args.startyear, args.endyear, args.percent_keep)
    except ValueError as e:
        print(f"Parameter validation error: {e}")
        sys.exit(1)

    # Parse and validate inputs
    syear, smonth, eyear, emonth = parse_dates(startyear, endyear)

    # Generate monthly date range from the start to the end
    time_bgn = np.datetime64("{}-{:02d}".format(syear, smonth), "M")
    time_end = np.datetime64("{}-{:02d}".format(eyear, emonth), "M")
    dtM = pd.to_datetime(
        np.arange(time_bgn, time_end + 1, dtype="datetime64[M]")
    )
    ntM = len(dtM)
    yy = dtM.year.values
    mm = dtM.month.values

    print(f"Generated {ntM} monthly time points from {time_bgn} to {time_end}")

    # Verify directories exist
    DIR_LIST = [BASE_DIR, WORK_DIR, RESOURCES_DIR, CORRELATIONS_DIR]

    for dr in DIR_LIST:
        os.makedirs(dr, exist_ok=True)
    print(len(DIR_LIST), "directories validated")

    # Set up the template output file
    ofile_tmpl = "{}_{}_corr_month_{}_{}_{}_{}_{}_{:03d}percent.nc"

    # Calculate correlations for ncvar variable wanted
    ncvar_list = [args.ncvar]
    for args.ncvar in ncvar_list:
        # Set up directories and output file parameters
        source_list = np.sort([args.source1, args.source2])
        source_comb = "{}_{}".format(source_list[0], source_list[1])
        odir = os.path.join(CORRELATIONS_DIR, source_comb)
        odir_filelist = os.listdir(odir)

        # Generate output filenaeme
        ofile = ofile_tmpl.format(
            args.ncvar,
            args.timeseries_type,
            source_list[0],
            source_list[1],
            "9km",
            startyear,
            endyear,
            int(args.percent_keep * 100),
        )

        # Skip if output file exists and overwrite is False
        if not args.overwrite and ofile in odir_filelist:
            print(f"{ofile} already exists")
            continue

    # Prepare ncgen input and output filenames
    now = datetime.now()
    ncgen_ofile_nc = f"ncgen_corr_ofile{now:%Y%m%d%H%M%S}.nc"
    ncgen_ifile_cdl = "correlations_{}_month_{}_{}_{}.cdl".format(
        args.ncvar, source_list[0], source_list[1], "9km"
    )

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

    # For ERDDAP function
    erddap_url = "http://localhost:8080/erddap"
    modis_id = "productivity_modis_aqua_monthly"
    noaa20_id = "productivity_viirs_noaa20_monthly"

    modis_ds = xr_open_ds(
        e_source=erddap_url, e_id=modis_id, var_name=args.ncvar
    )
    print("MODIS Dataset:", modis_ds)

    noaa20_ds = xr_open_ds(
        e_source=erddap_url, e_id=noaa20_id, var_name=args.ncvar
    )
    print("NOAA20 Dataset:", noaa20_ds)

    # Intialize an empty list to store datasets
    ds_list = []

    # ERDDAP Datasets
    datasets = [modis_ds, noaa20_ds]

    for ds in datasets:
        ds_dates_list = []  # To store individual monthly slices
        for i in range(ntM):
            # Select a single time slice
            start_date = f"{yy[i]}-{mm[i]:02d}-16"
            try:
                # Extract data for the specific time slice
                time_slice = ds.sel(time=start_date)
                ds_dates_list.append(time_slice)
            except KeyError as e:
                print(f"Time {start_date} not found in dataset: {e}")
                continue

        # Append the time series for the dataset
        ds_list.append(ds_dates_list)

    print(f"Loaded data for {len(ds_list)} datasets.")

    # Intialize matrices for correlations and p-values
    _, ny1, nx1 = ds_list[0][0].shape
    nt1 = len(ds_list[0])

    print(f"Temporal Dimension (nt1): {nt1}")
    print(f"Latitude Dimension (ny1): {ny1}, Longitude dimensions (nx1): {nx1}")

    # Intialize matrices with NaN values
    corr_mtrx = np.zeros([ny1, nx1]) * np.nan
    pval_mtrx = np.zeros([ny1, nx1]) * np.nan
    n_mtrx = np.zeros([ny1, nx1]) * np.nan

    # Reshape data into latitude blocks
    indx_block = reshape_data_block(
        data=ds_dates_list[0].values, block_size=args.ny_block
    )

    # Calculate the number of latitude blocks
    num_block = indx_block.shape[0]

    print(f"Latitude divided into {num_block} blocks of size {args.ny_block}")

    for i in range(num_block):
        print(f"Processing block {i + 1}/{num_block}")
        # Construct the two data matrix of shape [ntM X ny_block X nx1]
        injk = indx_block[i, :]
        data_block_mtrx = np.zeros([len(ds_list), ntM, args.ny_block, nx1])

        # Populate data_block_mtrx with data from each source and time step
        for j in range(ntM):
            for k in range(len(ds_list)):
                data_block_mtrx[k, j, :, :] = ds_list[k][j].values[0, injk, :]

        # Compute correlations for each latitude subset in the block
        for j in range(args.ny_block):
            data_mtrx = data_block_mtrx[:, :, j, :]

            # Only keep grid points that have non-missing numbers above 'percent_keep'
            ones_mtrx = data_mtrx / data_mtrx
            ones_mtrx_comb = np.sum(ones_mtrx, axis=0) / 2
            sum_one_mtrx = np.nansum(ones_mtrx_comb, axis=0)
            in_keep = np.where(sum_one_mtrx > ntM * args.percent_keep)[0]

            # Create data array with time as coordinate, useful for
            # Calculating anoms in xarray
            time_series_data = xr.DataArray(
                data_mtrx[:, :, in_keep],
                coords=[source_list, dtM.astype("datetime64[ns]"), in_keep],
                dims=["source", "time", "in_keep"],
            )

            # Correlations for either data or anom time series
            if args.timeseries_type == "anom":
                # Calculate monthly climatology
                monthly_climatology = time_series_data.groupby(
                    "time.month"
                ).mean("time")
                # Get anomalies by subtracting climatology from the data
                compare_data = (
                    time_series_data.groupby("time.month") - monthly_climatology
                )
            elif args.timeseries_type == "data":
                # Use raw data directly for correlation
                compare_data = time_series_data

            # Calculate Pearson correlation between sources
            correlation = xr.corr(
                compare_data.sel(source=source_list[0]),
                compare_data.sel(source=source_list[1]),
                dim="time",
            )

            # Use xskillscore to get p values
            p_value = xs.pearson_r_p_value(
                compare_data.sel(source=source_list[0]),
                compare_data.sel(source=source_list[1]),
                dim="time",
                skipna=True,
            )

            # Place corr, pval values for latitude subset in final global data matrix
            corr_mtrx[injk[j], in_keep] = correlation.data
            pval_mtrx[injk[j], in_keep] = p_value.data
            # Store the count of valid data points
            n_mtrx[injk[j], in_keep] = sum_one_mtrx[in_keep]

    print("Correlations and p-values computed.")
    print(
        f"Min: {np.nanmin(corr_mtrx)}, Max: {np.nanmax(corr_mtrx)}, Mean: {np.nanmean(corr_mtrx)}"
    )
    print(
        f"Min: {np.nanmin(pval_mtrx)}, Max: {np.nanmax(pval_mtrx)}, Mean: {np.nanmean(pval_mtrx)}"
    )
    print(
        f"Min: {np.nanmin(n_mtrx)}, Max: {np.nanmax(n_mtrx)}, Mean: {np.nanmean(n_mtrx)}"
    )

    # Mask invalid values in correlation, p-value, and count matrices
    corr_mtrx = ma.masked_invalid(corr_mtrx)
    pval_mtrx = ma.masked_invalid(pval_mtrx)
    n_mtrx = ma.masked_invalid(n_mtrx)

    # Labels and matrices to save in NetCDF file
    cpn_lbl = ["corr", "pval", "n"]
    cpn_list = [corr_mtrx, pval_mtrx, n_mtrx]

    # Save results to NetCDF file
    print("Saving the results to a NetCDF file")

    # Open temporary file and load data into it
    with netCDF4.Dataset(os.path.join(WORK_DIR, ncgen_ofile_nc), "a") as nc:
        # Set global attributes
        nc.acknowledgement = (
            "The project was supported by funding from the "
            "Portfolio Management Branch of NESDIS and NOAA CoastWatch."
        )
        nc.contributors = (
            "Dale Robinson, Isaac Shroeder, Ryan Vandermeulen, "
            "Jonathan Sherman, Jesse Espinoza, & Madison Richardson"
        )
        nc.title = (
            f"Pearson correlation coefficients and p-values of "
            f"monthly primary productivity fields between "
            f"{args.source1.upper()} and {args.source2.upper()}"
        )
        nc.summary = (
            f"Correlations between primary productivity or "
            f"PAR, chlorophyl and SST from {args.source1.upper()} "
            f"and {args.source2.upper()}. These are 9km products "
            f"generated from time series of monthly means. "
            f"See Melin et al 2017 for more details"
        )
        nc.source = f"{args.source1.upper()} and {args.source2.upper()}"
        nc.instrument = f"{args.source1.upper()} and {args.source2.upper()}"
        nc.id = f"correlations_{args.source1}_{args.source2}_monthly_9km"
        nc.platform = f"{args.source1.upper()} and {args.source2.upper()}"

        # Place corr, pval, n in nc
        for j in range(len(cpn_lbl)):
            nc["{}".format(cpn_lbl[j])][0, :, :] = cpn_list[j]

    # Compress and save the temporary file to the final file name
    myCmd = " ".join(
        [
            "nccopy",
            "-d6",
            os.path.join(WORK_DIR, ncgen_ofile_nc),
            os.path.join(odir, ofile),
        ]
    )
    print("nccopy", subprocess.call(myCmd, shell=True))
    print("Done with", ofile)

    # Clean up temporary file after saving the final output
    os.remove(os.path.join(WORK_DIR, ncgen_ofile_nc))
