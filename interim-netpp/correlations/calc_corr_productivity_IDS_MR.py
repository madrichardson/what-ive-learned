"""Create correlations and p values between two time series using primary
productivity monthly products.

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


import xarray as xr
import pandas as pd
import numpy as np
import os
import argparse
import netCDF4
import numpy.ma as ma
import sys
import xskillscore as xs
from datetime import datetime
import subprocess

BASE_DIR = '/Users/madisonrichardson/netpp/monthly'
BIN_DIR = '/Users/madisonrichardson/netpp/bin'
DATA_DIR = os.path.join(BASE_DIR, 'data', 'monthly')
WORK_DIR = os.path.join(BASE_DIR, 'data', 'work')
RESOURCES_DIR = os.path.join(BASE_DIR, 'resources')
CORRELATIONS_DIR = os.path.join(DATA_DIR, 'correlations')

# Functions for processing data


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


# Run main function


def main():
    """Create correlations on monthly means of primary productivity products.

    Generate Pearson correlation coefficients and p-values of monthly
    primary productivity fields made from satellite data using the method of
    Behrenfeld and Falkowski 1997.
    Pearson correlation to measure the linear correlation between two time
    series (ref).
    Users set the two satellite sources, resolution, start/stop dates, dataset,
    data or anomaly time series, and the lat storage size (this only is used
    for the calculation of the correlations).
    """
    # -----------------------------------------------------------------------------
    # Start: Parser
    # -----------------------------------------------------------------------------
    # Parse arguments
    print('start main')
    help_txt = {
        'describe': 'Correlation analysis by ' +
        'setting two satellite sources, resolution, start/stop date, dataset, '
        'time series type, and number of values in time series',
        'resolution': 'Set the resolution of the data',
        'start': 'Date farthest in past, format YYmm as an int',
        'end': 'Most recent date, format YYmm as an int',
        'source1': '[F]irst source satellite data: modis, nasa',
        'source2': '[S]econd source satellite data: modis, nasa',
        'dataset': 'Dataset in primary productivity netcdf: productivity',
        'timeseries_type': 'Time series using data or anomalies',
        'percent_keep': 'Correlations for concurrent grid points with enough '
        'values: range between 0. to 1.0',
        'overwrite': 'Set to overwrite existing output files.',
        'ny_block': 'Set number of lat-grid points to store in data matrix'
    }

    doc_formatter = argparse.RawDescriptionHelpFormatter
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=doc_formatter)

    parser.add_argument('-a', '--startyear',
                        help=help_txt['start'],
                        type=int,
                        required=True)
    parser.add_argument('-z', '--endyear',
                        help=help_txt['end'],
                        type=int,
                        required=True)
    parser.add_argument('-f', '--source1',
                        help=help_txt['source1'],
                        choices=["snpp", "nasa", "noaa20", "modis"],
                        required=True)
    parser.add_argument('-s', '--source2',
                        help=help_txt['source2'],
                        choices=["snpp", "nasa", "noaa20", "modis"],
                        required=True)
    parser.add_argument('-d', '--dataset',
                        help=help_txt['dataset'],
                        choices=['productivity'],
                        required=True)
    parser.add_argument('-t', '--timeseries_type',
                        help=help_txt['timeseries_type'],
                        choices=['data', 'anom'],
                        required=True)
    parser.add_argument('-p', '--percent_keep',
                        help=help_txt['percent_keep'],
                        type=float,
                        required=True)
    parser.add_argument('-o', '--overwrite',
                        action='store_true',
                        help=help_txt['overwrite'])
    parser.add_argument('-n', '--ny_block',
                        choices=[2, 10, 20, 30, 40],
                        type=int,
                        help=help_txt['ny_block'],
                        required=True)

    args = parser.parse_args()
    startyear = args.startyear
    endyear = args.endyear
    source1 = args.source1
    source2 = args.source2
    ncvar_wnt = args.dataset
    timeseries_type = args.timeseries_type
    percent_keep = args.percent_keep
    overwrite = args.overwrite
    ny_block = args.ny_block

    # Validate parameters
    try:
        validate_params(startyear, endyear, percent_keep)
    except ValueError as e:
        print(f"Parameter validation error: {e}")
        sys.exit(1)

    # Parse and validate inputs
    syear, smonth, eyear, emonth = parse_dates(startyear, endyear)

    # Generate monthly date range from the start to the end
    time_bgn = np.datetime64('{}-{:02d}'.format(syear, smonth), 'M')
    time_end = np.datetime64('{}-{:02d}'.format(eyear, emonth), 'M')
    dtM = pd.to_datetime(
        np.arange(time_bgn, time_end+1, dtype='datetime64[M]')
    )
    ntM = len(dtM)
    yy = dtM.year.values
    mm = dtM.month.values

    # Verify directories exist
    DIR_LIST = [
        BASE_DIR,
        BIN_DIR,
        DATA_DIR,
        WORK_DIR,
        RESOURCES_DIR,
        CORRELATIONS_DIR
    ]

    for dr in DIR_LIST:
        os.makedirs(dr, exist_ok=True)
    print(len(DIR_LIST), 'directories validated')

    # Set up templates for the input and output files
    ifile_tmpl = 'productivity_month_{}_{}_9km.nc'
    ofile_tmpl = '{}_{}_corr_month_{}_{}_{}_{}_{}_{:03d}percent.nc'

    # Calculate correlations for ncvar variable wanted
    ncvar_list = [ncvar_wnt]
    for ncvar in ncvar_list:
        # Set up directories and output file parameters
        source_list = np.sort([source1, source2])
        source_comb = '{}_{}'.format(source_list[0], source_list[1])
        odir = os.path.join(
            CORRELATIONS_DIR,
            source_comb
        )
        odir_filelist = os.listdir(odir)

        # Generate output filenaeme
        ofile = ofile_tmpl.format(
            ncvar,
            timeseries_type,
            source_list[0],
            source_list[1],
            '9km',
            startyear,
            endyear,
            int(percent_keep*100)
        )

        # Skip if output file exists and overwrite is False
        if not overwrite and ofile in odir_filelist:
            print(f'{ofile} already exists')
            continue

        # Prepare ncgen input and output filenames
        now = datetime.now()
        ncgen_ofile_nc = f"ncgen_corr_ofile{now:%Y%m%d%H%M%S}.nc"
        ncgen_ifile_cdl = 'correlations_{}_month_{}_{}_{}.cdl'.format(
            ncvar,
            source_list[0],
            source_list[1],
            '9km'
        )

        # Run ncgen command
        myCmd1 = ' '.join(['ncgen',
                           '-o',
                           os.path.join(WORK_DIR, ncgen_ofile_nc),
                           os.path.join(RESOURCES_DIR, ncgen_ifile_cdl)
                           ])
        print('ncgen', subprocess.call(myCmd1, shell=True))

        # Initialize directory list for each source
        idir_list = []
        for ii in range(len(source_list)):
            idir_list.append(os.path.join(BASE_DIR, source_list[ii]))

        # Open productivity monthly files and place in list
        ds_list = []
        for ii in range(len(source_list)):
            ds_dates_list = []
            for i in range(ntM):
                ifile = ifile_tmpl.format(
                    source_list[ii], str(yy[i]) + f'{mm[i]:02}', '9km'
                )
                if not os.path.isfile(os.path.join(idir_list[ii], ifile)):
                    print('file not found for', ifile)
                    continue
                ds_dates_list.append(
                    xr.open_dataset(os.path.join(idir_list[ii], ifile))
                )
            ds_list.append(ds_dates_list)

        # Intialize matrices for correlations and p-values
        nt1, ny1, nx1 = ds_dates_list[0][ncvar_list[0]].shape
        corr_mtrx = np.zeros([ny1, nx1])*np.nan
        pval_mtrx = np.zeros([ny1, nx1])*np.nan
        n_mtrx = np.zeros([ny1, nx1])*np.nan

        # Reshape data into latitude blocks
        indx_block = reshape_data_block(
            data=ds_dates_list[0][ncvar].values,
            block_size=ny_block
        )
        num_block = indx_block.shape[0]

        for i in range(num_block):
            print(i)
            # Construct the two data matrix of shape [ntM X ny_block X nx1]
            injk = indx_block[i, :]
            data_block_mtrx = np.zeros([len(ds_list), ntM, ny_block, nx1])

            # Populate data_block_mtrx with data from each source and time step
            for j in range(ntM):
                for k in range(len(ds_list)):
                    data_block_mtrx[k, j, :, :] = (
                        ds_list[k][j][ncvar][0, injk, :].data
                    )

            # Compute correlations for each latitude subset in the block
            for j in range(ny_block):
                data_mtrx = data_block_mtrx[:, :, j, :]
                # Only keep grid points that have non-missing numbers
                # above 'percent_keep'
                ones_mtrx = data_mtrx/data_mtrx
                ones_mtrx12 = np.sum(ones_mtrx, axis=0)/2
                sum_one_mtrx = np.nansum(ones_mtrx12, axis=0)
                in_keep = np.where(sum_one_mtrx > ntM*percent_keep)[0]

                # Create data array with time as coordinate, useful for
                # Calculating anoms in xarray
                time_series_data = xr.DataArray(
                    data_mtrx[:, :, in_keep],
                    coords=[
                        source_list,
                        dtM.astype('datetime64[ns]'),
                        in_keep
                    ],
                    dims=['source', 'time', 'in_keep'])

                # Correlations for either data or anom time series
                if timeseries_type == 'anom':
                    # Calculate monthly climatology
                    monthly_climatology = (
                        time_series_data.groupby('time.month').mean('time')
                    )
                    # Get anomalies by subtracting climatology from the data
                    compare_data = (
                        time_series_data.groupby('time.month') - monthly_climatology
                    )
                elif timeseries_type == 'data':
                    # Use raw data directly for correlation
                    compare_data = time_series_data

                # Calculate Pearson correlation between sources
                correlation = xr.corr(
                    compare_data.sel(source=source_list[0]),
                    compare_data.sel(source=source_list[1]),
                    dim='time'
                )

                # Use xskillscore to get p values
                p_value = xs.pearson_r_p_value(
                    compare_data.sel(source=source_list[0]),
                    compare_data.sel(source=source_list[1]),
                    dim='time',
                    skipna=True
                )

                # Place corr, pval values for latitude subset in final
                # global data matrix
                corr_mtrx[injk[j], in_keep] = correlation.data
                pval_mtrx[injk[j], in_keep] = p_value.data
                # Store the count of valid data points
                n_mtrx[injk[j], in_keep] = sum_one_mtrx[in_keep]

        # Mask invalid values in correlation, p-value,
        # and count matrices (not sure if need)
        corr_mtrx = ma.masked_invalid(corr_mtrx)
        pval_mtrx = ma.masked_invalid(pval_mtrx)
        n_mtrx = ma.masked_invalid(n_mtrx)

        # Labels and matrices to save in NetCDF file
        cpn_lbl = ['corr', 'pval', 'n']
        cpn_list = [corr_mtrx, pval_mtrx, n_mtrx]

        # Open temporary file and load data into it
        with netCDF4.Dataset(
            os.path.join(WORK_DIR, ncgen_ofile_nc), 'a'
        ) as nc:
            # Set global attribute
            nc.title = (
                f"Pearson correlation coefficients and p-values of "
                f"monthly primary productivity fields between "
                f"{source1.upper()} and {source2.upper()}"
                )
            nc.summary = (
                f"Correlations between primary productivity or "
                f"PAR, chlorophyl and SST from {source1.upper()} "
                f"and {source2.upper()}. These are 9km products "
                f"generated from time series of monthly means. "
                f"See Melin et al 2017 for more details")
            nc.source = f"{source1.upper()} and {source2.upper()}"
            nc.instrument = f"{source1.upper()} and {source2.upper()}"
            nc.id = f"correlations_{source1}_{source2}_monthly_9km"
            nc.platform = f"{source1.upper()} and {source2.upper()}"

            # Place corr, pval, n in nc
            for j in range(len(cpn_lbl)):
                nc['{}'.format(cpn_lbl[j])][0, :, :] = cpn_list[j]

        # Compress and save the temporary file to the final file name
        myCmd = ' '.join(['nccopy',
                          '-d6',
                          os.path.join(WORK_DIR, ncgen_ofile_nc),
                          os.path.join(odir, ofile)
                          ])
        print('nccopy', subprocess.call(myCmd, shell=True))
        print('Done with', ofile)

        # Clean up temporary file after saving the final output
        os.remove(os.path.join(WORK_DIR, ncgen_ofile_nc))


if __name__ == "__main__":
    main()
