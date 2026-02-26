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
        'setting two satellite sources, resolution, start/stop date, dataset, time series type, and number of values in time series',
        'resolution': 'Set the resolution of the data',
        'start': 'Date farthest in past, format YYmm as an int',
        'end': 'Most recent date, format YYmm as an int',
        'source1': '[F]irst source satellite data: modis, nasa',
        'source2': '[S]econd source satellite data: modis, nasa',
        'dataset': 'Dataset in primary productivity netcdf: productivity',
        'timeseries_type': 'Time series using data or anomalies',
        'percent_keep': 'Correlations for concurrent grid points with enough values: range between 0. to 1.0',
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
    prcnt_keep = args.percent_keep
    overwrite = args.overwrite
    ny_block = args.ny_block
    # -----------------------------------------------------------------------------
    # END: Parser
    # -----------------------------------------------------------------------------

    # .1) test the input args
    if startyear > endyear:
        print('start_date, -a, must be earlier than end_date, -z')
        sys.exit()

    if prcnt_keep < 0.0 or prcnt_keep > 1.0:
        print('percent_keep, -p, must be a float between 0 and 1')
        sys.exit()

    # .2) setup stuff
    # templates for the output file
    ifile_tmpl = 'productivity_month_{}_{}_9km.nc'
    ofile_tmpl = '{}_{}_corr_month_{}_{}_{}_{}_{}_{:03d}percent.nc'

    # calculate correlations for ncvar variable wanted
    ncvar_list = [ncvar_wnt]
    for ncvar in ncvar_list:
        # .3) only calculate the correlations if not already done
        # list files in output directory
        source_list = np.sort([source1, source2])
        source12 = '{}_{}'.format(source_list[0], source_list[1])

        odir = os.path.join(DATA_DIR,
                            'correlations',
                            source12,
                            )
        odir_filelist = os.listdir(odir)

        # generate the name of output file
        idir_list = []
        for ii in range(len(source_list)):
            idir_list.append(os.path.join(BASE_DIR, source_list[ii]))

        ofile = ofile_tmpl.format(ncvar, timeseries_type, source_list[0], source_list[1], '9km', startyear, endyear, int(prcnt_keep*100))

        # .4) input and output files for ncgen
        ncgen_ofile_nc_tmpl = 'ncgen_corr_ofile{}.nc'
        now = datetime.now()
        ncgen_str_time = '{0:%Y%m%d%H%M%S}'.format(now)
        ncgen_ofile_nc = ncgen_ofile_nc_tmpl.format(ncgen_str_time)

        # create temperary file to accept output data from a .cdl file
        # cdl is unique for each ncvar
        ncgen_ifile_cdl = 'correlations_{}_month_{}_{}_{}.cdl'.format(ncvar, source_list[0], source_list[1], '9km')

        myCmd1 = ' '.join(['ncgen',
                           '-o',
                           os.path.join(WORK_DIR, ncgen_ofile_nc),
                           os.path.join(RESOURCES_DIR, ncgen_ifile_cdl)
                           ])
        print('ncgen', subprocess.call(myCmd1, shell=True))

        # .5) calculate the correlations
        # first check to see if file exists, if overwrite=False then skip the calculation
        if not overwrite:
            if ofile in odir_filelist:
                print(ofile, 'already exists')
                continue

        # get start/end year and month
        yy1 = int(np.floor(startyear/100))
        mm1 = startyear - yy1*100
        yy2 = int(np.floor(endyear/100))
        mm2 = endyear - yy2*100

        time_bgn = np.datetime64('{}-{:02d}'.format(yy1, mm1), 'M')
        time_end = np.datetime64('{}-{:02d}'.format(yy2, mm2), 'M')
        dtM = pd.to_datetime(np.arange(time_bgn, time_end+1, dtype='datetime64[M]'))
        ntM = len(dtM)

        yy = dtM.year.values
        mm = dtM.month.values

        # open productivity monthly files and place in list
        ds_list = []
        for ii in range(len(source_list)):
            ds_dates_list = []
            for i in range(ntM):
                ifile = ifile_tmpl.format(source_list[ii], str(yy[i]) + f'{mm[i]:02}', '9km')
                if not os.path.isfile(os.path.join(idir_list[ii], ifile)):
                    print('file not found for', ifile)
                    continue
                ds_dates_list.append(xr.open_dataset(os.path.join(idir_list[ii], ifile)))
            ds_list.append(ds_dates_list)

        # size of data for calculating correlations, since they have same size get
        # the size from source1
        nt1, ny1, nx1 = ds_dates_list[0][ncvar_list[0]].shape
        ntM_list = len(ds_list[0])

        # initialize correlations and p values matrix
        corr_mtrx = np.zeros([ny1, nx1])*np.nan
        pval_mtrx = np.zeros([ny1, nx1])*np.nan
        n_mtrx = np.zeros([ny1, nx1])*np.nan

        # correlations need to done on two matrix with shape [time X lat X lon]
        # the global domain and all dates is too large
        # break the storage into latitude "blocks" of size "ny_block"
        ny_index = np.arange(0, ny1)
        num_block = int(ny1/ny_block)
        indx_block = np.reshape(ny_index, [num_block, ny_block])

        for i in range(num_block):
            print(i)
            # construct the two data matrix of shape [ntM X ny_block X nx1]
            injk = indx_block[i, :]
            data_block_mtrx = np.zeros([len(ds_list), ntM_list, ny_block, nx1])
            for j in range(ntM_list):
                for k in range(len(ds_list)):
                    data_block_mtrx[k, j, :, :] = ds_list[k][j][ncvar][0, injk, :].data

            # find the correlation
            for j in range(ny_block):
                data_mtrx = data_block_mtrx[:, :, j, :]
                # only keep grid points that have non-missing numbers above 'prcnt_keep'
                ones_mtrx = data_mtrx/data_mtrx
                ones_mtrx12 = np.sum(ones_mtrx, axis=0)/2
                sum_one_mtrx = np.nansum(ones_mtrx12, axis=0)
                in_keep = np.where(sum_one_mtrx > ntM_list*prcnt_keep)[0]

                # create datarray with time as coordinate, useful for
                # calculating anoms in xarray
                da1_time = xr.DataArray(
                    data_mtrx[:, :, in_keep],
                    coords=[source_list, dtM.astype('datetime64[ns]'), in_keep],
                    dims=['source', 'time', 'in_keep'])

                # correlaitons for either data or anom time series
                if timeseries_type == 'anom':
                    # anoms are created from subtracting long-term clim from data
                    da1_clim = da1_time.groupby('time.month').mean('time')
                    da1_compare = da1_time.groupby('time.month') - da1_clim
                    # da1_compare = da1_anom
                elif timeseries_type == 'data':
                    da1_compare = da1_time

                # pearson correlation
                da1_corr = xr.corr(da1_compare.sel(source=source_list[0]), da1_compare.sel(source=source_list[1]), dim='time')

                # use xskillscore to get p values
                da1_pval = xs.pearson_r_p_value(da1_compare.sel(source=source_list[0]), da1_compare.sel(source=source_list[1]), dim='time', skipna=True)

                # place corr, pval values for latitude subset in final global data matrix
                corr_mtrx[injk[j], in_keep] = da1_corr.data
                pval_mtrx[injk[j], in_keep] = da1_pval.data
                # n_mtrx[i, in_keep] = da1_n.data
                n_mtrx[injk[j], in_keep] = sum_one_mtrx[in_keep]

        # create masked array (not sure if this necessary)
        corr_mtrx = ma.masked_invalid(corr_mtrx)
        pval_mtrx = ma.masked_invalid(pval_mtrx)
        n_mtrx = ma.masked_invalid(n_mtrx)

        # save beta, pval and n
        cpn_lbl = ['corr', 'pval', 'n']
        cpn_list = [corr_mtrx, pval_mtrx, n_mtrx]

        # Open temporary file and load data into it
        with netCDF4.Dataset(os.path.join(WORK_DIR, ncgen_ofile_nc), 'a') as nc:
            # place corr, pval, n in nc
            for j in range(len(cpn_lbl)):
                nc['{}'.format(cpn_lbl[j])][0, :, :] = cpn_list[j]

        # save the temporary file to the final file name
        myCmd = ' '.join(['nccopy',
                          '-d6',
                          os.path.join(WORK_DIR, ncgen_ofile_nc),
                          os.path.join(odir, ofile)
                          ])
        print('nccopy', subprocess.call(myCmd, shell=True))
        print('Done with', ofile)

        # Clean up directories
        os.remove(os.path.join(WORK_DIR, ncgen_ofile_nc))


if __name__ == "__main__":
    main()
