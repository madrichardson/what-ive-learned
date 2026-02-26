"""Create trends and p values for time series using netPP monthly products.

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

import xarray as xr
import pandas as pd
import numpy as np
import os
import argparse
import netCDF4
import numpy.ma as ma
import sys
from scipy import stats
from datetime import datetime
from dateutil.parser import parse
import subprocess

BASE_DIR = "/Users/madisonrichardson/netpp"
BIN_DIR = os.path.join(BASE_DIR, "bin")
DATA_DIR = os.path.join(BASE_DIR, "data", "monthly")  # /snpp/
IDATA_DIR = os.path.join(DATA_DIR, "results")
ODATA_DIR = os.path.join(DATA_DIR, "stats")
WORK_DIR = os.path.join(DATA_DIR, "work")
RESOURCES_DIR = os.path.join(BASE_DIR, "resources")


def main():
    """Create trends on monthly means of primary productivity products.

    Generate trends and significance of monthly primary productivity fields
    made from satellite data using the method of Behrenfeld and Falkowski 1997.
    Trend analysis described in Melin et al 2017.
    Users set the satellite source, resolution, start/stop dates, dataset, data
    or anomaly time series.
    """
    # -----------------------------------------------------------------------------
    # Start: Parser
    # -----------------------------------------------------------------------------
    # Parse arguments
    print("start main")
    help_txt = {
        "describe": (
            "Trend analysis by "
            "setting satellite source, resolution, start/stop date, "
            "dataset, time series type, and number of values in time series"
            ),
        "resolution": "Set the resolution of the data",
        "start": "Date farthest in past, format YYmm as an int",
        "end": "Most recent date, format YYmm as an int",
        "source": "Source satellite data: snpp, noaa20, modis",
        "dataset": "Dataset in primary productivity netcdf: productivity",
        "timeseries_type": "Time series using data or anomalies",
        "percent_keep": ("Trends for grid points with enough values, "
                         "range 0. to 1.0"),
        "overwrite": "Set to overwrite existing output files.",
    }

    doc_formatter = argparse.RawDescriptionHelpFormatter
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=doc_formatter)

    parser.add_argument(
                        "-a",
                        "--startyear",
                        help=help_txt["start"],
                        required=True
                        )
    parser.add_argument(
                        "-z",
                        "--endyear",
                        help=help_txt["end"],
                        required=True
                        )
    parser.add_argument(
                        "-s",
                        "--source",
                        help=help_txt["source"],
                        choices=["snpp", "nasa", "noaa20", "modis"],
                        required=True,
                        )
    #parser.add_argument(
                        #"-d",
                        #"--dataset",
                        #help=help_txt["dataset"],
                        #choices=["productivity"],
                        #required=True,
                        #)
    parser.add_argument(
                        "-t",
                        "--timeseries_type",
                        help=help_txt["timeseries_type"],
                        choices=["data", "anom"],
                        required=True,
                        )
    parser.add_argument(
                        "-p",
                        "--percent_keep",
                        help=help_txt["percent_keep"],
                        type=float,
                        required=True
                        )
    parser.add_argument(
                        "-o", "--overwrite",
                        action="store_true",
                        help=help_txt["overwrite"]
                        )

    args = parser.parse_args()
    # startyear = args.startyear
    # start_date = parse(args.startyear)
    if '-' in args.startyear:
        start_date = parse(args.startyear + '-01')
    else:
        start_date = parse(args.startyear + '01')
        
    if '-' in args.startyear:
        end_date = parse(args.endyear + '-01')
    else:
        end_date = parse(args.endyear + '01')

    # endyear = args.endyear
    print('end_date', end_date)
    #end_date = parse(args.endyear)
    source = args.source
    ncvar_wnt = "productivity"
    timeseries_type = args.timeseries_type
    prcnt_keep = args.percent_keep
    overwrite = args.overwrite
    # -----------------------------------------------------------------------------
    # END: Parser
    # -----------------------------------------------------------------------------

    # .1) test the input args
    if start_date > end_date:
        print("start_date, -a, must be earlier than end_date, -z")
        sys.exit()

    if prcnt_keep < 0.0 or prcnt_keep > 1.0:
        print("percent_keep, -p, must be a float between 0 and 1")
        sys.exit()

    # .2) setup stuff
    # templates for the output file
    ifile_tmpl = "productivity_month_{}_{}_{}.nc"
    ofile_tmpl = "{}_{}_trend_month_{}_{}_{}_{}_{:03d}percent.nc"

    # calculate trends for ncvar variable wanted
    ncvar_list = [ncvar_wnt]
    for ncvar in ncvar_list:
        # .3) only calculate the trends if not already done
        # list files in output directory
        #odir = os.path.join(
            #ODATA_DIR,
            #source
        #)
        odir_filelist = os.listdir(ODATA_DIR)

        # generate the name of output file
        ofile = ofile_tmpl.format(
            ncvar,
            timeseries_type,
            source,
            "9km",
            int(start_date.year),  # startyear
            int(end_date.year),  # endyear
            int(prcnt_keep * 100),
            )

        idir = os.path.join(IDATA_DIR, source)

        # .4) input and output files for ncgen
        ncgen_ofile_nc_tmpl = "ncgen_trend_ofile{}.nc"
        now = datetime.now()
        ncgen_str_time = "{0:%Y%m%d%H%M%S}".format(now)
        ncgen_ofile_nc = ncgen_ofile_nc_tmpl.format(ncgen_str_time)

        # create temporary file to accept output data from a .cdl file
        # cdl is unique for each ncvar
        ncgen_ifile_cdl = "trends_{}_month_{}_{}.cdl".format(ncvar,
                                                             source,
                                                             "9km")

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


        # .5) calculate the trends
        # check if file exists, if overwrite=False then skip the calculation.
        if not overwrite:
            if ofile in odir_filelist:
                print(ofile, "already exists")
                continue

        # get start/end year and month
        # yy1 = int(np.floor(startyear / 100))
        # mm1 = startyear - yy1 * 100
        # yy2 = int(np.floor(endyear / 100))
        # mm2 = endyear - yy2 * 100

        # time_bgn = np.datetime64("{}-{:02d}".format(yy1, mm1), "M")
        # time_end = np.datetime64("{}-{:02d}".format(yy2, mm2), "M")
        # dtM = pd.to_datetime(np.arange(start_date,  # time_bgn
                                         # end_date.replace(year=end_date.year + 1),  # time_end + 1
                                         # dtype="datetime64[M]")
                               #)
        
        dates_obj = pd.date_range(start_date, end_date, freq='MS')
        
        date_len = len(dates_obj)
        #yy = dates_obj.year.values
        #mm = dates_obj.month.values

        # create year variable for use as a coordinate for dataarrays
        # used in finding the trends
        # year1 = yy + (mm - 1) / 12
        yr_fractional = dates_obj.year + (dates_obj.month - 1) / 12

        # open productivity monthly files and place in list
        ds_list = []
        for i in range(date_len):
            #ifile = ifile_tmpl.format(source,
                                      #dates_obj[0].strftime("%Y%m"),  # str(yy[i]) + f"{mm[i]:02}"
                                      #"9km")
            # DHR I think the code above has an error:
            # dates_obj[0] should be dates_obj[i]
            ifile = ifile_tmpl.format(source,
                                      dates_obj[i].strftime("%Y%m"),  # str(yy[i]) + f"{mm[i]:02}"
                                      "9km")
            if not os.path.isfile(os.path.join(idir, ifile)):
                print("file not found for", ifile)
                continue
            ds_list.append(xr.open_dataset(os.path.join(idir, ifile)))

        # size of data for calculating trends
        _, lat_len, lon_len = ds_list[0][ncvar_list[0]].shape
        ds_len = len(ds_list)

        # initialize matrices to hold trend p, beta, and n values
        beta_mtrx = np.zeros([lat_len, lon_len]) * np.nan
        pval_mtrx = np.zeros([lat_len, lon_len]) * np.nan
        n_mtrx = np.zeros([lat_len, lon_len]) * np.nan

        # subset the data matrix by a single latitude
        for i in range(lat_len):
            print(
                "{}: {} ({}) filter by # of dates above {} percent".format(
                    i, lat_len, ncvar, prcnt_keep
                )
            )

            # for the given latitude place each series into the data matrix
            # DHR: This section determines if a series has enough data to use it in calculations
            #   points to each monthly file 
            #   For each latitude, pulls all the values at each longitude
            #   the result is a 2D array with 120 lists of length = len(longitude)  
            #   i.e shape of (120, len(lon)) 
            data_mtrx = np.zeros([ds_len, lon_len])
            for j in range(ds_len):
                data_mtrx[j, :] = ds_list[j][ncvar][0, i, :].data

            # keep grid points that have non-missing numbers above 'prcnt_keep'  
            # DHR: if data_mtrx item has a value > 0 it will = 1, others will be nan
            ones_mtrx =  data_mtrx / data_mtrx
            # DHR: sum total along time axis
            sum_one_mtrx = np.nansum(ones_mtrx, axis=0)
            # DHR: get lon indices where % of valid values over time is above threshold
            in_keep = np.where(sum_one_mtrx > ds_len * prcnt_keep)[0]

            # create datarray with time as coordinate, useful for
            # calculating anoms in xarray
            # DHR Saves values at indices of lon where total in above > % limit
            #   associated with time
            da1_time = xr.DataArray(
                data_mtrx[:, in_keep],
                coords=[dates_obj.astype("datetime64[ns]"), in_keep],
                dims=["time", "in_keep"],
            )

            # trends for either data or anom time series, for either
            # create new dataarray with year as coordinate, this will
            # result in the trends have units of "(data units)/year"
            if timeseries_type == "anom":
                # anoms are created from subtracting long-term clim from data
                da1_clim = da1_time.groupby("time.month").mean("time")
                da1_anom = da1_time.groupby("time.month") - da1_clim
                da1_compare = xr.DataArray(
                                           da1_anom.data,
                                           coords=[yr_fractional, in_keep],
                                           dims=["yr_fractional", "in_keep"]
                                           )
            elif timeseries_type == "data":  # DHR: change to else?
                da1_compare = xr.DataArray(
                                           da1_time.data,
                                           coords=[yr_fractional, in_keep],
                                           dims=["yr_fractional", "in_keep"]
                                           )

            # DHR start up here -> create a time index, used for calculating sb
            # DHR Make an array of all year_months, with
            yr_fractional_mtrx = (np.expand_dims(yr_fractional, axis=1)
                                  * np.ones([1, len(in_keep)])
                                  )

            # create time dataarray
            da1_tt = xr.DataArray(
                                  yr_fractional_mtrx,
                                  coords=[yr_fractional, in_keep],
                                  dims=["yr_fractional", "in_keep"]
                                  )

            # mask time variable with missing from da1_compare
            # Keeps values >0, < nan
            da1_tt = da1_tt.where(da1_compare.data > 0, np.nan)

            # create n dataarray
            da1_n = (da1_compare / da1_compare).sum("yr_fractional")

            # get the slope and intercept
            ds1_coeffs = da1_compare.polyfit("yr_fractional", deg=1, full=True)

            ds1_pred = xr.polyval(da1_compare["yr_fractional"], ds1_coeffs)

            # create the beta dataarray
            da1_beta = ds1_coeffs["polyfit_coefficients"].sel(degree=1)

            # create sb dataarray
            # DHR: look in paper Eq. 1
            da1_sb = np.sqrt(
                (1 / (da1_n - 2))
                * (
                    np.square(da1_compare.std("yr_fractional")) / np.square(da1_tt.std("yr_fractional"))
                    - np.square(da1_beta)
                )
            )

            # create t stat
            da1_tstat = da1_beta / da1_sb

            # create degrees of freedom
            da1_df = da1_n - 2

            # create pval dataarray
            da1_pval = (
                xr.apply_ufunc(stats.distributions.t.sf,
                               np.abs(da1_tstat),
                               da1_df
                               ) * 2
            )
            # place beta, pval values for latitude subset in final global data matrix
            beta_mtrx[i, in_keep] = da1_beta.data
            pval_mtrx[i, in_keep] = da1_pval.data
            n_mtrx[i, in_keep] = da1_n.data

        # create masked array (not sure if this necessary)
        beta_mtrx = ma.masked_invalid(beta_mtrx)
        pval_mtrx = ma.masked_invalid(pval_mtrx)
  
        # save beta, pval and n
        bp_lbl = ["beta", "pval", "n"]
        bp_list = [beta_mtrx, pval_mtrx, n_mtrx]

        # Open temporary file and load data into it
        with netCDF4.Dataset(os.path.join(WORK_DIR,
                                          ncgen_ofile_nc
                                          ), "a") as nc:
            # place both beta and pval in nc
            for j in range(len(bp_lbl)):
                nc["{}".format(bp_lbl[j])][0, :, :] = bp_list[j]

        # save the temporary file to the final file name
        myCmd = " ".join(
            [
                "nccopy",
                "-d6",
                os.path.join(WORK_DIR, ncgen_ofile_nc),
                os.path.join(ODATA_DIR, ofile),
            ]
        )
        print("nccopy", subprocess.call(myCmd, shell=True))
        print("Done with", ofile)

        # Clean up directories
        os.remove(os.path.join(WORK_DIR, ncgen_ofile_nc))


if __name__ == "__main__":
    main()
