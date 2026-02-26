"Contingency Tables for Legacy NetPP Product MODIS and Interim NetPP product VIIRS"

"Using a contingency matrix, we can evaluate the agreement "
"between the legacy NetPP product, MODIS, and interim NetPP "
"product, VIIRS, globally. This analysis highlights areas of "
"agreement or disagreement, as well as regions with "
"positive and negative trends."

"Cohen's Kappa provides a summary metric on a scale of 0 to 1, "
"where 1 indicates perfect agreement and 0 indicates "
"no agreement."

# Import packages
import argparse
import geopandas
import regionmask
import os
import xarray as xr
import pandas as pd
import numpy as np
from shapely.geometry import box
from dateutil.parser import parse

# Create global variables
BASE_DIR = "/Users/madisonrichardson/netpp/monthly"
BIN_DIR = os.path.join(BASE_DIR, "bin")
DATA_DIR = os.path.join(BASE_DIR, "data", "monthly")
WORK_DIR = os.path.join(DATA_DIR, "work")
RESOURCES_DIR = os.path.join('/Users/madisonrichardson/netpp/resources')

# Template for ifile
ifile_tmpl = '{}_{}_trend_month_{}_{}_{}_{}_{:03d}percent.nc'

# Define functions
# Function for loading trend data


def load_trends(
        source_list,
        DATA_DIR, ifile_tmpl,
        ncvar,
        timeseries_type,
        startyear,
        endyear,
        prcnt_keep,
        rgn_wnt
):
    """
    Load trends data for the specified sources
    and extract lat/lon ranges.

    Args:
        source_list (list): List of data sources to load.
        DATA_DIR (str): Base directory for the data.
        ifile_tmpl (str): Template for the input file name.
        ncvar (str): Variable name used in the filename template
        (e.g., 'productivity').
        timeseries_type (str): Type of time series (e.g., 'data' or 'anom').
        startyear (str): Start year for the data range.
        endyear (str): End year for the data range.
        prcnt_keep (float): Percentage of data to retain expressed as a
        decimal (e.g., 0.5).
        rgn_wnt (list): Bounding region as [[min_lon, max_lon], [min_lat,
        max_lat]].

    Returns:
        tuple: A list of loaded datsets (loaded_datasets), longitude data
        (region_lons), and latitude data (region_lats).
    """
    loaded_datasets = []
    for source in source_list:
        # COnstruct input directory and file path
        idir = os.path.join(
            DATA_DIR,
            'trends',
            source,
        )
        ifile = ifile_tmpl.format(
            ncvar,
            timeseries_type,
            source,
            '9km',
            startyear,
            endyear,
            int(prcnt_keep*100)
        )

        # Load dataset using xarray and append list
        loaded_data = xr.open_dataset(os.path.join(idir, ifile))
        loaded_datasets.append(loaded_data)

    # Extract lat and lon for region of interest
    region_lons = loaded_data['longitude'].sel(
        longitude=slice(rgn_wnt[0][0], rgn_wnt[0][1])
    ).data
    region_lats = loaded_data['latitude'].sel(
        latitude=slice(rgn_wnt[1][1], rgn_wnt[1][0])
    ).data

    return loaded_datasets, region_lons, region_lats


# Function for finding common grid points among the sources


def get_common_points(source_list, beta_data_list, pval_data_list):
    """
    Extract common grid points with valid data between
    two sources. Reshapes 2D beta and p-value arrays
    into 1D vectors for comparison. Only grid points
    present in both datasets are included in the output.

    Args:
        source_list (list): List of source names (e.g., 'modis' or 'noaa20')
        beta_data_list (list): List of arrays containing beta values for each
        source dataset.
        pval_data_list (list): List of arrays containing p-values for each
        source dataset.

    Returns:
        tuple:
            beta_common (list): A list containing beta value for common grid
            points from both datasets.
            pval_common (list): A list containing p-values for common grid
            points from both datasets.
            num_common (int): The total number of common grid points with
            valid data.
    """
    # Find the common grid points with trends
    beta_re_list = []  # reshaped beta values
    pval_re_list = []  # reshaped p-values
    ind_re_list = []  # indices of non-missing values

    for i in range(len(source_list)):
        # Reshape 2D arrays into 1D vectors
        ny, nx = beta_data_list[i].shape
        beta_re = np.reshape(beta_data_list[i].data, ny*nx)
        pval_re = np.reshape(pval_data_list[i].data, ny*nx)

        # Find indices of finites values (non-NaN)
        ind_re = np.isfinite(beta_re).nonzero()[0]

        # Replace zeros with small values to avoid issues
        trend_indices_1 = np.where(beta_re == 0.0)
        beta_re[trend_indices_1] = 0.00001
        pval_re[trend_indices_1] = 0.00001

        # Append results to lists
        beta_re_list.append(beta_re)
        pval_re_list.append(pval_re)
        ind_re_list.append(ind_re)

    # Find common indices between the two datasets
    common_indices_12 = np.isin(
        ind_re_list[0],
        ind_re_list[1]
    )  # common in source1 and source2
    common_indices_21 = np.isin(
        ind_re_list[1],
        ind_re_list[0]
    )  # common in source2 and source1

    # Extract beta and p-value data for the common grid points
    beta_common = [
        beta_re_list[0][ind_re_list[0]][common_indices_12],
        beta_re_list[1][ind_re_list[1]][common_indices_21]
    ]
    pval_common = [
        pval_re_list[0][ind_re_list[0]][common_indices_12],
        pval_re_list[1][ind_re_list[1]][common_indices_21]
    ]

    # Number of common points
    num_common = len(beta_common[0])

    return beta_common, pval_common, num_common


# Function for 2x2 Contingency table
def make_2x2_contingency_table(beta_common, num_common, source_list):
    """
    Create a 2x2 contingency table comparing
    trends between two sources.

    Args:
        beta_common (list): A list containing two arrays of beta values for
        common grid points.
        num_common (int): Total number of common grid points.
        source_list (list): List of source names (e.g., 'modis', 'noaa20')

    Returns:
        pd.Dataframe: A Pandas DataFrame representing the 2x2 contingency
        table with percentages.
    """
    # Define the labels for a table
    operator_lbl = [r'$\beta>=0$', r'$\beta<0$']
    col_labels_2x2 = []
    row_label_2x2 = []
    for i in range(2):
        col_labels_2x2.append(
            '{} {}'.format(source_list[0].upper(),
                           operator_lbl[i])
        )
        row_label_2x2.append(
            '{} {}'.format(source_list[1].upper(),
                           operator_lbl[i])
        )

    # Define trend sign combinations for the contingency table
    trend_sign_combos = [[1, 1], [-1, 1], [1, -1], [-1, -1]]
    # Intialize the table as a vector
    table22_vec = np.zeros(len(trend_sign_combos))

    # Intialize the table as a vector
    trend_sign_combos = [[1, 1], [-1, 1], [1, -1], [-1, -1]]
    table22_vec = np.zeros(len(trend_sign_combos))
    for i in range(len(trend_sign_combos)):
        # Identify points with the specified trends
        sign_combo = trend_sign_combos[i]
        trend_indices_1 = np.where(sign_combo[0]*beta_common[0] > 0)[0]
        trend_indices_2 = np.where(sign_combo[1]*beta_common[1] > 0)[0]
        common_indices_sig1 = np.isin(trend_indices_1, trend_indices_2)
        num1 = len(common_indices_sig1.nonzero()[0])
        table22_vec[i] = 100*num1/num_common

    # Convert vector to 2x2 matrix
    contingency_table_2x2 = table22_vec.reshape(2, 2)

    # Create a Pandas DataFrame from the 2x2 table
    contingency_df_2x2 = pd.DataFrame(
        contingency_table_2x2, columns=col_labels_2x2
    )
    contingency_df_2x2 = contingency_df_2x2.set_index(
        np.array(row_label_2x2)
    )

    return (
        contingency_table_2x2,
        contingency_df_2x2
    )


# Function for 3x3 Contingency table


def make_3x3_contingency_table(
        beta_common,
        pval_common,
        num_common,
        alpha,
        source_list
):
    """
    Create a 3x3 contingency table comparing trends and
    significance between two sources.

    Args:
        beta_common (list): A list containing two arrays of beta values for
        common grid points.
        pval_common (list): A list containing two arrays of p-values for
        common grid points.
        num_common (int): Total number of common grid points.
        alpha (float): Significance threshold for p-values.
        source_list (list): List of source names (e.g., 'modis', 'noaa20')

    Returns:
        pd.Dataframe: A Pandas DataFrame representing the 3x3 contingency
        table with percentages.
    """

    # Initialize the 3x3 table
    contingency_table_3x3 = np.zeros([3, 3])

    # Define and create labels for the 3x3 matrix
    operator_lbl = ['n.s.', r'$\beta>=0$', r'$\beta<0$']
    col_labels_3x3 = []
    row_label_3x3 = []
    for i in range(3):
        col_labels_3x3.append(
            '{} {}'.format(source_list[0].upper(), operator_lbl[i])
        )
        row_label_3x3.append(
            '{} {}'.format(source_list[1].upper(), operator_lbl[i])
        )

    # Calculate values for the first row/column in non-significant cases (ns)
    ns_indices_1 = np.where(pval_common[0] >= alpha)[0]
    ns_indices_2 = np.where(pval_common[1] >= alpha)[0]
    common_indices_ns = np.isin(ns_indices_1, ns_indices_2).nonzero()[0]
    contingency_table_3x3[0, 0] = 100*len(common_indices_ns)/num_common

    # Calculate when beta1 is ns and beta2 is sig (pos/neg trends)
    # Identify indices in the beta2 where trends are sig
    sig_indices_2 = np.where(pval_common[1][ns_indices_1] < alpha)[0]

    # Extract sig beta values from beta2
    beta2_sig = beta_common[1][ns_indices_1][sig_indices_2]

    # Separate sig beta2 values into pos and neg trends
    in2_sig_pos = np.where(beta2_sig >= 0)[0]
    in2_sig_neg = np.where(beta2_sig < 0)[0]

    # Calculate the percentage of positive and negative trends for beta2
    contingency_table_3x3[1, 0] = 100*len(in2_sig_pos)/num_common  # % of pos trends when beta1 is ns
    contingency_table_3x3[2, 0] = 100*len(in2_sig_neg)/num_common  # % of neg trends when beta1 is ns

    # Calculate the percentage of positive and negative trends for beta2
    contingency_table_3x3[1, 0] = 100*len(in2_sig_pos)/num_common  # % of pos trends when beta1 is ns
    contingency_table_3x3[2, 0] = 100*len(in2_sig_neg)/num_common  # % of neg trends when beta1 is ns

    # Calculate when beta2 is ns and beta1 is sig (pos/neg trends)
    # Identify indices in the beta1 where trends are sig
    sig_indices_1 = np.where(pval_common[0][ns_indices_2] < alpha)[0]

    # Extract the sig beta1 values
    beta1_sig = beta_common[0][ns_indices_2][sig_indices_1]

    # Separate sig beta1 values into pos/neg trends
    trend_indices_2_sig_pos = np.where(beta1_sig >= 0)[0]
    trend_indices_2_sig_neg = np.where(beta1_sig < 0)[0]

    # Calculate the percentage of pos trends for beta1
    contingency_table_3x3[0, 1] = 100*len(trend_indices_2_sig_pos)/num_common  # % of pos trends when beta2 is ns
    contingency_table_3x3[0, 2] = 100*len(trend_indices_2_sig_neg)/num_common  # % of neg trends when beta2 is ns

    # Calculate when both beta1 and beta2 have sig trends
    # Identify sig trends for both datasets
    sig_indices_1 = np.where(pval_common[0] < alpha)[0]
    sig_indices_2 = np.where(pval_common[1] < alpha)[0]

    # Find common sig trends for both datsets
    common_indices_sig1 = np.isin(sig_indices_1, sig_indices_2).nonzero()[0]
    common_indices_sig2 = np.isin(sig_indices_2, sig_indices_1).nonzero()[0]

    # Extract beta values for common sig points
    reshaped_data = [
        beta_common[0][sig_indices_1][common_indices_sig1],
        beta_common[1][sig_indices_2][common_indices_sig2]
    ]

    # Define sign combinations for trends in both datasets
    trend_sign_combos = [[1, 1], [-1, 1], [1, -1], [-1, -1]]

    # Initialize a vector to store the percentages for each combo
    table33_vec = np.zeros(len(trend_sign_combos))

    # Calculate the percentage of grid points for each trend combo
    for i in range(len(trend_sign_combos)):
        sign_combo = trend_sign_combos[i]
        trend_indices_1 = np.where(sign_combo[0]*reshaped_data[0] > 0)[0]
        trend_indices_2 = np.where(sign_combo[1]*reshaped_data[1] > 0)[0]
        common_indices_sig1 = np.isin(trend_indices_1, trend_indices_2)
        num1 = len(common_indices_sig1.nonzero()[0])
        table33_vec[i] = 100*num1/num_common

    # Populate the lower-right 4 cells
    contingency_table_3x3[1:, 1:] = np.reshape(table33_vec, [2, 2])

    # Create a Pandas DataFrame from the 3x3 table
    contingency_df_3x3 = pd.DataFrame(
        contingency_table_3x3,
        columns=col_labels_3x3
    )
    contingency_df_3x3 = contingency_df_3x3.set_index(np.array(row_label_3x3))
    contingency_df_3x3 = np.round(contingency_df_3x3*100)/100

    return (
        contingency_table_3x3,
        contingency_df_3x3
    )


# Function for finding Cohen's Kappa


def find_kappa(contingency_table):
    """
    Calculate Cohen's Kappa for a given contingency table.

    Args:
        contingency_table (numpy.ndarray): A square matrix.
        representing the contingency table (e.g., 2x2 or 3x3).

    Returns:
        float: Cohen's Kappa statistic.
    """

    # Calculate row and column sums
    row_sum = np.sum(contingency_table, axis=1)
    clmn_sum = np.sum(contingency_table, axis=0)

    # Calculate the trace
    trc_sum = np.trace(contingency_table)

    # Calculate the overall sum (total number of grid points)
    overall_sum = np.sum(row_sum)

    # Calculate expected freqs for agreement by chance
    expected_freqs = row_sum * clmn_sum / overall_sum
    expected_freqs_sum = np.sum(expected_freqs)

    # Compute kappa
    kappa = (trc_sum - expected_freqs_sum) / (overall_sum - expected_freqs_sum)

    return kappa


# Run main function


def main():
    """
    Runs the main function.

    Raises:
        ValueError: ProvCode does not exist in the shapefile.
    """

    doc_formatter = argparse.RawDescriptionHelpFormatter
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=doc_formatter
    )

    parser.add_argument("-a", "--start",
                        type=int,
                        required=True,
                        help="Start date in YYYY format")
    parser.add_argument("-z",
                        "--end",
                        type=int,
                        required=True,
                        help="End date in YYYY format")
    parser.add_argument("-p",
                        "--provcode",
                        type=str,
                        required=True,
                        help="Longhurst Province code or 'box' for custom region")
    parser.add_argument("--lat_rgn",
                        nargs=2,
                        type=float,
                        metavar=("LON_MIN", "LON_MAX"),
                        help="Longitude range as two values: min max")
    parser.add_argument("--lon_rgn",
                        nargs=2,
                        type=float,
                        metavar=("LAT_MIN", "LAT_MAX"),
                        help="Latitude range as two values: min max")
    parser.add_argument("-s1",
                        "--source1",
                        type=str,
                        required=True,
                        choices=["modis", "noaa20"],
                        help="First data source")
    parser.add_argument("-s2",
                        "--source2",
                        type=str,
                        required=True,
                        choices=["modis", "noaa20"],
                        help="Second data source")
    parser.add_argument("-pk",
                        "--prcnt_keep",
                        type=float,
                        required=True,
                        default=0.5,
                        help="Percentage of grid points to keep")
    parser.add_argument("-n",
                        "--ncvar",
                        type=str,
                        required=True,
                        default="productivity",
                        help="NetCDF variable to analyze")
    parser.add_argument("-t",
                        "--timeseries_type",
                        type=str,
                        required=True,
                        choices=["data", "anom"],
                        help="Type of timeseries")
    parser.add_argument("--alpha",
                        type=float,
                        required=True,
                        default=0.05,
                        help="Significance threshold for p-values")

    args = parser.parse_args()

    # Set input variables
    ProvCode = args.provcode
    lat_rgn = args.lat_rgn
    lon_rgn = args.lon_rgn
    startyear = args.start
    endyear = args.end
    source1 = args.source1
    source2 = args.source2
    prcnt_keep = args.prcnt_keep
    ncvar = args.ncvar
    timeseries_type = args.timeseries_type
    alpha = args.alpha

    # Load Longhurst Provinces
    # Determine region of interest
    if ProvCode == 'box':
        # If 'box' is specified, use predefined longitude and latitude ranges
        rgn_wnt = [lon_rgn, lat_rgn]
    else:
        # Load the Longhurst Provinces shape files into a geopandas dataframe
        shape_path = os.path.join(
            RESOURCES_DIR,
            'Longhurst',
            'Longhurst_world_v4_2010.shp'
        )
        shapefiles = geopandas.read_file(shape_path)
        # Display the first 8 rows of the shapefile for reference
        shapefiles.head(8)

        # Ensure ProvCode exists in the shapefile
        if ProvCode not in shapefiles["ProvCode"].values:
            raise ValueError(
                f"ProvCode '{ProvCode}' does not exist in the shapefile."
            )

        # Locate the row with the ProvCode code
        prov_wnt = shapefiles.loc[shapefiles["ProvCode"] == ProvCode]

        # Find the coordinates of the bounding box
        # The bounding box is the smallest rectangle that will completely 
        # enclose the province.
        # We will use the bounding box coordinates to subset the satellite data
        gs_bnds = prov_wnt.bounds

        # Save bounding boc coordinates as the region of interest
        rgn_wnt = [
            [gs_bnds.minx.item(), gs_bnds.maxx.item()],
            [gs_bnds.miny.item(), gs_bnds.maxy.item()]
        ]

    # Load trend datasets
    source_list = [source1, source2]
    print(f"Loading trend datasets for sources: {source_list}")
    loaded_datasets, region_lons, region_lats = load_trends(
        source_list=source_list,
        DATA_DIR=DATA_DIR,
        ifile_tmpl=ifile_tmpl,
        ncvar=ncvar,
        timeseries_type=timeseries_type,
        startyear=startyear,
        endyear=endyear,
        prcnt_keep=prcnt_keep,
        rgn_wnt=rgn_wnt
    )

    # Get region from shapefile
    # ProvCode='box', use geopandas to create a new GeoDataFrame 
    # Provcode='Longhurt', region from the province shape file
    # Create the region from the shape file
    if ProvCode == 'box':
        # Create a rectangular box region if 'box' is specified
        box_rgn = box(
            region_lons[0],
            region_lats[0],
            region_lons[-1],
            region_lats[-1]
        )
        newdata = geopandas.GeoDataFrame(
            index=[0], crs='epsg:4326', geometry=[box_rgn]
        )
        region = regionmask.from_geopandas(newdata)
    else:
        # Use the geometry of the selected Longhurst Province as the region
        region = regionmask.from_geopandas(prov_wnt)

    # Mask the data based on region
    # Mask the trend data for both source types
    # Apply the mask to the DataArray using the 'where' function.
    # The 'where' function sets any gridpoints outside the mask to a NaN value.
    beta_data_list = []
    pval_data_list = []
    for i in range(len(source_list)):
        # Extract trend (beta) and significance (pval) from the dataset
        beta_data = (
            loaded_datasets[i]['beta']
            .squeeze()
            .sel(latitude=region_lats)
            .sel(longitude=region_lons)
        )
        pval_data = (
            loaded_datasets[i]['pval']
            .squeeze()
            .sel(latitude=region_lats)
            .sel(longitude=region_lons)
        )

        # Create region masked using specified geometry
        mask = region.mask(beta_data.longitude, beta_data.latitude)

        # Apply mask the the satellite data
        beta_data_mask = beta_data.where(mask == region.numbers[0])
        beta_data_list.append(beta_data_mask)
        pval_data_mask = pval_data.where(mask == region.numbers[0])
        pval_data_list.append(pval_data_mask)

    print("Region masking completed for all datasets.")

    # Find common grid points that have trend values
    beta_common, pval_common, num_common = get_common_points(
        source_list,
        beta_data_list,
        pval_data_list
    )

    print(f"Found {num_common} common grid points between sources.")

    print('Province: {}\n'.format(ProvCode))
    print(
        'Region: {:6.2f} to {:6.2f}, {:6.2f} to {:6.2f}\n'.format(
            rgn_wnt[0][0], rgn_wnt[0][1], rgn_wnt[1][0], rgn_wnt[1][1]
            )
    )

    # Create 2x2 contingency table
    contingency_table_2x2, contingency_df_2x2 = make_2x2_contingency_table(
       beta_common, num_common, source_list
       )
    print('Contingency table 2 x 2:\n')
    print(contingency_df_2x2)

    # Create 3x3 contingency table
    contingency_table_3x3, contingency_df_3x3 = make_3x3_contingency_table(
        beta_common,
        pval_common,
        num_common,
        alpha,
        source_list
    )
    print('Contigency table 3 x 3:\n')
    print(contingency_df_3x3)

    # Find Cohen's Kappa for 2x2 and 3x3 Contingency Tables
    kappa_2x2 = find_kappa(contingency_table_2x2)
    kappa_3x3 = find_kappa(contingency_table_3x3)
    print(f"2x2 Kappa: {kappa_2x2:.2f}")
    print(f"3x3 Kappa: {kappa_3x3:.2f}")


if __name__ == '__main__':
    main()
