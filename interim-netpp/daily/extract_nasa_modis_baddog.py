"Extract MODIS data from NASA to make 8-day composite NetPP"

# Import packages
import argparse
import os
import sys
import xarray as xr
from datetime import datetime, timedelta

# Select days of the year with 8 day data

# NASA makes an 8 day every 8th day beginning with Jan.1, or day of the year 1.
# So, for example, the next day of the year will be 9.

# The opendap url to Jan. 1, 2023 is this:
# https://oceandata.sci.gsfc.nasa.gov/opendap/VIIRS/L3SMI/2023/0101/SNPP_VIIRS.20230101_20230108.L3m.8D.CHL.chlor_a.4km.nc

# http://oceandata.sci.gsfc.nasa.gov/opendap/MODISA/L3SMI/2023/0101/AQUA_MODIS.20230101_20230108.L3m.8D.CHL.chlor_a.4km.nc

# Breaking that down into parts you get:
# https://oceandata.sci.gsfc.nasa.gov/opendap/VIIRS/L3SMI/
# 2023/  # year
# 0101/  # month/day
# AQUA_MODIS.20230101_20230108.L3m.8D.CHL.chlor_a.9km.nc   # file name

# And the file name pattern is:
# `sensor`.`startdate_enddate`.L3m.8D.CHL.chlor_a.9km.nc'
# Where:
# sensor: 'SNPP_VIIRS'
# startdate = yyyy.mm.dd
# enddate = yyyy.mm.dd  # startdate plus 7 days

# So, for every 8 days you have to generate the correct part
# of the url, i.e year and month/day
# Plus the file name

# Use the same logic to get the PAR and SST
# Everything in the URL is the same except the file name
# AQUA_MODIS.20230101_20230108.L3m.8D.PAR.par.9km.nc
# AQUA_MODIS.20230101_20230108.L3m.8D.SST.sst.9km.nc

# Create global variables
ROOT_DIR = ("/home/madison/projects/netpp")
OUTPUT_DIR_T = os.path.join(ROOT_DIR, "data/{}/8_day")

# Define file name templates for each product
TEMPLATES = {
    "CHL": "{}.{}_{}.L3m.8D.CHL.chlor_a.9km.nc",
    "PAR": "{}.{}_{}.L3m.8D.PAR.par.9km.nc",
    "SST": "{}.{}_{}.L3m.8D.SST.sst.9km.nc",
}

# Run main function


def main():

    help_txt = {
        "describe": (
            "Download 8-day composite NASA data from ",
            "MODIS or VIIRS by specifying sensor, start ",
            "data, and data product type."
        ),
        "start_date": "Start date in the format YYYY-MM_DD, e.g 2023-04-07",
        "end_date": "End date in the format YYYY-MM-DD",
        "sensor": "Sensor name (AQUA_MODIS or SNPP_VIIRS)",
        "product": "Data product to download (CHL, PAR, SST, or ALL for all products)",
        "platform": "Platform identifier in URL (default: MODISA)",
        "overwrite": "Set to overwrite a file that exists"
    }

    doc_formatter = argparse.RawDescriptionHelpFormatter
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=doc_formatter)

    parser.add_argument(
                        "-a",
                        "--start_date",
                        help=help_txt["start_date"],
                        required=True
                        )
    parser.add_argument(
                        "-z",
                        "--end_date",
                        help=help_txt["end_date"],
                        required=True
                        )
    parser.add_argument(
                        "-s",
                        "--sensor",
                        help=help_txt["sensor"],
                        choices=["AQUA_MODIS", "SNPP_VIIRS"],
                        required=True,
                        )
    parser.add_argument(
                        "-p",
                        "--product",
                        help=help_txt["product"],
                        choices=["CHL", "PAR", "SST", "ALL"],
                        default="ALL",
                        required=True,
                        )
    parser.add_argument(
                        "--platform",
                        help=help_txt["platform"],
                        choices=["MODISA", "VIIRS"],
                        required=True,
                        )
    parser.add_argument("-o", "--overwrite", required=False,
                        action='store_true',
                        help=help_txt["overwrite"])

    args = parser.parse_args()

    # Determine which products to process
    if args.product == "ALL":
        products = TEMPLATES.keys()
    else:
        products = [args.product]

    # Create dynamic directories and verify
    OUTPUT_DIR = OUTPUT_DIR_T.format(args.sensor.lower())

    DIR_LIST = [
        ROOT_DIR,
        OUTPUT_DIR
    ]

    for dr in DIR_LIST:
        os.makedirs(dr, exist_ok=True)
    print(len(DIR_LIST), 'directories validated')

    # Convert start and end dates to datetime objects
    start_date_obj = datetime.strptime(args.start_date, "%Y-%m-%d")
    end_date_obj = datetime.strptime(args.end_date, "%Y-%m-%d")

    # Check if start date is after end date
    if start_date_obj > end_date_obj:
        print("Error: Start date cannot be later than end date")
        sys.exit(1)

    # Loop through each product and generate its specific URL
    for product in products:
        # Create subdirectory for each product within OUTPUT_DIR
        product_dir = os.path.join(OUTPUT_DIR, product.lower())
        os.makedirs(product_dir, exist_ok=True)

        # Loop over the date range in 8-day increments
        current_date = start_date_obj
        while current_date <= end_date_obj:
            interval_end_date = current_date + timedelta(days=7)

        # Generate file name based on the sensor, product, and dates
            template = TEMPLATES[product]
            file_name = template.format(
                args.sensor,
                current_date.strftime('%Y%m%d'),
                interval_end_date.strftime('%Y%m%d')
            )

            # Generate the full URL
            url = '/'.join([
                'http://oceandata.sci.gsfc.nasa.gov/opendap',
                args.platform,
                'L3SMI',
                current_date.strftime('%Y'),
                current_date.strftime('%m%d'),
                file_name
            ])

            # Open the dataset and print summary information
            try:
                ds = xr.open_dataset(url)
                print(f"Dataset opened successfully for {file_name}")

                # Save the dataset to a NetCDF file
                output_path = os.path.join(product_dir, file_name)
                print(f"Saving dataset to {output_path}")
                ds.to_netcdf(output_path)

            except Exception as e:
                print(
                    f"Dataset for {file_name} not available, moving to"
                    f"next period. Error: {e}"
                )

            # Move to the next 8-day interval
            current_date += timedelta(days=8)


if __name__ == "__main__":
    main()
