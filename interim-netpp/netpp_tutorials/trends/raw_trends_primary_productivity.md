---
date: "2025-11-17"
output: html_document
editor_options: 
  markdown: 
    wrap: 72
---

# Calculating Trends and P-Values for Legacy and Interim NetPP Products

> Update: November 2025

## Objectives

This script generates trends for time series of monthly means of primary
productivity described by Behrenfeld and Falkowski 1997. The trend
analysis follows methods outlined in Melin et al 2017, see the section
2.3 "Trend estimates and comparison of trends". Trends can also be
calculated for dataset of PAR, chlorophyll, SST, etc.

------------------------------------------------------------------------

## Datasets Overview

We will be computing pixel-by-pixel trends for the VIIRS-NOAA20 product.
Here is the dataset we will be using and the other primary productivity
products, including VIIRS-SNPP and MODIS-Aqua:

1.  **Primary Productivity, MODIS Aqua, Science Quality, Global, 9km,
    2013-2022, (Monthly Composite)**

-   Distributed via the West Coast Node ERDDAP dataset at the following
    link:
    <https://coastwatch.pfeg.noaa.gov/wcn/erddap/griddap/productivity_modis_aqua_monthly_9km.graph>

2.  **Primary Productivity, VIIRS SNPP, Science Quality, Global, 9km,
    2013-2022 (Monthly Composite)**

-   Distributed via the West Coast Node ERDDAP dataset at the following
    link:
    <https://coastwatch.pfeg.noaa.gov/wcn/erddap/griddap/productivity_viirs_snpp_monthly_9km.graph>

3.  **Primary Productivity, VIIRS NOAA20, Science Quality, Global, 9km,
    2018-2022, (Monthly Composite)**

-   Distributed via the West Coast Node ERDDAP dataset at the following
    link:
    <https://coastwatch.pfeg.noaa.gov/wcn/erddap/griddap/productivity_viirs_noaa20_monthly_9km.graph>

## Resource requirements

-   **R Studio** with the modules included within the *Install and Load
    Required Packages* section below

-   **Internet connection**

## Tutorial for this notebook

For this tutorial, we will be using the **VIIRS-NOAA20** dataset to
compute monthly primary productivity trends at 9km resolution for the
West Coast of the US.

## Install and Load Required Packages

```{r}
pkges = installed.packages()[,"Package"]
# Function to check if pkgs are installed, install missing pkgs, and load
pkgTest <- function(x)
{
  if (!require(x,character.only = TRUE))
  {
    install.packages(x,dep=TRUE,repos='http://cran.us.r-project.org')
    if(!require(x,character.only = TRUE)) stop(x, " :Package not found")
  }
}

# create list of required packages
list.of.packages <- c("ncdf4", "rerddap","plotdap", "parsedate", 
                      "sp", "ggplot2", "RColorBrewer", "sf", 
                      "reshape2", "maps", "mapdata", 
                      "jsonlite", "rerddapXtracto", "dplyr",
                      "lubridate", "rlang")

# Run install and load function
for (pk in list.of.packages) {
  pkgTest(pk)
}

# create list of installed packages
pkges = installed.packages()[,"Package"]
```

## Create Global Variables

Global variables are used to set up the directory paths required for
processing trends.

-   **BASE_DIR** is the base directory defined for the project, and the
    root for all other subdirectories.
-   **ODATA_DIR** is where the output trends file will be located.

```{r}
ROOT_DIR <- "/Users/madisonrichardson/netpp"
ODATA_DIR <- file.path(ROOT_DIR, "data", "trends")

dir_list <- c(ROOT_DIR, ODATA_DIR)

# Create if missing
invisible(lapply(dir_list, dir.create, showWarnings = FALSE, recursive = TRUE))

cat(length(dir_list), "directories validated\n")

# 2 directories validated
```


## Make Useful Functions

```{r}
#' Write linear trend fields (β, p, n) to an ERDDAP-style NetCDF file
#'
#' @description
#' Create a CF/ERDDAP-style NetCDF file containing 3-D fields of
#' \eqn{\beta} (slope), \eqn{p} (p-value), and \eqn{n} (sample size)
#' on a regular lon/lat/time grid, using metadata patterned after the
#' `productivity_anom_trend_month_nasa_9km_201301_202212_050percent` CDL.
#'
#' The function:
#' \enumerate{
#'   \item Validates that the input arrays share consistent dimensions.
#'   \item Derives spatial bounds and grid resolution from the latitude
#'         and longitude vectors.
#'   \item Constructs an output filename that encodes the lon/lat bounds.
#'   \item Defines CF-style dimensions (`time`, `latitude`, `longitude`)
#'         and `beta`, `pval`, `n` variables on \code{[lon, lat, time]}.
#'   \item Writes the data to disk and attaches ERDDAP-like coordinate and
#'         global metadata consistent with the trends CDL template.
#' }
#'
#' @param beta_array Numeric 3-D array of trend slopes (β) with dimensions
#'   \code{[time, lat, lon]}.
#' @param pval_array Numeric 3-D array of p-values with dimensions
#'   \code{[time, lat, lon]}.
#' @param n_array Numeric 3-D array of sample sizes (number of months used
#'   for the trend) with dimensions \code{[time, lat, lon]}.
#' @param time_index Numeric vector of length \code{dim(beta_array)[1]}
#'   giving the time coordinates. In the CDL template, \code{time} has
#'   units \code{"count"} and represents month index.
#' @param lat_vals Numeric vector of length \code{dim(beta_array)[2]}
#'   giving latitude coordinates in degrees north.
#' @param lon_vals Numeric vector of length \code{dim(beta_array)[3]}
#'   giving longitude coordinates in degrees east.
#' @param out_dir Character path to the directory where the NetCDF file will
#'   be written. The directory will be created if it does not exist.
#' @param overwrite Logical. If \code{FALSE} (default), the function stops
#'   with an error if a file with the same name already exists. If
#'   \code{TRUE}, an existing file is deleted and recreated.
#' @param source Character. Short source/platform code used in `id`,
#'   `instrument`, `platform`, and `source` global attributes
#'   (e.g. `"nasa"`, `"modis"`, `"viirs"`). Default `"nasa"`.
#' @param platform_name Character. Human-readable platform name used in
#'   `title` and `summary` (e.g. `"MODIS-Aqua"`, `"VIIRS-SNPP"`).
#'
#' @return
#' Invisibly returns the full file path (\code{character}) of the created
#' NetCDF file. A message is printed indicating where the file was written.
#'
#' @export
write_trends_netcdf <- function(beta_array,
                                pval_array,
                                n_array,
                                time_index,
                                lat_vals,
                                lon_vals,
                                out_dir,
                                overwrite = FALSE,
                                source = "nasa",
                                platform_name = "Aqua") {

  # Ensure ncdf4 is available
  if (!requireNamespace("ncdf4", quietly = TRUE)) {
    stop("Package 'ncdf4' is required but not installed.")
  }

  # ------------------------- Basic checks ------------------------------------
  if (!is.numeric(beta_array) || !is.numeric(pval_array) || !is.numeric(n_array)) {
    stop("'beta_array', 'pval_array', and 'n_array' must be numeric.")
  }

  if (length(dim(beta_array)) != 3L ||
      length(dim(pval_array)) != 3L ||
      length(dim(n_array))    != 3L) {
    stop("All input arrays must be 3-D [time, lat, lon].")
  }

  if (!identical(dim(beta_array), dim(pval_array)) ||
      !identical(dim(beta_array), dim(n_array))) {
    stop("All input arrays must have identical dimensions [time, lat, lon].")
  }

  n_time <- length(time_index)
  n_lat  <- length(lat_vals)
  n_lon  <- length(lon_vals)

  if (!identical(dim(beta_array), c(n_time, n_lat, n_lon))) {
    stop(
      "dim(beta_array) must be [length(time_index), length(lat_vals), length(lon_vals)].\n",
      "Found: ", paste(dim(beta_array), collapse = " x "), " vs ",
      n_time, " x ", n_lat, " x ", n_lon
    )
  }

  # --------------------- Bounds & filename -----------------------------------
  lat_min <- min(lat_vals); lat_max <- max(lat_vals)
  lon_min <- min(lon_vals); lon_max <- max(lon_vals)

  fmt <- function(x) sprintf("%.2f", x)

  nc_filename <- sprintf(
    "trends_%s_monthly_9km_lon%s_to_%s_lat%s_to_%s.nc",
    tolower(source),
    fmt(lon_min), fmt(lon_max),
    fmt(lat_min), fmt(lat_max)
  )

  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }
  nc_file_path <- file.path(out_dir, nc_filename)

  # --------------------- Define dimensions -----------------------------------
  dim_time <- ncdf4::ncdim_def(
    name  = "time",
    units = "count",          # from CDL
    vals  = time_index,
    unlim = TRUE,
    create_dimvar = TRUE
  )

  dim_lat <- ncdf4::ncdim_def(
    name  = "latitude",
    units = "degrees_north",
    vals  = lat_vals,
    create_dimvar = TRUE
  )

  dim_lon <- ncdf4::ncdim_def(
    name  = "longitude",
    units = "degrees_east",
    vals  = lon_vals,
    create_dimvar = TRUE
  )

  # --------------------- Define variables ------------------------------------
  # In NetCDF, variables are [lon, lat, time] → R dim order list(lon, lat, time)
  var_beta <- ncdf4::ncvar_def(
    name     = "beta",
    units    = "(mg C m-2 d-1) year-1",
    dim      = list(dim_lon, dim_lat, dim_time),
    longname = "Slopes of Linear Regression",
    missval  = -999.0,
    prec     = "float"
  )

  var_pval <- ncdf4::ncvar_def(
    name     = "pval",
    units    = "",
    dim      = list(dim_lon, dim_lat, dim_time),
    longname = "Level of Significance for the Slopes of Linear Regression",
    missval  = -999.0,
    prec     = "float"
  )

  var_n <- ncdf4::ncvar_def(
    name     = "n",
    units    = "",
    dim      = list(dim_lon, dim_lat, dim_time),
    longname = "Number of months in time series used for trend analysis.",
    missval  = -999.0,
    prec     = "float"
  )

  # --------------------- Create file & write data ----------------------------
  if (file.exists(nc_file_path)) {
    if (!overwrite) {
      stop(
        "NetCDF file already exists at:\n  ", nc_file_path,
        "\nSet overwrite = TRUE or delete the file manually if you want to recreate it."
      )
    } else {
      message("File exists, overwriting: ", nc_file_path)
      unlink(nc_file_path)
    }
  }

  nc <- ncdf4::nc_create(
    nc_file_path,
    vars = list(var_beta, var_pval, var_n)
  )

  # Arrays are [time, lat, lon]; NetCDF expects [lon, lat, time]
  beta_for_nc <- aperm(beta_array, perm = c(3, 2, 1))
  pval_for_nc <- aperm(pval_array, perm = c(3, 2, 1))
  n_for_nc    <- aperm(n_array,    perm = c(3, 2, 1))

  ncdf4::ncvar_put(nc, "beta", beta_for_nc)
  ncdf4::ncvar_put(nc, "pval", pval_for_nc)
  ncdf4::ncvar_put(nc, "n",    n_for_nc)

  # Coordinate variables
  ncdf4::ncvar_put(nc, "time",      time_index)
  ncdf4::ncvar_put(nc, "latitude",  lat_vals)
  ncdf4::ncvar_put(nc, "longitude", lon_vals)

  # --------------------- Coordinate attributes -------------------------------
  now <- Sys.time()

  # time
  ncdf4::ncatt_put(nc, "time", "_CoordinateAxisType", "Time",                 prec = "text")
  ncdf4::ncatt_put(nc, "time", "axis",                "T",                    prec = "text")
  ncdf4::ncatt_put(nc, "time", "calendar",            "gregorian",            prec = "text")
  ncdf4::ncatt_put(nc, "time", "ioos_category",       "Time",                 prec = "text")
  ncdf4::ncatt_put(nc, "time", "long_name",           "Month of Climatology", prec = "text")
  ncdf4::ncatt_put(nc, "time", "standard_name",       "time",                 prec = "text")

  # latitude
  ncdf4::ncatt_put(nc, "latitude", "_CoordinateAxisType", "Lat",          prec = "text")
  ncdf4::ncatt_put(nc, "latitude", "axis",                "Y",            prec = "text")
  ncdf4::ncatt_put(nc, "latitude", "ioos_category",       "Location",     prec = "text")
  ncdf4::ncatt_put(nc, "latitude", "long_name",           "Latitude",     prec = "text")
  ncdf4::ncatt_put(nc, "latitude", "standard_name",       "latitude",     prec = "text")
  ncdf4::ncatt_put(nc, "latitude", "units",               "degrees_north",prec = "text")
  ncdf4::ncatt_put(nc, "latitude", "valid_min",           -90.0,          prec = "float")
  ncdf4::ncatt_put(nc, "latitude", "valid_max",           90.0,           prec = "float")

  # longitude
  ncdf4::ncatt_put(nc, "longitude", "_CoordinateAxisType", "Lon",        prec = "text")
  ncdf4::ncatt_put(nc, "longitude", "axis",                "X",          prec = "text")
  ncdf4::ncatt_put(nc, "longitude", "ioos_category",       "Location",   prec = "text")
  ncdf4::ncatt_put(nc, "longitude", "long_name",           "Longitude",  prec = "text")
  ncdf4::ncatt_put(nc, "longitude", "source_name",         "cols",       prec = "text")
  ncdf4::ncatt_put(nc, "longitude", "standard_name",       "longitude",  prec = "text")
  ncdf4::ncatt_put(nc, "longitude", "units",               "degrees_east",prec = "text")
  ncdf4::ncatt_put(
    nc, "longitude", "actual_range",
    c(min(lon_vals), max(lon_vals)), prec = "double"
  )

  # ---------------------- Data variable attributes --------------------------
  # beta attributes
  ncdf4::ncatt_put(nc, "beta", "_FillValue",            -999.0,                        prec = "float")
  ncdf4::ncatt_put(nc, "beta", "colorBarMaximum",       1100.0,                        prec = "float")
  ncdf4::ncatt_put(nc, "beta", "colorBarMinimum",       -1100.0,                       prec = "float")
  ncdf4::ncatt_put(nc, "beta", "coverage_content_type", "modelResult",                 prec = "text")
  ncdf4::ncatt_put(nc, "beta", "ioos_category",         "Other",                       prec = "text")
  ncdf4::ncatt_put(nc, "beta", "long_name",             "Slopes of Linear Regression", prec = "text")
  ncdf4::ncatt_put(nc, "beta", "missing_value",         -999.0,                        prec = "float")
  ncdf4::ncatt_put(nc, "beta", "source",
                   "NOAA CoastWatch West Coast Node",                                 prec = "text")
  ncdf4::ncatt_put(nc, "beta", "standard_name", "beta",                               prec = "text")
  ncdf4::ncatt_put(nc, "beta", "units", "(mg C m-2 d-1) year-1",                      prec = "text")

  # pval attributes
  ncdf4::ncatt_put(nc, "pval", "_FillValue",            -999.0,                        prec = "float")
  ncdf4::ncatt_put(nc, "pval", "colorBarMaximum",       1.0,                           prec = "float")
  ncdf4::ncatt_put(nc, "pval", "colorBarMinimum",       0.0,                           prec = "float")
  ncdf4::ncatt_put(nc, "pval", "coverage_content_type", "modelResult",                 prec = "text")
  ncdf4::ncatt_put(nc, "pval", "ioos_category",         "Other",                       prec = "text")
  ncdf4::ncatt_put(nc, "pval", "long_name",
                   "Level of Significance for the Slopes of Linear Regression",
                   prec = "text")
  ncdf4::ncatt_put(nc, "pval", "missing_value", -999.0,                                prec = "float")
  ncdf4::ncatt_put(nc, "pval", "source",
                   "NOAA CoastWatch West Coast Node",                                 prec = "text")
  ncdf4::ncatt_put(nc, "pval", "standard_name", "p",                                   prec = "text")
  ncdf4::ncatt_put(nc, "pval", "units", "",                                            prec = "text")

  # n attributes
  ncdf4::ncatt_put(nc, "n", "_FillValue",            -999.0,                           prec = "float")
  ncdf4::ncatt_put(nc, "n", "colorBarMaximum",       120.0,                            prec = "float")
  ncdf4::ncatt_put(nc, "n", "colorBarMinimum",       0.0,                              prec = "float")
  ncdf4::ncatt_put(nc, "n", "coverage_content_type", "modelResult",                    prec = "text")
  ncdf4::ncatt_put(nc, "n", "ioos_category",         "Other",                          prec = "text")
  ncdf4::ncatt_put(nc, "n", "long_name",
                   "Number of months in time series used for trend analysis.",
                   prec = "text")
  ncdf4::ncatt_put(nc, "n", "missing_value", -999.0,                                   prec = "float")
  ncdf4::ncatt_put(nc, "n", "source",
                   "NOAA CoastWatch West Coast Node",                                 prec = "text")
  ncdf4::ncatt_put(nc, "n", "standard_name", "n",                                      prec = "text")
  ncdf4::ncatt_put(nc, "n", "units", "",                                              prec = "text")

  # ---------------------- Global attributes ----------------------------------
  lat_step <- if (length(lat_vals) > 1L) abs(diff(lat_vals)[1]) else NA_real_
  lon_step <- if (length(lon_vals) > 1L) abs(diff(lon_vals)[1]) else NA_real_

  src_upper <- toupper(source)

  # From your Python snippet
  ncdf4::ncatt_put(
    nc, 0, "acknowledgement",
    "The project was supported by funding from the Portfolio Management Branch of NESDIS and NOAA CoastWatch.",
    prec = "text"
  )
  ncdf4::ncatt_put(
    nc, 0, "contributors",
    "Dale Robinson, Isaac Shroeder, Ryan Vandermeulen, Jonathan Sherman, Jesse Espinoza, & Madison Richardson",
    prec = "text"
  )
  ncdf4::ncatt_put(
    nc, 0, "date_created",
    format(as.POSIXct(now, tz = "UTC"), "%Y-%m-%dT%H:%M:%SZ"),
    prec = "text"
  )
  ncdf4::ncatt_put(nc, 0, "instrument", src_upper, prec = "text")
  ncdf4::ncatt_put(
    nc, 0, "id",
    sprintf("trends_%s_monthly_9km", tolower(source)),
    prec = "text"
  )
  ncdf4::ncatt_put(nc, 0, "platform", src_upper, prec = "text")
  ncdf4::ncatt_put(nc, 0, "source",   src_upper, prec = "text")

  ncdf4::ncatt_put(
    nc, 0, "title",
    sprintf(
      "Trend coefficients and p values for monthly primary productivity from %s globally at 9km resolution",
      platform_name
    ),
    prec = "text"
  )

  ncdf4::ncatt_put(
    nc, 0, "summary",
    sprintf(
      "Trends between primary productivity or PAR, chlorophyll and SST from %s. These are 9km products generated from time series of monthly means. Trends are slopes of linear regressions and the level of significance of the trend are computed from a t-test. See Melin et al 2017 for more details.",
      platform_name
    ),
    prec = "text"
  )

  # Keep other CF/ACDD globals
  ncdf4::ncatt_put(nc, 0, "cdm_data_type",       "Grid",                         prec = "text")
  ncdf4::ncatt_put(nc, 0, "Conventions",        "CF-1.6, COARDS, ACDD-1.3",      prec = "text")
  ncdf4::ncatt_put(nc, 0, "creator_name",       "NOAA CoastWatch West Coast Node", prec = "text")
  ncdf4::ncatt_put(nc, 0, "creator_type",       "institution",                   prec = "text")
  ncdf4::ncatt_put(nc, 0, "creator_url",        "https://coastwatch.pfeg.noaa.gov/", prec = "text")
  ncdf4::ncatt_put(nc, 0, "map_projection",     "geographic",                    prec = "text")
  ncdf4::ncatt_put(nc, 0, "time_coverage_resolution", "PD1",                     prec = "text")
  ncdf4::ncatt_put(nc, 0, "project",
                   "NOAA CoastWatch Interim Primary Productivity",               prec = "text")
  ncdf4::ncatt_put(nc, 0, "latitude_units",     "degrees_north",                 prec = "text")
  ncdf4::ncatt_put(nc, 0, "longitude_units",    "degrees_east",                  prec = "text")

  # Geospatial coverage
  ncdf4::ncatt_put(nc, 0, "northernmost_latitude", lat_max,    prec = "float")
  ncdf4::ncatt_put(nc, 0, "southernmost_latitude", lat_min,    prec = "float")
  ncdf4::ncatt_put(nc, 0, "westernmost_longitude", lon_min,    prec = "float")
  ncdf4::ncatt_put(nc, 0, "easternmost_longitude", lon_max,    prec = "float")
  ncdf4::ncatt_put(nc, 0, "geospatial_lat_max",    lat_max,    prec = "float")
  ncdf4::ncatt_put(nc, 0, "geospatial_lat_min",    lat_min,    prec = "float")
  ncdf4::ncatt_put(nc, 0, "geospatial_lon_max",    lon_max,    prec = "float")
  ncdf4::ncatt_put(nc, 0, "geospatial_lon_min",    lon_min,    prec = "float")

  if (!is.na(lat_step)) {
    ncdf4::ncatt_put(nc, 0, "latitude_step",             lat_step, prec = "float")
    ncdf4::ncatt_put(nc, 0, "geospatial_lat_resolution", lat_step, prec = "float")
  }
  if (!is.na(lon_step)) {
    ncdf4::ncatt_put(nc, 0, "longitude_step",            lon_step, prec = "float")
    ncdf4::ncatt_put(nc, 0, "geospatial_lon_resolution", lon_step, prec = "float")
  }

  ncdf4::ncatt_put(nc, 0, "geospatial_lat_units", "degrees_north", prec = "text")
  ncdf4::ncatt_put(nc, 0, "geospatial_lon_units", "degrees_east",  prec = "text")
  ncdf4::ncatt_put(nc, 0, "spatialResolution",    "9.28 km",       prec = "text")

  ncdf4::ncatt_put(nc, 0, "naming_authority",     "gov.noaa.pfeg.coastwatch", prec = "text")
  ncdf4::ncatt_put(
    nc, 0, "keywords",
    paste(
      "EARTH SCIENCE > BIOSPHERE > ECOLOGICAL DYNAMICS > ECOSYSTEM FUNCTIONS > PRIMARY PRODUCTION,",
      "Earth Science > Oceans > Ocean Temperature > Sea Surface Temperature,",
      "EARTH SCIENCE > OCEANS > OCEAN CHEMISTRY > CHLOROPHYLL,",
      "EARTH SCIENCE > OCEANS > OCEAN OPTICS > PHOTOSYNTHETICALLY ACTIVE RADIATION,",
      "EARTH SCIENCE > BIOSPHERE > ECOSYSTEMS > AQUATIC ECOSYSTEMS > PLANKTON > PHYTOPLANKTON,",
      "primary production, primary productivity, production, productivity,",
      "CoastWatch, NOAA, West Coast Node, nesdis, noaa, aqua, ocean, oceans,",
      "satellite, sea_surface_temperature, sst, par, temperature, NASA"
    ),
    prec = "text"
  )
  ncdf4::ncatt_put(nc, 0, "keywords_vocabulary", "GCMD Science Keywords",  prec = "text")
  ncdf4::ncatt_put(
    nc, 0, "proj4_string",
    "+proj=eqc +lat_ts=0 +lat_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs +lon_0=0.000000",
    prec = "text"
  )
  ncdf4::ncatt_put(
    nc, 0, "publisher_name",
    "NOAA/NMFS/SWFSC/ERD, NOAA/NESDIS/STAR/CoastWatch/West Coast Node",
    prec = "text"
  )
  ncdf4::ncatt_put(nc, 0, "publisher_type", "group",              prec = "text")
  ncdf4::ncatt_put(nc, 0, "publisher_url",  "https://coastwatch.pfeg.noaa.gov", prec = "text")
  ncdf4::ncatt_put(
    nc, 0, "references",
    "Behrenfeld and Falkowski (1997) https://doi.org/10.4319/lo.1997.42.1.0001",
    prec = "text"
  )
  ncdf4::ncatt_put(nc, 0, "standard_name_vocabulary",
                   "CF Standard Name Table v70",            prec = "text")

  # Close and finish
  ncdf4::nc_close(nc)

  message("Wrote ERDDAP-like trends NetCDF to: ", nc_file_path)

  invisible(nc_file_path)
}

```

## Set Up Parameters

Initialize the parameters required for the trends calculation:

-   **start_date** and **end_date**: Define the date range for the
    analysis (format: YYYY-MM-DDTXX:XX:XXZ). Monthly composites are
    centered on the 16th of each month.

-   **sensor**: Specify the satellite sensor for the trends analysis
    ('modis', 'snpp', or 'noaa20').

-   **ncvar**: The name of the variable being analyzed (we are using
    'productivity').

-   **overwrite**: If 'True", existing output files will be replaced.

-   **percent_keep**: The minimum proportion of non-missing data
    required for trends calculation. We will be using 50%.

```{r}
start_date = "2018-01-16T12:00:00Z"
end_date = "2022-12-16T12:00:00Z"

# Define sensor (either "snpp", "noaa20", or "modis")
sensor = "noaa20"

# Define variable
ncvar = "productivity"

overwrite = FALSE

percent_keep = 0.5

```

## Select Satellite Dataset from ERDDAP

We will be using the VIIRS-NOAA20 netpp dataset from the West Coast Node
ERDDAP Server. The dataset ID is:
**productivity_viirs_noaa20_monthly_9km**. We will use the info function
from the **rerddap** package to first obtain information about the
dataset of interest, then we will import the data.

```{r}
# Set ERDDAP URL
erddap_url = "https://coastwatch.pfeg.noaa.gov/wcn/erddap/"

# Set dataset ID
noaa20_id = "productivity_viirs_noaa20_monthly_9km"

# Get dataset info
noaa20_dataInfo <- rerddap::info(noaa20_id, url=erddap_url)  
print(noaa20_dataInfo)

# <ERDDAP info> productivity_viirs_noaa20_monthly_9km 
# Base URL: https://coastwatch.pfeg.noaa.gov/wcn/erddap 
# Dataset Type: griddap 
# Dimensions (range):  
#     time: (2018-01-16T12:00:00Z, 2022-12-16T12:00:00Z) 
#     latitude: (-89.95834, 89.95834) 
#     longitude: (-179.9583, 179.9583) 
# Variables:  
#     productivity: 
#         Units: mg C m-2 day-1 
```

## Subset Dataset

The dataset is too large to use the entire globe, so we are subsetting
to use the West Coast of the US.

```{r}
# Define West Coast latitude and longitudes
wc_lon <- c(-135, -115)
wc_lat <- c(30, 50)

# Subset noaa20 dataset
noaa20_ds <- griddap(
  datasetx = noaa20_id,
  url = erddap_url,
  time = c(start_date, end_date),
  longitude = wc_lon,
  latitude = wc_lat,
  fields = ncvar
)

print(noaa20_ds)

# longitude latitude time productivity
#-134.9583	50.04166	2018-01-16T12:00:00Z	245.31010	
#-134.8750	50.04166	2018-01-16T12:00:00Z	238.91531	
#-134.7917	50.04166	2018-01-16T12:00:00Z	235.14050	
#-134.7083	50.04166	2018-01-16T12:00:00Z	237.93399	
#-134.6250	50.04166	2018-01-16T12:00:00Z	241.79921	
#-134.5417	50.04166	2018-01-16T12:00:00Z	211.73715	
#-134.4583	50.04166	2018-01-16T12:00:00Z	186.05602	
#-134.3750	50.04166	2018-01-16T12:00:00Z	188.03001	
#-134.2917	50.04166	2018-01-16T12:00:00Z	207.34081	
#-134.2083	50.04166	2018-01-16T12:00:00Z	205.43228	
#...
#1-10 of 3,484,860 rows

```

## Convert griddap() output to data frame and standardize time

The `griddap()` function returns its results in a list structure where
the actual data are stored in the `$data` element.

```{r}
noaa20_df <- noaa20_ds$data

```

## Filter valid pixels based on `percent_keep`

Before computing trends, we need to ensure that each spatial pixel has
enough valid observations to support a meaningful regression. Pixels
with too many missing values can distort the analysis or produce
unstable estimates. To address this, we keep only those pixels that
contain at least `percent_keep` (here, 50%) non-missing observations
across the full time series.

```{r}
message(sprintf("Filtering valid data points using percent_keep = %g", percent_keep))

# Number of time steps
n_time <- dim(noaa20_df)[3]

# Count non-NA values along the time dimension for each [lon, lat] pixel
valid_counts <- apply(!is.na(noaa20_df), c(1, 2), sum)  # dims: [lon, lat]

# Build a mask of pixels that meet the percent_keep threshold
valid_pixels <- valid_counts >= (percent_keep * n_time)  # dims: [lon, lat]

# Expand that mask across time so we can apply it to the full 3D array
valid_mask_3d <- array(valid_pixels, dim = dim(noaa20_df))  # dims: [lon, lat, time]

# Apply mask: keep only pixels that are valid in all time slices; others → NA
da_filtered <- noaa20_df
da_filtered[!valid_mask_3d] <- NA_real_

# Filtering valid data points using percent_keep = 0.5
```

## Count Valid Time Points

To understand the temporal coverage of each spatial pixel, we count how
many valid (non-missing) observations remain after applying the
percent-keep mask. This helps confirm whether each pixel still contains
usable information for the trends analysis.

```{r}
lon_col  <- "longitude"
lat_col  <- "latitude"
time_col <- "time"

# Count valid (non-NA) time points per pixel
n_df <- da_filtered %>%
  group_by(.data[[lon_col]], .data[[lat_col]]) %>%
  summarise(
    n = sum(!is.na(.data[[ncvar]])),
    .groups = "drop"
  )

# Apply the "n > 0" filter, like n.where(n > 0)
n_df <- n_df %>%
  mutate(n = ifelse(n > 0, n, NA_integer_))

# Inspect unique n values
unique_n <- sort(unique(n_df$n))
message("Unique n values after filtering: ", paste(unique_n, collapse = ", "))

# Unique n values after filtering: 1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 14, 15, 17, 18, 19, 22, 23, 24, 26, 27, 28, 29, # 30, 31, 33, 34, 35, 37, 38, 39, 40, 41, 42, 43, 44, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60

```

## Visualize the number of time points

We are checking whether most pixels are well-sampled through time.

```{r}
n_df %>%
  filter(!is.na(n)) %>%
  ggplot(aes(x = n)) +
  geom_histogram(bins = 20) + 
  labs(
    title = "Histogram of valid time points per pixel",
    x = "Number of valid time points",
    y = "Count of pixels"
  )

```
![](images/hist_of_n.png)

## Convert to fractional years

We convert each timestamp into a fractional year (e.g., 2018.04 instead
of “2018-02-15”). Using a continuous time variable ensures that the
regression slope represents change per year and generally yields more
interpretable and numerically stable trend estimates. Because the
fractional-year scale is larger than a simple month index, the resulting
slope values also become larger and easier to compare across sensors.

```{r}
# Build lookup of unique times → month + fractional year
time_tbl <- da_filtered %>%
  distinct(time) %>%              # one row per unique time
  arrange(time) %>%
  mutate(
    time_posix = ymd_hms(time, tz = "UTC"),   # parse "2018-01-16T12:00:00Z"
    year       = year(time_posix),
    doy        = yday(time_posix),
    years      = year + (doy - 1) / 365.25,   # fractional year
    month      = month(time_posix)            # 1–12
  )

# Quick check that years are fractional and vary
head(time_tbl$years, 10)
length(unique(time_tbl$years))    # should be 60

# [1] 2018.041 2018.126 2018.203 2018.287 2018.370 2018.454 2018.537 2018.621 2018.706 2018.789
# [1] 60
```
## Compute Anomalies
We remove the seasonal cycle by computing a monthly climatology for each
pixel and subtracting it from the observed productivity values. The
resulting anomalies highlight departures from typical monthly conditions
and provide a cleaner signal for trend analysis. We defined ncvar as
"productivity" earlier.
```{r}
da_anom <- da_filtered %>%
  # Attach month + fractional years to every row via join
  left_join(time_tbl %>% select(time, years, month),
            by = "time") %>%
  # Compute monthly climatology per pixel (lon, lat, month)
  group_by(longitude, latitude, month) %>%
  mutate(
    monthly_climatology = mean(.data[[ncvar]], na.rm = TRUE),
    anomaly             = .data[[ncvar]] - monthly_climatology
  ) %>%
  ungroup()

```

## Check that the years are fractional

```{r}
# Check unique fractional years actually look fractional
print(head(unique(da_anom$years), 10))

# [1] 2018.041 2018.126 2018.203 2018.287 2018.370 2018.454 2018.537 2018.621 2018.706 2018.789
```

## Fit a Linear Regression Model

For each pixel, we fit a simple linear regression of anomaly versus
fractional year to estimate the long-term trend. This produces per-pixel
slope, variance, and standard-error values that summarize how
productivity is changing over time.

```{r}
# Global variance of X (years) – same for all pixels
years_vec  <- sort(unique(da_anom$years))
sigma_X_sq <- var(years_vec, na.rm = TRUE)

# Per-pixel slope, Y variance, n, and SE 
trend_results <- da_anom %>%
  group_by(longitude, latitude) %>%
  summarise(
    n = sum(!is.na(anomaly)),
    
    # slope from simple linear regression anomaly ~ years
    slope = if (n > 1) {
      fit <- lm(anomaly ~ years)
      coef(fit)[["years"]]
    } else {
      NA_real_
    },
    
    # per-pixel variance of Y (anomaly) over time
    sigma_Y_sq = if (n > 1) {
      var(anomaly, na.rm = TRUE)
    } else {
      NA_real_
    },
    
    # standard error of the slope
    se = if (n > 2 && is.finite(sigma_Y_sq)) {
      sqrt((1 / (n - 2)) * ((sigma_Y_sq / sigma_X_sq) - slope^2))
    } else {
      NA_real_
    },
    
    .groups = "drop"
  )
  
# Summary stats for slope (ignoring NAs)
slope_valid <- trend_results$slope[is.finite(trend_results$slope)]

min_slope  <- min(slope_valid)
max_slope  <- max(slope_valid)
mean_slope <- mean(slope_valid)

message(sprintf(
  "Minimum slope: %g, Maximum slope: %g, Mean slope: %g",
  min_slope, max_slope, mean_slope
))
  
# Minimum slope: -4896.52, Maximum slope: 961.197, Mean slope: 4.80867
```

## Compute the t-statistic

We calculate the t-statistic for each pixel by dividing the slope by its
standard error, providing a measure of how confidently the trend differs
from zero. Insufficiently sampled pixels are masked to ensure only
meaningful t-values are retained.

```{r}
# Compute t = slope / se, masking invalid cases
trend_results$t_stat <- with(trend_results, {
  t_val <- slope / se
  # Mask non-finite or undefined (n <= 2) values
  t_val[!is.finite(t_val) | n <= 2] <- NA_real_
  t_val
})

# Summary stats for t (ignoring NAs)
t_valid <- trend_results$t_stat[is.finite(trend_results$t_stat)]

min_t  <- min(t_valid)
max_t  <- max(t_valid)
mean_t <- mean(t_valid)

message(sprintf(
  "Minimum t-stat: %g, Maximum t-stat: %g, Mean t-stat: %g",
  min_t, max_t, mean_t
))

# Minimum t-stat: -4.42051, Maximum t-stat: 7.3612, Mean t-stat: 0.888995
```

## Compute the p-value

We convert each pixel’s t-statistic into a two-tailed p-value to assess
whether the slope is statistically different from zero. Pixels with
invalid statistics are masked to ensure only reliable p-values are
included.

```{r}
# Two-tailed p-values for slope (H0: slope = 0)
trend_results$p_value <- with(trend_results, {
  p_val <- 2 * pt(-abs(t_stat), df = n - 2)
  p_val[!is.finite(p_val) | n <= 2] <- NA_real_
  p_val
})

# Summary stats for p (ignoring NAs)
p_valid <- trend_results$p_value[is.finite(trend_results$p_value)]

min_p  <- min(p_valid)
max_p  <- max(p_valid)
mean_p <- mean(p_valid)

message(sprintf(
  "Minimum p-value: %g, Maximum p-value: %g, Mean p-value: %g",
  min_p, max_p, mean_p
))

# Minimum p-value: 7.24423e-10, Maximum p-value: 0.999983, Mean p-value: 0.322973
```

## Build 3-D arrays and write NetCDF

We organize the per-pixel trend results (slope, p-value, and sample
size) into 3-D [time, lat, lon] arrays so they match the structure
expected by ERDDAP-style NetCDF outputs. These arrays are then written
to disk using `write_trends_netcdf()`, producing a properly formatted
file for analysis or distribution.

```{r}
# Coordinate vectors from the anomaly data (West Coast subset)
lon_vals <- sort(unique(da_anom$longitude))
lat_vals <- sort(unique(da_anom$latitude))

n_lon <- length(lon_vals)
n_lat <- length(lat_vals)

# For this product we’re saving ONE trend per pixel over the full period
n_time     <- 1L
time_index <- 1L

# Initialize [time, lat, lon] arrays for beta, p-value, and n
beta_array <- array(NA_real_, dim = c(n_time, n_lat, n_lon))
pval_array <- array(NA_real_, dim = c(n_time, n_lat, n_lon))
n_array    <- array(NA_real_, dim = c(n_time, n_lat, n_lon))

# Map each (lon, lat) row in trend_results into the arrays trend_results has columns: longitude, latitude, slope, p_value, n, ...

lon_idx <- match(trend_results$longitude, lon_vals)
lat_idx <- match(trend_results$latitude, lat_vals)

for (k in seq_len(nrow(trend_results))) {
  i_lon <- lon_idx[k]
  i_lat <- lat_idx[k]

  if (!is.na(i_lon) && !is.na(i_lat)) {
    beta_array[1, i_lat, i_lon] <- trend_results$slope[k]
    pval_array[1, i_lat, i_lon] <- trend_results$p_value[k]
    n_array[1, i_lat, i_lon]    <- trend_results$n[k]
  }
}

beta_array[!is.finite(beta_array)] <- -999.0
pval_array[!is.finite(pval_array)] <- -999.0
n_array[!is.finite(n_array)]       <- -999.0

# Write the NetCDF file `sensor` is already defined earlier 
out_file <- write_trends_netcdf(
  beta_array    = beta_array,
  pval_array    = pval_array,
  n_array       = n_array,
  time_index    = time_index,
  lat_vals      = lat_vals,
  lon_vals      = lon_vals,
  out_dir       = ODATA_DIR,
  overwrite     = TRUE,
  source        = sensor,          # e.g. "noaa20"
  platform_name = "VIIRS-NOAA20"   # used in title/summary
)

# Wrote ERDDAP-like trends NetCDF to:  #/Users/madisonrichardson/netpp/data/trends/trends_noaa20_monthly_9km_lon-134.96_to_-114.96_lat30.04_to_50.04.nc

```
