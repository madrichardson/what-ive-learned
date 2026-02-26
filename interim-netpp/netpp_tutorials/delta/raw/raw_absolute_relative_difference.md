---
output: html_document
date: "2025-11-17"
editor_options: 
  markdown: 
    wrap: 72
---

# Absolute Relative Difference Analysis of the Interim and Legacy Primary Productivity Products

> Update: November 2025

## Objectives

Calculate the mean absolute relative difference $\Delta^{\text{netPP}}$
of the interim VIIRS netPP ($\text{netPP}_{\text{VIIRS}}$) and legacy
MODIS netPP ($\text{netPP}_{\text{MODIS}}$) products for each month from
the timeseries of the user-specified sensors.

We are using this statistic to track the similarities in netPP values
between interim VIIRS and legacy MODIS datasets to provide validation
that the interim netPP products can be reliably used for continuity in
long-term productivity analyses.

For **VIIRS-SNPP**, the result will be a **120-month (10-year)** mean
monthly timeseries of $\Delta^{\text{netPP}}$.

For **VIIRS-NOAA20**, the result will be a **60-month (5-year)** mean
monthly timeseries of $\Delta^{\text{netPP}}$.

The absolute relative difference $\Delta^{\text{netPP}}$ is calculated
for each pixel as follows:

$$
\Delta^{\text{netPP}} = \frac{\text{netPP}_{\text{VIIRS}} - \text{netPP}_{\text{MODIS}}}{\text{netPP}_{\text{MODIS}}}
$$

Where:

-   The difference between VIIRS and MODIS values for that pixel is
    divided by MODIS values.

------------------------------------------------------------------------

## Datasets Overview

We will be creating the pixel-by-pixel absolute relative difference
$\Delta^{\text{netPP}}$ for each month at a 9km resolution across the
VIIRS-NOAA20 and MODIS-Aqua datasets. Here are the available datasets
including VIIRS-SNPP:

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

For this tutorial, we will be using the **VIIRS-NOAA20** and
**MODIS-Aqua** datasets to generate a timeseries of monthly mean values
of $\Delta^{\text{netPP}}$.

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
-   **ODATA_DIR** is where the output delta file will be located.

```{r}
ROOT_DIR <- "/Users/madisonrichardson/netpp"
ODATA_DIR <- file.path(ROOT_DIR, "data", "matrix", "delta")

dir_list <- c(ROOT_DIR, ODATA_DIR)

# Create if missing
invisible(lapply(dir_list, dir.create, showWarnings = FALSE, recursive = TRUE))

cat(length(dir_list), "directories validated\n")

# 2 directories validated
```

## Make Some Useful Functions

```{r}
#' Center timestamps on the 16th (and optionally at midnight)
#'
#' @description
#' For monthly composites, standardize timestamps so they fall on the **16th**
#' of each month. Designed for data frames returned by `rerddap::griddap(...)$data`.
#' Optionally zeroes the time-of-day to `00:00:00` (recommended for joins).
#'
#' @param x Either:
#'   - a `data.frame`/`tibble` with a date-time column, or
#'   - a `POSIXct`/`Date` vector.
#' @param time_col Character. Name of the time column if `x` is a data frame.
#'   Ignored otherwise. Default `"time"`.
#' @param normalize_midnight Logical. If `TRUE`, set hour/min/sec to `00:00:00`.
#'   Default `TRUE` so datasets with different hours (e.g., 08Z vs 12Z) align.
#' @param quiet Logical. If `FALSE`, prints a short message. Default `TRUE`.
#'
#' @returns
#' Object of the same class as `x`, with timestamps on the 16th (and, by default,
#' at midnight). If input is a `Date`, it is returned as `Date`; if `POSIXct`,
#' returned as `POSIXct`. For data frames, only `time_col` is modified.
#'

ensure_time_centered_on_16th <- function(x,
                                         time_col = "time",
                                         normalize_midnight = TRUE,
                                         quiet = TRUE) {
  # Vector method --------------------------------------------------------------
  if (inherits(x, c("POSIXct", "POSIXlt", "Date"))) {
    # Track type to return in same class
    is_date <- inherits(x, "Date")
    tz_attr <- if (!is_date) attr(x, "tzone") else "UTC"
    
    # Work in POSIXlt for component edits
    lt <- if (is_date) as.POSIXlt(x) else as.POSIXlt(x, tz = tz_attr)
    all16 <- all(as.integer(format(lt, "%d")) == 16L, na.rm = TRUE)
    
    if (!all16 && !quiet) message("Replacing day with 16...")
    lt$mday <- 16L
    
    if (normalize_midnight) {
      lt$hour <- 0L; lt$min <- 0L; lt$sec <- 0
    }
    
    if (is_date) {
      return(as.Date(lt))
    } else {
      return(as.POSIXct(lt, tz = tz_attr))
    }
  }
  
  # Data-frame method ----------------------------------------------------------
  if (is.data.frame(x)) {
    if (!time_col %in% names(x)) stop("Column '", time_col, "' not found.")
    t <- x[[time_col]]
    
    if (!inherits(t, c("POSIXct", "POSIXlt", "Date"))) {
      stop("Column '", time_col, "' must be POSIXct/POSIXlt or Date.")
    }
    
    # Preserve original class for the column
    is_date <- inherits(t, "Date")
    tz_attr <- if (!is_date) attr(t, "tzone") else "UTC"
    
    lt <- if (is_date) as.POSIXlt(t) else as.POSIXlt(t, tz = tz_attr)
    all16 <- all(as.integer(format(lt, "%d")) == 16L, na.rm = TRUE)
    
    if (!all16 && !quiet) message("Replacing day with 16...")
    lt$mday <- 16L
    
    if (normalize_midnight) {
      lt$hour <- 0L; lt$min <- 0L; lt$sec <- 0
    }
    
    x[[time_col]] <- if (is_date) as.Date(lt) else as.POSIXct(lt, tz = tz_attr)
    return(x)
  }
  
  stop("Unsupported input: pass a POSIXct/Date vector or a data.frame with '", time_col, "'.")
}

```

```{r}

#' Relative difference (Δ) between two numeric datasets
#'
#' @title Calculate Δ = (A - B) / B with NA masking
#'
#' @description
#' Computes the **relative difference** between a *minuend* dataset
#' `A` and a *subtrahend* dataset `B`, defined as
#' \deqn{\delta = \frac{A - B}{B}.}
#' This mirrors the Python implementation:
#' \code{delta = (minuend_data - subtrahend_data) / subtrahend_data}.
#' Missing values are respected via `NA`; wherever either input is `NA`, the
#' result is `NA`. Division-by-zero (i.e., where \eqn{B = 0}) is handled according
#' to `zero_denominator`.
#'
#' @details
#' - *Minuend* is the quantity **from which** another is subtracted (A).
#' - *Subtrahend* is the quantity **that is subtracted** (B) **and** the
#'   denominator for the relative difference.
#' - Inputs must be numeric vectors, matrices, or arrays of identical shape.
#' - If you are starting from masked arrays (Python `numpy.ma` semantics),
#'   convert masked elements to `NA_real_` before calling this function.
#' - For gridded data classes that support vectorized arithmetic (e.g., **stars**,
#'   **terra**), you can usually pass them directly if they share geometry and
#'   alignment; otherwise extract numeric arrays first and then re-wrap.
#'
#' @param minuend_data Numeric vector/matrix/array. The dataset **A** (minuend).
#' @param subtrahend_data Numeric vector/matrix/array. The dataset **B** (subtrahend),
#'   used as the denominator. Must be the same shape as `minuend_data`.
#' @param zero_denominator Character scalar controlling how to treat elements
#'   where \eqn{B = 0}. One of:
#'   - `"NA"` (default): set Δ to `NA_real_` at those locations.
#'   - `"Inf"`: keep machine `Inf/-Inf` results.
#'   - `"Zero"`: force Delta to `0` at those locations.
#'
#' @returns
#' A numeric object of the same type and dimensions as the inputs (vector/matrix/array)
#' containing Δ, with `NA` wherever inputs were `NA` (or where division-by-zero was
#' mapped to `NA`).
#'
calculate_relative_diff <- function(minuend_data,
                                    subtrahend_data,
                                    zero_denominator = c("NA", "Inf", "Zero")) {

  zero_denominator <- match.arg(zero_denominator)

  if (!is.numeric(minuend_data) || !is.numeric(subtrahend_data)) {
    stop("Both 'minuend_data' and 'subtrahend_data' must be numeric.")
  }
  if (!identical(dim(minuend_data), dim(subtrahend_data))) {
    # For vectors, dim() is NULL; identical(NULL, NULL) is TRUE.
    if (!(is.null(dim(minuend_data)) && is.null(dim(subtrahend_data)) &&
          length(minuend_data) == length(subtrahend_data))) {
      stop("Inputs must have identical shapes (same length/dimensions).")
    }
  }

  # Element-wise numerator and denominator
  denom <- subtrahend_data
  num   <- minuend_data - subtrahend_data

  # Relative difference Δ = (A - B) / B
  rel_diff <- num / denom

  # Handle exact-zero denominators in a controlled way
  denom_zero_idx <- !is.na(denom) & (denom == 0)

  if (any(denom_zero_idx, na.rm = TRUE)) {
    if (zero_denominator == "NA") {
      rel_diff[denom_zero_idx] <- NA_real_
    } else if (zero_denominator == "Zero") {
      rel_diff[denom_zero_idx] <- 0
    } else {
      # "Inf": leave as produced by num/denom (Inf/-Inf or NaN if num also 0)
      invisible(NULL)
    }
  }

  # Ensure NA propagation where either input is NA
  na_mask <- is.na(minuend_data) | is.na(subtrahend_data)
  if (any(na_mask)) rel_diff[na_mask] <- NA_real_

  return(rel_diff)
}


```

```{r}

#' Write Δ (relative difference) to an ERDDAP-style NetCDF file
#'
#' @description
#' Create a CF/ERDDAP-style NetCDF file containing a 3-D field of relative
#' difference (Δ, “delta”) on a regular lon/lat/time grid. The function:
#' \enumerate{
#'   \item Validates the input array dimensions against the coordinate vectors.
#'   \item Derives spatial bounds from the latitude/longitude vectors.
#'   \item Constructs an output filename that encodes sensor names and
#'         lat/lon bounds.
#'   \item Defines CF-style dimensions (`time`, `latitude`, `longitude`) and a
#'         `delta` variable on \code{[lon, lat, time]}.
#'   \item Writes the data to disk and attaches ERDDAP-like coordinate and
#'         global metadata consistent with the delta CDL template.
#' }
#'
#' @details
#' The relative difference is defined as:
#' \deqn{\Delta = \frac{A - B}{B},}
#' where \eqn{A} is the minuend (e.g., VIIRS netPP) and \eqn{B} is the
#' subtrahend (e.g., MODIS netPP).
#'
#' The output filename has the form:
#' \preformatted{
#'   netpp_delta_<sensor1>_<sensor2>_lon<min>_to_<max>_lat<min>_to_<max>.nc
#' }
#' where the lat/lon bounds are taken from \code{lat_vals} and \code{lon_vals}
#' and formatted to one decimal place.
#'
#' Internally, the \code{delta_array} is transposed from
#' \code{[time, lat, lon]} to \code{[lon, lat, time]} before writing, so that
#' the resulting NetCDF is consistent with ERDDAP’s \code{[longitude, latitude,
#' time]} ordering. The \code{time} dimension is created as an unlimited
#' dimension to allow appending in other workflows.
#'
#' The NetCDF file includes:
#' \itemize{
#'   \item CF/ERDDAP-style attributes on \code{time}, \code{latitude},
#'         \code{longitude}, and \code{delta}, following your delta CDL.
#'   \item Geospatial coverage metadata (\code{geospatial_*} and
#'         \code{southernmost_latitude} etc.).
#'   \item A project/title/summary block tailored to the interim
#'         primary productivity use case and the Δ (delta) definition.
#' }
#'
#' @param delta_array Numeric 3-D array of Δ (delta) values with dimensions
#'   \code{[time, lat, lon]}. This is usually the result of computing pixel-wise
#'   relative differences between two gridded products and stacking over time.
#' @param time_secs Numeric vector of length \code{dim(delta_array)[1]} giving
#'   time coordinates in \strong{seconds since 1970-01-01T00:00:00Z}.
#' @param lat_vals Numeric vector of length \code{dim(delta_array)[2]} giving
#'   latitude coordinates in degrees north.
#' @param lon_vals Numeric vector of length \code{dim(delta_array)[3]} giving
#'   longitude coordinates in degrees east (e.g., \code{[-180, 180]} or
#'   \code{[0, 360]}).
#' @param sensor1 Character scalar giving the name of the VIIRS (or other
#'   interim) sensor, e.g. \code{"noaa20"} or \code{"snpp"}. Used in the
#'   filename and title as the minuend dataset in Δ = sensor1 – sensor2,
#'   divided by sensor2.
#' @param sensor2 Character scalar giving the name of the legacy/reference
#'   sensor, default \code{"modis"}. Used in the filename and title as the
#'   subtrahend dataset in Δ = sensor1 – sensor2.
#' @param out_dir Character path to the directory where the NetCDF file will
#'   be written. The directory will be created if it does not exist.
#' @param region_label Character label describing the spatial subset, used
#'   in the global \code{title} and \code{summary} attributes
#'   (e.g. \code{"West Coast subset"}).
#' @param overwrite Logical. If \code{FALSE} (default), the function stops
#'   with an error if a file with the same name already exists. If
#'   \code{TRUE}, an existing file is deleted and recreated.
#'
#' @return
#' Invisibly returns the full file path (\code{character}) of the created
#' NetCDF file. A message is printed to the console indicating where the
#' file was written.
#'
#' @export
write_delta_netcdf <- function(delta_array,
                               time_secs,
                               lat_vals,
                               lon_vals,
                               sensor1,
                               sensor2 = "modis",
                               out_dir,
                               region_label = "West Coast subset",
                               overwrite = FALSE) {
  # Ensure ncdf4 is available
  if (!requireNamespace("ncdf4", quietly = TRUE)) {
    stop("Package 'ncdf4' is required but not installed.")
  }

  # Basic checks --------------------------------------------------------------
  if (!is.numeric(delta_array)) stop("'delta_array' must be numeric.")
  if (length(dim(delta_array)) != 3L) {
    stop("'delta_array' must be a 3-D array [time, lat, lon].")
  }

  n_time <- length(time_secs)
  n_lat  <- length(lat_vals)
  n_lon  <- length(lon_vals)

  if (!identical(dim(delta_array), c(n_time, n_lat, n_lon))) {
    stop(
      "dim(delta_array) must be [length(time_secs), length(lat_vals), length(lon_vals)].\n",
      "Found: ", paste(dim(delta_array), collapse = " x "), " vs ",
      n_time, " x ", n_lat, " x ", n_lon
    )
  }

  # Derive spatial bounds for filename & metadata -----------------------------
  lat_min <- min(lat_vals); lat_max <- max(lat_vals)
  lon_min <- min(lon_vals); lon_max <- max(lon_vals)

  # Compact numeric formatting for filename
  fmt <- function(x) sprintf("%.1f", x)

  nc_filename <- sprintf(
    "netpp_delta_%s_%s_lon%s_to_%s_lat%s_to_%s.nc",
    tolower(sensor1),
    tolower(sensor2),
    fmt(lon_min), fmt(lon_max),
    fmt(lat_min), fmt(lat_max)
  )

  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }
  nc_file_path <- file.path(out_dir, nc_filename)

  # Dimensions and variable definition ----------------------------------------
  dim_time <- ncdf4::ncdim_def(
    name  = "time",
    units = "seconds since 1970-01-01T00:00:00Z",
    vals  = time_secs,
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

  # delta stored as [lon, lat, time] in NetCDF → dims in R order [lon, lat, time]
  var_delta <- ncdf4::ncvar_def(
    name     = "delta",
    units    = "1",  # unitless
    dim      = list(dim_lon, dim_lat, dim_time),
    longname = "Summation for Calculating Relative Difference",
    missval  = -999.0,
    prec     = "float"
  )

  # ------------------------------------------------------------
  # Create NetCDF file and write data
  # ------------------------------------------------------------
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

  nc <- ncdf4::nc_create(nc_file_path, vars = list(var_delta))

  # delta_array is [time, lat, lon]; NetCDF expects [lon, lat, time]
  delta_for_nc <- aperm(delta_array, perm = c(3, 2, 1))
  ncdf4::ncvar_put(nc, "delta", delta_for_nc)

  # Coordinate variables (dimvars auto-created by ncdim_def)
  ncdf4::ncvar_put(nc, "time",      time_secs)
  ncdf4::ncvar_put(nc, "latitude",  lat_vals)
  ncdf4::ncvar_put(nc, "longitude", lon_vals)

  # ------------------------------------------------------------
  # Add ERDDAP-style attributes (matching the delta CDL)
  # ------------------------------------------------------------
  now <- Sys.time()

  # time attributes
  ncdf4::ncatt_put(nc, "time", "_CoordinateAxisType", "Time",      prec = "text")
  ncdf4::ncatt_put(nc, "time", "axis",                "T",         prec = "text")
  ncdf4::ncatt_put(nc, "time", "calendar",            "gregorian", prec = "text")
  ncdf4::ncatt_put(nc, "time", "ioos_category",       "Time",      prec = "text")
  ncdf4::ncatt_put(nc, "time", "long_name",           "Time",      prec = "text")
  ncdf4::ncatt_put(nc, "time", "standard_name",       "time",      prec = "text")
  ncdf4::ncatt_put(nc, "time", "time_origin",
                   "01-JAN-1970 00:00:00",           prec = "text")

  # latitude attributes
  ncdf4::ncatt_put(nc, "latitude", "_CoordinateAxisType", "Lat",       prec = "text")
  ncdf4::ncatt_put(nc, "latitude", "axis",                "Y",         prec = "text")
  ncdf4::ncatt_put(nc, "latitude", "ioos_category",       "Location",  prec = "text")
  ncdf4::ncatt_put(nc, "latitude", "long_name",           "Latitude",  prec = "text")
  ncdf4::ncatt_put(nc, "latitude", "standard_name",       "latitude",  prec = "text")
  ncdf4::ncatt_put(nc, "latitude", "units",               "degrees_north", prec = "text")
  ncdf4::ncatt_put(nc, "latitude", "valid_min",           -90.0,       prec = "float")
  ncdf4::ncatt_put(nc, "latitude", "valid_max",           90.0,        prec = "float")

  # longitude attributes
  ncdf4::ncatt_put(nc, "longitude", "_CoordinateAxisType", "Lon",       prec = "text")
  ncdf4::ncatt_put(nc, "longitude", "axis",                "X",         prec = "text")
  ncdf4::ncatt_put(nc, "longitude", "ioos_category",       "Location",  prec = "text")
  ncdf4::ncatt_put(nc, "longitude", "long_name",           "Longitude", prec = "text")
  ncdf4::ncatt_put(nc, "longitude", "source_name",         "cols",      prec = "text")
  ncdf4::ncatt_put(nc, "longitude", "standard_name",       "longitude", prec = "text")
  ncdf4::ncatt_put(nc, "longitude", "units",               "degrees_east", prec = "text")
  ncdf4::ncatt_put(nc, "longitude", "valid_min",           -180.0,      prec = "float")
  ncdf4::ncatt_put(nc, "longitude", "valid_max",           180.0,       prec = "float")

  # delta variable attributes
  ncdf4::ncatt_put(nc, "delta", "ioos_category", "Other",    prec = "text")
  ncdf4::ncatt_put(nc, "delta", "long_name",     "Summation for Calculating Relative Difference",
                   prec = "text")
  ncdf4::ncatt_put(nc, "delta", "missing_value", -999.0,     prec = "float")
  ncdf4::ncatt_put(nc, "delta", "units",         "1",        prec = "text")

  # Global attributes ---------------------------------------------------------
  lat_step <- if (length(lat_vals) > 1L) abs(diff(lat_vals)[1]) else NA_real_
  lon_step <- if (length(lon_vals) > 1L) abs(diff(lon_vals)[1]) else NA_real_

  ncdf4::ncatt_put(nc, 0, "cdm_data_type",            "Grid",        prec = "text")
  ncdf4::ncatt_put(nc, 0, "map_projection",           "geographic",  prec = "text")
  ncdf4::ncatt_put(nc, 0, "time_coverage_resolution", "PDM",         prec = "text")

  ncdf4::ncatt_put(
    nc, 0, "project",
    "NOAA CoastWatch Interim Primary Productivity",
    prec = "text"
  )

  ncdf4::ncatt_put(
    nc, 0, "title",
    sprintf(
      "Pixel by pixel Delta, VIIRS %s minus %s Aqua",
      toupper(sensor1), toupper(sensor2)
    ),
    prec = "text"
  )

  ncdf4::ncatt_put(
    nc, 0, "summary",
    sprintf(
      paste0(
        "The relative difference (Delta) between primary productivity (netPP) ",
        "calculated using VIIRS %s data and %s Aqua data. For each pixel, delta is ",
        "calculated as netPP(%s) - netPP(%s) divided by netPP(%s). ",
        "Primary productivity was calculated as described by Behrenfeld and Falkowski 1997. ",
        "The data are at 9km resolution. Input data for primary productivity were obtained ",
        "from NASA and included chlorophyll_a, sea surface temperature, and ",
        "photosynthetically active radiation from either MODIS Aqua or VIIRS %s."
      ),
      toupper(sensor1), toupper(sensor2),
      tolower(sensor1), tolower(sensor2), tolower(sensor2),
      toupper(sensor1)
    ),
    prec = "text"
  )

  ncdf4::ncatt_put(
    nc, 0, "date_created",
    format(as.POSIXct(now, tz = "UTC"), "%Y-%m-%dT%H:%M:%SZ"),
    prec = "text"
  )

  # Global coordinate/coverage metadata
  ncdf4::ncatt_put(nc, 0, "latitude_units",    "degrees_north", prec = "text")
  ncdf4::ncatt_put(nc, 0, "longitude_units",   "degrees_east",  prec = "text")

  # Use actual subset for these, but structure matches CDL
  ncdf4::ncatt_put(nc, 0, "northernmost_latitude",   lat_max,  prec = "float")
  ncdf4::ncatt_put(nc, 0, "southernmost_latitude",   lat_min,  prec = "float")
  ncdf4::ncatt_put(nc, 0, "westernmost_longitude",   lon_min,  prec = "float")
  ncdf4::ncatt_put(nc, 0, "easternmost_longitude",   lon_max,  prec = "float")
  ncdf4::ncatt_put(nc, 0, "geospatial_lat_max",      lat_max,  prec = "float")
  ncdf4::ncatt_put(nc, 0, "geospatial_lat_min",      lat_min,  prec = "float")
  ncdf4::ncatt_put(nc, 0, "geospatial_lon_max",      lon_max,  prec = "float")
  ncdf4::ncatt_put(nc, 0, "geospatial_lon_min",      lon_min,  prec = "float")

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

  # From CDL: spatialResolution + proj4_string
  ncdf4::ncatt_put(nc, 0, "spatialResolution", "9.28 km", prec = "text")
  ncdf4::ncatt_put(
    nc, 0, "proj4_string",
    "+proj=eqc +lat_ts=0 +lat_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs +lon_0=0.000000",
    prec = "text"
  )

  ncdf4::ncatt_put(nc, 0, "creator_name",
                   "NOAA CoastWatch West Coast Node", prec = "text")
  ncdf4::ncatt_put(nc, 0, "creator_url",
                   "https://coastwatch.pfeg.noaa.gov/", prec = "text")

  ncdf4::nc_close(nc)

  message("Wrote ERDDAP-like Delta NetCDF to: ", nc_file_path)

  invisible(nc_file_path)
}

```

## Set Up Parameters

Initialize the parameters required for the absolute difference
calculation:

-   **start_date** and **end_date**: Define the date range for the
    analysis (format: YYYY-MM-DDTXX:XX:XXZ). Monthly composites are
    centered on the 16th of each month.

-   **sensor1** and **sensor2**: Specify the two satellite sensors for
    unbiased relative difference analysis ('modis', 'snpp', or
    'noaa20').

-   **ncvar**: The name of the variable being analyzed (we are using
    'productivity').

-   **overwrite**: If 'True", existing output files will be replaced.

```{r}
start_date_noaa20 = "2018-01-16T12:00:00Z"
end_date_noaa20 = "2022-12-16T12:00:00Z"

start_date_modis = "2018-01-16T08:00:00Z"
end_date_modis = "2022-12-16T08:00:00Z"

# Define sensor (either "snpp", "noaa20", or "modis")
sensor1 = "noaa20"
sensor2 = "modis"

# Valid sensors list
valid_sensors <- c("modis", "snpp", "noaa20")

# Check both sensors are valid
for (s in c(sensor1, sensor2)) {
  if (!s %in% valid_sensors) {
    stop(sprintf("The sensor '%s' is not a valid option. Choose from: %s",
                 s, paste(valid_sensors, collapse = ", ")))
  }
}

cat("Sensors validated:", sensor1, "and", sensor2, "\n")

ncvar = "productivity"
overwrite = FALSE

# Sensors validated: noaa20 and modis 
```

## Select Satellite Dataset from ERDDAP

We will be using the VIIRS-NOAA20 - MODIS-AQUA netpp datasets from the
West Coast Node ERDDAP Server. The dataset IDs are:
**productivity_viirs_noaa20_monthly_9km** and
**productivity_modis_aqua_monthly_9km**. We will use the info function
from the **rerddap** package to first obtain information about the
dataset of interest, then we will import the data.

```{r}
# Set ERDDAP URL
erddap_url = "https://coastwatch.pfeg.noaa.gov/wcn/erddap/"

# Set dataset IDs
noaa20_id = "productivity_viirs_noaa20_monthly_9km"
modis_id = "productivity_modis_aqua_monthly_9km"


# Get datasets info
noaa20_dataInfo <- rerddap::info(noaa20_id, url=erddap_url)  
print(noaa20_dataInfo)

modis_dataInfo <- rerddap::info(modis_id, url=erddap_url)  
print(modis_dataInfo)

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
#<ERDDAP info> productivity_modis_aqua_monthly_9km 
# Base URL: https://coastwatch.pfeg.noaa.gov/wcn/erddap 
# Dataset Type: griddap 
# Dimensions (range):  
#     time: (2013-01-16T08:00:00Z, 2022-12-16T08:00:00Z) 
#     latitude: (-89.95834, 89.95834) 
#     longitude: (-179.9583, 179.9583) 
# Variables:  
#     productivity: 
#         Units: mg C m-2 d-1 
```

## Subset Datasets

The datasets are too large to use the entire globe, so we are subsetting
to use the West Coast of the US.

```{r}
# Define West Coast latitude and longitudes
wc_lon <- c(-135, -115)
wc_lat <- c(30, 50)

# Subset noaa20 dataset
noaa20_ds <- griddap(
  datasetx = noaa20_id,
  url = erddap_url,
  time = c(start_date_noaa20, end_date_noaa20),
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

# Subset modis dataset
modis_ds <- griddap(
  datasetx = modis_id,
  url = erddap_url,
  time = c(start_date_modis, end_date_modis),
  longitude = wc_lon,
  latitude = wc_lat,
  fields = ncvar
)

print(modis_ds)

# longitude latitude time productivity
#-134.9583	50.04166	2018-01-16T08:00:00Z	NA	
#-134.8750	50.04166	2018-01-16T08:00:00Z	190.3552	
#-134.7917	50.04166	2018-01-16T08:00:00Z	230.7118	
#-134.7083	50.04166	2018-01-16T08:00:00Z	183.1599	
#-134.6250	50.04166	2018-01-16T08:00:00Z	183.2100	
#-134.5417	50.04166	2018-01-16T08:00:00Z	192.2304	
#-134.4583	50.04166	2018-01-16T08:00:00Z	182.2389	
#-134.3750	50.04166	2018-01-16T08:00:00Z	178.6511	
#-134.2917	50.04166	2018-01-16T08:00:00Z	183.1943	
#-134.2083	50.04166	2018-01-16T08:00:00Z	192.7871
```

## Convert griddap() output to data frames and standardize time

The `griddap()` function returns its results in a list structure where
the actual data are stored in the `$data` element. Here, we extract
those data tables for VIIRS-NOAA20 and MODIS-Aqua and convert the `time`
column to `POSIXct` in UTC so the two products use a consistent time
format before we align and compare them.

```{r}
noaa20_df <- noaa20_ds$data
modis_df <- modis_ds$data


# Make sure time is POSIXct before calling the function
noaa20_df$time <- as.POSIXct(noaa20_df$time, tz = "UTC")
modis_df$time <- as.POSIXct(modis_df$time, tz = "UTC")

```

## Normalize timestamps to monthly centers (16th of each month)

The MODIS and VIIRS monthly composites do not always share identical
timestamp formats, some products use 08:00Z, others 12:00Z, and
occasionally different days within the month. To ensure a clean
one-to-one match between the two datasets, we standardize all timestamps
so they fall on the **16th of each month at 00:00 UTC**, which is the
convention used in ERDDAP monthly products.

After centering the timestamps, we create a `date_key` column that
stores only the date (YYYY-MM-16). This key is used to merge MODIS and
VIIRS pixels by their shared month.

```{r}
noaa20_df <- ensure_time_centered_on_16th(noaa20_df, time_col = "time")
modis_df <- ensure_time_centered_on_16th(modis_df, time_col = "time")

# Make a date key (now all should be YYYY-mm-16)
modis_df   <- modis_df   %>% mutate(date_key = as.Date(time))
noaa20_df <- noaa20_df %>% mutate(date_key = as.Date(time))

```

## Identify months with data from both sensors

Before computing $\Delta$, we restrict the analysis to months where
**both** MODIS-Aqua and VIIRS-NOAA20 have valid data. We do this by
intersecting the `date_key` columns from each dataset (each key
represents a monthly composite centered on the 16th). The resulting
`common_dates` vector defines the shared time axis for all subsequent
$\Delta$ calculations.

We also include a quick check that `ncvar` is a single character string
(for example `"productivity"`), ensuring that the code references a
single, well-defined variable from each data frame.

```{r}
# Build the common date set
common_dates <- sort(intersect(unique(modis_df$date_key), unique(noaa20_df$date_key)))

message("Found ", length(common_dates), " matching dates")
stopifnot(length(common_dates) > 0)   # fail fast if nothing matches

# ncvar should be a single string with the variable name, e.g. "productivity"
stopifnot(is.character(ncvar), length(ncvar) == 1)

# Found 60 matching dates
```

## Calculate monthly pixel-wise $\Delta$ (VIIRS – MODIS)

For each month where both MODIS-Aqua and VIIRS-NOAA20 have data, we
compute the absolute relative difference delta at every grid cell.

1.  We subset each dataset to a single month (`dt`) and keep only
    latitude, longitude, and the productivity variable.
2.  We join the MODIS and VIIRS slices on their shared
    `(latitude, longitude)` coordinates so that each row represents a
    single pixel observed by both sensors.
3.  We convert the productivity values to numeric (defensive step in
    case they arrive as character/factor).
4.  We compute $\Delta$ for each pixel using: $$
      \Delta^{\text{netPP}} = \frac{\text{netPP}_{\text{VIIRS}} -     \text{netPP}_{\text{MODIS}}}{\text{netPP}_{\text{MODIS}}}
      $$
5.  We reshape the $\Delta$ values into a 2-D matrix with rows ordered
    by latitude (south → north) and columns by longitude (west → east).
6.  For each month, we store the Unix timestamp, $\Delta$ matrix, and
    associated lat/lon grid in a list. This list becomes the basis for
    building the NetCDF file and subsequent timeseries and map
    visualizations.

```{r}
delta_results <- list()

for (kk in seq_along(common_dates)) {
  dt <- common_dates[kk]
  
  # Select coords + the variable, then rename value columns to keep them distinct
  s_slice <- modis_df  |>
    dplyr::filter(date_key == dt) |>
    dplyr::select(latitude, longitude, !!all_of(ncvar)) |>
    dplyr::rename(value_modis = !!sym(ncvar))
  
  n_slice <- noaa20_df |>
    dplyr::filter(date_key == dt) |>
    dplyr::select(latitude, longitude, !!all_of(ncvar)) |>
    dplyr::rename(value_noaa20 = !!sym(ncvar))
  
  # Intersection on (lat, lon)
  paired <- dplyr::inner_join(s_slice, n_slice, by = c("latitude", "longitude"))
  if (nrow(paired) == 0) {
    message("No overlapping cells on ", dt, "; skipping")
    next
  }
  
  # Ensure numeric (in case the field came in as character/factor)
  paired <- paired |>
    mutate(
      value_modis   = as.numeric(value_modis),
      value_noaa20 = as.numeric(value_noaa20)
    )
  
  # Delta = (A - B) / B 
  paired <- paired |>
    mutate(delta = calculate_relative_diff(value_noaa20, value_modis, zero_denominator = "NA"))
  
  # Build Delta matrix: rows = latitude (S→N), cols = longitude (W→E)
  lon_vals <- sort(unique(paired$longitude))
  lat_vals <- sort(unique(paired$latitude))
  delta_mat  <- matrix(NA_real_, nrow = length(lat_vals), ncol = length(lon_vals))
  
  r_idx <- match(paired$latitude,  lat_vals)
  c_idx <- match(paired$longitude, lon_vals)
  delta_mat[cbind(r_idx, c_idx)] <- paired$delta
  
  # Unix timestamp @ 16th 00:00:00 UTC (already normalized earlier)
  ts_utc <- as.numeric(as.POSIXct(dt, tz = "UTC"))
  
  delta_results[[length(delta_results) + 1L]] <- list(
    timestamp = ts_utc,
    delta     = delta_mat,
    lon       = lon_vals,
    lat       = lat_vals
  )
  
  message("Appended Delta for ", format(as.Date(dt, origin="1970-01-01"), "%Y-%m-%d"))
}

# Appended Delta for 2018-01-16
#Appended Delta for 2018-02-16
#Appended Delta for 2018-03-16
#Appended Delta for 2018-04-16
#Appended Delta for 2018-05-16
#Appended Delta for 2018-06-16
#Appended Delta for 2018-07-16
#Appended Delta for 2018-08-16
#Appended Delta for 2018-09-16
#Appended Delta for 2018-10-16
#Appended Delta for 2018-11-16
#Appended Delta for 2018-12-16
#Appended Delta for 2019-01-16
#Appended Delta for 2019-02-16
#Appended Delta for 2019-03-16
#Appended Delta for 2019-04-16
#Appended Delta for 2019-05-16
#Appended Delta for 2019-06-16
#Appended Delta for 2019-07-16
#Appended Delta for 2019-08-16
#Appended Delta for 2019-09-16
#Appended Delta for 2019-10-16
#Appended Delta for 2019-11-16
#Appended Delta for 2019-12-16
#Appended Delta for 2020-01-16
#Appended Delta for 2020-02-16
#Appended Delta for 2020-03-16
#Appended Delta for 2020-04-16
#Appended Delta for 2020-05-16
#Appended Delta for 2020-06-16
#Appended Delta for 2020-07-16
#Appended Delta for 2020-08-16
#Appended Delta for 2020-09-16
#Appended Delta for 2020-10-16
#Appended Delta for 2020-11-16
#Appended Delta for 2020-12-16
#Appended Delta for 2021-01-16
#Appended Delta for 2021-02-16
#Appended Delta for 2021-03-16
#Appended Delta for 2021-04-16
#Appended Delta for 2021-05-16
#Appended Delta for 2021-06-16
#Appended Delta for 2021-07-16
#Appended Delta for 2021-08-16
#Appended Delta for 2021-09-16
#Appended Delta for 2021-10-16
#Appended Delta for 2021-11-16
#Appended Delta for 2021-12-16
#Appended Delta for 2022-01-16
#Appended Delta for 2022-02-16
#Appended Delta for 2022-03-16
#Appended Delta for 2022-04-16
#Appended Delta for 2022-05-16
#Appended Delta for 2022-06-16
#Appended Delta for 2022-07-16
#Appended Delta for 2022-08-16
#Appended Delta for 2022-09-16
#Appended Delta for 2022-10-16
#Appended Delta for 2022-11-16
#Appended Delta for 2022-12-16
```

## Build a 3D $\Delta$ array for NetCDF output

At this stage, `delta_results` is a list where each element contains the
$\Delta$ field for a single month:

-   `delta_results[[k]]$delta` – 2D matrix of $\Delta$ values on a
    lat–lon grid
-   `delta_results[[k]]$lat` – latitude vector for that month
-   `delta_results[[k]]$lon` – longitude vector for that month
-   `delta_results[[k]]$timestamp` – Unix time stamp for the monthly
    center

To write an ERDDAP-style NetCDF file, we first need to:

1.  Extract the **time axis** as a numeric vector (seconds since
    1970-01-01).
2.  Confirm that all monthly $\Delta$ fields share the **same latitude
    and longitude grid** (so they can be stacked into a single cube).
3.  Assemble a 3D array `delta_array[time, lat, lon]` by stacking the
    monthly $\Delta$ matrices along the time dimension.

This 3D array, together with `time_secs`, `lat_vals`, and `lon_vals`,
forms the core data structure that we pass into the NetCDF writer.

```{r}
n_time <- length(delta_results)
stopifnot(n_time > 0L)

time_secs <- vapply(
  delta_results,
  function(x) as.numeric(x$timestamp),
  numeric(1)
)

# Take grid directly from first delta_results entry
lat_vals <- delta_results[[1]]$lat
lon_vals <- delta_results[[1]]$lon

n_lat <- length(lat_vals)
n_lon <- length(lon_vals)

delta_array <- array(NA_real_, dim = c(n_time, n_lat, n_lon))

for (ti in seq_len(n_time)) {
  slab     <- delta_results[[ti]]$delta   # matrix [lat, lon]
  slab_lat <- delta_results[[ti]]$lat
  slab_lon <- delta_results[[ti]]$lon
  
  # Ensure grid is consistent across all time slices
  if (!identical(slab_lat, lat_vals) || !identical(slab_lon, lon_vals)) {
    stop("delta_results[[", ti, "]] lat/lon grid does not match the first slice.")
  }
  if (!identical(dim(slab), c(n_lat, n_lon))) {
    stop("delta_results[[", ti, "]]$delta has dim ",
         paste(dim(slab), collapse = "x"),
         " but expected ", n_lat, "x", n_lon)
  }
  
  # Directly stack on [time, lat, lon]
  delta_array[ti, , ] <- slab
}

```

## Save results to a NetCDF file

Now that we have a complete 3D $\Delta$ cube
(`delta_array[time, lat, lon]`) and the corresponding coordinate
vectors, we write the results to disk as a CF/ERDDAP-style NetCDF file.
The helper function `write_delta_netcdf()`:

-   Validates that the array and coordinate dimensions are consistent.
-   Derives the latitude/longitude bounds from `lat_vals` and
    `lon_vals`.
-   Builds a descriptive filename that encodes the sensors and spatial
    subset.
-   Creates a NetCDF file with `time`, `latitude`, `longitude`, and
    `delta` variables in ERDDAP-friendly order.
-   Attaches global metadata describing the project, region, and
    geospatial coverage.

The returned `out_path` points to the NetCDF file on disk, which we can
then use for plotting, sharing, or publishing to ERDDAP.

```{r}
out_path <- write_delta_netcdf(
  delta_array   = delta_array,          # [time, lat, lon]
  time_secs   = time_secs,          # numeric seconds since 1970-01-01
  lat_vals    = lat_vals,
  lon_vals    = lon_vals,
  sensor1     = sensor1,            # "noaa20"
  sensor2     = sensor2,            # "modis"
  out_dir     = ODATA_DIR,
  region_label = "West Coast subset",
  overwrite   = TRUE                # or FALSE if you want safety
)

# Wrote ERDDAP-like Delta NetCDF to: /Users/madisonrichardson/netpp/data/matrix/delta/netpp_delta_noaa20_modis_lon-134.9_to_-115.1_lat30.1_to_49.9.nc
```

## Read $\Delta$ NetCDF and reconstruct the time axis

To generate summary statistics and visualizations, we first reopen the
NetCDF file we just wrote and pull out the core pieces we need:

-   3D $\Delta$ array on the longitude, latitude, time grid
-   1D longitude and latitude coordinate vectors
-   Raw time coordinate

The `time` variable is stored as an offset from 1970-01-01T00:00:00Z
(either in days or seconds, depending on the writer). We convert this
numeric time axis back into calendar dates so that we can group $\Delta$
by month and build timeseries plots.

```{r}
delta_nc <- nc_open(out_path)

# delta is stored as [lon, lat, time] based on your ncvar_def
delta_arr  <- ncvar_get(delta_nc, "delta")        # 3D array
lon_vals <- ncvar_get(delta_nc, "longitude")  # 1D
lat_vals <- ncvar_get(delta_nc, "latitude")   # 1D
time_raw <- ncvar_get(delta_nc, "time")       # seconds since 1970-01-01T00:00:00Z

nc_close(delta_nc)

# Heuristic: treat values < 1e6 as days since 1970, otherwise seconds
if (max(time_raw, na.rm = TRUE) < 1e6) {
  # time_raw is in days since 1970-01-01
  date_vec <- as.Date(time_raw, origin = "1970-01-01")
} else {
  # time_raw is in seconds since 1970-01-01
  date_vec <- as.Date(
    as.POSIXct(time_raw, origin = "1970-01-01", tz = "UTC")
  )
}

```

## Convert the $\Delta$ cube to a monthly mean time series

The NetCDF file stores $\Delta$ as a full 3D field over longitude,
latitude, and time. To build a 1D time series of spatially averaged
$\Delta$, we:

1.  Expand the lon–lat–time grid into a long-format data frame where
    each row represents a single pixel at a single month.
2.  Flatten the 3D $\Delta$ array into a vector and attach it to this
    grid in the same order.
3.  Attach the corresponding calendar date for each time index.
4.  Group by year and month and compute the **mean** $\Delta$ across all
    pixels for each month, producing a single $\Delta$ value per month.
5.  Construct a plotting date (using the 16th of each month) so the time
    series can be easily visualized and compared to other monthly
    products.

```{r}
# Expand grid over lon, lat, and *index* of time, then attach date
ntime <- length(date_vec)

delta_df <- expand.grid(
  longitude = lon_vals,
  latitude  = lat_vals,
  t_index   = seq_len(ntime)
)

# Flatten delta_arr in a matching order: [lon, lat, time]
delta_df$delta <- as.vector(delta_arr)

# Attach the Date for each t_index
delta_df$date <- date_vec[delta_df$t_index]

# At this point, delta_df$date should *not* be 1970-01-01 everywhere:
# head(delta_df$date)

monthly_data <- delta_df %>%
  mutate(
    year  = year(date),
    month = month(date, label = TRUE)
  ) %>%
  group_by(year, month) %>%
  summarise(
    monthly_mean = mean(delta, na.rm = TRUE),
    .groups = "drop"
  )

# Convert month abbrev → numeric
monthly_data$month_num <- match(monthly_data$month, month.abb)

# Build a plotting date (16th of each month)
monthly_data$date <- as.Date(
  sprintf("%04d-%02d-16", as.integer(monthly_data$year), monthly_data$month_num)
)

```

## Visualize time series of spatially averaged $\Delta$

Finally, we visualize the monthly mean $\Delta$ values as a time series.
Each point represents the **spatial average of** $\Delta$ over the West
Coast subset for a given month, so values near zero indicate good
agreement between VIIRS and MODIS netPP, while larger positive or
negative values indicate months where one sensor consistently reports
higher or lower productivity than the other.

We also add horizontal dashed lines at ±0.05 to highlight a ±5%
difference band. Points that fall within this band correspond to months
where the two products generally agree within about 5%, whereas values
outside this range indicate larger differences.

```{r}
ggplot(monthly_data, aes(x = date, y = monthly_mean)) +
  geom_point(color = "black", size = 2) +
  geom_line(color = "black", linewidth = 1) +
  geom_hline(yintercept = 0.05,  linetype = "dashed", color = "red") +
  geom_hline(yintercept = -0.05, linetype = "dashed", color = "blue") +
  scale_x_date(
    date_breaks = "3 month",        # tick every month
    date_labels = "%b %Y"           # e.g., "Jan 2018"
  ) +
  labs(
    title = expression(bar(Delta)[month] * " Timeseries: Absolute Relative Difference NOAA-20 - MODIS (West Coast)"),
    x = "Month",
    y = "Absolute Relative Difference (Delta)"
  ) +
  ylim(-0.1, 0.1) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)  # tilt labels so they fit
  )


```

![](images/raw_delta_WC.png)

## Visualize the spatial variability of mean $\Delta$

In addition to the time series, it is useful to see **where** VIIRS and
MODIS tend to differ the most. Here we collapse the full $\Delta$ cube
across time and compute the **long-term mean** $\Delta$ at each grid
cell:

-   Positive values (red) indicate locations where VIIRS netPP is, on
    average, higher than MODIS.
-   Negative values (blue) indicate locations where VIIRS netPP is, on
    average, lower than MODIS.
-   Values near zero (white) correspond to areas of good long-term
    agreement.

We then plot this mean $\Delta$ field over the West Coast subset. The
color scale is clipped to the range $$-0.2, 0.2$$ to highlight
differences within ±20%, and the coastline is added for geographic
context.

```{r}
mean_map_df <- delta_df %>%
  group_by(longitude, latitude) %>%
  summarise(
    delta_mean = mean(delta, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(mean_map_df, aes(x = longitude, y = latitude, fill = delta_mean)) +
  geom_raster() +
  borders("world",
          xlim = c(-135, -115),
          ylim = c(30, 50),
          colour = "black",
          size = 0.4) +
  coord_quickmap(xlim = c(-135, -115), ylim = c(30, 50)) +
  scale_fill_gradient2(
    name = expression(bar(Delta)),
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(-0.2, 0.2),
    na.value = "grey80"
  ) +
  labs(
    title = expression("Spatial Variability " * bar(Delta) * " (NOAA-20 – MODIS), West Coast"),
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal()

```

![](images/raw_delta_WC_map.png)
