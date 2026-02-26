---
title: "Contingency Tables: Assessing the Similarity Between Two Primary Productivity (NetPP) Products"
output: html_document
date: "2025-09-12"
---

## Objectives

In this tutorial, we will construct and analyze contingency matrices that compare trends from the legacy primary productivity product, MODIS-Aqua, and an interim product, VIIRS-SNPP or VIIRS-NOAA20.

We have already calculated the long-term linear trend for NetPP using linear regression analysis. Key details of the analysis include:

-   **Pixel-by-Pixel Analysis:** Linear regression was applied over time to each pixel for both legacy and interim NetPP products.

-   **Output Results:** Pixel-by-pixel values were produced for:

    -   The slope of the linear trend.

    -   The number of observations used to calculate the trend (n).

    -   The p-values for statistical significance.

    -   Only pixels with a time series containing 50% of the maximum n-values were included in the analysis.

For this tutorial, we will be comparing MODIS-AQUA and VIIRS-SNPP.

### Constructing a Contingency Matrix

Using the results of the regression analysis, a contingency matrix is constructed to compare the trends between the legacy and interim product. Pixels are categorized by the matrix into three groups based on the sign of the slope:

1.  **Positive Trend:** Increasing trend (positive slope).

2.  **Negative Trend:** Decreasing trend (negative slope).

3.  **No significant Trend:** Slope not significantly different than zero.

The contingency matrix shows the percentage of a pixels where the two NetPP products agree or disagree on the sign of the slope.

### Quantifying Agreement with Cohen's Kappa

The contingency matrix is used to compute Cohen's Kappa, a statistical measure that quantifies the level of agreement between the two NetPP products. Cohen's Kappa accounts for the agreement occurring by chance and provides a standardized metric of similarity.

## Datasets Overview

1.  **Trend coefficients and p values for monthly primary productivity from MODIS-AQUA globally at 9km resolution**

-   Distributed via the West Coast Node ERDDAP dataset at the following link: \> <http://localhost:8080/erddap/griddap/trends_modis_monthly_9km.graph>

2.  **Trend coefficients and p values for monthly primary productivity from VIIRS-SNPP globally at 9km resolution**

-   Distributed via the West Coast Node ERDDAP dataset at the following link: \> <http://localhost:8080/erddap/griddap/trends_snpp_monthly_9km.graph>

3.  **Trend coefficients and p values for monthly primary productivity from VIIRS-NOAA20 globally at 9km resolution**

-   Distributed via the West Coast Node ERDDAP dataset at the following link: \> <http://localhost:8080/erddap/griddap/trends_noaa20_monthly_9km.graph>

## Shapefiles

#### Longhurst Marine Provinces

The dataset represents the division of the world oceans into provinces as defined by Longhurst (1995; 1998; 2006). This division has been based on the prevailing role of physical forcing as a regulator of phytoplankton distribution. The Longhurst Marine Provinces dataset is available online (<https://www.marineregions.org/downloads.php>) and within the shapes folder associated with this repository.

![](images/longhurst.png)

**For our example we will use the shapefile for the "California Upwelling Coastal Province" (ProvCode: CCAL) within the Longhurst Marine Provinces classification**.

## Resource requirements

-   **R Studio** with the modules included within the *Install and Load Required Packages* section below

-   **Shapefile** of your area of interest

    -   If you don't have shapefile, we will include some workarounds in the notebook.

-   **Internet connection**

## Install and Load Required Packages

```{r setup, message=FALSE, warning=FALSE}
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
                      "lubridate", "tidyr", "psych", "gridExtra",
                      "grid")

# Run install and load function
for (pk in list.of.packages) {
  pkgTest(pk)
}

# create list of installed packages
pkges = installed.packages()[,"Package"]


```

### Define Extraction Method: Bounding Box or Longhurst Province

For this tutorial, we are setting **use_bbox = FALSE** because we will be looking at the CCAL province. If you do not have the Longhurst Province shapefiles, set **use_bbox = TRUE** and manually define the bounding box of interest.

```{r define-method, message=FALSE, warning=FALSE}
# User option: Use bounding box or province-based extraction
use_bbox <- FALSE  

# Path to shapefile
shapefile_path <- "/Users/madisonrichardson/netpp/resources/Longhurst/Longhurst_world_v4_2010.shp"

# Read shapefile
shapes <- read_sf(dsn = shapefile_path, layer = "Longhurst_world_v4_2010")

# Example List of all the province names
shapes$ProvCode

if (!use_bbox) {
  # Set Province Code
  ProvCode <- "CCAL"
  
  # Extract the province region
  selected_region <- shapes[shapes$ProvCode == ProvCode,]
  
  # Get bounding box of the province
  bbox <- st_bbox(selected_region)
  lon_min <- bbox["xmin"]
  lon_max <- bbox["xmax"]
  lat_min <- bbox["ymin"]
  lat_max <- bbox["ymax"]
  
  # Extract longitude & latitude for polygon
  longitude <- st_coordinates(selected_region)[,1]
  latitude  <- st_coordinates(selected_region)[,2]
  
} else {
  # Manually set bounding box
  lon_min <- -128.0
  lon_max <- -124.0
  lat_min <-  42.0
  lat_max <-  46.0
  
  longitude <- c(lon_min, lon_max, lon_max, lon_min, lon_min)
  latitude  <- c(lat_min, lat_min, lat_max, lat_max, lat_min)
}

# Print bounding box
print(paste("Bounding Box:", lon_min, lon_max, lat_min, lat_max))

```

## Select Satellite Dataset from ERDDAP

We will be using two datasets, the MODIS-AQUA trends and VIIRS-NOAA20 trends dataset from the West Coast Node ERDDAP Server. The MODIS-AQUA dataset ID is: **trends_modis_monthly_9km**. The VIIRS-SNPP dataset ID is: **trends_snpp_monthly_9km**. We will use the info function from the **rerddap** package to first obtain information about the dataset of interest, then we will import the data.

```{r sat-data, message=FALSE, warning=FALSE}
# Set ERDDAP URL
erddap_url = "http://localhost:8080/erddap"

# Get dataset info
dataInfo1 <- rerddap::info('trends_modis_monthly_9km', url=erddap_url)
print(dataInfo1)

dataInfo2 <- rerddap::info('trends_snpp_monthly_9km', url=erddap_url)
print(dataInfo2)

```

## Fetching Data from ERDDAP

Using the function **griddap**, we fetch gridded data from an ERDDAP server

```{r extract, message=FALSE, warning=FALSE}
# Set ERDDAP datasets IDs
dataset_ID1 <- "trends_modis_monthly_9km"
dataset_ID2 <- "trends_snpp_monthly_9km"

# Retrieve MODIS data
dataset1 <- griddap(
    dataset_ID1,
    url = erddap_url,
    latitude = c(lat_min, lat_max),
    longitude = c(lon_min, lon_max),
    fields = c("beta", "pval")
)$data

# Retrieve SNPP data
dataset2 <- griddap(
    dataset_ID2,
    url = erddap_url,
    latitude = c(lat_min, lat_max),
    longitude = c(lon_min, lon_max),
    fields = c("beta", "pval")
)$data

# View the data
print(dataset1)
print(dataset2)

```

## Extract Trends (beta) and Statistical Significance (pval)

Each dataset contains gridded values for `beta` and `pval`, which are store in separate objects for further comparison.

```{r beta-pval, message=FALSE, warning=FALSE}
# Extract beta and pval from MODIS dataset
beta1 <- dataset1$beta
pval1 <- dataset1$pval

# Extract beta and pval from SNPP dataset
beta2 <- dataset2$beta
pval2 <- dataset2$pval

```

## Finding Common Grid Points for Comparison

To compare trends between MODIS and SNPP datasets, we need to identify grid points that have valid values in both datasets. The function `find_common_points` filters out missing (NA) values and extracts the corresponding `beta` and `pval` values from both datasets.

```{r common-points, message=FALSE, warning=FALSE}
# Function to find common grid points
find_common_points <- function(beta1, beta2, pval1, pval2) {
  mask <- !is.na(beta1) & !is.na(beta2)
  beta_common <- list(beta1[mask], beta2[mask])
  pval_common <- list(pval1[mask], pval2[mask])
  num_common <- sum(mask)
  return(list(beta_common = beta_common, pval_common = pval_common, num_common = num_common))
}

# Compute common grid points
common_data <- find_common_points(beta1, beta2, pval1, pval2)

```

## Creating a 2x2 Contingency Table

The 2x2 contingency table categorizes and compares trend directions (β values) from the legacy (MODIS) and interim (SNPP) datasets at common grid points. The table categorizes pixels based on whether both datasets agree on the trend direction being either **positive** ($\beta>=0$) or **negative** ($\beta<0$).

The function `create_2x2_contingency_table` constructs a 2x2 matrix that classifies each grid point based on positive or negative trends in each of the datasets, and computes the common grid points in each category.

```{r 2x2, message=FALSE, warning=FALSE}
# Function to create 2x2 contingency table
create_2x2_contingency_table <- function(beta_common, num_common) {
  table <- matrix(0, 2, 2)
  rownames(table) <- c("MODIS β >= 0", "MODIS β < 0")
  colnames(table) <- c("SNPP β >= 0", "SNPP β < 0")
  
  table[1, 1] <- sum(beta_common[[1]] >= 0 & beta_common[[2]] >= 0) / num_common * 100
  table[1, 2] <- sum(beta_common[[1]] >= 0 & beta_common[[2]] < 0) / num_common * 100
  table[2, 1] <- sum(beta_common[[1]] < 0 & beta_common[[2]] >= 0) / num_common * 100
  table[2, 2] <- sum(beta_common[[1]] < 0 & beta_common[[2]] < 0) / num_common * 100
  
  return(table)
}

# Create 2x2 contingency table
contingency_table_2x2 <- create_2x2_contingency_table(common_data$beta_common, common_data$num_common)

# Print 2x2 contingency table
print(contingency_table_2x2)

```

## Creating a 3x3 Contingency Table

The 3x3 contingency table categorizes and compares trend directions (β values) from the legacy (MODIS) and interim (SNPP) datasets at common grid points and incorporates statistical significance (p-values) comparing agreement in non-significant trends.

The function `create_3x3_contingency_table` constructs a 3x3 matrix that categorizes each grid point based on trend significance and direction in each of the datasets and compute the percentage of common grid points in each category.

```{r 3x3, message=FALSE, warning=FALSE}
# Function to create 3x3 contingency table
create_3x3_contingency_table <- function(beta_common, pval_common, num_common, alpha = 0.05) {
  table <- matrix(0, 3, 3)
  rownames(table) <- c("MODIS n.s.", "MODIS β >= 0", "MODIS β < 0")
  colnames(table) <- c("SNPP n.s.", "SNPP β >= 0", "SNPP β < 0")
  
  ns_indices_1 <- which(pval_common[[1]] >= alpha)
  ns_indices_2 <- which(pval_common[[2]] >= alpha)
  
  table[1, 1] <- sum(ns_indices_1 %in% ns_indices_2) / num_common * 100
  table[2, 1] <- sum(ns_indices_1 %in% which(pval_common[[2]] < alpha & beta_common[[2]] >= 0)) / num_common * 100
  table[3, 1] <- sum(ns_indices_1 %in% which(pval_common[[2]] < alpha & beta_common[[2]] < 0)) / num_common * 100
  
  table[1, 2] <- sum(ns_indices_2 %in% which(pval_common[[1]] < alpha & beta_common[[1]] >= 0)) / num_common * 100
  table[1, 3] <- sum(ns_indices_2 %in% which(pval_common[[1]] < alpha & beta_common[[1]] < 0)) / num_common * 100
  
  sig_indices_1 <- which(pval_common[[1]] < alpha)
  sig_indices_2 <- which(pval_common[[2]] < alpha)
  common_indices_sig1 <- sig_indices_1[sig_indices_1 %in% sig_indices_2]
  common_indices_sig2 <- sig_indices_2[sig_indices_2 %in% sig_indices_1]
  
  reshaped_data <- list(
    beta_common[[1]][common_indices_sig1],
    beta_common[[2]][common_indices_sig2]
  )
  
  table[2, 2] <- sum(reshaped_data[[1]] >= 0 & reshaped_data[[2]] >= 0) / num_common * 100
  table[2, 3] <- sum(reshaped_data[[1]] >= 0 & reshaped_data[[2]] < 0) / num_common * 100
  table[3, 2] <- sum(reshaped_data[[1]] < 0 & reshaped_data[[2]] >= 0) / num_common * 100
  table[3, 3] <- sum(reshaped_data[[1]] < 0 & reshaped_data[[2]] < 0) / num_common * 100
  
  return(table)
}

# Create 3x3 contingency table
contingency_table_3x3 <- create_3x3_contingency_table(common_data$beta_common, common_data$pval_common, common_data$num_common)

# Print 3x3 contingency table
print(contingency_table_3x3)

```

## Computing Cohen's Kappa

Cohen's Kappa quantifies the level of agreement between trend classifications in the MODIS and SNPP dataset.

Kappa values range from 0 to 1 where:

-   1.0 = Perfect Agreement

-   0.0 = No Agreement

The function `compute_kappa` calculates Cohen's Kappa for both the 2x2 and 3x3 contingency tables.

```{r kappa, message=FALSE, warning=FALSE}
# Function to compute Cohen's Kappa
compute_kappa <- function(contingency_table) {
  kappa_result <- cohen.kappa(contingency_table)$kappa
  return(kappa_result)
}

# Compute Cohen's Kappa for both contingency tables
kappa_2x2 <- compute_kappa(contingency_table_2x2)
kappa_3x3 <- compute_kappa(contingency_table_3x3)

```

## Visualizing Contingency Tables with Cohen's Kappa

The function `plot_contingency_table` generates a contingency table with Cohen's Kappa displayed below the plot. The contingency table represents the percentage of grid points in each classification category, while Cohen's Kappa quantifies the agreement between the two datasets.

```{r visualize, message=FALSE, warning=FALSE}
# Function to plot contingency table with Cohen's Kappa below
plot_contingency_table <- function(table, title, kappa_value) {
  df <- as.data.frame(as.table(table))
  
  # Create the contingency table plot
  p <- ggplot(df, aes(Var1, Var2, fill = Freq)) +
    geom_tile() +
    geom_text(aes(label = round(Freq, 1)), color = "white", size = 6) +
    scale_fill_gradient(low = "lightblue", high = "darkblue") +
    labs(title = title, fill = "Percentage", x = "", y = "") +
    theme_minimal()
  
  # Create a separate text grob for Kappa
  kappa_text <- grid::textGrob(
    label = paste("Cohen's Kappa:", round(kappa_value, 2)),
    gp = grid::gpar(fontsize = 12, fontface = "bold")
  )
  
  # Arrange the plot and Kappa text below
  full_plot <- gridExtra::grid.arrange(p, kappa_text, ncol = 1, heights = c(4, 0.5))
  
  return(full_plot)
}

```

## Plot 2x2 Contingency Table

```{r 2x2_plot, message=FALSE, warning=FALSE}
# Display and save 2x2 contingency table plot with Cohen's Kappa
p_2x2 <- plot_contingency_table(contingency_table_2x2, "2x2 Contingency Table for CCAL", kappa_2x2)

```

### Results for the 2x2 Contingency Table

-   42.8% of the grid points had a positive trend in both MODIS and SNPP.
-   39.8% of the grid points had a negative trend in both MODIS and SNPP.
-   15.3% of the grid points had a negative trend in SNPP but a positive trend in MODIS.
-   2.1% of the grid points had a negative trend in MODIS but a positive trend in SNPP.
-   Cohen's Kappa = 0.66 indicating substantial agreement between MODIS and SNPP datasets.

## Plot 3x3 Contingency Table

```{r 3x3_plot, message=FALSE, warning=FALSE}
# Display and save 3x3 contingency table plot with Cohen's Kappa
p_3x3 <- plot_contingency_table(contingency_table_3x3, "3x3 Contingency Table for CCAL", kappa_3x3)

```

### Results for the 3x3 Contingency Table

-   83.6% of the grid points were non-significant in both MODIS and SNPP.
-   5.1% of the grid points had a significant positive trend in both datasets.
-   3.3% of the grid points had a significant negative trend in both datasets.
-   4.1% of the grid points had a significant positive trend in SNPP but were non-significant in MODIS.
-   2.4% of the grid points had a significant negative trend in MODIS but were non-significant in SNPP.
-   Cohen's Kappa = 0.64 indicating substantial agreement between MODIS and SNPP datasets.
