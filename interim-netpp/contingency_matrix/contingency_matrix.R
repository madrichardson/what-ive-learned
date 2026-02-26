# Contingency Matrix Tutorial 3/7/25

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

# Set ERDDAP URL
erddap_url = "http://localhost:8080/erddap"

# Get dataset info
dataInfo1 <- rerddap::info('trends_modis_monthly_9km', url=erddap_url)
print(dataInfo1)

dataInfo2 <- rerddap::info('trends_snpp_monthly_9km', url=erddap_url)
print(dataInfo2)

# Fetch data from ERDDAP
fetch_data <- function(dataset_id) {
  data <- griddap(
    dataset_id,
    url = erddap_url,
    latitude = c(lat_min, lat_max),
    longitude = c(lon_min, lon_max),
    fields = c("beta", "pval")
  )
  return(data$data)
}

dataset1 <- fetch_data("trends_modis_monthly_9km")
dataset2 <- fetch_data("trends_snpp_monthly_9km")

# Extract beta and pval
beta1 <- dataset1$beta
pval1 <- dataset1$pval
beta2 <- dataset2$beta
pval2 <- dataset2$pval

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

# Function to compute Cohen's Kappa
compute_kappa <- function(contingency_table) {
  kappa_result <- cohen.kappa(contingency_table)$kappa
  return(kappa_result)
}

# Compute Cohen's Kappa for both contingency tables
kappa_2x2 <- compute_kappa(contingency_table_2x2)
kappa_3x3 <- compute_kappa(contingency_table_3x3)

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

# Display and save 2x2 contingency table plot with Cohen's Kappa
p_2x2 <- plot_contingency_table(contingency_table_2x2, "2x2 Contingency Table for CCAL", kappa_2x2)

# Display and save 3x3 contingency table plot with Cohen's Kappa
p_3x3 <- plot_contingency_table(contingency_table_3x3, "3x3 Contingency Table for CCAL", kappa_3x3)
