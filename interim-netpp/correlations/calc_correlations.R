#### Correlations Tutorial 2/27/25

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
                      "jsonlite", "rerddapXtracto")

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
dataInfo <- rerddap::info('correlations_modis_snpp_monthly_9km', url=erddap_url)  
print(dataInfo$alldata$time)
# Set extraction parameters
parameter <- 'corr'
trends_time <- c("2018-01-16T12:00:00Z")

# ✅ Extract **full bounding box data** (like Python's unmasked bbox plot)
bbox_data <- rxtracto_3D(
  dataInfo,
  parameter = parameter,
  xcoord = c(lon_min, lon_max),
  ycoord = c(lat_min, lat_max),
  tcoord = c("2018-01-16 12:00:00", "2018-01-16 12:00:00")
)

# ✅ Extract **province-masked data** (like Python's masked province plot)
prov_data <- rxtractogon(
  dataInfo,
  parameter = parameter,
  xcoord = longitude,
  ycoord = latitude,
  tcoord = trends_time
)

# ✅ Plot both datasets side by side, like Python
par(mfrow = c(1, 2))  # Set side-by-side layout

# Set minimum
prov_data$corr[prov_data$corr < -1] <- -1
bbox_data$corr[bbox_data$corr < -1] <- -1
# Set maximum
prov_data$corr[prov_data$corr > 1] <- 1
bbox_data$corr[bbox_data$corr > 1] <- 1

custom_colors <- colorRampPalette(c("red", "white", "blue"))(100)

plotBBox(bbox_data,
         plotColor = custom_colors,
         maxpixels = 1000000,
         name = 'Correlation coefficients'
)

plotBBox(prov_data,
         plotColor = custom_colors,
         maxpixels = 1000000,
         name = 'Correlation coefficients'
)



































