library(terra)
library(raster)
library(future)
library(future.apply)
library(doParallel)
library(zoom)
library(tidyverse)
library(janitor)
library(grid)
library(gridExtra)
library(circular)
source("../mineral_research/functions.R")





area_750_1250 <- (rast("../AVIRIS/2023/AV320231005t185153/area_750_1250.tif"))
area_1300_1600 <- (rast("../AVIRIS/2023/AV320231005t185153/area_1300_1600.tif"))
area_1600_2200 <- (rast("../AVIRIS/2023/AV320231005t185153/area_1600_2200.tif"))
area_1750_2100 <- (rast("../AVIRIS/2023/AV320231005t185153/area_1750_2100.tif"))
area_2050_2250 <- (rast("../AVIRIS/2023/AV320231005t185153/area_2050_2250.tif"))
area_2150_2250 <- (rast("../AVIRIS/2023/AV320231005t185153/area_2150_2250.tif"))
area_2275_2375 <- (rast("../AVIRIS/2023/AV320231005t185153/area_2275_2375.tif"))
area_2250_2450 <- (rast("../AVIRIS/2023/AV320231005t185153/area_2250_2450.tif"))
area_2200_2400 <- (rast("../AVIRIS/2023/AV320231005t185153/area_2200_2400.tif"))

#-----------------------------------------------------------

# #reading in 153 dataset
# raster_153 <- stack("../AVIRIS/2023/AV320231005t185153/L2A_OE/rectified.tif")
# 
# #removing layers with max value of -0.01
# raster_153_smaller <- (raster_153[[c(1:129,142:190,214:nlayers(raster_153))]])
# 
# #reading in 211 dataset
# raster_211 <- stack("../AVIRIS/2023/AV320231005t190211/L2A_OE/rectified.tif")
# 
# #removing layers with max value of -0.01
# raster_211_smaller <- (raster_211[[c(1:129,142:190,214:nlayers(raster_211))]])
# 
# raster_153_smaller <- rast(raster_153_smaller)
# raster_211_smaller <- rast(raster_211_smaller)
# 
# # Define your desired extent
desired_extent <- ext(402000, 404500, 4246000, 4248000)
# 
# ext(desired_extent)
# # Time the cropping operation
# timing_result <- system.time({
#   cropped_153 <- crop(raster_153_smaller, desired_extent)
# })
# 
# #reading in 2023 AVIRIS
cropped_153 <- (rast("../AVIRIS/2023/AV320231005t185153/L2A_OE/cropped_1.tif"))
# # 
cropped_153 <- crop(cropped_153,desired_extent)
# # 
# # #performing PCA
#PCA <- prcomp(cropped_153)


# # Extract scores and loadings from PCA_result
#scores <- predict(cropped_153,PCA, index=1:5)
# 
# #cropping again
# scores <- crop(scores,desired_extent)
scores <- rast("../AVIRIS/2023/AV320231005t185153/scores.tif")

#reading in 2019 AVIRIS
aviris_2019 <- rast("../AVIRIS/2019/rectified.tif")

#cropping 2019 AVIRIS
aviris_2019 <- crop(aviris_2019,desired_extent)

#performing PCA
PCA2 <- prcomp(aviris_2019)
plot(PCA2)
#turning into rasater of PCA
scores2 <- predict(aviris_2019,PCA2,index=1:5)

# Get the CRS of scores2
crs_scores2 <- crs(scores2)

#reading in geologic map
malmsten_cropped <- rast("../Quads/Malmsten_peak_quad/cropped_malmsten.tif")

#reading in shape file of geologic units
shape <- vect("../Quads/Malmsten_peak_quad/units_malmsten.shp")

#projecting shape units to the same reference system
shape <- project(shape,crs_scores2)

#cropping shapefile
shape_cropped <- crop(shape,desired_extent)

# Specify the directory containing the .tif files
directory <- "../landsat/LC09_L2SP_038033_20231019_20231020_02_T1/LC09_L2SP_038033_20231019_20231020_02_T1/"

# Get a list of all .tif files in the directory
tif_files <- list.files(directory, pattern = "\\.TIF$", full.names = TRUE)

tif_files <- rast(tif_files)

tif_files <- crop(tif_files,desired_extent)

PCA3 <- prcomp(tif_files)

scores3 <- predict(tif_files,PCA3)

# Define custom colors for each unit
unit_colors <- c(Qa = "cornsilk1", Ql = "burlywood1", Qp = "lightgreen", Tdk = "hotpink", Tlv = "red4", Tn = "red", Tda = "purple",Tla = "darkcyan")

# Plot the RGB image
plotRGB(scores,1,2,3, axes=TRUE, stretch="lin", main="2023 AVIRIS (3m resolution) PCA RGB")
plotRGB(scores2,1,2,3, axes=TRUE, stretch="hist",main="2019 AVIRIS (14m resolution) PCA RGB")
plot(malmsten_cropped, main = "Geologic map")
plot(shape_cropped, col = unit_colors[shape_cropped$unit],add=T)
legend("top", legend = names(unit_colors), fill = unit_colors, cex = 0.5,inset = c(0, 0.3))

#reading in spectral mineral signatures
list <-
  list.files(path ="../usgs_spectral_library/usgs_splib07 (1)/ASCIIdata/ASCIIdata_splib07b_cvAVIRISc2014/ChapterM_Minerals/",
             full.names = TRUE)

#list applying as csv and into a data frame
minerals <- list %>%
  lapply(read.csv) %>%
  as.data.frame() %>%
  clean_names() %>%
  mutate(Bands = row_number())

#removing everything before first underscore
colnames(minerals) = sub(".*?_", "", colnames(minerals))

#removing everything before first underscore
colnames(minerals) = sub(".*?_", "", colnames(minerals))

#removing everything before first underscore
colnames(minerals) = sub(".*?_", "", colnames(minerals))

#keeping only after first two underscores
colnames(minerals) <- str_extract(colnames(minerals), "[^_]*_[^_]*[^_]*")

#renaming blank column to band
colnames(minerals)[1277] = "Band"

# Generate the vector
vector <- seq(from = 400, by = 9.375, length.out = 224)

#adding wavelength as a column
minerals <- minerals %>% 
  mutate(wavelength = vector)

#------------------------------------
# Identify columns with negative values
neg_cols <- sapply(minerals, function(x) any(x < 0))

# Subset the data, excluding columns with negative values
data_subset <- minerals[, !neg_cols]

#alunite spectral signature
data_subset %>% 
  ggplot(aes(x=wavelength,y=data_subset$'222_alunite')) +
  geom_line(color="red") +
  ylim(0,1) +
  labs(y="Reflectance",
       x="Wavelength") +
  theme_bw()

# Define custom breaks for the x-axis
custom_breaks <- seq(400, 2500, by = 100) 

#chlorite spectral signature
data_subset %>% 
  ggplot(aes(x=wavelength)) +
  geom_line(color="blue",aes(y=data_subset$'805_chlorite'),size=1) +
  geom_line(color="blue",aes(y=data_subset$'1140_epidote'),size=1)+
  geom_line(color="green",aes(y=data_subset$'2206_muscovite'),size=1)+
  geom_line(color="green",aes(y=data_subset$'1768_kaolinite'),size=1)+
  geom_line(color="red",aes(y=data_subset$'1412_halloysite'),size=1)+
  geom_line(color="red",aes(y=data_subset$'1031_dickite'),size=1)+
  geom_line(color="orange",aes(y=data_subset$'2788_quartz'),size=1)+
  geom_line(color="orange",aes(y=data_subset$'742_chalcedony'),size=1)+
  geom_line(color="purple",aes(y=data_subset$'2736_pyrophyllite'),size=1)+
  geom_line(color="purple",aes(y=data_subset$'1020_diaspore'),size=1)+
  geom_line(color="yellow",aes(y=data_subset$'581_biotite'),size=1)+
  geom_line(color="black",aes(y=data_subset$'90_albite'),size=1)+
  geom_line(color="black",aes(y=data_subset$'2341_nepheline'),size=1)+
  ylim(0,1) +
  xlim(420,2500) +
  labs(y="Reflectance",
       x="Wavelength (nm)") +
  scale_x_continuous(breaks = custom_breaks) +
  theme_bw()

data_subset %>%
  select(matches("scolecite|wavelength")) %>% 
  pivot_longer(cols = ends_with("scolecite")) %>% 
  ggplot(aes(x=wavelength,y=value,col = name)) +
  geom_line() +
  labs(x="Wavelength (nm)",
       y="Reflectance") +
  theme_bw()

data_subset %>%
  select(matches("albite|wavelength")) %>% 
  pivot_longer(cols = ends_with("albite")) %>%
  filter(name %in% unique(head(name, 10))) %>%  # Keep data for the first 8 unique chlorites
  ggplot(aes(x = wavelength, y = value, col = name)) +
  geom_line() +
  labs(x = "Wavelength (nm)",
       y = "Reflectance",
       col = NULL) +
  theme_bw()

#kaolinite 1300-1600 and 2100-2250

data_subset %>%
  select(matches("chlorite|wavelength")) %>% 
  pivot_longer(cols = ends_with("_chlorite")) %>% 
  ggplot(aes(x=wavelength,y=value,col = name)) +
  geom_line() +
  labs(x="Wavelength (nm)",
       y="Reflectance") +
  theme_bw()

#absorption zones
# black: 1900 to 1950
# blue: 2150 to 2400
# red: 2100 to 2250
# Orange: 1900 to 2050
# yellow: 2300 to 2400
# green:2100 to 2250
# purple:1375 to 1450

# Assuming your data frame is called data_subset
data_subset_long <- pivot_longer(data_subset, 
                                 cols = -c(Band, wavelength),
                                 names_to = "Sample_Mineral", 
                                 values_to = "Reflectance")

# Extract sample and mineral names
data_subset_long <- separate(data_subset_long, Sample_Mineral, into = c("Sample", "Mineral"), sep = "_")

# Convert wavelength column to numeric
data_subset_long$wavelength <- as.numeric(data_subset_long$wavelength)

# Rename the 'wavelength' column to 'Wavelength'
data_subset_long <- rename(data_subset_long, Wavelength = wavelength)

visualize_area_df_mineral(data_subset_long, mineral_name = "halloysite",wavelength1 = 1300,wavelength2 = 1600)

# Call the function with your data
results <- visualize_area_all_minerals(data_subset_long, wavelength1 = 750, wavelength2 = 1250)

results <- arrange(results,desc(Area))

# Assuming your data frame is named 'df'
mean_area_by_mineral <- results %>%
  group_by(Mineral) %>%
  summarize(mean_area = mean(Area),
            sd_area = sd(Area))

# Print the result
print(arrange(mean_area_by_mineral,desc(mean_area)),n =210)
slope
#-----------------------------------------------------------

dem <- rast("../DEMs/dry_creek/output_USGS1m.tif")

# # Define your desired extent
desired_extent <- ext(402000, 404500, 4246000, 4248000)

dem <- crop(dem,desired_extent)

slope <- terrain(dem, v = "slope", unit = "degrees")
aspect <- terrain(dem, v = "aspect", unit = "degrees")
plot(shade(slope*pi/180, aspect*pi/180, angle = 10, direction = 270), col = grey(0:100/100), legend = F)

plot(area_1300_1600, breaks = c(17-3.26,17,17+3.26), col = c("red","red"), add = T, legend = F) #analcime, halloysite, saponite, heulandite
plot(area_2150_2250, breaks = c(5.47-1.06,5.47,5.47+1.06), col = c("orange","orange"), add = T, legend = F) #dickite nacrite, paragonite, kaolinite
plot(area_750_1250, breaks = c(37-8.34,37,37+8.34), col = c("green","green"), add = T, legend = F) #bronzite, jarosite, hypersthene, rhodonite
plot(area_1750_2100, breaks = c(22.3-1.45,22.3,22.3+1.45), col = c("blue","blue"), add = T, legend = F) #analcime, saponite, laumontite, clinoptilolite
plot(area_2050_2250, breaks = c(8.64-2.68,8.64,8.64+2.68), col = c("purple","purple"), add = T, legend = F) #zunyite, dickite, alunite, pyrophyllite


#-----------------------------------------------------------

# Creating a directional filter for emphasizing horizontal lineaments
kernel_horizontal <- matrix(c(-1, -1, -1, 2, 2, 2, -1, -1, -1), nrow=3, byrow=TRUE)
response_horizontal <- focal(dem, w=kernel_horizontal, fun=sum)

# Define a vertical Sobel filter
kernel_vertical <- matrix(c(-1, 0, 1,
                            -2, 0, 2,
                            -1, 0, 1), byrow=TRUE, nrow=3)

# Apply the filter
vertical_response <- focal(dem, w=kernel_vertical, fun=sum, na.rm=TRUE)

# Calculate a threshold, for example using a high percentile
threshold_vertical <- quantile(values(vertical_response, na.rm=TRUE), probs=0.95)

# Create a binary raster where high values are considered potential lineaments
vertical_lineaments <- vertical_response > threshold_vertical

# Thresholding the response to get significant features
threshold_horizontal <- quantile(na.omit(values(response_horizontal)), probs=0.95)
lineaments_horizontal <- response_horizontal > threshold_horizontal

#-----------------------------------------------------------

# Calculate gradients (Sobel filters)
gx <- focal(dem, matrix(c(-1, 0, 1, -2, 0, 2, -1, 0, 1), nrow = 3), fun = sum)
gy <- focal(dem, matrix(c(1, 2, 1, 0, 0, 0, -1, -2, -1), nrow = 3), fun = sum)

edges = sqrt(gx^2 + gy^2)
# Apply threshold to gradients to identify significant features
threshold <- 20 # Adjust threshold as needed
gx[edges <= threshold] <- NA
gy[edges <= threshold] <- NA
plot(shade(slope*pi/180, aspect*pi/180, angle = 10, direction = 270), col = grey(0:100/100), legend = F)
plot(edges > threshold, col = c("transparent","red"), add = T)

# Convert edge data to slope angles (assuming they represent the orientations)
lineament_angles <- atan2(na.omit(values(gy)), na.omit(values(gx))) * 180 / pi

# Adjust angles to be in the range of 0 to 360 degrees
lineament_angles <- (lineament_angles + 360) %% 360

# Convert angles to radians for circular statistics
lineament_angles_rad <- lineament_angles * pi / 180

# Restrict angles to the range 0-180 and 180-360
lineament_angles[lineament_angles > 180] <- lineament_angles[lineament_angles > 180] - 180

# Create circular data object
angles_circ <- circular(lineament_angles, units = "degrees")

# Define the breaks for grouping the angles
breaks <- seq(0, 360, by = 360 / 24)

# Use the cut function to group the angles into bins
angle_bins <- cut(lineament_angles, breaks = breaks, include.lowest = TRUE)

# Calculate the frequency of each angle bin
angle_counts <- table(angle_bins)

# Calculate the maximum count to use as the reference for scaling
max_count <- max(angle_counts)

# Scale the lengths of the bars based on the frequency of each angle bin
scaled_lengths <- angle_counts / max_count

# Plot the rose diagram with adjusted lengths
rose.diag(angles_circ, bins = 24, prop = 2, col = "red")
rose.diag(angles_circ, bins = 24, col = "red",prop = 2, zero =pi, add = T, axes = F)


# Plot the lineaments
plot(shade(slope*pi/180, aspect*pi/180, angle = 10, direction = 270), col = grey(0:100/100), legend = F)
plot(vertical_lineaments, col = c("transparent","red"), legend = F, add = T)
plot(lineaments_horizontal, col = c("transparent","lightgreen"),add = T, legend = F)


#-----------------------------------------------------------
# Define colors and associated minerals
colors <- c("red", "orange", "green", "blue", "purple")
minerals <- list(
  c("Analcime", "Halloysite", "Saponite", "Heulandite"),
  c("Dickite", "Nacrite", "Paragonite", "Kaolinite"),
  c("Bronzite", "Jarosite", "Hypersthene", "Rhodonite"),
  c("Analcime", "Saponite", "Laumontite", "Clinoptilolite"),
  c("Zunyite", "Dickite", "Alunite", "Pyrophyllite")
)

plot(shade(slope*pi/180, aspect*pi/180, angle = 10, direction = 270), col = grey(0:100/100), legend = F)

# Add other plots with 'add = TRUE' to overlay them
plot(area_1300_1600, breaks = c(17-3.26,17,17+3.26), col = c("red","red"), legend = FALSE, add = T)
plot(area_2150_2250, breaks = c(5.47-1.06,5.47,5.47+1.06), col = c("orange","orange"), add = TRUE, legend = FALSE)
plot(area_750_1250, breaks = c(37-8.34,37,37+8.34), col = c("green","green"), add = TRUE, legend = FALSE)
plot(area_1750_2100, breaks = c(22.3-1.45,22.3,22.3+1.45), col = c("blue","blue"), add = TRUE, legend = FALSE)
plot(area_2050_2250, breaks = c(8.64-2.68,8.64,8.64+2.68), col = c("purple","purple"), add = TRUE, legend = FALSE)

# Plot colored rectangles with mineral labels
legend_x <- 402000
legend_y <- 4247800
legend_rect_size <- 80
legend_text_x <- legend_x + legend_rect_size + 50
for (i in seq_along(minerals)) {
  legend_labels <- minerals[[i]]
  legend_color <- colors[i]
  legend_rect_y <- legend_y - i * legend_rect_size
  rect(legend_x, legend_rect_y, legend_x + legend_rect_size, legend_rect_y - legend_rect_size, col = legend_color)
  text(legend_text_x, legend_rect_y - legend_rect_size / 2, paste(legend_labels, collapse = ", "), adj = 0)
}



##### reading in sample location sites
locat <- read.csv("./Locations_minerals.csv")

locat_v <- vect(locat, geom = c("Longitude", "Latitude"), crs = "EPSG:4326")

plot(shade(slope*pi/180, aspect*pi/180, angle = 10, direction = 270), col = grey(0:100/100), legend = F, mar = 3)


plot(project(locat_v,slope), "Title", add = T, col = c("red","orange","green","blue","purple","hotpink"))

cell_graph(image = cropped_153,wavelength1 = 1300,wavelength2 = 1600)
