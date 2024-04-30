library(terra)


dem <- rast("../DEMs/dry_creek/output_USGS1m.tif")

# # Define your desired extent
desired_extent <- ext(402000, 404500, 4246000, 4248000)

dem <- crop(dem,desired_extent)

slope <- terrain(dem, v = "slope", unit = "degrees")
aspect <- terrain(dem, v = "aspect", unit = "degrees")
TPI <- terrain(dem, "TPI")
TRI <- terrain(dem, "TRI")
TRIriley <- terrain(dem, "TRIriley")
TRIrmsd <- terrain(dem, "TRIrmsd")
roughness <- terrain(dem, "roughness")
flowdir <- terrain(dem, "flowdir")

plot(slope)
plot(aspect)
plot(TPI)
plot(TRI)
plot(TRIriley)
plot(TRIrmsd)
plot(roughness)
plot(flowdir)

?terrain
plot(shade(slope*pi/180, aspect*pi/180, angle = 10, direction = 270), col = grey(0:100/100))



#Smoothing dem

smoothed_dem <- focal(dem, w = matrix(1, nrow = 3, ncol = 3)/9, fun = mean)

