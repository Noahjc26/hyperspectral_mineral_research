library(tidyverse)
library(terra)
library(raster)
library(foreach)
library(doParallel)
library(parallel)
library(snow)

calc_area_aviris <- function(image, wavelength1, wavelength2, cores) {
  # Register parallel backend
  registerDoParallel(cores = cores)
  
  # Define the loop function to calculate area
  calculate_area <- function(i) {
    cell_value <- as.data.frame(raster::extract(image, i))
    
    reflectance_values <- cell_value %>%
      pivot_longer(cols = starts_with("c"),
                   names_to = "wavelength",
                   values_to = "value") %>%
      mutate(wavelength = gsub("^channel_\\d+\\.\\.(\\d+\\.?\\d*)\\.nanometers\\.$", "\\1", wavelength),
             wavelength = as.numeric(wavelength))
    
    if (any(is.na(reflectance_values))) {
      area <- NA
    } else {
      index1 <- which.min(abs(reflectance_values$wavelength - wavelength1))
      index2 <- which.min(abs(reflectance_values$wavelength - wavelength2))
      
      subset_reflectance <- reflectance_values[index1:index2, ]
      
      value1 <- subset_reflectance$value[which.min(abs(subset_reflectance$wavelength - wavelength1))]
      value2 <- subset_reflectance$value[which.min(abs(subset_reflectance$wavelength - wavelength2))]
      
      slope <- (value2 - value1) / (wavelength2 - wavelength1)
      intercept <- value1 - slope * wavelength1
      
      area <- 0
      for (j in 1:(nrow(subset_reflectance) - 1)) {
        if (subset_reflectance$value[j] < slope * subset_reflectance$wavelength[j] + intercept &
            subset_reflectance$value[j + 1] < slope * subset_reflectance$wavelength[j + 1] + intercept) {
          segment_area <- ((subset_reflectance$wavelength[j + 1] - subset_reflectance$wavelength[j]) * 
                             (subset_reflectance$value[j] + subset_reflectance$value[j + 1]) / 2)
          area <- area + segment_area
        }
      }
    }
    print(area)
  }
  
  # Apply the loop function in parallel and collect the results
  results <- foreach(i = 1:ncell(image), .combine = 'rbind') %dopar% {
    calculate_area(i)
  }
  
  # Extract the cell numbers from the results
  cell_numbers <- unique(as.vector(results))
  
  # Create a copy of the raster with the same extent and CRS
  area_raster <- raster(image)
  values(area_raster) <- NA
  
  # Replace cell values with the calculated area values
  values(area_raster) <- results
  
  return(area_raster)
}


######################################
######################################
######################################
######################################
######################################
######################################
######################################
######################################

cell_graph <- function(image, wavelength1, wavelength2) {
  
  # Use a mouse click to select a pixel
  clicked_pixel <- locator(n = 1, type = "p")
  
  image <- stack(image)
  
  # Extract x and y coordinates from the clicked_pixel list
  clicked_x <- clicked_pixel$x
  clicked_y <- clicked_pixel$y
  
  # Convert the clicked pixel coordinates to cell indices
  clicked_cell <- cellFromXY(image, c(clicked_x, clicked_y))
  
  reflectance <- as.data.frame(terra::extract(image, clicked_cell))
  
  # Pivot the data frame longer
  reflectance <- reflectance %>%
    pivot_longer(cols = starts_with("c"),
                 names_to = "wavelength",
                 values_to = "value") %>%
    mutate(wavelength = gsub("^channel_\\d+\\.\\.(\\d+\\.?\\d*)\\.nanometers\\.$", "\\1", wavelength),
           wavelength = as.numeric(wavelength))
  
  # Find the indices of the closest wavelengths
  index1 <- which.min(abs(reflectance$wavelength - wavelength1))
  index2 <- which.min(abs(reflectance$wavelength - wavelength2))
  
  # Extract reflectance values within the desired wavelength range
  subset_reflectance <- subset(reflectance, wavelength >= wavelength1 & wavelength <= wavelength2)
  
  # Extract reflectance values at the selected wavelengths
  value1 <- subset_reflectance$value[which.min(abs(subset_reflectance$wavelength - wavelength1))]
  value2 <- subset_reflectance$value[which.min(abs(subset_reflectance$wavelength - wavelength2))]
  
  # Create a new data frame for the two points
  points_df <- data.frame(
    wavelength = c(wavelength1, wavelength2),
    value = c(value1, value2)
  )
  
  # Calculate the slope and intercept of the line between the two points
  slope <- (value2 - value1) / (wavelength2 - wavelength1)
  intercept <- value1 - slope * wavelength1
  
  # Initialize area variable
  area <- 0
  
  # Loop through each pair of consecutive points within the desired wavelength range
  for (i in 1:(nrow(subset_reflectance) - 1)) {
    # Check if both points are below the new line
    if (subset_reflectance$value[i] < slope * subset_reflectance$wavelength[i] + intercept &
        subset_reflectance$value[i + 1] < slope * subset_reflectance$wavelength[i + 1] + intercept) {
      # Calculate trapezoidal area for this segment and add it to the total area
      segment_area <- ((subset_reflectance$wavelength[i + 1] - subset_reflectance$wavelength[i]) * 
                         (subset_reflectance$value[i] + subset_reflectance$value[i + 1]) / 2)
      area <- area + segment_area
    }
  }
  
  # Print the area
  print(area)
  
  # Create a ggplot object with line, points, shaded area, and dashed lines
  p <- ggplot(reflectance, aes(x = wavelength, y = value)) +
    geom_line() +
    # geom_point(data = points_df, aes(x = wavelength, y = value), color = "red", size = 3) +
    # geom_ribbon(data = subset(subset_reflectance, value < (slope * wavelength + intercept)),
    #             aes(x = wavelength, ymin = value, ymax = slope * wavelength + intercept),
    #             fill = "gray", alpha = 0.5) +
    # geom_abline(intercept = intercept, slope = slope, linetype = "dashed", color = "blue") +
    theme_bw()
  
  return(p)
}

######################################
######################################
######################################
######################################
######################################
######################################
######################################
######################################

visualize_area_df_mineral <- function(data_frame, mineral_name, wavelength1, wavelength2) {
  # Filter data for the specified mineral
  mineral_data <- subset(data_frame, Mineral == mineral_name & 
                           Wavelength >= wavelength1 & Wavelength <= wavelength2)
  
  # Calculate slope and intercept for each sample
  mineral_data <- mineral_data %>%
    group_by(Sample) %>%
    mutate(
      slope = (last(Reflectance) - first(Reflectance)) / (last(Wavelength) - first(Wavelength)),
      intercept = first(Reflectance) - slope * first(Wavelength)
    )
  
  # Calculate area beneath the line for each sample
  mineral_data <- mineral_data %>%
    mutate(
      area = ifelse(Reflectance < slope * Wavelength + intercept,
                    abs((slope * Wavelength + intercept - Reflectance) * (Wavelength - lag(Wavelength, default = first(Wavelength)))) / 2,
                    0)
    )
  
  # Plot reflectance spectra and shaded areas for each sample
  p_list <- list()
  areas_df <- data.frame(Sample = character(), Area = numeric())
  for (sample_id in unique(mineral_data$Sample)) {
    sample_data <- filter(mineral_data, Sample == sample_id)
    p <- ggplot(sample_data, aes(x = Wavelength, y = Reflectance)) +
      geom_line() +
      geom_ribbon(data = subset(sample_data, Reflectance < (slope * Wavelength + intercept)),
                  aes(ymin = Reflectance, ymax = slope * Wavelength + intercept),
                  fill = "gray", alpha = 0.5) +
      geom_abline(aes(intercept = intercept, slope = slope), linetype = "dashed", color = "blue") +
      geom_point() +
      geom_text(data = subset(sample_data, Reflectance < (slope * Wavelength + intercept)),
                aes(x = max(Wavelength), y = max(Reflectance), label = paste("Area:", round(sum(area), 2))), 
                hjust = 1, vjust = 1) +
      labs(title = paste("Reflectance Spectra and Area Calculation for", mineral_name, "-", sample_id),
           x = "Wavelength",
           y = "Reflectance") +
      theme_bw()
    
    p_list[[sample_id]] <- p
    
    # Calculate total area for the sample and add to areas_df
    total_area <- sum(sample_data$area)
    areas_df <- rbind(areas_df, data.frame(Sample = sample_id, Area = total_area))
  }
  
  return(list(plots = p_list, areas_df = areas_df))
}


visualize_area_all_minerals <- function(data_frame, wavelength1, wavelength2) {
  mineral_names <- unique(data_frame$Mineral)
  all_areas <- data.frame(Mineral = character(), Sample = character(), Area = numeric())
  
  for (mineral_name in mineral_names) {
    mineral_data <- filter(data_frame, Mineral == mineral_name & Wavelength >= wavelength1 & Wavelength <= wavelength2)
    
    if (nrow(mineral_data) == 0) {
      cat("No data found for mineral:", mineral_name, "\n")
      next
    }
    
    mineral_data <- mineral_data %>%
      group_by(Sample) %>%
      mutate(
        slope = (last(Reflectance) - first(Reflectance)) / (last(Wavelength) - first(Wavelength)),
        intercept = first(Reflectance) - slope * first(Wavelength)
      ) %>%
      mutate(
        area = ifelse(Reflectance < slope * Wavelength + intercept,
                      abs((slope * Wavelength + intercept - Reflectance) * (Wavelength - lag(Wavelength, default = first(Wavelength)))) / 2,
                      0)
      )
    
    areas_df <- data.frame(
      Mineral = mineral_name,
      Sample = unique(mineral_data$Sample),
      Area = sapply(unique(mineral_data$Sample), function(sample_id) {
        sum(filter(mineral_data, Sample == sample_id)$area)
      })
    )
    
    all_areas <- rbind(all_areas, areas_df)
  }
  
  return(all_areas)
}
