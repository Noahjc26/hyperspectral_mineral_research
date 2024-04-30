library(tidyverse)
library(terra)


files <- list.files("./spectrometer_data/2024_Apr_22/", pattern = ".sed", full.names = T)

for (i in 1:length(files)){
  df = read.table(files[i], skip = 27, sep= "")
  colnames(df) = c("Wavelength", "Reflectance")
  plot(Reflectance~Wavelength, df, type="l", lwd = 1.5)
}

