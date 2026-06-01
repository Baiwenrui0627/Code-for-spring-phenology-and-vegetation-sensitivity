## Overview
This R code repository identifies flash drought events, extracts vegetation phenology, and quantifies the vegetation sensitivities to flash drought events, using NDVI/LAI/GOSIF/VOD as proxies. 

## 1. System Requirements

* **Operating System:** Tested on Windows 11.
* **Software:** R version 4.3.2 or higher.
* **Dependencies:** All required R packages for data processing are explicitly listed and loaded at the beginning of each corresponding R script.
* **Hardware Requirements:** No non-standard hardware is required.

## 2. Installation Guide

To set up the environment, please ensure you have R installed. You can install any missing packages using the standard R `install.packages()` command based on the libraries called at the top of each script.
* **Typical Install Time:** Less than 5 minutes on a normal desktop computer.

## 3. Demo 

* **Sample Data:** The full study area covers the Northern Hemisphere north of 23.5°N (665 latitudinal bands at 0.1° intervals). For demonstration and testing purposes, we provide a sample dataset representing the **400th latitude band**. 
* **Instructions to Run:** Download the repository, navigate to the folder, and run the scripts using the provided sample dataset.
* **Expected Output:** Successfully running the demo will output a final file containing the relative decline (%) of vegetation proxies under flash droughts.
* **Expected Run Time:** The demo dataset takes approximately 5-10 minutes to run on a normal desktop computer.

## 4. Instructions for Use

### 4.1 `flash_drought_identification.R`
* **Function:** Identifies flash drought timings within growing seasons. Only flash droughts with evapotranspiration <= 0 are retained.
* **Data Source:** 
  * Soil moisture data from GLEAM v4.2a (https://www.gleam.eu/). 
  * Soil moisture data derived from ERA5-Land (https://cds.climate.copernicus.eu/datasets/reanalysis-era5-land?tab=overview).
* **Temporal Coverage:** 2001–2021 (1980–2021)
* **Spatial Resolution:** 0.1° × 0.1°
* **Preprocessing:** Daily soil moisture data are averaged into pentad (5-day) means and transformed into percentiles.

### 4.2 `phenology_extraction.R`
* **Function:** Extracts spring and autumn phenology.
* **Data Source:** GLASS NDVI data (https://www.glass.hku.hk/download.html).
* **Temporal Coverage:** 2001–2021
* **Spatial Resolution:** 0.1° × 0.1°

### 4.3 `vegetation_sensitivity_to_flash_drought.R`
* **Function:** Quantifies vegetation sensitivity to flash drought, represented as the relative decline (%) of vegetation proxies under drought.
* **Data Source:** 
  * GLASS NDVI and LAI (https://www.glass.hku.hk/download.html). 
  * GOSIF (http://data.globalecology.unh.edu/data/GOSIF_v2/). 
  * LPDRV3 VOD data (http://files.ntsg.umt.edu/data/LPDR_v3/GeoTif).
* **Temporal Coverage:** 2001–2021
* **Spatial Resolution:** 0.1° × 0.1°
* **Preprocessing:** All data are averaged into pentad (5-day) means. The vegetation proxy data were detrended to obtain detrended data, which were then deseasonalized into z-scores.
