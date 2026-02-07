1.	Overview
   This R code identifies flash drought events, extracts vegetation phenology, and quantifies the vegetation sensitivities to flash drought events, using NDVI/LAI/GOSIF/VOD as proxies.
2.	System Requirements
   •Software: R V4.3.2
3.	Code instruction 
3.1 flash_drought_identification.R
  	•Data Source: The soil moisture data from GLEAMv4.1a are available at https://www.gleam.eu/. The soil moisture data derived from ERA5-Land are available through https://cds.climate.copernicus.eu/datasets/reanalysis-era5-land?tab=overview.
  	•Temporal Coverage: 2001–2021 •Spatial Resolution: 0.1° × 0.1°
  	•Preprocessing: Daily soil moisture data are averaged into pentad (5-day) means and transformed into percentiles.
  	•Output: Identified flash drought timings within growing seasons. Only flash droughts with evapotranspiration<0 were retained. 
3.2 phenology_extraction.R
  	•Data Source: GLASS NDVI data are obtained from https://www.glass.hku.hk/download.html.
  	•Temporal Coverage: 2001–2021
  	•Spatial Resolution: 0.1° × 0.1°.
  	•Output: Thresholds for extracting spring phenology and autumn phenology. Results of spring phenology and autumn phenology.
3.3 vegetation_sensitivity_to_flash_drought.R
  	•Data Source: Flash drought data was extracted from flash_drought_identification.R. Vegetation proxies include NDVI, LAI, GOSIF, and VOD. GLASS NDVI and LAI are obtained from https://www.glass.hku.hk/download.html. GOSIF is available at http://data.globalecology.unh.edu/data/GOSIF_v2/. LPDRV3 VOD data can be downloaded at http://files.ntsg.umt.edu/data/LPDR_v3/GeoTif. •Temporal Coverage: 2001–2021
  	•Spatial Resolution: 0.1° × 0.1°.
  	•Preprocessing: All data are averaged into pentad (5-day) means. The vegetation proxy data were detrended to obtain detrended data, which were then deseasonalized into z-scores.
  	•Output: Vegetation sensitivity to flash drought, represented as the relative decline (%) of vegetation proxies under drought.
