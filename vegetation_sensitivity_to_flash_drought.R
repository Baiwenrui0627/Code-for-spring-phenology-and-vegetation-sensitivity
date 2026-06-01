
################ Calculate resilience-related parameters under FD using 0.1° NDVI data ################
library(data.table)

# impact_time_point,threshold_zscore,impact_min_zscore,actual_value,potential_value ####
resilience<-function(data_col,mean_sd_col,detrend_col,drought_col){
  # data_col: NDVI_norm_anomaly
  # mean_sd_col: multi-year mean and standard deviation of detrended NDVI
  # detrend_col: detrended NDVI
  if(length(na.omit(drought_col))>1 & length(na.omit(data_col))>0){
    drought_col<-na.omit(drought_col)
    onset_start<-drought_col[1:(length(drought_col)/6)]
    impact_time_point<-c()
    threshold_zscore<-c()
    impact_min_zscore<-c()
    actual_value<-c()
    potential_value<-c()
    
    for(k in 1:length(onset_start)){ 
      end_idx <- min(onset_start[k] + 17, length(data_col))
      norm_anomaly_search <- data_col[(onset_start[k]):end_idx]
      anomaly_min_time<-(which(norm_anomaly_search==min(norm_anomaly_search))[1])+onset_start[k]-1 #anomaly_min_time的实际时间
      threshold_anomaly<-data_col[onset_start[k]-1]
      data_decline<-threshold_anomaly-data_col[anomaly_min_time] 
      if(data_decline>0){ 
        impact_time_point[k]<-anomaly_min_time
        threshold_zscore[k]<-threshold_anomaly
        impact_min_zscore[k]<-data_col[anomaly_min_time]
        actual_value[k]<-detrend_col[anomaly_min_time]
        sd_value<-mean_sd_col[((anomaly_min_time-1)%%73)+1+73]
        mean_value<-mean_sd_col[((anomaly_min_time-1)%%73+1)]
        potential_value[k]<-threshold_anomaly*sd_value+mean_value
      }else{
        impact_time_point[k]<-(-999)
        threshold_zscore[k]<-(-999)
        impact_min_zscore[k]<-(-999)
        actual_value[k]<-(-999)
        potential_value[k]<-(-999)
      }
    }
    return(c(impact_time_point,threshold_zscore,impact_min_zscore,actual_value,potential_value,rep(NA,500-length(impact_time_point)*5)))
  }else{return(rep(NA,500))}
}

for(i in 1:665) {
  ndvi_file1 <- paste0("NDVI_norm_anomaly_0.1_", sprintf("%03d", i), ".csv")
  ndvi_file2 <- paste0("NDVI_pentad_means_sds_0.1_", sprintf("%03d", i), ".csv")
  ndvi_file3 <- paste0("detrend_NDVI_", sprintf("%03d", i), ".csv")
  if(file.exists(ndvi_file1) && file.exists(ndvi_file2) && file.exists(ndvi_file3)) {
    data.csv    <- as.data.frame(fread(ndvi_file1))
    mean.sd.csv <- as.data.frame(fread(ndvi_file2))
    detrend.csv <- as.data.frame(fread(ndvi_file3))
      FD_file <- paste0("FD_select_", sprintf("%03d", i), ".csv")
      if(file.exists(FD_file)) {
        flashdrought.csv <- as.data.frame(fread(FD_file))
        final.data <- as.data.frame(mapply(resilience, data.csv, mean.sd.csv, detrend.csv, flashdrought.csv, SIMPLIFY = TRUE))
        out_file <- paste0("resilience_", sprintf("%03d", i), ".csv")
        fwrite(final.data, out_file, row.names = FALSE, na = "NA", quote = FALSE)
      }
    Sys.sleep(0.5)
  }
  if (i %% 50 == 0) gc() 
}

############ Convert data to data.table format and calculate vegetation sensitivity ##############

cal_sos_zscore <- function(phe_data, all_col) {
  real_all_col <- all_col[!is.na(all_col)]
  len <- length(real_all_col) / 5
  if (len == 0 || sum(!is.na(phe_data)) == 0) {
    return(data.table(event_ID = integer(0),sos_zscore = numeric(0),impact_time_point = numeric(0),impact_min_zscore = numeric(0),actual_value = numeric(0),potential_value = numeric(0)))
  }
  impact_time_point <- as.numeric(real_all_col[1:len])
  impact_min_zscore <- as.numeric(real_all_col[(len*2+1):(len*3)])
  actual_value      <- as.numeric(real_all_col[(len*3+1):(len*4)])
  potential_value   <- as.numeric(real_all_col[(len*4+1):(len*5)])
  
  sos_zscore <- rep(NA, len)
  sos <- phe_data[1:21]
  sos_number <- which(is.na(sos))
  sos[is.na(sos)] <- round(median(sos, na.rm = TRUE), digits = 0)
  sd_sos <- sd(sos, na.rm = TRUE)
  if (!is.na(sd_sos) && sd_sos != 0) {
    zscore <- round((sos - mean(sos, na.rm = TRUE)) / sd_sos, digits = 3)
    zscore[sos_number] <- NA
    valid_idx <- which(!is.na(impact_time_point) & impact_time_point > 0) #剔除NA和-999
    if (length(valid_idx) > 0) {
      event_years <- ceiling(impact_time_point[valid_idx] / 73)
      sos_zscore[valid_idx] <- zscore[event_years]
    }
  }
  return(data.table(event_ID = 1:len,sos_zscore = sos_zscore,impact_time_point = impact_time_point,impact_min_zscore = impact_min_zscore,actual_value = actual_value,potential_value = potential_value))
}

final_out_file <- paste0("RF_combined.csv")
for(i in 1:665) {
  resilience_file <- paste0("resilience_", sprintf("%03d", i), ".csv")
  phe_file <- paste0("phenology_select_th0.2_", sprintf("%03d", i), ".csv")
  if(file.exists(resilience_file) && file.exists(phe_file)) {
    resilience.csv <- as.data.frame(fread(resilience_file))
    phe.csv <- as.data.frame(fread(phe_file))
    list_of_tables <- mapply(cal_sos_zscore, phe.csv, resilience.csv, SIMPLIFY = FALSE)
    long_dt <- rbindlist(list_of_tables, idcol = "col")
    long_dt[, col := as.integer(gsub("\\D", "", col))]
    long_dt[, row := i]
    long_dt <- long_dt[!is.na(impact_time_point) & impact_time_point != -999 & !is.na(sos_zscore)]
    setcolorder(long_dt, c("row", "col", "event_ID", "impact_time_point", "sos_zscore", "impact_min_zscore", "actual_value", "potential_value"))
    if (nrow(long_dt) > 0) {
      fwrite(long_dt, final_out_file, append = TRUE, row.names = FALSE, na = "NA", quote = FALSE)
    }
    Sys.sleep(0.5)
    print(paste0("finish: ", i))
    if (i %% 30 == 0) gc()
  }
}

# Calculate relative_decline(vegetation sensitivity) while excluding potential_value < 0.1 (NDVI) to avoid outliers

dt <- fread(final_out_file)
result_final <- dt[potential_value >= 0.1 & actual_value>0]
result_final[, relative_decline := (abs(actual_value - potential_value) / potential_value) * 100]
head(result_final)
fwrite(result_final, final_out_file, row.names = FALSE, na = "NA", quote = FALSE)
