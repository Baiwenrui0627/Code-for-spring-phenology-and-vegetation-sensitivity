################## flash_drought_identification ###################
library(data.table)

# Convert GLEAM 1980-2021 pentad top-1m soil moisture data to percentiles ####

for(i in c(1:665)){
  file <- paste0("pentad_GLEAM_sm_top1m_", sprintf("%03d", i), ".csv")
  # file is 1980-2021 pentad top-1m soil moisture data for a specific latitude
  
  if(file.exists(file)){
    dt <- fread(file, nrows = 3066)
    dt[, pentad := rep(1:73, length.out = .N)]
    cols_to_calc <- setdiff(names(dt), "pentad") 
    dt[, (cols_to_calc) := lapply(.SD, function(x) { 
      valid_n <- sum(!is.na(x))
      if (valid_n == 0) return(rep(NA, length(x))) 
      percentile <- (rank(x, na.last = "keep", ties.method = "max") / valid_n) * 100 
      return(round(percentile, digits = 4))
    }), by = pentad, .SDcols = cols_to_calc]
    dt[, pentad := NULL]
    out_file <- paste0("top1m_pentad_percentile_", sprintf("%03d", i), ".csv")
    fwrite(dt, out_file, row.names = FALSE, na = "NA", quote = FALSE)
    Sys.sleep(2)
  }
  print(paste0("finish: ", i))
  if(i %% 20 == 0) { gc() }
}

# Extract FD using 2001-2021 soil moisture percentiles from the 1980-2021 GLEAM dataset ####

find_FD_events <- function(soilm_data, high_threshold = 40, low_threshold = 20, min_rate = 5, number = 900) {
  if (sum(!is.na(soilm_data))<1000) {return(rep(NA, number))}
  n <- length(soilm_data)
  max_events <- 150 
  c_start_idx <- numeric(max_events)
  c_first_drop_idx <- numeric(max_events)
  c_drop_end_idx <- numeric(max_events)
  c_recovery_idx <- numeric(max_events)
  c_total_duration <- numeric(max_events)
  c_rate <- numeric(max_events)
  
  event_count <- 0
  state <- 0 
  i <- 2 
  while (i <= n) {
    if (is.na(soilm_data[i])) {
      state <- 0
      i <- i + 1
      next
    }
    if (state == 0) { 
      if (!is.na(soilm_data[i-1]) && soilm_data[i-1] >= high_threshold && soilm_data[i] < high_threshold) {
        temp_start_idx <- i
        temp_start_value <- soilm_data[i-1] #i-1
        state <- 1
      }
    } 
    else if (state == 1) { 
      if (soilm_data[i] >= high_threshold) {state <- 0} 
      else if (soilm_data[i] < low_threshold) {
        temp_first_drop_idx <- i
        temp_lowest_idx <- i
        temp_lowest_value <- soilm_data[i]
        state <- 2
      }
    }
    else if (state == 2) {
      if (soilm_data[i] < temp_lowest_value) {
        temp_lowest_idx <- i
        temp_lowest_value <- soilm_data[i]
      } else if (soilm_data[i] > temp_lowest_value) {
        rate <- (temp_start_value - temp_lowest_value) / (temp_lowest_idx - temp_start_idx + 1)
        if (rate > min_rate) {
          state <- 3 
          i <- i - 1 
        } else {state <- 0} 
      }
    }
    else if (state == 3) { 
      if (soilm_data[i] >= low_threshold) {
        temp_recovery_idx <- i
        duration <- temp_recovery_idx - temp_start_idx 
        if (duration >= 4 && duration <= 12) { 
          event_count <- event_count + 1
          c_start_idx[event_count] <- temp_start_idx 
          c_first_drop_idx[event_count] <- temp_first_drop_idx 
          c_drop_end_idx[event_count] <- temp_lowest_idx 
          c_recovery_idx[event_count] <- temp_recovery_idx 
          c_total_duration[event_count] <- duration 
          c_rate[event_count] <- rate 
        }
        state <- 0 
      }
    }
    i <- i + 1
  }
  if (event_count > 0) {
    res_vec <- c(c_start_idx[1:event_count],c_first_drop_idx[1:event_count],c_drop_end_idx[1:event_count],c_recovery_idx[1:event_count],c_total_duration[1:event_count],c_rate[1:event_count])
    res_vec <- c(res_vec, rep(NA, number - length(res_vec)))
    return(res_vec)
  } else {return(rep(NA, number))}
}

for(i in c(1:665)){
  file<-paste0("top1m_pentad_percentile_", sprintf("%03d", i), ".csv")
  # file is soil moisture percentiles data for a specific latitude
  if(file.exists(file)){
    data.percentile<-as.data.frame(fread(file))[1534:3066,] # select 2001-2021
    final.data <- as.data.frame(mapply(find_FD_events, data.percentile, SIMPLIFY = T))
    fwrite(final.data,paste0("flashdrought_",sprintf("%03d",i),".csv"),row.names = F, na = "NA", quote = F)
  }
  print(paste0("finish+",i))
}

######## Filter FD events within the growing season & ET <= 0 during the FD period ########

DOY<-rep(c(1:73),21)
FD_select <- function(drought_data, phe_data, ET_col) {
  sos <- phe_data[1:21]
  eos <- phe_data[22:42]
  if(length(na.omit(sos)) > 10 && length(na.omit(eos)) > 10 && length(na.omit(drought_data)) > 0 && length(na.omit(ET_col)) > 0) {
    drought_data <- na.omit(drought_data)
    len <- length(drought_data) / 6
    FD_start <- drought_data[1:len]
    FD_end <- (drought_data[(len * 3 + 1):(len * 4)])-1
    sos[is.na(sos)] <- round(median(sos[!is.na(sos)]), digits = 0)
    eos[is.na(eos)] <- round(median(eos[!is.na(eos)]), digits = 0)
    sos_pentad <- 73 * (0:20) + ceiling(sos / 5)
    eos_pentad <- 73 * (0:20) + ceiling(eos / 5)
    number1 <- which(sapply(FD_start, function(x) any(x >= sos_pentad & x <= eos_pentad)))
    number2 <- which(sapply(FD_end, function(x) any(x >= sos_pentad & x <= eos_pentad)))
    number <- intersect(number1, number2)
    if(length(number) == 0) {return(rep(NA, 600))}
    long_term_daily_means <- tapply(ET_col, DOY, mean, na.rm = TRUE)
    means <- long_term_daily_means[as.character(DOY)]
    anomaly <- round((ET_col - means), digits = 4)
    names(anomaly) <- NULL
    ET_anomaly <- numeric(max(number, 0)) 
    ET_anomaly[] <- NA 
    for(k in number) {ET_anomaly[k] <- round(mean(anomaly[FD_start[k]:FD_end[k]], na.rm = TRUE), digits = 4)}
    number3 <- which(ET_anomaly <= 0)
    select_number <- intersect(number, number3)
    if(length(select_number) > 0) {
      number_final <- c(select_number, select_number + len, select_number + len * 2, select_number + len * 3, select_number + len * 4, select_number + len * 5)
      result <- c(drought_data[number_final], rep(NA, 600 - length(number_final)))
    } else {result <- rep(NA, 600)}
    return(result)
  } else {return(rep(NA, 600))}
}

for(i in 1:665) {
  FD_file <- paste0("flashdrought_", sprintf("%03d", i), ".csv")
  phe_file <- paste0("phenology_select_th0.2_", sprintf("%03d", i), ".csv")
  ET_file <- paste0("pentad_E_", sprintf("%03d", i), ".csv")
  if(file.exists(FD_file) && file.exists(phe_file) && file.exists(ET_file)) {
    flashdrought <- as.data.frame(fread(FD_file))
    data.phe <- as.data.frame(fread(phe_file))
    data.ET <- as.data.frame(fread(ET_file))
    final.data <- as.data.frame(mapply(FD_select, flashdrought, data.phe, data.ET, SIMPLIFY = TRUE))
    out_file <- paste0("FD_select_", sprintf("%03d", i), ".csv")
    fwrite(final.data, out_file, row.names = FALSE, na = "NA", quote = FALSE)
    Sys.sleep(1)
  }
}
gc() 
