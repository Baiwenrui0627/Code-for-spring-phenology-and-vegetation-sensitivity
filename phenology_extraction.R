################ Extract spring and autumn phenology using 0.1° NDVI data and different thresholds ##########
library(data.table)

# HANTS curve fitting ####
range_matrix<-as.matrix(raster("final_range_res0.1.tif"))
time<-seq(1,361,8)
FFT_smooth<-function(data,K=8){
  N <- length(data)
  z<-fft(data) 
  mod_z <- Mod(z)
  threshold <- sort(mod_z, decreasing = TRUE)[K] 
  z[mod_z < threshold] <- 0 
  y <- Re(fft(z, inverse = TRUE)) / N 
  return(y)
}
HANTS<-function(ndvi_col,range){
  if(length(na.omit(range))==1){
    result.new<-numeric(7665)
    for(year in 1:21){
      ndvi.raw.data<-ndvi_col[c(((year-1)*46+1):(year*46))]
      x<-time
      xout<-c(1:365)
      chazhi <- approx(x, ndvi.raw.data, xout, method = "linear", rule = 2)
      yout<-chazhi$y
      result<-FFT_smooth(yout,K=8)
      result_x<-result[x]
      diff<-result_x-ndvi.raw.data
      while (any(diff >= 0.05)) {
        omit <- which(diff >= 0.05)
        x <- x[-omit]
        if (length(x) < 30) { break }
        ndvi.raw.data <- ndvi.raw.data[-omit]
        chazhi_new <- approx(x, ndvi.raw.data, xout, method = "linear", rule = 2)
        yout_new <- chazhi_new$y
        result <- FFT_smooth(yout_new, K = 8)
        result_x <- result[x]
        diff <- result_x - ndvi.raw.data
      }
      result.new[((year - 1) * 365 + 1):(year * 365)] <- result
    }
    return(round(result.new,digits=5))
  }else{return(c(rep(NA,7665)))}
}
for(i in c(1:665)){
  if(length(na.omit(range_matrix[i,]))>0){
    data.ndvi<-as.data.frame(fread(paste0("NDVI_0.1_raw_",sprintf("%03d",i),".csv")))
    data.range<-as.data.frame(t(range_matrix[i,]))
    final.data <- as.data.frame(mapply(HANTS, data.ndvi, data.range, SIMPLIFY = T))
    Sys.sleep(1)
    fwrite(final.data,paste0("HANTS_curve_res0.1_",sprintf("%03d",i),".csv"),row.names=F,na = "NA",quote = FALSE)
  }
  print(paste0("finish+",i))
  gc()
}

# Extract spring and autumn phenology using different thresholds ####
phenology_extract <- function(curve_col, threshold = 0.2) {
  if (sum(!is.na(curve_col)) < 7000) {  
    return(rep(NA, 65)) 
  } # Returns 2 thresholds + 21*SOS + 21*POS + 21*EOS = 65 NAs
  curve_mat <- matrix(curve_col, nrow = 365) 
  mean_curve <- rowMeans(curve_mat, na.rm = TRUE) 
  max_meandata <- max(mean_curve, na.rm = TRUE) 
  max_time <- round(median(which(mean_curve == max_meandata)),digits=0)
  
  left_min <- min(mean_curve[1:max_time], na.rm = TRUE)
  left_yuzhi <- threshold * (max_meandata - left_min) + left_min
  right_min <- min(mean_curve[max_time:365], na.rm = TRUE)
  right_yuzhi <- threshold * (max_meandata - right_min) + right_min
  
  result_SOS <- numeric(21)
  result_POS <- numeric(21)
  result_EOS <- numeric(21)
  for (year in 1:21) {
    yearly_curve <- curve_mat[, year] 
    POS <- round(median(which(yearly_curve == max(yearly_curve, na.rm = TRUE))),digits=0) 
    result_POS[year] <- POS
    # --- spring phenology ---
    left_part <- yearly_curve[1:POS]
    min_left_idx <- which.min(left_part)
    search_left <- left_part[min_left_idx:length(left_part)] 
    cross_up <- which(search_left[-length(search_left)] < left_yuzhi & search_left[-1] >= left_yuzhi) 
    if (length(cross_up) > 0) {
      result_SOS[year] <- min_left_idx - 1 + cross_up[1] 
    } else { result_SOS[year] <- NA }
    # --- autumn phenology ---
    right_part <- yearly_curve[POS:365]
    min_right_idx <- which.min(right_part)
    search_right <- right_part[1:min_right_idx]
    cross_down <- which(search_right[-length(search_right)] >= right_yuzhi & search_right[-1] < right_yuzhi) 
    if (length(cross_down) > 0) {
      result_EOS[year] <- (POS - 1) + cross_down[length(cross_down)]
    } else {result_EOS[year] <- NA}
  }
  return(c(left_yuzhi, right_yuzhi, result_SOS, result_POS, result_EOS))
}
thresholds <- c(0.2, 0.3, 0.4)
for (i in c(1:665)) {
  file <- paste0("HANTS_curve_res0.1_", sprintf("%03d", i), ".csv")
  if (file.exists(file)) {
    curve.df<-as.data.frame(fread(file))
    for (th in thresholds) {
      final.data <- as.data.frame(mapply(phenology_extract, curve.df, MoreArgs = list(threshold = th), SIMPLIFY = TRUE))
      Sys.sleep(1)
      out_file <- paste0("HANTS_phenology_th", th, "_", sprintf("%03d", i), ".csv")
      fwrite(final.data, out_file,row.names=F,na = "NA",quote = FALSE)
    }
  }
  gc()
  print(paste0("finish: ", i))
}

# Exclude spring and autumn phenology data exceeding 2.5 times the Median Absolute Deviation (MAD) ########
fun_mad <- function(data){
  mad.data <- mad(data, na.rm = TRUE)
  median.data <- median(data, na.rm = TRUE)
  if (is.na(mad.data) || mad.data == 0) return(data)
  data[abs(data - median.data) > 2.5 * mad.data] <- NA
  return(data)
}
thresholds <- c(0.2, 0.3, 0.4)
for (th in thresholds) {
  th_str <- paste0("th", th)
  for(i in 1:665){
    file <- paste0("HANTS_phenology_", th_str, "_", sprintf("%03d", i), ".csv")
    if(file.exists(file)){
      data.all <- fread(file)
      data.SOS <- data.all[3:23, ]
      data.EOS <- data.all[45:65, ]
      data.SOS_clean <- data.SOS[, lapply(.SD, fun_mad)]
      data.EOS_clean <- data.EOS[, lapply(.SD, fun_mad)]
      out_file <- paste0("phenology_select_", th_str, "_",sprintf("%03d", i), ".csv")
      fwrite(rbind(data.SOS_clean, data.EOS_clean), out_file, row.names = FALSE, na = "NA", quote = FALSE)
      Sys.sleep(0.5)
    }
  }
  print(paste0("========== Threshold ", th_str, " Finished! =========="))
  gc() 
}
