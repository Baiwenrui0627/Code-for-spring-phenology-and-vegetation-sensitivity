########## calculation_of_vegetation_sensitivity_to_flash_drought ###############

#### transfer NDVI data to 5day resolution ####
years<-2001:2021
days_of_year <- sprintf("%03d",seq(1,361,8))
dates <- c()
for (year in years) {
  for (day in days_of_year) {
    date <- as.Date(paste(year, day, sep = "-"), format = "%Y-%j")
    dates <- append(dates, date)
  }
}
dates <- as.numeric(dates)
full_dates <- seq(from = as.Date(paste(min(years), "-01-01", sep = "")),
                  to = as.Date(paste(max(years), "-12-31", sep = "")),
                  by = "day")
leap_years <- seq(2004,max(years), by = 4)
last_day_of_leap_years <- as.Date(paste(leap_years, "-12-31", sep = ""))
full_dates <- full_dates[!full_dates %in% last_day_of_leap_years] 
full_dates <- as.numeric(full_dates)
group <- ceiling(c(1:7665)/ 5)

five_days_mean<-function(de_data){
  if(length(na.omit(de_data))==966){
    daily <- approx(x = dates, y = de_data, xout = full_dates, method = "linear", rule = 2)$y
    fiveday_avg<-tapply(daily, group, mean, na.rm = TRUE)
    return(fiveday_avg)
  }else{return(rep(NA,1533))}
}

for(i in c(1:665)){
  file<-paste0("path_to_NDVI_8_day_row_data.csv")
  if(file.exists(file)){
    NDVI.data <- read.csv(file)
    final.data <- as.data.frame(mapply(five_days_mean, NDVI.data, SIMPLIFY = T))    
    fwrite(final.data,"path_to_NDVI_5_day_row_data.csv",row.names=F)
  }
  print(paste0("finish+",i))
}


#### detrend and deseason NDVI ####
years <- 2001:2021
days_of_year <- sprintf("%03d",seq(3,365,5))
dates <- c()
for (year in years) {
  for (day in days_of_year) {
    date <- as.Date(paste(year, day, sep = "-"), format = "%Y-%j")
    dates <- append(dates, date)
  }
}
DOY <- as.numeric(format(dates, "%j"))
X <- c(1:1533)
deseason_detrend_linear <- function(data_col) {
  if(length(na.omit(data_col))==1533){
    model <- lm(data_col ~ X) 
    linear_trend<-(model$coefficients[2])*(X[seq(37,1533,73)])
    final_linear_trend<-rep(linear_trend,each=73)
    detrended_data <- data_col - final_linear_trend
    return(round(detrended_data,digits = 6))
  }else{return(rep(NA,1533))}
}

for(i in c(1:664)){
  file<-paste0("path_to_NDVI_5_day_row_data.csv")
  if(file.exists(file)){
    NDVI.data <- read.csv(file)
    final.data <- as.data.frame(mapply(deseason_detrend_linear, NDVI.data, SIMPLIFY = T))    
    fwrite(final.data,"path_to_deseason_detrend_NDVI_5_day_row_data.csv",row.names=F)
  }
  print(paste0("finish+",i))
}


#### transfer to normalized NDVI anomaly ####
DOY<-rep(c(1:73),21)
cal_norm_anomaly<-function(data_col){
  if(length(na.omit(data_col))>0){
    long_term_daily_means <- tapply(data_col, DOY, mean, na.rm = TRUE) 
    long_term_daily_sd <- tapply(data_col, DOY, sd, na.rm = TRUE)
    means<-long_term_daily_means[as.character(DOY)]
    sds<-long_term_daily_sd[as.character(DOY)]
    z_score_data <- round((data_col - means)/sds,digits=6)
    names(z_score_data)<-NULL
    names(long_term_daily_means)<-NULL
    names(long_term_daily_sd)<-NULL
    return(c(long_term_daily_means,long_term_daily_sd,z_score_data))
  }else{return(rep(NA,(1533+73*2)))}
}

for(i in c(1:664)){
  file<-paste0("path_to_deseason_detrend_NDVI_5_day_row_data.csv")
  if(file.exists(file)){
    data.csv <- read.csv(file)
    final.data <- as.data.frame(mapply(cal_norm_anomaly, data.csv, SIMPLIFY = T))
    fwrite(final.data[c(1:146),],"path_to_NDVI_pentad_means_sds_row_data.csv",row.names=F)
    fwrite(final.data[c(147:1679),],"path_to_NDVI_pentad_norm_anomaly_row_data.csv",row.names=F)
  }
  print(paste0("finish+",i))
}


##### calculate (1)total_actual_loss,(2)Impact,(3)Recovery Rate,(4)impact_time_point,(5)recovery_time_point,(6)threshold_norm_anomaly ######
# These are process variables that may be used in subsequent calculations of vegetation sensitivity to flash drought in our research.
# The LAI in the following function refer to vegetation indicators (NDVI, GOSIF, VOD, LAI).

resilience<-function(data_col,mean_sd_col,detrend_col,drought_col,phe_data){
  # data_col:NDVI_norm_anomaly
  # mean_sd_col:NDVI means and sds
  # detrend_col:actual NDVI
  if(length(na.omit(drought_col))>1 & length(na.omit(data_col))>0 & length(na.omit(phe_data))>0){
    sos<-phe_data[1:21]
    eos<-phe_data[22:42]
    sos[is.na(sos)]<-round(median(sos[!is.na(sos)]),digits=0)
    eos[is.na(eos)]<-round(median(eos[!is.na(eos)]),digits=0)
    sos_pentad<-73*c(0:20)+ceiling(sos/5)
    eos_pentad<-73*c(0:20)+ceiling(eos/5)
    drought_col<-na.omit(drought_col)
    onset_start<-drought_col[1:(length(drought_col)/5)]
    
    total_actual_LAI_loss<-c()
    Impact<-c()
    Rate<-c()
    impact_time_point<-c()
    recovery_time_point<-c()
    threshold_norm_anomaly<-c()
    
    for(k in 1:length(onset_start)){ 
      year<-ceiling(onset_start[k]/73)
      if((onset_start[k]+17)<=eos_pentad[year]){
        norm_anomaly_search<-data_col[(onset_start[k]):(onset_start[k]+17)]
        anomaly_min_time<-(which(norm_anomaly_search==min(norm_anomaly_search))[1])+onset_start[k]-1 
        threshold_anomaly<-data_col[onset_start[k]-1]
        data_decline<-threshold_anomaly-data_col[anomaly_min_time] 
        if(data_decline>0){ 
          impact_time_point[k]<-anomaly_min_time
          Impact[k]<-data_decline
          threshold_norm_anomaly[k]<-threshold_anomaly
          r_threshold<-0.95*data_decline+data_col[anomaly_min_time] 
          data_select<-data_col[anomaly_min_time:(eos_pentad[year])]
          indices <- (which(data_select > r_threshold & 
                              c(FALSE, data_select[-length(data_select)] > r_threshold) & 
                              c(FALSE, FALSE, data_select[c(-length(data_select),-(length(data_select)-1))] > r_threshold))-2)+anomaly_min_time-1
          if(length(indices)>0){
            recovery_time_point[k]<-indices[1]
            Rate[k]<-(data_decline/(indices[1]-anomaly_min_time))
            actual_LAI_threshold <- (threshold_anomaly*mean_sd_col[(c((onset_start[k]-1):indices[1])%%73)+73])+mean_sd_col[c((onset_start[k]-1):indices[1])%%73]
            actual_LAI <- detrend_col[c((onset_start[k]-1):indices[1])]
            total_actual_LAI_loss[k]<-sum((actual_LAI_threshold-actual_LAI),na.rm=T)
          }else{
            actual_LAI_threshold <- (threshold_anomaly*mean_sd_col[(c((onset_start[k]-1):eos_pentad[year])%%73)+73])+mean_sd_col[c((onset_start[k]-1):eos_pentad[year])%%73]
            actual_LAI <- detrend_col[c((onset_start[k]-1):eos_pentad[year])]
            total_actual_LAI_loss[k]<-sum((actual_LAI_threshold-actual_LAI),na.rm=T)
            Rate[k]<-"LAInotrecovery_inthisGS"
            recovery_time_point[k]<-"LAInotrecovery_inthisGS"
          }
        }else{
          total_actual_LAI_loss[k]<-"LAInotdecline"
          Impact[k]<-"LAInotdecline"
          Rate[k]<-"LAInotdecline"
          impact_time_point[k]<-"LAInotdecline"
          recovery_time_point[k]<-"LAInotdecline"
          threshold_norm_anomaly[k]<-"LAInotdecline"
        }
      }else{
        norm_anomaly_search<-data_col[(onset_start[k]):eos_pentad[year]]
        anomaly_min_time<-(which(norm_anomaly_search==min(norm_anomaly_search))[1])+onset_start[k]-1 
        threshold_anomaly<-data_col[onset_start[k]-1]
        data_decline<-threshold_anomaly-data_col[anomaly_min_time] 
        if(data_decline>0){
          if(anomaly_min_time<=(eos_pentad[year]-2)){
            impact_time_point[k]<-anomaly_min_time
            Impact[k]<-data_decline
            threshold_norm_anomaly[k]<-threshold_anomaly
            r_threshold<-0.95*data_decline+data_col[anomaly_min_time] 
            data_select<-data_col[anomaly_min_time:(eos_pentad[year])]
            indices <- (which(data_select > r_threshold & 
                                c(FALSE, data_select[-length(data_select)] > r_threshold) & 
                                c(FALSE, FALSE, data_select[c(-length(data_select),-(length(data_select)-1))] > r_threshold))-2)+anomaly_min_time-1
            if(length(indices)>0){
              recovery_time_point[k]<-indices[1]
              Rate[k]<-(data_decline/(indices[1]-anomaly_min_time))
              actual_LAI_threshold <- (threshold_anomaly*mean_sd_col[(c((onset_start[k]-1):indices[1])%%73)+73])+mean_sd_col[c((onset_start[k]-1):indices[1])%%73]
              actual_LAI <- detrend_col[c((onset_start[k]-1):indices[1])]
              total_actual_LAI_loss[k]<-sum((actual_LAI_threshold-actual_LAI),na.rm=T)
            }else{
              actual_LAI_threshold <- (threshold_anomaly*mean_sd_col[(c((onset_start[k]-1):eos_pentad[year])%%73)+73])+mean_sd_col[c((onset_start[k]-1):eos_pentad[year])%%73]
              actual_LAI <- detrend_col[c((onset_start[k]-1):eos_pentad[year])]
              total_actual_LAI_loss[k]<-sum((actual_LAI_threshold-actual_LAI),na.rm=T)
              Rate[k]<-"LAInotrecovery_inthisGS"
              recovery_time_point[k]<-"LAInotrecovery_inthisGS"
            }
          }else{
            impact_time_point[k]<-"LAInotrecovery_inthisGS"
            Impact[k]<-"LAInotrecovery_inthisGS"
            threshold_norm_anomaly[k]<-threshold_anomaly
            actual_LAI_threshold <- (threshold_anomaly*mean_sd_col[(c((onset_start[k]-1):eos_pentad[year])%%73)+73])+mean_sd_col[c((onset_start[k]-1):eos_pentad[year])%%73]
            actual_LAI <- detrend_col[c((onset_start[k]-1):eos_pentad[year])]
            total_actual_LAI_loss[k]<-sum((actual_LAI_threshold-actual_LAI),na.rm=T)
            Rate[k]<-"LAInotrecovery_inthisGS"
            recovery_time_point[k]<-"LAInotrecovery_inthisGS"
          }
        }else{
          total_actual_LAI_loss[k]<-"LAInotdecline"
          Impact[k]<-"LAInotdecline"
          Rate[k]<-"LAInotdecline"
          impact_time_point[k]<-"LAInotdecline"
          recovery_time_point[k]<-"LAInotdecline"
          threshold_norm_anomaly[k]<-"LAInotdecline"
        }
      } 
    }
    return(c(total_actual_LAI_loss,Impact,Rate,impact_time_point,recovery_time_point,threshold_norm_anomaly,rep(NA,500-length(Impact)*6)))
  }else{return(rep(NA,500))}
}

for(i in c(1:665)){
  file<-"path_to_flash_drought"
  if(file.exists(file)){
    flashdrought.csv <- read.csv(file)
    data.csv<-read.csv("path_to_NDVI_norm_anomaly")
    mean.sd.csv<-read.csv("path_to_NDVI_means_and_sds")
    detrend.csv<-read.csv("path_to_detrend_NDVI")
    data.phe<-read.csv("path_to_spring_and_autumn_phenology")
    final.data <- as.data.frame(mapply(resilience, data.csv, mean.sd.csv, detrend.csv, flashdrought.csv, data.phe, SIMPLIFY = T))
    write.csv(final.data,"output_path",row.names=F)
  }
  print(paste0("finish+",i))
}


######### vegetation sensitivity to flash drought is indicated by relative decline (%) ########
cal_percent_decline<-function(actual_col,mean_sd_col,all_col){
  # actual_col:actual NDVI
  # mean_sd_col:NDVI means and sds
  # all_col: output data from the above "resilience" function
  all_col<-na.omit(all_col)
  if(length(all_col)>0 & length(na.omit(actual_col))>0){
    len<-length(all_col)/6
    impact<-as.numeric(all_col[(len+1):(len*2)])
    number<-which(!is.na(impact))
    if(length(number)>0){
      impact_time_point<-as.numeric((all_col[(len*3+1):(len*4)])[number])
      threshold_norm_anomaly<-as.numeric((all_col[(len*5+1):(len*6)])[number])
      
      actual_data<-actual_col[impact_time_point]
      potential_data<-c()
      for(m in 1:length(threshold_norm_anomaly)){
        potential_data[m]<-(threshold_norm_anomaly[m]*mean_sd_col[((impact_time_point[m]-1)%%73)+1+73])+mean_sd_col[((impact_time_point[m]-1)%%73+1)]
      }
      percent_decline<-((actual_data-potential_data)/potential_data)*100 # change to %
      percent_decline[which(potential_data<0.1)]<-NA
      return(c(percent_decline,rep(NA,40-length(percent_decline))))
    }else{return(rep(NA,40))}
  }else{return(rep(NA,40))}
}

for(i in c(1:665)){
  file<-"path_to_output_data_in_the_above_resilience_function"
  if(file.exists(file)){
    all.csv <- read.csv(file)
    actual.csv<-read.csv("path_to_detrend_NDVI")
    mean.sd.csv<-read.csv("path_to_NDVI_means_and_sds")
    final.data <- as.data.frame(mapply(cal_percent_decline, actual.csv, mean.sd.csv, all.csv, SIMPLIFY = T))
    write.csv(final.data,"path_to_output_relative_decline",row.names=F)
  }
  print(paste0("finish+",i))
}


