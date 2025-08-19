############# Flash drought identification ##################

# raw_flash_drought ####
is_non_increasing <- function(x) {
  for (k in seq_len(length(x) - 1)) {
    if (x[k] < x[k + 1]) return(FALSE)
  }
  return(TRUE)
}
flash_drought<-function(percentile){
  if(length(na.omit(percentile))==1533){
    
    point1<-which(percentile<20)
    point2<-point1-1
    point2[point2==0]<-NA
    point2<-point2[!is.na(point2)]
    point3<-point2[which(percentile[point2]>20)]+1
    
    point4<-which(percentile>40)
    point5<-point4+1
    point5[point5==1534]<-NA
    point5<-point5[!is.na(point5)]
    point6<-point5[which(percentile[point5]<40)]
    
    point7<-which(percentile>20)
    point8<-point7-1
    point8[point8==0]<-NA
    point8<-point8[!is.na(point8)]
    point9<-point8[which(percentile[point8]<20)]

    if(length(point3)>0 & length(point6)>0){
      valid_point6 <- integer(0)
      valid_point3 <- integer(0)
      if(length(point6)>1){
        for (i in 1:(length(point6)-1)) {
          possible_point3 <- min(point3[point3 >= point6[i]])
          if(is.numeric(possible_point3) & !is.infinite(possible_point3) & possible_point3<point6[(i+1)]){
            sequence <- percentile[point6[i]:possible_point3]
            rate <- (percentile[point6[i]-1]-percentile[possible_point3])/(possible_point3-(point6[i]-1))
            if (is_non_increasing(sequence) & rate>5) {
              valid_point6 <- c(valid_point6, point6[i])
              valid_point3 <- c(valid_point3, possible_point3)
            }
          }
        }
        for (i in length(point6)) {
          possible_point3 <- min(point3[point3 >= point6[i]])
          if(is.numeric(possible_point3) & !is.infinite(possible_point3)){
            sequence <- percentile[point6[i]:possible_point3]
            rate <- (percentile[point6[i]-1]-percentile[possible_point3])/(possible_point3-(point6[i]-1))
            if (is_non_increasing(sequence) & rate>5) {
              valid_point6 <- c(valid_point6, point6[i])
              valid_point3 <- c(valid_point3, possible_point3)
            }
          }
        }
      }else{
        for (i in length(point6)) {
          possible_point3 <- min(point3[point3 >= point6[i]])
          if(is.numeric(possible_point3) & !is.infinite(possible_point3)){
            sequence <- percentile[point6[i]:possible_point3]
            rate <- (percentile[point6[i]-1]-percentile[possible_point3])/(possible_point3-(point6[i]-1))
            if (is_non_increasing(sequence) & rate>5) {
              valid_point6 <- c(valid_point6, point6[i])
              valid_point3 <- c(valid_point3, possible_point3)
            }
          }
        }
      }
    }else{return(rep(NA,800))}
 
    if(length(point9)>0 & length(valid_point3)>0){
      valid_point9 <- integer(0)
      valid_point3_new <- integer(0)
      valid_point6_new <- integer(0)
      if(length(valid_point3)>1){
        for (k in 1:(length(valid_point3)-1)){
          possible_point9 <- min(point9[point9 >= valid_point3[k]])
          if(is.numeric(possible_point9) & !is.infinite(possible_point9) & possible_point9<valid_point3[(k+1)]){
            valid_point9 <- c(valid_point9, possible_point9)
            valid_point3_new <- c(valid_point3_new, valid_point3[k])
            valid_point6_new <- c(valid_point6_new, valid_point6[k])
          }
        }
        for (k in length(valid_point3)){
          possible_point9 <- min(point9[point9 >= valid_point3[k]])
          if(is.numeric(possible_point9) & !is.infinite(possible_point9)){
            valid_point9 <- c(valid_point9, possible_point9)
            valid_point3_new <- c(valid_point3_new, valid_point3[k])
            valid_point6_new <- c(valid_point6_new, valid_point6[k])
          }
        }
      }else{
        for (k in length(valid_point3)){
          possible_point9 <- min(point9[point9 >= valid_point3[k]])
          if(is.numeric(possible_point9) & !is.infinite(possible_point9)){
            valid_point9 <- c(valid_point9, possible_point9)
            valid_point3_new <- c(valid_point3_new, valid_point3[k])
            valid_point6_new <- c(valid_point6_new, valid_point6[k])
          }
        }
      }
    }else{return(rep(NA,800))}

    if(length(valid_point3_new)>0){
      onset_start <- integer(0)
      first_below_20 <- integer(0)
      onset_end <- integer(0)
      min_soilm <- integer(0)
      drought_end <- integer(0)
      for(m in 1:length(valid_point3_new)){
        point_min_sw<-which.min(percentile[valid_point3_new[m]:valid_point9[m]])+valid_point3_new[m]-1
        if(point_min_sw==valid_point3_new[m]){
          onset_start <- c(onset_start,valid_point6_new[m])
          first_below_20 <- c(first_below_20,valid_point3_new[m])
          onset_end <- c(onset_end,valid_point3_new[m])
          min_soilm <- c(min_soilm,valid_point3_new[m])
          drought_end <- c(drought_end,valid_point9[m])
        }else{
          speed<-percentile[valid_point3_new[m]:(point_min_sw-1)]-percentile[(valid_point3_new[m]+1):point_min_sw]
          if(speed[1]<5){
            onset_start <- c(onset_start,valid_point6_new[m])
            first_below_20 <- c(first_below_20,valid_point3_new[m])
            onset_end <- c(onset_end,valid_point3_new[m])
            min_soilm <- c(min_soilm,point_min_sw)
            drought_end <- c(drought_end,valid_point9[m])
          }else{
            onset_start <- c(onset_start,valid_point6_new[m])
            first_below_20 <- c(first_below_20,valid_point3_new[m])
            min_soilm <- c(min_soilm,point_min_sw)
            drought_end <- c(drought_end,valid_point9[m])
            onset_end <- c(onset_end,(min(which(cumsum(speed<5)==0))+valid_point3_new[m]))
          }
        }
      }
    }else{return(rep(NA,800))}
    number<-which((drought_end-onset_start+1)>=4)
    if(length(number)>0){
      return(c(onset_start[number],first_below_20[number],onset_end[number],min_soilm[number],drought_end[number],rep(NA,800-length(onset_start[number])*5)))
    }else{return(rep(NA,800))}
  }else{return(rep(NA,800))}
}

for(i in c(1:665)){
  file<-paste0("path_to_soil_moisture_percentile_row_data.csv")
  if(file.exists(file)){
    data.percentile<-read.csv(file)
    final.data <- as.data.frame(mapply(flash_drought, data.percentile, SIMPLIFY = T))
    write.csv(final.data,"path_to_flash_drought_row_data.csv",row.names=F)
  }
  print(paste0("finish+",i))
}

# flash_drought_within_growing_season ####
between_GS<-function(drought_data,phe_data){
  sos<-phe_data[1:21]
  eos<-phe_data[22:42]
  if(length(na.omit(sos))>14 & length(na.omit(eos))>14 & length(na.omit(drought_data))>0){
    drought_data<-na.omit(drought_data)
    onset_start<-drought_data[1:(length(drought_data)/5)]
    sos[is.na(sos)]<-round(median(sos[!is.na(sos)]),digits=0)
    eos[is.na(eos)]<-round(median(eos[!is.na(eos)]),digits=0)
    sos_pentad<-73*c(0:20)+ceiling(sos/5)
    eos_pentad<-73*c(0:20)+ceiling(eos/5)
    number<-which(sapply(onset_start,function(x) any(x>=sos_pentad & x<=eos_pentad))==T)
    all.number<-c(number,(number+length(onset_start)),(number+length(onset_start)*2),
                  (number+length(onset_start)*3),(number+length(onset_start)*4))
    drought_inGS<-drought_data[all.number]
    return(c(drought_inGS,rep(NA,(500-length(drought_inGS)))))
  }else{return(rep(NA,500))}
}

# flash_drought_with_ET<0 ####
DOY<-rep(c(1:73),21)
ET_anomaly_lessthan_0<-function(drought_col,ET_col){
  if(length(na.omit(drought_col))>0 & length(na.omit(ET_col))>0){
    drought_col<-na.omit(drought_col)
    len<-length(drought_col)/5
    onset_start<-drought_col[1:len]
    drought_r_end<-drought_col[(len*4+1):(len*5)]
    
    long_term_daily_means <- tapply(ET_col, DOY, mean, na.rm = TRUE) 
    means<-long_term_daily_means[as.character(DOY)]
    anomaly <- round((ET_col - means),digits=4)
    names(anomaly)<-NULL
    ET_anomaly<-c()
    for(k in 1:len){
      ET_anomaly[k]<-round(mean(anomaly[onset_start[k]:drought_r_end[k]],na.rm=T),digits=4)
    }
    number<-which(ET_anomaly<0)
    if(length(number)>0){
      number_final<-c(number,number+len,number+len*2,number+len*3,number+len*4)
      result<-c(drought_col[number_final],rep(NA,400-length(number_final)))
    }else{result<-rep(NA,400)}
    return(result)
  }else{return(rep(NA,400))}
}