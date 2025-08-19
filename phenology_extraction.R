######## phenology extraction ##########

##### HANTS #####
range_matrix<-as.matrix(raster("path_to_study_range"))
time<-seq(1,361,8)
FFT_smooth<-function(data,K=9){
  z<-fft(data) 
  order<-order(Mod(z),decreasing = T) 
  y<-rep(0,times=length(data))
  for(h in 1:length(data)){
    for(k in order[1:K]){
      y[h]=y[h]+z[k]*exp(2*pi*1i*(k-1)*(h-1)/length(data))
    }
  }
  return(Re(y)/length(data))
}
HANTS<-function(ndvi_col,range){
  if(length(na.omit(range))==1){
    result.new<-c()
    for(year in 1:21){
      ndvi.raw.data<-ndvi_col[c(((year-1)*46+1):(year*46))]
      x<-time
      xout<-c(1:361)
      chazhi<-approx(x,ndvi.raw.data,xout,method = "linear")
      yout<-chazhi[[2]]
      result<-FFT_smooth(yout,K=9)
      result_x<-result[x]
      diff<-result_x-ndvi.raw.data
      while(length(diff[diff>=0.05])>0){
        omit<-which(diff>=0.05)
        x<-x[-omit]
        if(length(x)<30){ break }
        ndvi.raw.data<-ndvi.raw.data[-omit]
        xout<-c(min(x):max(x))
        chazhi_new<-approx(x,ndvi.raw.data,xout,method = "linear")
        yout_new<-chazhi_new[[2]]
        result<-FFT_smooth(yout_new,K=9)
        result_x<-result[x]
        diff<-result_x-ndvi.raw.data
      }
      result.new[c(((year-1)*365+1):(year*365))]<-c(rep(NA,(min(xout)-1)),result,rep(NA,(365-max(xout))))
    }
    return(result.new)
  }else{return(c(rep(NA,7665)))}
}
for(i in c(1:665)){
  if(length(na.omit(range_matrix[i,]))>0){
    data.ndvi<-read.csv("path_to_NDVI_row_data.csv")[-1,]
    data.range<-as.data.frame(t(range_matrix[i,]))
    final.data <- as.data.frame(mapply(HANTS, data.ndvi,data.range, SIMPLIFY = T))
    write.csv(final.data,"path_to_saved_row_HANTS_curve.csv",row.names=F)
  }
  print(paste0("finish+",i))
}

##### Extract spring phenology and autumn phenology #####

range_matrix<-as.matrix(raster("path_to_study_range"))
DOY<-c(1:365)
row_indices <- split(1:7665, rep(1:365, length.out=7665))
phenology_extract<-function(curve_col,range){
  if(length(na.omit(range))==1){
    result<-c()
    for (time_point in 1:365) { 
      time_rows <- row_indices[[time_point]]
      result[time_point]<-mean(curve_col[time_rows],na.rm=T)
    } 
    mean.df<-na.omit(data.frame(DOY,result))
    max.meandata<-max(mean.df$result,na.rm=T)
    max.time<-round(median(mean.df$DOY[which(mean.df$result==max.meandata)]),digits=0)
    left.curve<-mean.df[c(1:(max.time-min(mean.df$DOY))),]
    left.min.meandata<-min(left.curve$result,na.rm=T)
    left.yuzhi<-(0.2*(max.meandata-left.min.meandata)+left.min.meandata)
    right.curve<-mean.df[c((max.time-min(mean.df$DOY)+2):length(mean.df$DOY)),]
    right.min.meandata<-min(right.curve$result,na.rm=T)
    right.yuzhi<-(0.2*(max.meandata-right.min.meandata)+right.min.meandata)
    result.SOS<-c()
    result.POS<-c()
    result.EOS<-c()
    for(year in 1:21){
      curve.data<-curve_col[c(((year-1)*365+1):(year*365))]
      curve.df<-na.omit(data.frame(DOY,curve.data))
      POS<-round(median(curve.df$DOY[which(curve.df$curve.data==max(curve.df$curve.data,na.rm=T))]),digits=0)
      result.POS[year]<-POS
      curve1<-curve.df[c(1:(POS-min(curve.df$DOY))),]
      curve1<-curve1[c(which(curve1$curve.data==min(curve1$curve.data,na.rm=T)):length(curve1$DOY)),]
      if(length(which((curve1$curve.data[-nrow(curve1)]-left.yuzhi) * (curve1$curve.data[-1]-left.yuzhi) < 0))>0){
        crossing_points <- min(which((curve1$curve.data[-nrow(curve1)]-left.yuzhi) * (curve1$curve.data[-1]-left.yuzhi) < 0),na.rm=T)
        SOS<-curve1$DOY[which((curve1$curve.data-left.yuzhi)==min(curve1$curve.data[crossing_points:(crossing_points+1)]-left.yuzhi,na.rm=T))]
      }else{SOS<-(-1)}
      result.SOS[year]<-SOS
      curve2<-curve.df[c((POS-min(curve.df$DOY)+2):length(curve.df$DOY)),]
      curve2<-curve2[c(1:which(curve2$curve.data==min(curve2$curve.data,na.rm=T))),]
      if(length(which((curve2$curve.data[-nrow(curve2)]-right.yuzhi) * (curve2$curve.data[-1]-right.yuzhi) < 0))>0){
        crossing_points2 <- max(which((curve2$curve.data[-nrow(curve2)]-right.yuzhi) * (curve2$curve.data[-1]-right.yuzhi) < 0),na.rm=T)
        EOS<-curve2$DOY[which((curve2$curve.data-right.yuzhi)==min(curve2$curve.data[crossing_points2:(crossing_points2+1)]-right.yuzhi,na.rm=T))]
      }else{EOS<-(-1)}
      result.EOS[year]<-EOS
    }
    return(c(left.yuzhi,right.yuzhi,result.SOS,result.POS,result.EOS))
  }else{return(rep(NA,105))}
}

for(i in c(1:665)){
  file<-paste0("path_to_saved_row_HANTS_curve.csv")
  if(file.exists(file) & length(na.omit(range_matrix[i,]))>0){
    data.range<-as.data.frame(t(range_matrix[i,]))
    curve.csv<-read.csv(file)
    final.data <- as.data.frame(mapply(phenology_extract,curve.csv,data.range, SIMPLIFY = T))
    write.csv(final.data,"path_to_saved_phenology.csv",row.names=F)
  }
  print(paste0("finish+",i))
}


##### Remove outliers #####

fun_mad<-function(data){
  mad.data<-mad(data,na.rm=T)
  median.data<-median(data,na.rm=T)
  data[data>(median.data+2.5*mad.data)|data<(median.data-2.5*mad.data)]<-NA
  return(data)
}

