################ Figure3 #############
library(ggplot2)
library(patchwork)
library(export)
library(cowplot)
library(pdp)

library(data.table)
library(ranger)       
library(dplyr)        
library(parallel)
library(doParallel)
library(fastshap)

library(sf)            
library(spatialsample) 
library(rsample)

############### RF #####################
all_rf_weighted <- ranger(formula = relative_decline ~ sos_zscore+MAT+MAP+Aridity_index+Sand_fraction+Clay_fraction+Biodiversity+Ra_lag+NPP_lag+LAI_pre+
                            Biome+root_depth+temp_pre+temp_lag+VPD_pre+VPD_lag+SSRD_pre+SSRD_lag+soilm_pre+drought_timing+onset_speed+pentads_below_20th+soilm_diff, 
                          data = RF_df_clean,verbose = TRUE,importance = "permutation",max.depth = 25,min.node.size=12,num.trees = 500,mtry=6,
                          num.threads=10,respect.unordered.factors = "order",case.weights = RF_df_clean$area_weight,seed = 42)

saveRDS(all_rf_weighted, file = "RF_model_all_GLEAM.rds",compress = TRUE)
saveRDS(RF_df_clean, file = "RF_data_clean_all_GLEAM.rds",compress = TRUE)

############### SHAP  ##############
all_rf_model <- readRDS("RF_model_all_GLEAM.rds")
RF_df_clean <- readRDS("RF_data_clean_all_GLEAM.rds")
X_train <- RF_df_clean %>% select(-row,-col,-lon,-lat,-relative_decline,-area_weight) %>% as.data.frame()
set.seed(2026)
pred_fun <- function(object, newdata) {
  require(ranger) 
  predict(object, data = newdata, num.threads = 1)$predictions
}
registerDoParallel(cores = 8)
start_time <- Sys.time()
shap_values <- explain(all_rf_model,X = X_train,newdata = X_train,nsim = 30,pred_wrapper = pred_fun,parallel = TRUE)
end_time <- Sys.time()
saveRDS(shap_values, file = "fastshap_all_GLEAM.rds", compress = T)

############## pdp #####################
RF_df_clean <- readRDS("RF_data_clean_all_GLEAM.rds")
cl <- makeCluster(8)
registerDoParallel(cl)
pred_wrapper <- function(object, newdata) {
  pred <- predict(object, data = newdata)$predictions
  return(mean(pred, na.rm = TRUE))
}
feature_all<-colnames(RF_df_clean)
for(k in 1:length(feature_all)){
  pdp_data <- partial(all_rf_model,pred.var = feature_all[k], train = RF_df_clean,paropts = list(.packages = "ranger"),
                      pred.fun = pred_wrapper,grid.resolution=100,parallel = TRUE)
  write.csv(pdp_data,paste0(feature_all[k],"_pdp_data.csv"))
}
stopCluster(cl)
registerDoSEQ()

############### spatial cv #############
RF_df_clean <- readRDS("work6/Github代码/202605_新版/Figure_data_and_code/Figure3/RF_data_clean_all_GLEAM.rds")
rf_sf <- st_as_sf(RF_df_clean, coords = c("lon", "lat"), crs = 4326)
set.seed(1234) 
spatial_folds <- spatial_block_cv(data = rf_sf,v = 5,method = "random",cellsize = 5,square = TRUE)
saveRDS(spatial_folds, file = "spatial_folds_5cv_cellsize5.rds",compress = TRUE)

cv_RMSE_results <- numeric(5)
cv_TestRsquare_results <- numeric(5)
cv_OOBRsquare_results <- numeric(5)
cv_OOB_RMSE_results <- numeric(5)
for (i in 1:5) {
  split_i <- spatial_folds$splits[[i]]
  train_data <- training(split_i) %>% st_drop_geometry() 
  test_data  <- testing(split_i) %>% st_drop_geometry()
  
  temp_model <- ranger(formula = relative_decline ~ sos_zscore+MAT+MAP+Aridity_index+Sand_fraction+Clay_fraction+Biodiversity+Ra_lag+NPP_lag+LAI_pre+
                         Biome+root_depth+temp_pre+temp_lag+VPD_pre+VPD_lag+SSRD_pre+SSRD_lag+soilm_pre+drought_timing+onset_speed+pentads_below_20th+soilm_diff,
                       verbose = TRUE,data = train_data,num.trees = 500,max.depth = 25,min.node.size=12,mtry=6,num.threads=10,respect.unordered.factors = "order",
                       case.weights = train_data$area_weight,seed = 42)
  
  preds <- predict(temp_model, data = test_data)$predictions
  cv_OOBRsquare_results[i] <- temp_model$r.squared
  cv_OOB_RMSE_results[i] <- sqrt(temp_model$prediction.error)
  cv_RMSE_results[i] <- sqrt(mean((test_data$relative_decline - preds)^2))
  cv_TestRsquare_results[i] <- cor(test_data$relative_decline, preds)^2
}
result_df<-data.frame(cv_RMSE_results,cv_TestRsquare_results,cv_OOBRsquare_results,cv_OOB_RMSE_results)
write.csv(result_df,"result_df.csv")

################ plot #####################
## importance plot ####
X_train<-RF_df_clean
shap_values<-readRDS("fastshap_all_GLEAM.rds")
shap_values<-as.data.frame(shap_values)
shap_importance <- sort(colMeans(abs(shap_values)), decreasing = TRUE)
imp_df <- data.frame(Feature = names(shap_importance),Importance = shap_importance)
# 绘制水平条形图
color_mapping <- c("LAI_pre"="#A7D398","soilm_diff"="#74A3D4","drought_timing"="#DD847E","temp_pre"="#74A3D4","sos_zscore"="#A7D398","NPP_lag"="#A7D398",
                   "pentads_below_20th"="#DD847E","VPD_lag"="#74A3D4","temp_lag"="#74A3D4","Aridity_index"="#E9D389","MAT"="#E9D389",
                   "VPD_pre"="#74A3D4","Biome"="#A7D398","MAP"="#E9D389","onset_speed"="#DD847E","soilm_pre"="#74A3D4","Biodiversity"="#A7D398",
                   "root_depth"="#A7D398","Clay_fraction"="#969696","Sand_fraction"="#969696","Ra_lag"="#A7D398","SSRD_pre"="#74A3D4","SSRD_lag"="#74A3D4")
imp_df$Color <- ifelse(imp_df$Feature %in% names(color_mapping),color_mapping[imp_df$Feature],"#969696")
name_mapping <- c("LAI_pre"=expression(LAI[pre]),"soilm_diff"=expression(Soilm[diff]),"drought_timing"="Drought timing","temp_pre"=expression(Temp[pre]),
                  "sos_zscore"="GUD","pentads_below_20th"=expression(Pentads["below 20th"]),"VPD_lag"=expression(VPD[lag]),"temp_lag"=expression(Temp[lag]),"Ra_lag"=expression(Ra[lag]),
                  "Aridity_index"="Aridity index","MAT"="MAT","VPD_pre"=expression(VPD[pre]),"Biome"="Biome","MAP"="MAP","onset_speed"="Onset speed","NPP_lag"=expression(NPP[lag]),
                  "soilm_pre"=expression(Soilm[pre]),"Biodiversity"="Biodiversity","root_depth"="Root depth","Clay_fraction"="Clay fraction","Sand_fraction"="Sand fraction",
                  "SSRD_pre"=expression(SSRD[pre]), "SSRD_lag"=expression(SSRD[lag]))
plot_importance<-ggplot(imp_df, aes(x = reorder(Feature, Importance),y = Importance,fill = Color)) +
  geom_bar(stat = "identity", width = 0.75) +
  scale_fill_identity() + 
  scale_y_continuous(limits=c(0,3.5),breaks=seq(0,3,1),labels=seq(0,3,1),expand = c(0, 0))+
  coord_flip() +scale_x_discrete(labels = name_mapping) +  
  labs(title = NULL, x = NULL, y = "Mean |SHAP|") +theme_bw() +
  theme(legend.position = "none",axis.text = element_text(size = 8),axis.title = element_text(size = 8),
        plot.margin = margin(t = 5, r=5, b = 5, l = 0,unit="pt"),panel.grid = element_blank(),
        panel.border = element_rect(color="black",size = 0.7, linetype="solid"),
        axis.line = element_blank(),plot.background = element_blank())

# drought_timing ####
feature<-"drought_timing"
plot_data <- data.frame(FeatureValue = X_train$drought_timing,SHAP = shap_values$drought_timing)
breaks <- seq(0, 70, by = 7)
plot_data$group <- cut(plot_data$FeatureValue,breaks = breaks,include.lowest = TRUE)
group_midpoints <- round(breaks[-length(breaks)] + 3.5, 2)
group_info <- data.frame(group = levels(plot_data$group),group_mid = group_midpoints,color = colorRampPalette(c("#B7D4EA", "#084A91"))(length(group_midpoints)))
plot_data <- merge(plot_data, group_info[, c("group", "group_mid")], by = "group")

drought_timing_plot<-ggplot(plot_data, aes(x = FeatureValue, y = SHAP)) +
  geom_point(color="darkgrey",fill="grey",alpha = 0.3, size = 1) +
  geom_boxplot(aes(x = group_mid, group = group, fill = group),  
               width = 0.38*70/5,outlier.shape = NA,show.legend = FALSE,size=0.3) +
  scale_fill_manual(values = group_info$color) +  
  geom_hline(yintercept = 0,linetype = "dashed",color = "black",size = 0.4) +
  scale_y_continuous(limits=c(-15,37),breaks=seq(-10,30,10),labels=seq(-10,30,10),expand = c(0, 0))+
  scale_x_continuous(limits=c(0,70),breaks=seq(10,70,20),labels=seq(10,70,20))+
  labs(title = NULL, x = "Drought timing", y = "SHAP values") +theme_bw() +
  theme(legend.position = "none",axis.text = element_text(size = 8),axis.title = element_text(size = 8),panel.grid = element_blank(),
        panel.border = element_rect(color="black",size = 0.7, linetype="solid"),
        axis.line = element_blank(),plot.background = element_blank())

pdp_data<-read.csv(paste0(feature,"_pdp_data.csv"))
pdp_data$yhat<-(pdp_data$yhat)
RF_PDP_drought_timing<-ggplot(pdp_data, aes(x = drought_timing, y = yhat)) +
  geom_line(color = "red", size = 0.5) +
  geom_ribbon(aes(ymin = yhat - sd(yhat), ymax = yhat + sd(yhat)), fill = "lightgrey", alpha = 0.5) +  
  labs(y = "Sen (%)",title = NULL,x=NULL) +
  scale_x_continuous(limits=c(0,70),breaks=seq(10,70,20),labels=seq(10,70,20)) +
  scale_y_continuous(limits = c(3,26.5), breaks = seq(5,25,10),labels=seq(5,25,10)) +
  theme_bw() +
  theme(panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 0.7),
        axis.title = element_text(size = 8),axis.line = element_blank(),plot.background = element_blank(),
        axis.text = element_text(size = 8),axis.ticks.x = element_blank(),axis.text.x = element_blank())
combined_plot_drought_timing <- ggdraw() +
  draw_plot(drought_timing_plot) +  
  draw_plot(RF_PDP_drought_timing,x = 0.19,y = 1, width = 0.5, height = 0.35, hjust = 0,vjust = 1)  #left-top

# LAI_pre ####
feature<-"LAI_pre"
plot_data <- data.frame(FeatureValue = X_train$LAI_pre,SHAP = shap_values$LAI_pre)
breaks <- seq(-2.5, 2.5, by = 0.5)
plot_data$group <- cut(plot_data$FeatureValue,breaks = breaks,include.lowest = TRUE)
group_midpoints <- round(breaks[-length(breaks)] + 0.25, 2)
group_info <- data.frame(group = levels(plot_data$group),group_mid = group_midpoints,color = colorRampPalette(c("#B7D4EA", "#084A91"))(length(group_midpoints)))
plot_data <- merge(plot_data, group_info[, c("group", "group_mid")], by = "group")

LAI_pre_plot<-ggplot(plot_data, aes(x = FeatureValue, y = SHAP)) +
  geom_point(color="darkgrey",fill="grey",alpha = 0.3, size = 1) +
  geom_boxplot(aes(x = group_mid, group = group, fill = group),  
               width = 0.38,alpha = 0.8,outlier.shape = NA,show.legend = FALSE,size=0.3) +
  scale_fill_manual(values = group_info$color) +  
  geom_hline(yintercept = 0,linetype = "dashed",color = "black",size = 0.4) +
  scale_y_continuous(limits=c(-14,26),breaks=seq(-10,20,10),labels=seq(-10,20,10),expand = c(0, 0))+
  scale_x_continuous(limits=c(-2.5,2.5),breaks=seq(-2,2,1),labels=seq(-2,2,1))+
  labs(title = NULL, x = expression(LAI[pre]), y = "SHAP values") +theme_bw() +
  theme(legend.position = "none",axis.text = element_text(size = 8),axis.title = element_text(size = 8),panel.grid = element_blank(),
        panel.border = element_rect(color="black",size = 0.7, linetype="solid"),
        axis.line = element_blank(),plot.background = element_blank())

pdp_data<-read.csv(paste0(feature,"_pdp_data.csv"))
pdp_data$yhat<-(pdp_data$yhat)
RF_PDP_LAI_pre<-ggplot(pdp_data, aes(x = LAI_pre, y = yhat)) +
  geom_line(color = "red", size = 0.5) +
  geom_ribbon(aes(ymin = yhat - sd(yhat), ymax = yhat + sd(yhat)), fill = "lightgrey", alpha = 0.5) +  
  labs(y = "Sen (%)",title = NULL,x=NULL) +
  scale_x_continuous(limits = c(-2.5, 2.5), breaks = seq(-2, 2, 1),labels=seq(-2,2,1)) +
  scale_y_continuous(limits = c(4, 23), breaks = seq(10,20, 10),labels=seq(10,20, 10)) +
  theme_bw() +
  theme(panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 0.7),
        axis.title = element_text(size = 8),axis.line = element_blank(),plot.background = element_blank(),
        axis.text = element_text(size = 8),axis.ticks.x = element_blank(),axis.text.x = element_blank())
combined_plot_LAI_pre <- ggdraw() +
  draw_plot(LAI_pre_plot) +  
  draw_plot(RF_PDP_LAI_pre,x = 0.19,y = 1, width = 0.5, height = 0.35, hjust = 0, vjust =1)  

# Aridity_index ####
feature<-"Aridity_index"
plot_data <- data.frame(FeatureValue = X_train$Aridity_index,SHAP = shap_values$Aridity_index)
breaks <- seq(0, 3, by = 0.3)
plot_data$group <- cut(plot_data$FeatureValue,breaks = breaks,include.lowest = TRUE)
group_midpoints <- round(breaks[-length(breaks)] + 0.15, 2)
group_labels <- sprintf("(%.2f, %.2f]", breaks[-length(breaks)], breaks[-1])
group_info <- data.frame(group = levels(plot_data$group),group_mid = group_midpoints,color = colorRampPalette(c("#B7D4EA", "#084A91"))(length(group_midpoints)))
plot_data <- merge(plot_data, group_info[, c("group", "group_mid")], by = "group")

AI_plot<-ggplot(plot_data, aes(x = FeatureValue, y = SHAP)) +
  geom_point(color="darkgrey",fill="grey",alpha = 0.3, size = 1) +
  geom_boxplot(aes(x=group_mid,group = group, fill = group),  
               width = 0.38*3/5,alpha = 0.8,outlier.shape = NA,show.legend = FALSE,size=0.3) +
  scale_fill_manual(values = group_info$color) +  
  geom_hline(yintercept = 0,linetype = "dashed",color = "black",size = 0.4) +
  scale_y_continuous(limits=c(-6.5,16.5),breaks=seq(-5,15,5),labels=seq(-5,15,5),expand = c(0, 0))+
  scale_x_continuous(limits=c(0,3),breaks=seq(0,3,1),labels=seq(0,3,1))+
  labs(title = NULL, x = "Aridity index", y = "SHAP values") +theme_bw() +
  theme(legend.position = "none",axis.text = element_text(size = 8),axis.title = element_text(size = 8),panel.grid = element_blank(),
        panel.border = element_rect(color="black",size = 0.7, linetype="solid"),
        axis.line = element_blank(),plot.background = element_blank())

pdp_data<-read.csv("Aridity_index_pdp_data.csv")
pdp_data$yhat<-(pdp_data$yhat)
RF_PDP_AI<-ggplot(pdp_data, aes(x = Aridity_index, y = yhat)) +
  geom_line(color = "red", size = 0.5) +
  geom_ribbon(aes(ymin = yhat - sd(yhat), ymax = yhat + sd(yhat)), fill = "lightgrey", alpha = 0.5) + 
  labs(y = "Sen (%)",title = NULL,x=NULL) +
  scale_x_continuous(limits=c(0,3),breaks=seq(0,3,1),labels=seq(0,3,1)) +
  scale_y_continuous(limits = c(9.5, 14.5), breaks = seq(10, 14, 2),labels=seq(10, 14, 2)) +
  theme_bw() +
  theme(panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 0.7),
        axis.title = element_text(size = 8),axis.line = element_blank(),plot.background = element_blank(),
        axis.text = element_text(size = 8),axis.ticks.x = element_blank(),axis.text.x = element_blank())
combined_plot_AI <- ggdraw() +
  draw_plot(AI_plot) +  
  draw_plot(RF_PDP_AI,x = 1,y =1, 
            width = 0.5, height = 0.35,  hjust = 1, vjust = 1)  
# MAT ####
feature<-"MAT"
plot_data <- data.frame(FeatureValue = X_train$MAT,SHAP = shap_values$MAT)
breaks <- seq(-15, 25, by = 4)
plot_data$group <- cut(plot_data$FeatureValue,breaks = breaks,include.lowest = TRUE)
group_midpoints <- round(breaks[-length(breaks)] +2, 2)
group_info <- data.frame(group = levels(plot_data$group),group_mid = group_midpoints,color = colorRampPalette(c("#B7D4EA", "#084A91"))(length(group_midpoints)))
plot_data <- merge(plot_data, group_info[, c("group", "group_mid")], by = "group")

MAT_plot<-ggplot(plot_data, aes(x = FeatureValue, y = SHAP)) +
  geom_point(color="darkgrey",fill="grey",alpha = 0.3, size = 1) +
  geom_boxplot(aes(x = group_mid, group = group, fill = group),  
               width = 0.38*40/5,outlier.shape = NA,show.legend = FALSE,size=0.3) +
  scale_fill_manual(values = group_info$color) +  
  geom_hline(yintercept = 0,linetype = "dashed",color = "black",size = 0.4) +
  scale_y_continuous(limits=c(-14.9,22),breaks=seq(-10,20,10),labels=seq(-10,20,10),expand = c(0, 0))+
  scale_x_continuous(limits=c(-15,25),breaks=seq(-10,20,10),labels=seq(-10,20,10))+
  labs(title = NULL, x = "MAT (°C)", y = "SHAP values") +theme_bw() +
  theme(legend.position = "none",axis.text = element_text(size = 8,color = "black"),
        axis.title = element_text(size = 8,color = "black"),panel.grid = element_blank(),
        panel.border = element_rect(color="black",size = 0.7, linetype="solid"),
        axis.line = element_blank(),plot.background = element_blank())

pdp_data<-read.csv("MAT_pdp_data.csv")
pdp_data$yhat<-(pdp_data$yhat)
RF_PDP_MAT<-ggplot(pdp_data, aes(x = MAT, y = yhat)) +
  geom_line(color = "red", size = 0.5) +
  geom_ribbon(aes(ymin = yhat - sd(yhat), ymax = yhat + sd(yhat)), fill = "lightgrey", alpha = 0.5) + 
  labs(y = "Sen (%)",title = NULL,x=NULL) +
  scale_x_continuous(limits=c(-15,25),breaks=seq(-10,20,10),labels=seq(-10,20,10)) +
  scale_y_continuous(limits = c(8, 16.8), breaks = seq(8, 16,4),labels=seq(8, 16,4)) +
  theme_bw() +
  theme(panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 0.7),
        axis.title = element_text(size = 8,color = "black"),axis.line = element_blank(),plot.background = element_blank(),
        axis.text = element_text(size = 8,color = "black"),axis.ticks.x = element_blank(),axis.text.x = element_blank())
combined_plot_MAT <- ggdraw() +
  draw_plot(MAT_plot) +  
  draw_plot(RF_PDP_MAT,x = 1,y =1, width = 0.5, height = 0.35, hjust = 1, vjust = 1) 

# SSRD_lag ####
feature<-"SSRD_lag"
plot_data <- data.frame(FeatureValue = X_train$SSRD_lag,SHAP = shap_values$SSRD_lag)
breaks <- seq(-2.5, 2.5, by = 0.5)
plot_data$group <- cut(plot_data$FeatureValue,breaks = breaks,include.lowest = TRUE)
group_midpoints <- round(breaks[-length(breaks)] + 0.25, 2)
group_info <- data.frame(group = levels(plot_data$group),group_mid = group_midpoints,color = colorRampPalette(c("#B7D4EA", "#084A91"))(length(group_midpoints)))
plot_data <- merge(plot_data, group_info[, c("group", "group_mid")], by = "group")

SSRD_lag_plot<-ggplot(plot_data, aes(x = FeatureValue, y = SHAP)) +
  geom_point(color="darkgrey",fill="grey",alpha = 0.3, size = 1) +
  geom_boxplot(aes(x = group_mid, group = group, fill = group),  
               width = 0.38,alpha = 0.8,outlier.shape = NA,show.legend = FALSE,size=0.3) +
  scale_fill_manual(values = group_info$color) +  #
  geom_hline(yintercept = 0,linetype = "dashed",color = "black",size = 0.4) +
  scale_y_continuous(limits=c(-16,21),breaks=seq(-15,15,10),labels=seq(-15,15,10),expand = c(0, 0))+
  scale_x_continuous(limits=c(-2.5,2.5),breaks=seq(-2,2,1),labels=seq(-2,2,1))+
  labs(title = NULL, x = expression(SSRD[lag]), y = "SHAP values") +theme_bw() +
  theme(legend.position = "none",axis.text = element_text(size = 8),axis.title = element_text(size = 8),panel.grid = element_blank(),
        panel.border = element_rect(color="black",size = 0.7, linetype="solid"),
        axis.line = element_blank(),plot.background = element_blank())

pdp_data<-read.csv(paste0("SSRD_lag_pdp_data.csv"))
pdp_data$yhat<-(pdp_data$yhat)
RF_PDP_SSRD_lag<-ggplot(pdp_data, aes(x = SSRD_lag, y = yhat)) +
  geom_line(color = "red", size = 0.5) +
  geom_ribbon(aes(ymin = yhat - sd(yhat), ymax = yhat + sd(yhat)), fill = "lightgrey", alpha = 0.5) +  
  labs(y = "Sen (%)",title = NULL,x=NULL) +
  scale_x_continuous(limits=c(-2.5,2.5),breaks=seq(-2,2,1),labels=seq(-2,2,1)) +
  scale_y_continuous(limits = c(3.3, 14.2), breaks = seq(4,12, 4),labels=seq(4,12, 4)) +
  theme_bw() +
  theme(panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 0.7),
        axis.title = element_text(size = 8),axis.line = element_blank(),plot.background = element_blank(),
        axis.text = element_text(size = 8),axis.ticks.x = element_blank(),axis.text.x = element_blank())
combined_plot_SSRD_lag <- ggdraw() +
  draw_plot(SSRD_lag_plot) +  
  draw_plot(RF_PDP_SSRD_lag,x = 1,y = 1, width = 0.5, height = 0.35,hjust = 1, vjust =1)      

# GUD ####
feature<-"sos"
plot_data <- data.frame(FeatureValue = X_train$sos_zscore,SHAP = shap_values$sos_zscore)
breaks <- seq(-3, 3, by = 0.6)
plot_data$sos_group <- cut(plot_data$FeatureValue,breaks = breaks,include.lowest = TRUE)
group_midpoints <- round(breaks[-length(breaks)] + 0.3, 2)
group_labels <- sprintf("(%.2f, %.2f]", breaks[-length(breaks)], breaks[-1])
group_info <- data.frame(sos_group = levels(plot_data$sos_group),group_mid = group_midpoints,color = colorRampPalette(c("#B7D4EA", "#084A91"))(length(group_midpoints)))
plot_data <- merge(plot_data, group_info[, c("sos_group", "group_mid")], by = "sos_group")

sos_plot<-ggplot(plot_data, aes(x = FeatureValue, y = SHAP)) +
  geom_point(color="darkgrey",fill="grey",alpha = 0.3, size = 1) +
  geom_boxplot(aes(x = group_mid, group = sos_group, fill = sos_group),  
               width = 0.38*6/5,alpha = 0.8,outlier.shape = NA,show.legend = FALSE,size=0.3) +
  scale_fill_manual(values = group_info$color) +  
  geom_hline(yintercept = 0,linetype = "dashed",color = "black",size = 0.4) +
  scale_y_continuous(limits=c(-5.2,15),breaks=seq(-5,15,5),labels=seq(-5,15,5),expand = c(0, 0))+
  scale_x_continuous(limits=c(-3,3),breaks=seq(-3,3,1),labels=seq(-3,3,1))+
  labs(title = NULL, x = "GUD", y = "SHAP values") +theme_bw() +
  theme(legend.position = "none",axis.text = element_text(size = 8),axis.title = element_text(size = 8),panel.grid = element_blank(),
        panel.border = element_rect(color="black",size = 0.7, linetype="solid"),
        axis.line = element_blank(),plot.background = element_blank())

pdp_data<-read.csv("sos_zscore_pdp_data.csv")
pdp_data$yhat<-(pdp_data$yhat)
RF_PDP_sos<-ggplot(pdp_data, aes(x = sos_zscore, y = yhat)) +
  geom_line(color = "red", size = 0.5) +
  geom_ribbon(aes(ymin = yhat - sd(yhat), ymax = yhat + sd(yhat)), fill = "lightgrey", alpha = 0.5) + 
  labs(y = "Sen (%)",title = NULL,x=NULL) +
  scale_x_continuous(limits = c(-3, 3), breaks = seq(-3, 3, 1),labels=seq(-3,3,1)) +
  scale_y_continuous(limits = c(8.6, 14.5), breaks = seq(10, 14, 2),labels=seq(10, 14, 2)) +
  theme_bw() +
  theme(panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 0.7),
        axis.title = element_text(size = 8),axis.line = element_blank(),plot.background = element_blank(),
        axis.text = element_text(size = 8),axis.ticks.x = element_blank(),axis.text.x = element_blank())
combined_plot_sos <- ggdraw() +
  draw_plot(sos_plot) +  
  draw_plot(RF_PDP_sos,x = 1,y = 1, width = 0.5, height = 0.35, hjust = 1,vjust = 1)

# Biome ####
feature<-"Biome"
plot_data <- data.frame(FeatureValue = X_train$Biome,SHAP = shap_values$Biome)
plot_data$FeatureValue<-as.numeric(plot_data$FeatureValue)
breaks <- seq(0, 7, by = 1)
plot_data$group <- cut(plot_data$FeatureValue,breaks = breaks,include.lowest = TRUE)
group_midpoints <- round(breaks[-length(breaks)] + 1, 2)
group_info <- data.frame(group = levels(plot_data$group),group_mid = group_midpoints,color = colorRampPalette(c("#B7D4EA", "#084A91"))(7))
plot_data <- merge(plot_data, group_info[, c("group", "group_mid")], by = "group")
desired_order <- c(1, 2, 5, 6, 3, 4, 7)
plot_data$FeatureValue <- factor(plot_data$FeatureValue, levels = desired_order)

biome_plot<-ggplot(plot_data, aes(x = FeatureValue, y = SHAP)) +
  geom_point(color="darkgrey",fill="grey",alpha = 0.3, size = 1) +
  geom_boxplot(aes(group = group, fill = group),  
               width = 0.38*7/5,alpha = 0.8,outlier.shape = NA,show.legend = FALSE,size=0.3) +
  scale_fill_manual(values = group_info$color) +  
  geom_hline(yintercept = 0,linetype = "dashed",color = "black",size = 0.4) +
  scale_y_continuous(limits=c(-7,17),breaks=seq(-5,15,5),labels=seq(-5,15,5),expand = c(0, 0))+
  scale_x_discrete(labels = c("Tundra", "BF", "TBMF", "TCF", "TGA", "DGA", "MGA")) +
  labs(title = NULL, x = "Biome", y = "SHAP values") +theme_bw() +
  theme(legend.position = "none",axis.text = element_text(size = 8),axis.title = element_text(size = 8),panel.grid = element_blank(),
        panel.border = element_rect(color="black",size = 0.7, linetype="solid"),
        axis.line = element_blank(),plot.background = element_blank())

pdp_data<-read.csv(paste0(feature,"_pdp_data.csv"))
pdp_data$yhat<-(pdp_data$yhat)
pdp_data$biome<-c("Tundra", "BF", "TGA", "DGA", "TBMF", "TCF", "MGA")
desired_order <- c("Tundra", "BF", "TBMF", "TCF", "TGA", "DGA", "MGA")
pdp_data$biome <- factor(pdp_data$biome, levels = desired_order)
RF_PDP_biome<-ggplot(pdp_data, aes(x = biome, y = yhat)) +
  geom_col(fill = "red", alpha = 0.8,width=0.5) + labs(y = "Sen (%)",title = NULL,x=NULL) +
  scale_y_continuous(limits = c(0, 12.5), breaks = seq(0, 10, 5),labels=seq(0, 10, 5),expand=c(0,0)) +
  theme_bw() +
  theme(panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 0.7),
        axis.title = element_text(size = 8),axis.line = element_blank(),plot.background = element_blank(),
        axis.text = element_text(size = 8),axis.ticks.x = element_blank(),axis.text.x = element_blank())
combined_plot_biome <- ggdraw() +
  draw_plot(biome_plot) + 
  draw_plot(RF_PDP_biome,x = 0.16,y = 1, width = 0.5, height = 0.35,hjust = 0, vjust = 1)    

# biodiversity ####
feature<-"biodiversity"
plot_data <- data.frame(FeatureValue = X_train$Biodiversity,SHAP = shap_values$Biodiversity)
breaks <- seq(0, 2000, by = 200)
plot_data$group <- cut(plot_data$FeatureValue,breaks = breaks,include.lowest = TRUE)
group_midpoints <- round(breaks[-length(breaks)] + 150, 2)
group_labels <- sprintf("(%.2f, %.2f]", breaks[-length(breaks)], breaks[-1])
group_info <- data.frame(group = levels(plot_data$group),group_mid = group_midpoints,color = colorRampPalette(c("#B7D4EA", "#084A91"))(length(group_midpoints)))
plot_data <- merge(plot_data, group_info[, c("group", "group_mid")], by = "group")

biodi_plot<-ggplot(plot_data, aes(x = FeatureValue, y = SHAP)) +
  geom_point(color="darkgrey",fill="grey",alpha = 0.3, size = 1) +
  geom_boxplot(aes(x=group_mid,group = group, fill = group),  
               width = 0.38*2000/5,alpha = 0.8,outlier.shape = NA,show.legend = FALSE,size=0.3) +
  scale_fill_manual(values = group_info$color) +  
  geom_hline(yintercept = 0,linetype = "dashed",color = "black",size = 0.4) +
  scale_y_continuous(limits=c(-7,15),breaks=seq(-5,15,5),labels=seq(-5,15,5),expand = c(0, 0))+
  scale_x_continuous(limits=c(0,2000),breaks=seq(0,2000,500),labels=seq(0,2000,500))+
  labs(title = NULL, x = "Biodiversity", y = "SHAP values") +theme_bw() +
  theme(legend.position = "none",axis.text = element_text(size = 8),axis.title = element_text(size = 8),panel.grid = element_blank(),
        panel.border = element_rect(color="black",size = 0.7, linetype="solid"),
        axis.line = element_blank(),plot.background = element_blank())

pdp_data<-read.csv("Biodiversity_pdp_data.csv")
pdp_data$yhat<-(pdp_data$yhat)
RF_PDP_biodi<-ggplot(pdp_data, aes(x = Biodiversity, y = yhat)) +
  geom_line(color = "red", size = 0.5) +
  geom_ribbon(aes(ymin = yhat - sd(yhat), ymax = yhat + sd(yhat)), fill = "lightgrey", alpha = 0.5) + 
  labs(y = "Sen (%)",title = NULL,x=NULL) +
  scale_x_continuous(limits=c(0,2000),breaks=seq(0,2000,500),labels=seq(0,2000,500)) +
  scale_y_continuous(limits = c(14.7, 9.7), breaks = seq(10, 14, 2),labels=seq(10, 14, 2)) +
  theme_bw() +
  theme(panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 0.7),
        axis.title = element_text(size = 8),axis.line = element_blank(),plot.background = element_blank(),
        axis.text = element_text(size = 8),axis.ticks.x = element_blank(),axis.text.x = element_blank())
combined_plot_biodi <- ggdraw() +
  draw_plot(biodi_plot) +  
  draw_plot(RF_PDP_biodi,x = 1,y = 1, width = 0.5, height = 0.35, hjust = 1,vjust = 1) 

# MAP ####
feature<-"MAP"
plot_data <- data.frame(FeatureValue = X_train$MAP,SHAP = shap_values$MAP)
breaks <- seq(0, 1600, by = 160)
plot_data$group <- cut(plot_data$FeatureValue,breaks = breaks,include.lowest = TRUE)
group_midpoints <- round(breaks[-length(breaks)] +80, 2)
group_info <- data.frame(group = levels(plot_data$group),group_mid = group_midpoints,color = colorRampPalette(c("#B7D4EA", "#084A91"))(length(group_midpoints)))
plot_data <- merge(plot_data, group_info[, c("group", "group_mid")], by = "group")

MAP_plot<-ggplot(plot_data, aes(x = FeatureValue, y = SHAP)) +
  geom_point(color="darkgrey",fill="grey",alpha = 0.3, size = 1) +
  geom_boxplot(aes(x = group_mid, group = group, fill = group),  
               width = 0.38*1600/5,outlier.shape = NA,show.legend = FALSE,size=0.3) +
  scale_fill_manual(values = group_info$color) + 
  geom_hline(yintercept = 0,linetype = "dashed",color = "black",size = 0.4) +
  scale_y_continuous(limits=c(-5.2,12),breaks=seq(-5,10,5),labels=seq(-5,10,5),expand = c(0, 0))+
  scale_x_continuous(limits=c(0,1600),breaks=seq(0,1600,500),labels=seq(0,1600,500))+
  labs(title = NULL, x = "MAP (mm)", y = "SHAP values") +theme_bw() +
  theme(legend.position = "none",axis.text = element_text(size = 8,color = "black"),
        axis.title = element_text(size = 8,color = "black"),panel.grid = element_blank(),
        panel.border = element_rect(color="black",size = 0.7, linetype="solid"),
        axis.line = element_blank(),plot.background = element_blank())

pdp_data<-read.csv("MAP_pdp_data.csv")
pdp_data$yhat<-(pdp_data$yhat)
RF_PDP_MAP<-ggplot(pdp_data, aes(x = MAP, y = yhat)) +
  geom_line(color = "red", size = 0.5) +
  geom_ribbon(aes(ymin = yhat - sd(yhat), ymax = yhat + sd(yhat)), fill = "lightgrey", alpha = 0.5) +  
  labs(y = "Sen (%)",title = NULL,x=NULL) +
  scale_x_continuous(limits=c(0,1600),breaks=seq(0,1600,500),labels=seq(0,1600,500)) +
  scale_y_continuous(limits = c(10.6, 13.6), breaks = seq(11, 13,1),labels=seq(11, 13,1)) +
  theme_bw() +
  theme(panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 0.7),
        axis.title = element_text(size = 8,color = "black"),axis.line = element_blank(),plot.background = element_blank(),
        axis.text = element_text(size = 8,color = "black"),axis.ticks.x = element_blank(),axis.text.x = element_blank())
combined_plot_MAP <- ggdraw() +
  draw_plot(MAP_plot) +  
  draw_plot(RF_PDP_MAP,x = 1,y =1,width = 0.5, height = 0.35,hjust = 1, vjust = 1) 

# NPP_lag ####
feature<-"NPP_lag"
plot_data <- data.frame(FeatureValue = X_train$NPP_lag,SHAP = shap_values$NPP_lag)
breaks <- seq(-3, 2, by = 0.5)
plot_data$group <- cut(plot_data$FeatureValue,breaks = breaks,include.lowest = TRUE)
group_midpoints <- round(breaks[-length(breaks)] + 0.25, 2)
group_info <- data.frame(group = levels(plot_data$group),group_mid = group_midpoints,color = colorRampPalette(c("#B7D4EA", "#084A91"))(length(group_midpoints)))
plot_data <- merge(plot_data, group_info[, c("group", "group_mid")], by = "group")

NPP_lag_plot<-ggplot(plot_data, aes(x = FeatureValue, y = SHAP)) +
  geom_point(color="darkgrey",fill="grey",alpha = 0.3, size = 1) +
  geom_boxplot(aes(x = group_mid, group = group, fill = group),  
               width = 0.38,alpha = 0.8,outlier.shape = NA,show.legend = FALSE,size=0.3) +
  scale_fill_manual(values = group_info$color) +  
  geom_hline(yintercept = 0,linetype = "dashed",color = "black",size = 0.4) +
  scale_y_continuous(limits=c(-6.05,11),breaks=seq(-4,10,4),labels=seq(-4,10,4),expand = c(0, 0))+
  scale_x_continuous(limits=c(-3,2),breaks=seq(-3,2,1),labels=seq(-3,2,1))+
  labs(title = NULL, x = expression(NPP[lag]), y = "SHAP values") +theme_bw() +
  theme(legend.position = "none",axis.text = element_text(size = 8),axis.title = element_text(size = 8),panel.grid = element_blank(),
        panel.border = element_rect(color="black",size = 0.7, linetype="solid"),
        axis.line = element_blank(),plot.background = element_blank())

pdp_data<-read.csv(paste0("NPP_lag_pdp_data.csv"))
pdp_data$yhat<-(pdp_data$yhat)
RF_PDP_NPP_lag<-ggplot(pdp_data, aes(x = NPP_lag, y = yhat)) +
  geom_line(color = "red", size = 0.5) +
  geom_ribbon(aes(ymin = yhat - sd(yhat), ymax = yhat + sd(yhat)), fill = "lightgrey", alpha = 0.5) +  
  labs(y = "Sen (%)",title = NULL,x=NULL) +
  scale_x_continuous(limits=c(-3,2),breaks=seq(-3,2,1),labels=seq(-3,2,1)) +
  scale_y_continuous(limits = c(9.8, 13.2), breaks = seq(10,12, 2),labels=seq(10,12, 2)) +
  theme_bw() +
  theme(panel.grid = element_blank(),panel.border = element_rect(color = "black", size = 0.7),
        axis.title = element_text(size = 8),axis.line = element_blank(),plot.background = element_blank(),
        axis.text = element_text(size = 8),axis.ticks.x = element_blank(),axis.text.x = element_blank())
combined_plot_NPP_lag <- ggdraw() +
  draw_plot(NPP_lag_plot) +  
  draw_plot(RF_PDP_NPP_lag,x = 1,y = 1, width = 0.5, height = 0.35, hjust = 1, vjust =1)  


#### 输出 ####
final_plot<-plot_grid(plot_grid(plot_importance,combined_plot_SSRD_lag,combined_plot_biodi, ncol = 1,rel_heights=c(2,1,1)),
                      plot_grid(combined_plot_drought_timing,combined_plot_AI,combined_plot_sos,combined_plot_MAP, ncol = 1,rel_heights=c(1,1,1,1)),
                      plot_grid(combined_plot_LAI_pre,combined_plot_MAT,combined_plot_biome,combined_plot_NPP_lag,ncol = 1,rel_heights=c(1,1,1,1)),
                      nrow = 1, rel_widths = c(1, 1, 1)) 
graph2tif(final_plot,"Figure3_raw",width=7.3,height=8.2,dpi=500,tiffcompression="lzw")
