############ Figure1 ###############
library(terra)
library(data.table)
library(dplyr)
library(ggplot2)
library(patchwork)
library(export)
library(cowplot)
library(RColorBrewer)
library(raster) 

#### Calculate relative_decline/impact_min_zscore under early/late spring phenology within 10x10 pixels (1°) ####
safe_wilcox <- function(x, y) {
  x <- na.omit(x)
  y <- na.omit(y)
  if (length(x) > 1 && length(y) > 1) {
    res <- tryCatch(wilcox.test(x, y, exact = FALSE)$p.value, error = function(e) NA_real_)
    return(res)
  } else {return(NA_real_)}
}
process_to_tif <- function(input_csv, output_tif) {
  dt <- fread(input_csv)
  dt <- dt[sos_zscore <= -0.5 | sos_zscore >= 0.5]
  dt[, grid_row := ceiling(row / 10)]
  dt[, grid_col := ceiling(col / 10)]
  res_grid <- dt[, .(
    all_impact_min_zscore     = median(impact_min_zscore, na.rm = TRUE),
    earlier_impact_min_zscore = median(impact_min_zscore[sos_zscore <= -0.5], na.rm = TRUE),
    later_impact_min_zscore   = median(impact_min_zscore[sos_zscore >= 0.5], na.rm = TRUE),
    all_relative_decline      = median(relative_decline, na.rm = TRUE),
    earlier_relative_decline  = median(relative_decline[sos_zscore <= -0.5], na.rm = TRUE),
    later_relative_decline    = median(relative_decline[sos_zscore >= 0.5], na.rm = TRUE),
    pval_impact_min_zscore    = safe_wilcox(impact_min_zscore[sos_zscore <= -0.5], impact_min_zscore[sos_zscore >= 0.5]),
    pval_relative_decline     = safe_wilcox(relative_decline[sos_zscore <= -0.5], relative_decline[sos_zscore >= 0.5])
  ), by = .(grid_row, grid_col)]
  res_grid[, diff_impact_min_zscore := earlier_impact_min_zscore - later_impact_min_zscore]
  res_grid[, diff_relative_decline := earlier_relative_decline - later_relative_decline]
  var_names <- c("all_impact_min_zscore", "earlier_impact_min_zscore", "later_impact_min_zscore", "diff_impact_min_zscore", "pval_impact_min_zscore",
                 "all_relative_decline", "earlier_relative_decline", "later_relative_decline", "diff_relative_decline", "pval_relative_decline")
  r_out <- rast(nrows = 67, ncols = 360, nlyrs = length(var_names), ext = c(-180, 180, 23, 90)) 
  names(r_out) <- var_names
  for (v in var_names) {
    mat <- matrix(NA, nrow = 67, ncol = 360)
    mat[cbind(res_grid$grid_row, res_grid$grid_col)] <- res_grid[[v]]
    r_out[[v]] <- mat
    single_tif_name <- paste0(output_tif, "/", v, ".tif")
    writeRaster(r_out[[v]], single_tif_name, overwrite = TRUE)
  }
}
task<-list(in_csv = "RF_combined.csv", out_tif = "GLEAM")
process_to_tif(task$in_csv, task$out_tif)


pval<-rast("GLEAM/pval_impact_min_zscore.tif")
pval[pval>=0.05]<-NA
pval[pval<0.05]<-1
writeRaster(pval, "GLEAM/pval_sig_impact_min_zscore.tif")
pval<-rast("GLEAM/pval_relative_decline.tif")
pval[pval>=0.05]<-NA
pval[pval<0.05]<-1
writeRaster(pval, "GLEAM/pval_sig_relative_decline.tif")

diff_sig<-rast("GLEAM/diff_impact_min_zscore.tif")*rast("GLEAM/pval_sig_impact_min_zscore.tif")
writeRaster(diff_sig, "GLEAM/diff_sig_impact_min_zscore.tif")

diff_sig<-rast("GLEAM/diff_relative_decline.tif")*rast("GLEAM/pval_sig_relative_decline.tif")
writeRaster(diff_sig, "GLEAM/diff_sig_relative_decline.tif")

# Calculate proportions
diff_sig<-rast("GLEAM/diff_sig_relative_decline.tif")
diff<-rast("GLEAM/diff_relative_decline.tif")
cell_area <- cellSize(diff, unit = "km")
v_all <- values(diff, mat = FALSE)
v_sig <- values(diff_sig, mat = FALSE)
v_area <- values(cell_area, mat = FALSE)
valid_all_idx <- !is.na(v_all)
valid_sig_idx <- !is.na(v_sig)
total_valid_area <- sum(v_area[valid_all_idx], na.rm = TRUE)          
total_sig_area <- sum(v_area[valid_sig_idx], na.rm = TRUE)            
total_nonsig_area <- total_valid_area - total_sig_area                
sig_positive_area <- sum(v_area[which(v_sig > 0)], na.rm = TRUE)
sig_negative_area <- sum(v_area[which(v_sig < 0)], na.rm = TRUE)
pct_nonsig <- (total_nonsig_area / total_valid_area) * 100
pct_sig    <- (total_sig_area / total_valid_area) * 100
pct_pos_in_sig <- (sig_positive_area / total_valid_area) * 100
pct_neg_in_sig <- (sig_negative_area / total_valid_area) * 100

cat(sprintf("Total valid area                  : %.2f km²\n", total_valid_area))
cat(sprintf("Non-significant area proportion   : %.2f%%\n", pct_nonsig))
cat(sprintf("Significant area proportion       : %.2f%%\n", pct_sig))
cat(sprintf("Positive (impact increase) ratio  : %.2f%%\n", pct_pos_in_sig))
cat(sprintf("Negative (impact decrease) ratio  : %.2f%%\n", pct_neg_in_sig))

########### plot ##############

###### Frequency distribution histograms & latitudinal line plots for early vs. late spring phenology ########

rgb_data <- data.frame(r = c(226,167,85,37,5),g = c(237,208,158,107,48),b = c(243,228,201,175,97))
blue_red_color<-apply(rgb_data, 1, function(row) {rgb(row["r"], row["g"], row["b"], maxColorValue = 255)})
data_df<-as.data.frame(as(raster("earlier_relative_decline.tif"),"SpatialPixelsDataFrame"))
data_df1<-as.data.frame(as(raster("later_relative_decline.tif"),"SpatialPixelsDataFrame"))
seq1<-c(0,3,6,9,12,1000)
fre<-c()
fre1<-c()
for(u in 1:5){
  fre[u]<-(length(which(seq1[u]<data_df[,1] & data_df[,1]<=seq1[(u+1)]))/length(data_df[,1]))*100
  fre1[u]<-(length(which(seq1[u]<data_df1[,1] & data_df1[,1]<=seq1[(u+1)]))/length(data_df1[,1]))*100
}
Interval<-c("1","2","3","4","5")
df<-data.frame(Interval,fre,fre1)
df$Interval<-as.factor(df$Interval)
df$X<-c(1:5)

g1<-ggplot(data=df,aes(y =fre, x =X)) +
  geom_bar(aes(fill=Interval),stat = "identity",width=0.7)+
  scale_x_continuous(breaks=c(seq(0.5,4.5,1),5.5),labels=c(0,3,6,9,12,">12"))+
  scale_y_continuous(limits=c(0,47),expand = c(0,0))+
  theme_bw()+
  theme(legend.position = "none",axis.text=element_text(size=10,color="black"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        plot.margin = margin(t = 0, r=10, b = 1, l = 0),
        panel.border = element_rect(color="black",size = 1, linetype="solid"),
        axis.line = element_blank(),plot.background = element_blank())+
  labs(x = "Sensitivity (%)",y = "Frequency (%)",title = NULL)+
  scale_fill_manual(values = blue_red_color)
g2<-ggplot(data=df,aes(y =fre1, x =X)) +
  geom_bar(aes(fill=Interval),stat = "identity",width=0.7)+
  scale_x_continuous(breaks=c(seq(0.5,4.5,1),5.5),labels=c(0,3,6,9,12,">12"))+
  scale_y_continuous(limits=c(0,47),expand = c(0,0))+
  theme_bw()+
  theme(legend.position = "none",axis.text=element_text(size=10,color="black"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        plot.margin = margin(t = 0, r=10, b = 1, l = 0),
        panel.border = element_rect(color="black",size = 1, linetype="solid"),
        axis.line = element_blank(),plot.background = element_blank())+
  labs(x = "Sensitivity (%)",y = "Frequency (%)",title = NULL)+
  scale_fill_manual(values = blue_red_color)

matrix<-as.matrix(raster("earlier_relative_decline.tif"))
row_means <- rowMeans(matrix, na.rm = TRUE)
row_sd <- apply(matrix, 1, sd, na.rm = TRUE)
n_non_na <- apply(matrix, 1, function(x) sum(!is.na(x)))
row_se <- row_sd / sqrt(n_non_na)
impact_stats<-data.frame(mean_impact=row_means,SE=row_se,SD=row_sd,lat=seq(89.5,23,-1))
impact_stats<-na.omit(impact_stats)
matrix<-as.matrix(raster("later_relative_decline.tif"))
row_means <- rowMeans(matrix, na.rm = TRUE)
row_sd <- apply(matrix, 1, sd, na.rm = TRUE)
n_non_na <- apply(matrix, 1, function(x) sum(!is.na(x)))
row_se <- row_sd / sqrt(n_non_na)
impact_stats1<-data.frame(mean_impact=row_means,SE=row_se,SD=row_sd,lat=seq(89.5,23,-1))
impact_stats1<-na.omit(impact_stats1)

g3<-ggplot(impact_stats, aes(x = lat, y = mean_impact)) +
  geom_ribbon(aes(ymin = mean_impact - SE, ymax = mean_impact + SE),fill = "grey60", alpha = 0.5) +
  geom_line(color = "black", linewidth = 0.6) +
  theme_bw() + labs(x="Latitude (°N)",y="Sensitivity (%)")+  
  scale_x_continuous(position = "top",limits=c(23.4,77),breaks=seq(30,70,10),labels=seq(30,70,10))+
  scale_y_continuous(limits=c(4,15),breaks=seq(0,15,5),labels=seq(0,15,5))+
  coord_flip() +
  theme(axis.title.y.right = element_text(angle = 90, vjust = 0.5),plot.margin = margin(t = 0, r=0, b = 1, l = 5,unit="pt"),
        axis.text=element_text(size=10,color="black"),panel.grid = element_blank(),
        panel.border = element_rect(color="black",size = 1, linetype="solid"),
        axis.line = element_blank(),plot.background = element_blank())
figure2d<-g1+g3+plot_layout(widths = c(1,1))
graph2tif(figure2d,file="Figure_d",width=3.8,height=3.1,dpi=500,tiffcompression="lzw")

g4<-ggplot(impact_stats1, aes(x = lat, y = mean_impact)) +
  geom_ribbon(aes(ymin = mean_impact - SE, ymax = mean_impact + SE),fill = "grey60", alpha = 0.5) +
  geom_line(color = "black", linewidth = 0.6) +
  theme_bw() + labs(x="Latitude (°N)",y="Sensitivity (%)")+
  scale_x_continuous(position = "top",limits=c(23.4,77),breaks=seq(30,70,10),labels=seq(30,70,10))+
  scale_y_continuous(limits=c(4,15),breaks=seq(0,15,5),labels=seq(0,15,5))+
  coord_flip() +
  theme(axis.title.y.right = element_text(angle = 90, vjust = 0.5),plot.margin = margin(t = 0, r=0, b = 1, l = 5,unit="pt"),
        axis.text=element_text(size=10,color="black"),panel.grid = element_blank(),
        panel.border = element_rect(color="black",size = 1, linetype="solid"),
        axis.line = element_blank(),plot.background = element_blank())
figure2e<-g2+g4+plot_layout(widths = c(1,1))
graph2tif(figure2e,file="Figure_e",width=3.8,height=3.1,dpi=500,tiffcompression="lzw")

#### Calculate the proportion of relative_decline diff across different biomes ####

filelist<-list.files("BIOME_shp",full.names = T,pattern=".shp$")[c(1,2,6,7,8,9,14)]
all_data<-rast("diff_relative_decline.tif")
sig_data<-rast("diff_sig_relative_decline.tif")
result<-matrix(data=NA,nrow=7,ncol=5)
for(k in 1:7){
  biome_all<-mask(crop(all_data,vect(filelist[k])),vect(filelist[k]))
  biome_sig<-mask(crop(sig_data,vect(filelist[k])),vect(filelist[k]))
  result[k,1]<-as.numeric(length(which(na.omit(as.data.frame(biome_all))[,1]>0))/length(na.omit(as.data.frame(biome_all))[,1]))
  result[k,2]<-as.numeric(length(which(na.omit(as.data.frame(biome_all))[,1]<0))/length(na.omit(as.data.frame(biome_all))[,1]))
  result[k,3]<-as.numeric(length(which(na.omit(as.data.frame(biome_sig))[,1]>0))/length(na.omit(as.data.frame(biome_all))[,1]))
  result[k,4]<-as.numeric(length(which(na.omit(as.data.frame(biome_sig))[,1]<0))/length(na.omit(as.data.frame(biome_all))[,1]))
}
result[,1]<-result[,1]-result[,3]
result[,2]<-result[,2]-result[,4]
result[,5]<-c("BF","DGA","MGA","TBMF","TCF","TGA","Tundra")
plot_df<-data.frame(biome=rep(c("BF","DGA","MGA","TBMF","TCF","TGA","Tundra"),each=4),group=rep(c("positive","negative"),14),
                    value=c(result[1,1:4],result[2,1:4],result[3,1:4],result[4,1:4],result[5,1:4],result[6,1:4],result[7,1:4]),
                    significance=rep(c("non-significant","non-significant","significant","significant"),7))
plot_df$value<-as.numeric(plot_df$value)
plot_df$value[which(plot_df$group=="negative")]<-(-plot_df$value[which(plot_df$group=="negative")])
colors <- c("non-significant.positive" = "grey","significant.positive" = "red",
            "non-significant.negative" = "grey","significant.negative" = "blue")
plot_df$fill<-interaction(plot_df$significance, plot_df$group)
plot_df$biome<-factor(plot_df$biome,level=c("Tundra","BF","TBMF","TCF","TGA","DGA","MGA"))

plot_biome<-ggplot(plot_df, aes(x = biome, y = value, fill = fill)) +
  geom_col(position = "stack", width = 0.6) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.8) +
  scale_fill_manual(values = colors,labels=NULL) +
  labs(x = "Biome",y = "Frequency (%)") +
  scale_y_continuous(limits = c(-0.48, 0.8), breaks=seq(-0.4,0.8,0.2),labels=c("-40","-20","0","20","40","60","80")) +
  theme_bw()+
  theme(axis.text=element_text(size=10,color="black"),axis.title=element_text(size=10,color="black"),
        axis.text.x = element_text(angle = 45, hjust = 1),legend.position = "none",
        panel.grid = element_blank(),plot.margin = margin(t = 5, r = 5, b = 5, l = 5),
        panel.border = element_rect(color="black",size = 1, linetype="solid"))
graph2tif(plot_biome,file="Figure_f",width=3.6,height=3.1,dpi=500,tiffcompression="lzw")


