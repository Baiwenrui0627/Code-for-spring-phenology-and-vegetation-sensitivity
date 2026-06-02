################ Figure 4 ####################
library(raster) 
library(terra)
library(dplyr)
library(tidyr)

library(ggplot2)
library(cowplot)
library(patchwork)
library(export)
library(data.table)
library(ranger)

data_df <- fread("pred_relative_decline_final.csv")
data_df[, i := floor((row - 1) / 10) + 1]
data_df[, k := floor((col - 1) / 10) + 1]
safe_lm <- function(y, x) {
  valid <- complete.cases(y, x)
  y_v <- y[valid]
  x_v <- x[valid]
  if(length(y_v) > 5 && var(x_v) > 0) {
    fit <- lm(y_v ~ x_v)
    sum_fit <- summary(fit)
    slope <- sum_fit$coefficients[2, 1]
    r2 <- sum_fit$r.squared
    f_stat <- sum_fit$fstatistic
    if(!is.null(f_stat)) {
      pval <- pf(f_stat[1], f_stat[2], f_stat[3], lower.tail = FALSE)
    } else {pval <- NA_real_}
    return(list(trend = slope, pval = pval, r2 = r2))
  } else {
    return(list(trend = NA_real_, pval = NA_real_, r2 = NA_real_))
  }
}

res_dt <- data_df[, {
  m1 <- safe_lm(predictions_consider_sos, year)
  m2 <- safe_lm(predictions_not_consider_sos, year)
  m3 <- safe_lm(sos, year)
  .(
    trend_c = m1$trend, pval_c = m1$pval, r2_c = m1$r2,
    trend_notc = m2$trend, pval_notc = m2$pval, r2_notc = m2$r2,
    trend_sos = m3$trend, pval_sos = m3$pval, r2_sos = m3$r2
  )
}, by = .(i, k)]

extent_val <- ext(-180, 180, 23.5, 90)
crs_val <- "EPSG:4326"

save_results <- function(val_vector, prefix) {
  mat <- matrix(NA, nrow = 67, ncol = 360)
  mat[cbind(res_dt$i, res_dt$k)] <- val_vector
  r <- rast(mat, extent = extent_val, crs = crs_val)
  writeRaster(r, paste0(prefix, "_1degree.tif"), overwrite = TRUE)
  return(r)
}

r_trend_c <- save_results(res_dt$trend_c, "RD_c_sos_trend")
r_pval_c  <- save_results(res_dt$pval_c,  "RD_c_sos_pvalue")
r_r2_c    <- save_results(res_dt$r2_c,    "RD_c_sos_R2")

r_trend_notc <- save_results(res_dt$trend_notc, "RD_notc_sos_trend")
r_pval_notc  <- save_results(res_dt$pval_notc,  "RD_notc_sos_pvalue")
r_r2_notc    <- save_results(res_dt$r2_notc,    "RD_notc_sos_R2")

r_trend_sos <- save_results(res_dt$trend_sos, "sos_trend")
r_pval_sos  <- save_results(res_dt$pval_sos,  "sos_pvalue")
r_r2_sos    <- save_results(res_dt$r2_sos,    "sos_R2")


extract_sig <- function(r_pval, out_name) {
  r_sig <- ifel(r_pval < 0.05, 1, NA)
  writeRaster(r_sig, paste0(out_dir, out_name, ".tif"), overwrite = TRUE)
}
extract_sig(r_pval_c,    "RD_c_sos_pvalue_1degree_only_sig")
extract_sig(r_pval_notc, "RD_notc_sos_pvalue_1degree_only_sig")
extract_sig(r_pval_sos,  "sos_pvalue_1degree_only_sig")


sos_trend <- rast(paste0("sos_trend_1degree.tif"))
consider_trend <- rast(paste0("RD_c_sos_trend_1degree.tif"))
not_consider_trend <- rast(paste0("RD_notc_sos_trend_1degree.tif"))
sig_mask <- rast(paste0("sos_pvalue_1degree_only_sig.tif")) 

diff_trend <- not_consider_trend - consider_trend
writeRaster(diff_trend, paste0("trend_diff_1degree.tif"), overwrite = TRUE)
class_rast <- ifel(sos_trend <= 0 & diff_trend <= 0, 1,
                   ifel(sos_trend > 0 & diff_trend > 0, 2,
                        ifel(sos_trend < 0 & diff_trend > 0, 3,
                             ifel(sos_trend > 0 & diff_trend < 0, 4, NA))))
writeRaster(class_rast, paste0("classification_new.tif"), overwrite = TRUE)


sos_trend <- rast(paste0("sos_trend_1degree.tif"))
c_trend   <- rast(paste0("RD_c_sos_trend_1degree.tif"))
notc_trend<- rast(paste0( "RD_notc_sos_trend_1degree.tif"))
sos_sig  <- rast(paste0( "sos_pvalue_1degree_only_sig.tif"))
c_sig    <- rast(paste0( "RD_c_sos_pvalue_1degree_only_sig.tif"))
notc_sig <- rast(paste0( "RD_notc_sos_pvalue_1degree_only_sig.tif"))
diff_trend <- rast(paste0( "trend_diff_1degree.tif"))
area_km2 <- cellSize(sos_trend, unit = "km")
valid_area <- mask(area_km2, sos_trend) 
total_area <- global(valid_area, "sum", na.rm = TRUE)[1, 1]
calc_trend_stats <- function(trend_rast, sig_rast) {
  area_pos     <- global(valid_area * (trend_rast > 0), "sum", na.rm = TRUE)[1, 1]
  area_neg     <- global(valid_area * (trend_rast < 0), "sum", na.rm = TRUE)[1, 1]
  area_sig_pos <- global(valid_area * (trend_rast > 0 & sig_rast == 1), "sum", na.rm = TRUE)[1, 1]
  area_sig_neg <- global(valid_area * (trend_rast < 0 & sig_rast == 1), "sum", na.rm = TRUE)[1, 1]
  return(c(
    pos_ratio     = (area_pos / total_area) * 100,
    neg_ratio     = (area_neg / total_area) * 100,
    sig_pos_ratio = (area_sig_pos / total_area) * 100,
    sig_neg_ratio = (area_sig_neg / total_area) * 100
  ))
}

stats_sos  <- calc_trend_stats(sos_trend, sos_sig)
stats_c    <- calc_trend_stats(c_trend, c_sig)
stats_notc <- calc_trend_stats(notc_trend, notc_sig)

area_diff_pos <- global(valid_area * (diff_trend > 0), "sum", na.rm = TRUE)[1, 1]
area_diff_neg <- global(valid_area * (diff_trend < 0), "sum", na.rm = TRUE)[1, 1]

diff_pos_ratio <- (area_diff_pos / total_area) * 100
diff_neg_ratio <- (area_diff_neg / total_area) * 100


cat("\n================ Global Area-Weighted Statistical Report ================\n")
cat(sprintf("%-20s | %-8s | %-8s | %-12s | %-12s\n", 
            "Data Type", "Trend>0", "Trend<0", "Sig. & >0", "Sig. & <0"))
cat(strrep("-", 72), "\n")
cat(sprintf("%-20s | %6.2f %% | %6.2f %% | %10.2f %% | %10.2f %%\n", 
            "SOS Trend", stats_sos["pos_ratio"], stats_sos["neg_ratio"], 
            stats_sos["sig_pos_ratio"], stats_sos["sig_neg_ratio"]))

cat(sprintf("%-20s | %6.2f %% | %6.2f %% | %10.2f %% | %10.2f %%\n", 
            "Consider SOS Trend", stats_c["pos_ratio"], stats_c["neg_ratio"], 
            stats_c["sig_pos_ratio"], stats_c["sig_neg_ratio"]))

cat(sprintf("%-20s | %6.2f %% | %6.2f %% | %10.2f %% | %10.2f %%\n", 
            "Not Consider SOS", stats_notc["pos_ratio"], stats_notc["neg_ratio"], 
            stats_notc["sig_pos_ratio"], stats_notc["sig_neg_ratio"]))

cat("\n[ Trend Diff Statistics (Not Considering SOS - Considering SOS) ]\n")
cat(sprintf("Diff > 0 (Impact Overestimated) Area Proportion: %6.2f %%\n", diff_pos_ratio))
cat(sprintf("Diff < 0 (Impact Underestimated) Area Proportion: %6.2f %%\n", diff_neg_ratio))


class_rast <- rast(paste0("classification_new.tif"))
area_km2 <- cellSize(class_rast, unit = "km")
valid_area_class <- mask(area_km2, class_rast)
total_area_class <- global(valid_area_class, "sum", na.rm = TRUE)[1, 1]
area_c1 <- global(valid_area_class * (class_rast == 1), "sum", na.rm = TRUE)[1, 1]
area_c2 <- global(valid_area_class * (class_rast == 2), "sum", na.rm = TRUE)[1, 1]
area_c3 <- global(valid_area_class * (class_rast == 3), "sum", na.rm = TRUE)[1, 1]
area_c4 <- global(valid_area_class * (class_rast == 4), "sum", na.rm = TRUE)[1, 1]
r1 <- (area_c1 / total_area_class) * 100
r2 <- (area_c2 / total_area_class) * 100
r3 <- (area_c3 / total_area_class) * 100
r4 <- (area_c4 / total_area_class) * 100
cat("\n================ Four-Category Area-Weighted Statistics =================\n")
cat("[ Area and Proportion by Category ]\n")
cat(sprintf("Category 1 (SOS <= 0 & Diff <= 0) : %12.2f km² | %6.2f %%\n", area_c1, r1))
cat(sprintf("Category 2 (SOS > 0  & Diff > 0)  : %12.2f km² | %6.2f %%\n", area_c2, r2))
cat(sprintf("Category 3 (SOS < 0  & Diff > 0)  : %12.2f km² | %6.2f %%\n", area_c3, r3))
cat(sprintf("Category 4 (SOS > 0  & Diff < 0)  : %12.2f km² | %6.2f %%\n", area_c4, r4))
