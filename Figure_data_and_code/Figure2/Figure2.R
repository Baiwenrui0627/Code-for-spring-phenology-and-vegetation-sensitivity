##### Figure2 #####
library(dplyr)
library(ggplot2)
library(cowplot)
library(patchwork)
library(export)
library(tidyr)
library(lavaan)
library(data.table)
library(purrr)

positive_sig_df <- fread("positive_sig_df.csv")
data_df <- positive_sig_df  %>% as_tibble() %>% mutate(Biome = as.factor(Biome),lon = -180 + (col - 0.5) * 0.1,lat = 90 - (row - 0.5) * 0.1) 
data_df <- data_df %>% mutate(area_weight_raw = cos(lat * pi / 180),area_weight = area_weight_raw / mean(area_weight_raw, na.rm = TRUE))
data_df <- data_df[which(data_df$relative_decline<100),]

data_df <- data_df %>% mutate(block_row = ceiling(row / 10),block_col = ceiling(col / 10),block_id = paste(block_row, block_col, sep = "_"))
model <- 'LAI_pre ~ sos_zscore + temp_pre + soilm_pre
          Ra_lag ~ LAI_pre
          soilm_diff ~  LAI_pre + prec_lag + SSRD_lag
          NPP_lag ~ LAI_pre
          relative_decline ~  Ra_lag + soilm_diff + NPP_lag'

safe_fit_sem <- function(sub_df) {
  if(nrow(sub_df) < 20) return(NULL)
  fit <- tryCatch({sem(model, data = sub_df, warn = FALSE)}, error = function(e) NULL, warning = function(w) NULL)
  if(is.null(fit)) return(NULL)
  if(!lavInspect(fit, "converged")) return(NULL)
  ests <- standardizedSolution(fit) 
  paths <- ests %>% dplyr::filter(op == "~") %>% dplyr::select(lhs, op, rhs, est.std, pvalue) %>% dplyr::mutate(path_name = paste(lhs, "<-", rhs))
  mean_lat <- mean(sub_df$lat, na.rm = TRUE) 
  return(data.frame(path = paths$path_name,coef = paths$est.std, p_value = paths$pvalue,n_pixels = nrow(sub_df),block_lat = mean_lat))
}

block_results <- data_df %>% group_split(block_id) %>% map_dfr(safe_fit_sem, .progress = TRUE) 
final_sem_summary <- block_results %>%
  mutate(area_weight = cos(block_lat * pi / 180)) %>%
  group_by(path) %>%
  summarise(valid_blocks = n(), 
            weighted_mean_coef = weighted.mean(coef, w = area_weight, na.rm = TRUE),
            weighted_sd_coef = sqrt(sum(area_weight * (coef - weighted_mean_coef)^2, na.rm = TRUE) / sum(area_weight, na.rm = TRUE)),
            se_coef = weighted_sd_coef / sqrt(valid_blocks), 
            t_value = weighted_mean_coef / se_coef,          
            global_p_value = 2 * pt(-abs(t_value), df = valid_blocks - 1), 
            .groups = "drop")
write.csv(final_sem_summary,"GLEAM_sem_result_spatial_average_ALL_pvalue.csv",row.names = F)


bar_data1 <- data.frame(Category = c("Physiological (C-assimilation)", "Physiological (C-consumption)","Biogeophysical"),
                        Value = c(-0.501, 0.475, -0.024),Fill = c("#66C2A5", "#8DA0CB", "#FC8D62"))
bar_plot_negative <- ggplot(bar_data1, aes(x = Category, y = Value)) +
  geom_col(aes(fill = Fill), color = "black", size = 0.4, width = 0.4) +
  scale_fill_identity() + labs(x = NULL,y = "Pathway effect") + theme_bw() +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 15),expand = c(0, 0.4)) +
  scale_y_continuous(limits=c(-0.55,0.5),breaks=c(-0.4,-0.2,0,0.2,0.4),labels=c(-0.4,-0.2,0,0.2,0.4),expand = c(0.01,0.01))+
  geom_hline(yintercept=0,linetype="solid",size=0.6)+
  theme(legend.position = "none",axis.text=element_text(size=9,color="black"),
        panel.grid = element_blank(),axis.title=element_text(size=9,color="black"),
        plot.margin = margin(t = 1, r=8, b = 1, l = 1),
        panel.border = element_rect(color="black",size = 0.6, linetype="solid"),
        axis.line = element_blank(),plot.background = element_blank())
graph2tif(bar_plot_negative,file="Figure2b",width=3.2,height=1.9,dpi=600,tiffcompression="lzw")
