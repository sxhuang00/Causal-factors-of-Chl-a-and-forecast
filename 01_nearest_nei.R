####
#### NN interpolation (rEDM version)
####
####


# Load the necessary library
library(rEDM); packageVersion("rEDM") # 0.7.5
library(dplyr)

# ---------------------------------------------------------
# 1. Generate Synthetic Data
# ---------------------------------------------------------
set.seed(8181)
n_samples <- 101

# Generate random longitude and latitude coordinates
lon <- runif(n_samples, min = 130, max = 145)
lat <- runif(n_samples, min = 30, max = 45)

# Generate temperature based on spatial coordinates with some noise
# This simulates a spatial gradient: temperature decreases as latitude increases
temp <- 40 - 0.6 * lat - 0.2 * lon + rnorm(n_samples, mean = 0, sd = 0.5)

# Generate synthetic data
df_all <- data.frame(lon = lon, lat = lat, temp = temp)


# ---------------------------------------------------------
# 2. Perform LOOCV-like Prediction for All Points
# ---------------------------------------------------------
# By setting lib and pred to the same range [1, n_samples], 
# simplex() will systematically exclude the target point 
# from the neighborhood search (effectively LOOCV).

results <- block_lnlp(df_all, method = "simplex",
                      lib = c(1, n_samples), pred = c(1, n_samples),
                      num_neighbors = 1,
                      # If this is 1, then it is simply the "nearest neighbor" method the reviewer mentioned.
                      columns = c("lat", "lon"), target = "temp", tp = 0, 
                      stats_only = FALSE)


# ---------------------------------------------------------
# 3. Combine Results
# ---------------------------------------------------------
# The 'Predictions' column contains the LOOCV estimate for every row
df_all$predicted_temp <- results$model_output[[1]]$pred
head(df_all)

################################################################################
##-----------------ABOVE IS DEMO BELOW WE DO ONE DATE --------------------------
################################################################################
d2 <- read.csv("~/2024_HKredtide/data/2002_2022_modis_wholeHK_tryagain_with_colname.csv", stringsAsFactors = FALSE)
rows_to_remove <- c(3:12, 17, 19, 20, 22:29, 35, 38:45, 53:61, 68:78, 87:94, 102, 
                    106:108, 115:118, 121:124, 130:134,
                    2, 18, 21, 36, 37, 52, 85, 86, 101, 105)

missing_site <- c(2, 7,8, 14, 15, 22, 35, 36, 43, 46)

d2_o <- d2[, -rows_to_remove]
d2_oo <- d2_o[, -missing_site]



head(d2_oo)
class(d2_o)


d2_meta <- readRDS("~/2025hk_chl/data_modis/df_transposed.obj") # 


ter_id <- c(2:11, 16, 18, 19, 21:28, 34, 37:44, 52:60, 67:77, 86:93, 101, 
            105:107, 114:117, 120:123, 129:133)
missing_id <- c(1, 6,7, 13, 14, 21, 34, 35, 42, 45)

d2_meta_o <- d2_meta[-ter_id,]
d2_meta_oo <- d2_meta_o[-missing_id,]


d2_meta_oo$site_code

missing_rate <- colSums(is.na(d2_o[2:102]))/240.
missing <- data.frame(site_code=d2_meta_o$site_code, missing_rate)
missing$Latitude <- d2_meta_o$lat
missing$Longitude<- d2_meta_o$lon
View(missing)
write.csv(missing, file="~/2025hk_chl/uic_sur_output/missing.csv")

new_colname <- as.character(d2_o[1, 1])

new_data <- as.data.frame(t(d2_o[1, 2:102])) 

nn_org <- cbind(d2_meta_o, new_data)

colnames(nn_org)[ncol(nn_org)] <- new_colname



for (i in 2:240) {
  
  new_colname <- as.character(d2_o[i, 1])
  
  new_data <- as.data.frame(t(d2_o[i, 2:102])) 
  
  nn_org <- cbind(nn_org, new_data)
  
  colnames(nn_org)[ncol(nn_org)] <- new_colname

}

View(nn_org)




results <- block_lnlp(nn_org, method = "simplex",
                      lib = c(1, n_samples), pred = c(1, n_samples),
                      num_neighbors = 4,
                      # If this is 1, then it is simply the "nearest neighbor" method the reviewer mentioned.
                      columns = c("lat", "lon"), target = "200208", tp = 0, 
                      stats_only = FALSE)

nn_pred <- data.frame(d2_meta_o$site_code)


tss <- as.character(d2_o$year_month)

for (ts in tss) {
 
  
  results <- block_lnlp(nn_org, 
                        method = "simplex",
                        lib = c(1, n_samples), 
                        pred = c(1, n_samples),
                        num_neighbors = 4,
                        columns = c("lat", "lon"), 
                        target = ts, 
                        tp = 0,
                        stats_only = FALSE)
  
  nn_pred <- cbind(nn_pred, results$model_output[[1]]$pred)
  
}
colnames(nn_pred)[2:241] <- tss

View(nn_pred)

nn_pred <- t(nn_pred)
nn_pred <- nn_pred[-1, ]

colnames(nn_pred)[1:101] <- d2_meta_o$site_code


saveRDS(nn_pred, file = "~/2025hk_chl/data_epd/nn_pred")
nn_pred <- readRDS("~/2025hk_chl/data_epd/nn_pred")



nn_pred[is.na(nn_pred)] <- 0

nn_pred <- apply(nn_pred, c(1, 2), function(x) as.numeric(trimws(x)))
nn_pred <- as.data.frame(nn_pred)

View(nn_pred)
nn_pred_diff <- as.data.frame(lapply(nn_pred, diff)) 

nn_pred_diff$year_month <- d2_diff$year_month


View(nn_pred_diff)
saveRDS(nn_pred_diff, file = "~/2025hk_chl/data_epd/nn_pred_diff")



#############################################################################
#############################################################################
head(nn_pred_diff)

f_smap_modis_nn <- data.frame()
for (i in 1:101) {
  time_series <- data.frame(time = nn_pred_diff$year_month, 
                            value = nn_pred_diff[[i]])
  
  colnames(time_series) <- c("time", "value")
  time_series$value <- scale(time_series$value)
  
  E_test <- rEDM::simplex(time_series = time_series,
                          E = 1:16,
                          lib = c(1, NROW(time_series)), 
                          pred = c(1, NROW(time_series)))
  E_opt <- E_test[which.max(E_test$rho), "E"]
  
  smap <- s_map (time_series = time_series,
                 lib = c(1, NROW(time_series)), 
                 pred = c(1, NROW(time_series)),
                 E = E_opt,
                 theta = c(0, 0.001, 0.01, 0.1, 0.5, 1, 2, 4, 8))
  
  y <- time_series$value  # Extract the observed values
  var_y <- var(y, na.rm = TRUE) 
  smap$nmse <- (smap$rmse)^2 / as.vector(var_y)
  
  max_index <- which.min(smap$nmse)
  
  # Select the row corresponding to the maximum rho value
  smap_best <- smap[max_index, ]
  
  f_smap_modis_nn <- rbind(f_smap_modis_nn, smap_best)
}
View(f_smap_modis_nn)

saveRDS(f_smap_modis_nn, file = "~/2025hk_chl/uic_sur_output/f_smap_modis_nn")

f_smap_modis_lin <- readRDS("~/2025hk_chl/uic_sur_output/f_smap_modis_lin")
f_smap_modis <- readRDS("~/2024_HKredtide/29_F_SMAP_OUT/f_smap_modis")
ff_smap_modis_din <- readRDS("~/2025hk_chl/uic_sur_output/ff_smap_modis_din")


colnames(nn_pred_diff)
missing_id <- c(1, 6,7, 13, 14, 21, 34, 35, 42, 45)


lin_din <- data.frame(
  value = c(f_smap_modis$nmse[-missing_id],
            
            f_smap_modis_lin$nmse[-missing_id],
            ff_smap_modis_din$nmse[-missing_id],
            f_smap_modis_nn$nmse[-missing_id]),
  variable = rep(c("Zero-filling", 
                   "Linear interpolation", "DINEOF", "Nearest neighbor"), 
                 times = c(91,91,91,91
                 )
  ))

lin_din$variable <- factor(lin_din$variable, 
                           levels = c("Zero-filling", 
                                      "Linear interpolation",
                                      "DINEOF", "Nearest neighbor"))




fff1 <- ggplot(lin_din, aes(x = variable, y = value)) +
  geom_boxplot(outlier.shape = NA) +  
  geom_jitter(aes(color = variable), size = 2, alpha = 0.5, width = 0.2) +  # Jittered points
  scale_color_manual(values = c("darkblue", "brown4","seagreen","gold2"))  +
  labs( y = "NMSE", x = "") +
  labs(title = expression("Forecasting performance of different treatments to missing Chl-a values")) +
  ylim(0.4, 1) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),  
    axis.line = element_line(color = "black"), 
    axis.text.x = element_text(size = 14),  # Adjust x-axis title size
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 14), 
    legend.position = "none" 
  )  

ggsave("~/2025hk_chl/uic_sur_output/lindinnn.tif", 
       plot = fff1, width = 10, height = 5, dpi = 300) 


table_s2 <- data.frame(linear=f_smap_modis_lin$nmse,
                       zero_filling=f_smap_modis$nmse,
                       DINEOF=ff_smap_modis_din$nmse,
                       nearest_neighbor=f_smap_modis_nn$nmse,
                       d2_meta_o)
table_s2 <- table_s2[-c(1, 6,7, 13, 14, 21, 34, 35, 42, 45),]
write.csv(table_s2, file="~/2024_HKredtide/29_F_SMAP_OUT/table_s2.csv")










