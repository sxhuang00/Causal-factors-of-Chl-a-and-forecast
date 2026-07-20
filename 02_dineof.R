# 载入数据
d2 <- read.csv("~/2024_HKredtide/data/2002_2022_modis_wholeHK_tryagain_with_colname.csv", stringsAsFactors = FALSE)
rows_to_remove <- c(3:12, 17, 19, 20, 22:29, 35, 38:45, 53:61, 68:78, 87:94, 102, 
                    106:108, 115:118, 121:124, 130:134)

d2_o <- d2[, -rows_to_remove]
head(d2_o)
# 选第2到177列，转矩阵
#data_sub0 <- as.matrix(d2_o[, 2:102])
data_sub <- as.matrix(d2_o[, 2:102])

# 把空字符串替换成NA
data_sub[data_sub == ""] <- NA

# 确保为数值型
storage.mode(data_sub) <- "double"

# 用每列均值初始填充缺失值
col_mean <- colMeans(data_sub, na.rm=TRUE)
inds <- which(is.na(data_sub), arr.ind=TRUE)
for(i in seq_len(nrow(inds))) {
  data_sub[inds[i,1], inds[i,2]] <- col_mean[inds[i,2]]
}

# 设置参数
max_iter <- 1000 # Ushio modified this from 100 to 1000 (but will not change the result dramatically)
tol <- 1e-5
rank <- 10

# 迭代
for(iter in seq_len(max_iter)) {
  svd_res <- svd(data_sub)
  s <- svd_res$d
  U <- svd_res$u
  Vt <- t(svd_res$v)
  
  # 选秩
  rank_i <- rank
  
  # 构建对角矩阵
  S <- diag(s[1:rank_i])
  
  # 重构
  X_reconstructed <- U[, 1:rank_i] %*% S %*% Vt[1:rank_i, ]
  
  # 计算误差，只在缺失点处
  diff <- X_reconstructed[inds] - data_sub[inds]
  error <- sqrt(sum(diff^2)) / sqrt(sum(data_sub[inds]^2))
  
  # 用重构填缺失
  data_sub[inds] <- X_reconstructed[inds]
  
  cat(sprintf("Iteration %d: error=%.8f\n", iter, error))
  
  if(error < tol) {
    cat(sprintf("Converged after %d iterations\n", iter))
    break
  }
  if(iter == max_iter) {
    cat("Reached max iterations without convergence\n")
  }
}

# 替换回数据框
view(data_sub)

# 保存
write.csv(data_sub, "~/2024_HKredtide/data/2002_2022_modis_wholeHK_tryagain_with_colname_dineof.csv", row.names=FALSE)

d2_dineof <- read.csv("~/2024_HKredtide/data/2002_2022_modis_wholeHK_tryagain_with_colname_dineof.csv")
head(d2_dineof)
d2_dineof[d2_dineof < 0] <- 0

d2_dineof_diff <- as.data.frame(lapply(d2_dineof, diff))  # Exclude the first column if it's not numeric (e.g., year_month)
write.csv(d2_dineof_diff, "~/2024_HKredtide/data/2002_2022_modis_wholeHK_tryagain_with_colname_dineof_diff.csv", row.names=FALSE)
d2_dineof_diff <- read.csv("~/2024_HKredtide/data/2002_2022_modis_wholeHK_tryagain_with_colname_dineof_diff.csv")


View(d2_dineof_diff)
d2_dineof_diff <- d2_dineof_diff %>%
  mutate(year_month = d2_diff$year_month) %>%
  select(year_month, everything())
#############################################################################
#############################################################################
head(d2_dineof_diff)

ff_smap_modis_din <- data.frame()
for (i in 2:102) {
  time_series <- data.frame(time = d2_dineof_diff$year_month, 
                            value = d2_dineof_diff[[i]])
  
  
  
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
  
  ff_smap_modis_din <- rbind(ff_smap_modis_din, smap_best)
}
View(ff_smap_modis_din)

missing_id <- c(1, 6,7, 13, 14, 21, 34, 35, 42, 45)
ff_smap_modis_din <- ff_smap_modis_din[-missing_id,]

saveRDS(ff_smap_modis_din, file = "~/2025hk_chl/uic_sur_output/ff_smap_modis_din")
write.csv(ff_smap_modis_din, file = "~/2025hk_chl/uic_sur_output/ff_smap_modis_din.csv")

f_smap_modis_lin <- readRDS("~/2025hk_chl/uic_sur_output/f_smap_modis_lin")
f_smap_modis <- readRDS("~/2024_HKredtide/29_F_SMAP_OUT/f_smap_modis")


lin_din <- data.frame(
  value = c(f_smap_modis$nmse,
            
            f_smap_modis_lin$nmse,
            f_smap_modis_din$nmse),
  variable = rep(c("Zero-filling", 
                   "Linear interpolation", "DINEOF"), 
                 times = c(101,101,101
                          )
))

lin_din$variable <- factor(lin_din$variable, 
                                      levels = c("Zero-filling", 
                                                 "Linear interpolation",
                                                 "DINEOF"))




ff1 <- ggplot(lin_din, aes(x = variable, y = value)) +
  geom_boxplot(outlier.shape = NA) +  
  geom_jitter(aes(color = variable), size = 2, alpha = 0.5, width = 0.2) +  # Jittered points
  scale_color_manual(values = c("darkblue", "brown4","seagreen"))  +
  labs( y = "NMSE", x = "") +
  labs(title = expression("Forecasting performance of different treatment to missing Chl-a value")) +
  ylim(0.4, 1) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),  
    axis.line = element_line(color = "black"), 
    axis.text.x = element_text(size = 14),  # Adjust x-axis title size
    axis.text.y = element_text(size = 14),
    legend.position = "none" 
  )  

ggsave("~/2025hk_chl/uic_sur_output/lindin.tif", 
       plot = ff1, width = 10, height = 5, dpi = 300) 

###################################################################################
###----------here is figure combing original data and interpolation------------
###################################################################################
#northeast=s49, mid=s119, northwest=s46
#southeast=s161, southwest=s176
dates <- seq.Date(from = as.Date("2002-08-01"), to = as.Date("2012-07-31"), by = "month")
mid_int <- data.frame(dates, 
                      #ne_o=d2_o$X49, ne_i=d2_dineof$X49,
                     # nw_o=d2_o$X46, nw_i=d2_dineof$X46,
                      mid_o=d2_o$X119[1:120], mid_i=d2_dineof$X119[1:120])
                     # se_o=d2_o$X161, se_i=d2_dineof$X161)

fig_long <- mid_int %>%
  pivot_longer(cols = -dates, 
               names_to = c("location", "type"), 
               names_sep = "_", 
               values_to = "value") 


fig_long <- fig_long %>% filter(!is.na(value))

g_mid <- ggplot(fig_long, aes(x = dates, y = value, color = type)) +
  geom_line(data = filter(fig_long, type == "i"), size = 1) + 
  geom_point(data = filter(fig_long, type == "i"), size = 1.5) +  
  geom_point(data = filter(fig_long, type == "o"), size = 1.5, na.rm = TRUE) +  
  facet_wrap(~ location, scales = "free_y", ncol = 2) +  
  labs(title = "Middle site",
       x = "Date",
       y = expression("Chl-a Concentration (mg m"^-3*")"),
       color = "Data Type") +
  scale_color_manual(values = c("o" = "blue4", "i" = "red3"), 
                     labels = c("o" = "Original Data", "i" = "Interpolated Data"),
                     name = "Data Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) +
  theme(legend.position = c(0.85, 0.92),
        axis.text.y = element_text(size = 10),   
        legend.text = element_text(size = 12),             
        legend.title = element_text(size = 12),
        plot.title = element_text(size = 16),
    panel.grid = element_blank(),  
    strip.text = element_blank(),
    axis.line = element_line(color = "black"))  


nw_int <- data.frame(dates, 
                     nw_o=d2_o$X33[1:120], nw_i=d2_dineof$X33[1:120])
                      # nw_o=d2_o$X46, nw_i=d2_dineof$X46,
                    #  mid_o=d2_o$X119[1:120], mid_i=d2_dineof$X119[1:120])
# se_o=d2_o$X161, se_i=d2_dineof$X161)

nw_long <- nw_int %>%
  pivot_longer(cols = -dates, 
               names_to = c("location", "type"), 
               names_sep = "_", 
               values_to = "value") 


nw_long <- nw_long %>% filter(!is.na(value))

g_nw <- ggplot(nw_long, aes(x = dates, y = value, color = type)) +
  geom_line(data = filter(nw_long, type == "i"), size = 1) + 
  geom_point(data = filter(nw_long, type == "i"), size = 1.5) +  
  geom_point(data = filter(nw_long, type == "o"), size = 1.5, na.rm = TRUE) +  
  facet_wrap(~ location, scales = "free_y", ncol = 2) +  
  labs(title = "Western site",
       x = "Date",
       y = expression("Chl-a Concentration (mg m"^-3*")"),
       color = "Data Type") +
  scale_color_manual(values = c("o" = "blue4", "i" = "red3"), 
                     labels = c("o" = "Original Data", "i" = "Interpolated Data"),
                     name = "Data Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) +
  theme(legend.position = c(0.85, 0.92),
        axis.text.y = element_text(size = 10),   
        legend.text = element_text(size = 12),             
        legend.title = element_text(size = 12),
        plot.title = element_text(size = 16),
        panel.grid = element_blank(),  
        strip.text = element_blank(),
        axis.line = element_line(color = "black"))  




se_int <- data.frame(dates, 
                   #  nw_o=d2_o$X33, nw_i=d2_dineof$X33)
# nw_o=d2_o$X46, nw_i=d2_dineof$X46,
#  mid_o=d2_o$X119[1:120], mid_i=d2_dineof$X119[1:120])
se_o=d2_o$X45[1:120], se_i=d2_dineof$X45[1:120])

se_long <- se_int %>%
  pivot_longer(cols = -dates, 
               names_to = c("location", "type"), 
               names_sep = "_", 
               values_to = "value") 


se_long <- se_long %>% filter(!is.na(value))

g_se <- ggplot(se_long, aes(x = dates, y = value, color = type)) +
  geom_line(data = filter(se_long, type == "i"), size = 1) + 
  geom_point(data = filter(se_long, type == "i"), size = 1.5) +  
  geom_point(data = filter(se_long, type == "o"), size = 1.5, na.rm = TRUE) +  
  facet_wrap(~ location, scales = "free_y", ncol = 2) +  
  labs(title = "Eastern site",
       x = "Date",
       y = expression("Chl-a Concentration (mg m"^-3*")"),
       color = "Data Type") +
  scale_color_manual(values = c("o" = "blue4", "i" = "red3"), 
                     labels = c("o" = "Original Data", "i" = "Interpolated Data"),
                     name = "Data Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) +
  theme(legend.position = c(0.8, 0.92),
        axis.text.y = element_text(size = 10),   
        legend.text = element_text(size = 12),             
        legend.title = element_text(size = 12),
        plot.title = element_text(size = 16),
        panel.grid = element_blank(),  
        strip.text = element_blank(),
        axis.line = element_line(color = "black"))  


d2_int_rep <- d2_meta[c(33, 45, 119), ]
d2_int_rep$fake <- c(1,1,1)

g_rep <- ggmap(hk_detail) + 
  geom_point(data = d2_int_rep, aes(x = lon, y = lat, size = 2,
                                       color = as.numeric(fake))) +
  scale_color_gradient2(low = "black", mid = "black", high = "black", 
                        limits = c(0.5, 1), midpoint = 0.8) +
  guides(color = "none")  +  #
  labs(title = expression("Representative sites")) +
  theme(legend.position = "none")+
  
  labs(x = "Longitude", y = "Latitude")+
  theme(
        plot.title = element_text(size = 16))

tif_int_rep <- grid.arrange(g_mid,g_nw,g_se, ncol = 1) 
ggsave("2024_HKredtide/29_F_SMAP_OUT/tif_int_rep.tif", 
       plot = tif_int_rep, width = 10, height = 8, dpi = 300) # Adjust dimensions as needed

head(d2)
head(d2_meta)

mean(colSums(is.na(d2_oo[, 2:92]))/240)
d2_meta$missing_rate <- colSums(is.na(d2[, 2:177]))/240

write.csv(d2_meta, file = "2024_HKredtide/29_F_SMAP_OUT/missing.csv" )


