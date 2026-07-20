library(tidyverse); packageVersion("tidyverse") # 2.0.0
library(rEDM); packageVersion("rEDM") # 1.15.3
library(ggmap)
library(rUIC)
library(pbapply)
library(zoo)
library(ggsci) 
library(ggplot2)
#remotes::install_github("yutakaos/rUIC")
#remotes::install_github("ong8181/macam", build_vignettes = TRUE, force = TRUE)
library(macam); packageVersion("macam")
vignette("macam")


# Install the legacy version of rEDM
remotes::install_github("ha0ye/rEDM")
library(rEDM); packageVersion("rEDM") # 0.7.5





##------------------------------------------------------------------------
##first we do defalut uic using original chl ts
##----------------------------------------------------------------------------
d2 <- read_csv("~/2025hk_chl/data_modis/2002_2022_modis_wholeHK_tryagain_with_colname.csv")
d2$year_month <- ymd(paste0(d2$year_month, "01"))
d2[is.na(d2)] <- 0
head(d2) ###modis 
tail(d2)

d_SST <- read_csv("~/2025hk_chl/data_modis/MODIS_SST.csv")

head(d_SST)








d2_diff <- as.data.frame(lapply(d2[, -1], diff))  # Exclude the first column if it's not numeric (e.g., year_month)

# Create a new year_month column for the differenced data
year_month_diff <- d2$year_month[-1]  # Adjust year_month to match the new length

# Combine year_month with the difference data
d2_diff <- data.frame(year_month = year_month_diff, d2_diff)
colnames(d2_diff) <- c("year_month", names(d2)[-1]) 
head(d2_diff)
tail(d2_diff)
colnames(d2_diff)

rows_to_remove <- c(3:12, 17, 19, 20, 22:29, 35, 38:45, 53:61, 68:78, 87:94, 102, 
                    106:108, 115:118, 121:124, 130:134)

d2_diff_o <- d2_diff[, -rows_to_remove]
head(d2_diff_o)

d_SST_o <- d_SST[, -rows_to_remove]
head(d_SST_o)


##########################################################################
##-------------------------MDR MODIS-------------------------------------------------------------------
############################################################################
MDR_analysis <- function(i) {
  set.seed(123)
  # Create block data frame for the current i
  block <- data.frame(chl = d2_diff_o[, i] , 
                      SST = d_SST_o[2:240, i]%>% pull
                     )
                      
  
  block <- apply(block, 2, function(x) as.numeric(scale(x)))
  block <- as.data.frame(block)
  
  effect_var <- "chl"
  y <- block[[effect_var]]  # Extract the observed values
  var_y <- var(y) 
  
  results_list_i <- list()
  
  ##step 1
  # Calculate variance of the observed values
  
  E_range <- 0:16
  simp_res <- rUIC::simplex(block, lib_var = effect_var, E = E_range, 
                            tau = 1, tp = 1)
  Ex <- simp_res[which.min(simp_res$rmse),"E"]
  
  
  tp_range = -4:0
  silent = FALSE
  # Step 2: Perform UIC to detect causality
  uic_res <- uic_across(block, effect_var, E_range = E_range, 
                        tp_range = tp_range, silent = silent)
  
  
  # Step 3: Make block to calculate multiview distance
  block_mvd <- macam::make_block_mvd (block, uic_res, effect_var,
                                     E_effect_var = Ex, 
                                     include_var = "all_significant")
  #block_mvd <- uic_res %>% filter(pval <= 0.05)
  #block_mvd <- block_mvd %>% filter(E == Ex)  
  
  
  tp = 1
  n_ssr = 10000
  k = 100
  random_seed = 1234
  # Step. 4: Compute multiview distance
  multiview_dist <- compute_mvd(block_mvd, effect_var, E = Ex ,
                                make_block_method = "rEDM",
                                make_block_max_lag = Ex,
                                n_ssr = n_ssr, k = k,
                                tp = tp, random_seed = random_seed)
  
  regularized = TRUE
  theta_values <- c(0, 0.001, 0.01, 0.1, 0.5, 1, 2, 4, 8)
  lambda_values <- c(0, 0.001, 0.01, 0.1, 0.5, 1, 2)
  alpha = 0
  glmnet_parallel = TRUE
  save_smap_coefficients = TRUE
  lib = c(1, nrow(block_mvd))
  # Loop over each theta value
  for (theta in theta_values) {
    for (lambda in lambda_values) {
      # Perform the MDR S-map with the current theta and lambda values
      mdr_res <- s_map_mdr(block_mvd,
                           dist_w = multiview_dist,
                           lib = lib,
                           pred = lib,
                           tp = tp,
                           theta = theta,
                           regularized = regularized,
                           lambda = lambda,
                           alpha = alpha,
                           glmnet_parallel = glmnet_parallel,
                           save_smap_coefficients = save_smap_coefficients,
                           random_seed = random_seed)
      
      # Extract the statistics and add the current theta and lambda values
      mdr_stats <- mdr_res$stats
      mdr_stats$theta <- theta  # Add the theta value to the statistics
      mdr_stats$lambda <- lambda  # Add the lambda value to the statistics
      
      # Calculate NMSE based on RMSE
      rmse <- mdr_stats$rmse  # Assuming 'rmse' is a column in mdr_stats
      nmse <- (rmse^2) / var_y  # Calculate NMSE
      mdr_stats$NMSE <- nmse  # Add NMSE to the statistics
      
      # Store the result in the list
      results_list_i[[paste(theta, lambda, sep = "_")]] <- mdr_stats
    }
  }
  
  # Combine all results for the current i into a single data frame
  final_mdr_i <- do.call(rbind, results_list_i)
  
  # Find the row with the minimum NMSE for the current i
  min_nmse_index <- which.min(final_mdr_i$NMSE)
  return(final_mdr_i[min_nmse_index, ])  # Return the best result
}

# Use pblapply to run the analysis in parallel
F_MDR_MODIS_tp <- pblapply(2:102, MDR_analysis, cl = 72)

# Combine all best results into a single data frame
F_MDR_MODIS_tp <- do.call(rbind, F_MDR_MODIS_tp)





saveRDS(F_MDR_MODIS_tp, file = "~/2025hk_chl/uic_sur_output/F_MDR_MODIS_tp")
F_MDR_MODIS_tp <- readRDS("~/2025hk_chl/uic_sur_output/F_MDR_MODIS_tp")

d2_meta <- readRDS("~/2025hk_chl/data_modis/df_transposed.obj") # 

ter_id <- c(2:11, 16, 18, 19, 21:28, 34, 37:44, 52:60, 67:77, 86:93, 101, 105:107, 114:117, 120:123, 129:133)
d2_meta_o <- d2_meta[-ter_id,]



F_MDR_MODIS_tp$site_code <- d2_meta_o$site_code


F_MDR_MODIS_tp <- F_MDR_MODIS_tp %>%
  filter(!is.na(as.numeric(NMSE)) )

view(F_MDR_MODIS_tp)


f_smap_modis <- readRDS("~/2025hk_chl/uic_sur_output/f_smap_modis")
ff_smap_modis <- data.frame(site_code = d2_meta_o$site_code, 
                            F_SMAP = f_smap_modis$nmse)
f_smap_modis_o <- ff_smap_modis %>%
  filter(site_code %in% F_MDR_MODIS_tp$site_code)


mdr_step2 <- data.frame(f_smap_modis_o, F_MDR_MODIS_tp)
view(mdr_step2)


##now check mdr works on which sites
mdr_step2$mdr_better <- as.numeric(mdr_step2$F_SMAP) - as.numeric(mdr_step2$NMSE)


##>0 => mdr works


##now check whther ete of SST higher, mdr nmse lower
d2_meta_mdr <- d2_meta %>%
  filter(site_code %in% mdr_step2$site_code)

mdr_step2 <- data.frame(mdr_step2, d2_meta_mdr)
view(mdr_step2)



ggmap(hk_detail) + 
  geom_point(data = mdr_step2, aes(x = lon, y = lat, color = mdr_better)) +
  scale_color_gradient2(mid = "darkblue", high = "yellow3",
                        name = "mdr_better")



ggmap(hk_detail) + 
  geom_point(data = mdr_step2, aes(x = lon, y = lat, color = mdr_better)) +
  scale_color_gradientn(colors = c("darkblue", "yellow"), 
                        values = scales::rescale(c(-1, 0, 1)), 
                        name = "NMSE_difference") +
  theme_minimal()  # Optional: Add a theme for better aesthetics









mean(f_smap_modis_o$F_SMAP)
mean(as.numeric(F_MDR_MODIS_tp$NMSE))

boxplot(
  f_smap_modis_o$F_SMAP, 
  as.numeric(F_MDR_MODIS_tp$NMSE),
  names = c("SMAP","MDR"),
  #main = "Summary of causal effect",
  ylab = "NMSE",
  col = c("lightblue", "coral3")
)



mdr_modis_sina <- data.frame(
  value = c(f_smap_modis_o$F_SMAP, as.numeric(F_MDR_MODIS_tp$NMSE)),
  variable = rep(c("Chl-a", "Chl-a + SST"), 
                 times = c(length(f_smap_modis_o$F_SMAP),  
                           length(as.numeric(F_MDR_MODIS_tp$NMSE))))
)

mdr_modis_sina$variable <- factor(mdr_modis_sina$variable, levels = c("Chl-a", 
                                                                      "Chl-a + SST"))



# Create the sina plot
g26 <- ggplot(mdr_modis_sina, aes(x = variable, y = value)) +
  geom_boxplot(outlier.shape = NA) +  # Boxplot for the summary
  geom_jitter(aes(color = variable), size = 2, alpha = 0.5, width = 0.2) +  # Jittered points
  scale_color_manual(values = c("blue4", "brown"))  +
  labs(title = "(a)", y = "NMSE", x = "") +
  ylim(0.4, 1) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),  # Remove grid lines
    axis.line = element_line(color = "black"),  # Keep axis lines
    legend.position = "none" 
  )  # K

data <- data.frame(
  SMAP = f_smap_modis_o$F_SMAP,
  MDR = as.numeric(F_MDR_MODIS_tp$NMSE))









g25 <- ggplot(data, aes(x = SMAP, y = MDR)) +
  geom_point(shape = 1, size = 3) +  # Add hollow circles
  geom_abline(slope = 1, intercept = 0) +  # Add the reference line
  xlim(0.4, 1) +  # Set x-axis limits
  ylim(0.4, 1) +  # Set y-axis limits
  labs(
    x = "NMSE of model Chl-a ",  # Label for x-axis
    y = "NMSE of model Chl-a + SST ",  # Label for y-axis
    title = expression("(b) " )  # Title of the plot
  ) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.background = element_blank()) +  # Optional: Remove panel background
  theme(axis.line = element_line())  # Keep axis lines



word_mdrmodis <- grid.arrange(g26, g25, ncol = 2) # Change ncol to 1 for vertical arrangement

# Save the combined figure
ggsave("~/2025hk_chl/uic_sur_output//word_mdrmodis.jpg", 
       plot = word_mdrmodis, width = 10, height = 5, dpi = 300) # Adjust dimensions as needed






ETE_SST$site_code <- colnames(d2_o)[-1]
ETE_SST_o <- ETE_SST%>%
  filter(site_code %in% F_MDR_MODIS_o$site_code)

ETE_SST_o$comp_modis <- f_smap_modis_o$nmse/as.numeric(F_MDR_MODIS_o$NMSE)
view(ETE_SST_o)


####################################################################################
########################EPD---------------------------------------------------------
####################################################################################
d_chl <- read_csv("~/2025hk_chl/data_epd/EPD_WholeHK.csv")
# Delete the first row
d_chl <- d_chl[-1, ]

# Delete the last five rows to keep time series is 200208-202207
d_chl <- d_chl[1:(nrow(d_chl) - 5), ]
head(d_chl)
tail(d_chl)
d_chl[is.na(d_chl)] <- 0

d_chl_diff <- as.data.frame(lapply(d_chl[, -1], diff))  # Exclude the first column if it's not numeric (e.g., year_month)

# Create a new year_month column for the differenced data
year_month_diff <- d_chl$year_month[-1]  # Adjust year_month to match the new length

# Combine year_month with the difference data
d_chl_diff <- data.frame(year_month = year_month_diff, d_chl_diff)
colnames(d_chl_diff) <- c("year_month", names(d_chl)[-1]) 
head(d_chl_diff)
tail(d_chl_diff)
colnames(d_chl_diff)



d_pH <- read_csv("~/2025hk_chl/data_epd/EPD_pH.csv")
d_sal <- read_csv("~/2025hk_chl/data_epd/EPD_salinity.csv")
d_temp <- read_csv("~/2025hk_chl/data_epd/EPD_temp.csv")
d_TN <- read_csv("~/2025hk_chl/data_epd/EPD_TN.csv")
d_TP <- read_csv("~/2025hk_chl/data_epd/EPD_TP.csv")

d_tur <- read_csv("~/2025hk_chl/data_epd/EPD_tur.csv")
d_sil <- read_csv("~/2025hk_chl/data_epd/EPD_sil.csv")
d_NP <- data.frame(d_TN$year_month, d_TN[2:77] / d_TP[2:77])


for (i in 2:77) {
  d_NP[[i]][!is.finite(d_NP[[i]])] <- NA
}

MDR_EPD <- function(i) {
  set.seed(1234)
  # Create block data frame for the current i
  block <- data.frame(d_chl_diff[, i]
                 # temp = d_temp[2:240, i],
                #    sal = d_sal[2:240, i],
                 #    pH= d_pH[2:240, i])
                    np = d_NP[, i])
                    # TP= d_TP[2:240, i])
                    #  sil = d_sil[2:240, i])
                    #  tur = d_tur[2:240, i])
  
  block <- apply(block, 2, function(x) as.numeric(scale(x)))
  block <- as.data.frame(block)
  
  effect_var <- "chl"
  y <- block[[effect_var]]  # Extract the observed values
  var_y <- var(y) 
  
  results_list_i <- list()
  
  ##step 1
  # Calculate variance of the observed values
  
  E_range <- 0:16
  simp_res <- rUIC::simplex(block, lib_var = effect_var, E = E_range, 
                            tau = 1, tp = 1)
  Ex <- simp_res[which.min(simp_res$rmse),"E"] 
  
  
  tp_range = -4:0
  silent = FALSE
  # Step 2: Perform UIC to detect causality
  uic_res <- uic_across(block, effect_var, E_range = E_range, 
                        tp_range = tp_range, silent = silent)
  
  
  # Step 3: Make block to calculate multiview distance
  block_mvd <- macam::make_block_mvd(block, uic_res, effect_var, E_effect_var = Ex ,
                            include_var = "all_significant")
  
  
  
  tp = 1
  n_ssr = 10000
  k = 100
  random_seed = 1234
  # Step. 4: Compute multiview distance
  
  multiview_dist <- compute_mvd(block_mvd, effect_var, E = Ex , #Ex
                                make_block_method = "rEDM",
                                make_block_max_lag = 3, #Ex
                                n_ssr = n_ssr, k = k,
                                tp = tp, random_seed = random_seed,
                                distance_only = TRUE)
  
  regularized = TRUE
  theta_values <- c(0, 0.001, 0.01, 0.1, 0.5, 1, 2, 4, 8)
  lambda_values <- c(0, 0.001, 0.01, 0.1, 0.5, 1, 2)
  alpha = 0
  glmnet_parallel = TRUE
  save_smap_coefficients = TRUE
  lib = c(1, nrow(block_mvd))
  # Loop over each theta value
  for (theta in theta_values) {
    for (lambda in lambda_values) {
      # Perform the MDR S-map with the current theta and lambda values
      mdr_res <- s_map_mdr(block_mvd,
                           dist_w = multiview_dist,
                           lib = lib,
                           pred = lib,
                           tp = tp,
                           theta = theta,
                           regularized = regularized,
                           lambda = lambda,
                           alpha = alpha,
                           glmnet_parallel = glmnet_parallel,
                           save_smap_coefficients = save_smap_coefficients,
                           random_seed = random_seed)
      
      # Extract the statistics and add the current theta and lambda values
      mdr_stats <- mdr_res$stats
      mdr_stats$theta <- theta  # Add the theta value to the statistics
      mdr_stats$lambda <- lambda  # Add the lambda value to the statistics
      
      # Calculate NMSE based on RMSE
      rmse <- mdr_stats$rmse  # Assuming 'rmse' is a column in mdr_stats
      nmse <- (rmse^2) / var_y  # Calculate NMSE
      mdr_stats$NMSE <- nmse  # Add NMSE to the statistics
      
      # Store the result in the list
      results_list_i[[paste(theta, lambda, sep = "_")]] <- mdr_stats
    }
  }
  
  # Combine all results for the current i into a single data frame
  final_mdr_i <- do.call(rbind, results_list_i)
  
  # Find the row with the minimum NMSE for the current i
  min_nmse_index <- which.min(final_mdr_i$NMSE)
  return(final_mdr_i[min_nmse_index, ])  # Return the best result
}

# Use pblapply to run the analysis in parallel
f_mdr_epd_tp <- pblapply(2:77, MDR_EPD, cl = 90)

f_mdr_epd_tp_np <- do.call(rbind, f_mdr_epd_tp)
  view(f_mdr_epd_tp_np)

saveRDS(f_mdr_epd_tp_np_nodiff, file = "~/2025hk_chl/uic_sur_output/f_mdr_epd_tp_np")

f_mdr_epd_tp_sal <- readRDS("~/2025hk_chl/uic_sur_output/f_mdr_epd_tp_sal")
f_mdr_epd_tp_temp <- readRDS("~/2025hk_chl/uic_sur_output/f_mdr_epd_tp_temp")
f_mdr_epd_tp_np <- readRDS("~/2025hk_chl/uic_sur_output/f_mdr_epd_tp_np")
f_mdr_epd_tp_sil <- readRDS("~/2025hk_chl/uic_sur_output/f_mdr_epd_tp_sil")
f_mdr_epd_tp_tur <- readRDS("~/2025hk_chl/uic_sur_output/f_mdr_epd_tp_tur")
f_mdr_epd_tp_ph <- readRDS("~/2025hk_chl/uic_sur_output/f_mdr_epd_tp_ph")
f_mdr_epd_nei1 <- readRDS("~/2025hk_chl/uic_sur_output/f_mdr_epd_nei1")



f_mdr_epd_tp_sal <- data.frame(f_mdr_epd_tp_sal, 
                               site_code = colnames(d_chl_diff[2:77]))

f_mdr_epd_tp_temp <- data.frame(f_mdr_epd_tp_temp, 
                               site_code = colnames(d_chl_diff[2:77]))

f_mdr_epd_tp_np <- data.frame(f_mdr_epd_tp_np, 
                               site_code = colnames(d_chl_diff[2:77]))

f_mdr_epd_tp_sil <- data.frame(f_mdr_epd_tp_sil, 
                               site_code = colnames(d_chl_diff[2:77]))

f_mdr_epd_tp_tur <- data.frame(f_mdr_epd_tp_tur, 
                               site_code = colnames(d_chl_diff[2:77]))

f_mdr_epd_tp_ph <- data.frame(f_mdr_epd_tp_ph, 
                               site_code = colnames(d_chl_diff[2:77]))

f_mdr_epd_tp_np_nodiff <- data.frame(f_mdr_epd_tp_np_nodiff, 
                              site_code = colnames(d_chl_diff[2:77]))


f_mdr_epd_tp_sal <- f_mdr_epd_tp_sal %>%
  filter(!is.na(as.numeric(NMSE)) )

f_mdr_epd_tp_temp <- f_mdr_epd_tp_temp %>%
  filter(!is.na(as.numeric(NMSE)) )

f_mdr_epd_tp_np <- f_mdr_epd_tp_np %>%
  filter(!is.na(as.numeric(NMSE)) )



f_mdr_epd_tp_ph <- f_mdr_epd_tp_ph %>%
  filter(!is.na(as.numeric(NMSE)) )

f_mdr_epd_tp_sil <- f_mdr_epd_tp_sil %>%
  filter(!is.na(as.numeric(NMSE)) )


f_mdr_epd_tp_tur <- f_mdr_epd_tp_tur %>%
  filter(!is.na(as.numeric(NMSE)) )

f_mdr_epd_tp_np_nodiff <- f_mdr_epd_tp_np_nodiff %>%
  filter(!is.na(as.numeric(NMSE)) )


view(f_mdr_epd_tp_np_nodiff)

f_smap_epd <- readRDS("~/2025hk_chl/uic_sur_output//f_smap_epd")
ff_smap_epd <- data.frame(site_code = colnames(d_chl)[2:77], 
                          F_SMAP = f_smap_epd$nmse)
f_smap_epd_tp_tempsal <- ff_smap_epd %>%
  filter(site_code %in% f_mdr_epd_tp_tempsal$site_code)

mean(ff_smap_epd$F_SMAP[-14]) #0.6601727
mean(as.numeric(f_mdr_epd_tp_sal$NMSE)) #0.665965
mean(as.numeric(f_mdr_epd_tp_temp$NMSE)) #0.662852
mean(as.numeric(f_mdr_epd_tp_np$NMSE)) #0.707885
mean(as.numeric(f_mdr_epd_tp_ph$NMSE)) #0.0.7082127
mean(as.numeric(f_mdr_epd_tp_sil$NMSE)) #0.0.6902192
mean(as.numeric(f_mdr_epd_tp_tur$NMSE)) #0.656304


##NOW WE SEE WHICH SITES WORKS USING MDR
comparingEPD <- data.frame(site_code = colnames(d_chl)[2:77], 
                           F_SMAP = f_smap_epd$nmse,
                           sal= f_mdr_epd_tp_sal$NMSE,
                           temp= f_mdr_epd_tp_temp$NMSE,
                           ph= f_mdr_epd_tp_ph$NMSE,
                           sil= f_mdr_epd_tp_sil$NMSE,
                           tur= f_mdr_epd_tp_tur$NMSE,
                           np= f_mdr_epd_tp_np$NMSE)



comparingEPD <- comparingEPD %>%
  mutate(across(-site_code, ~ as.numeric(as.character(.))))

comparingEPD$comp_sal <- ifelse(is.na(comparingEPD$F_SMAP) | is.na(comparingEPD$sal), 
                                NA, comparingEPD$F_SMAP - comparingEPD$sal)

comparingEPD$comp_temp <- ifelse(is.na(comparingEPD$F_SMAP) | is.na(comparingEPD$temp), 
                                NA, comparingEPD$F_SMAP - comparingEPD$temp)


comparingEPD$comp_ph <- ifelse(is.na(comparingEPD$F_SMAP) | is.na(comparingEPD$ph), 
                                 NA, comparingEPD$F_SMAP - comparingEPD$ph)

comparingEPD$comp_sil <- ifelse(is.na(comparingEPD$F_SMAP) | is.na(comparingEPD$sil), 
                               NA, comparingEPD$F_SMAP - comparingEPD$sil)

comparingEPD$comp_tur <- ifelse(is.na(comparingEPD$F_SMAP) | is.na(comparingEPD$tur), 
                                NA, comparingEPD$F_SMAP - comparingEPD$tur)


comparingEPD$comp_np <- ifelse(is.na(comparingEPD$F_SMAP) | is.na(comparingEPD$np), 
                                NA, comparingEPD$F_SMAP - comparingEPD$np)


comparingEPD$nei1 <- f_mdr_epd_nei1$NMSE
comparingEPD$comp_nei1 <- ifelse(is.na(comparingEPD$F_SMAP) | is.na(comparingEPD$nei1), 
                               NA, comparingEPD$F_SMAP - comparingEPD$nei1)




view(comparingEPD)

mdr_epd_sina <- data.frame(
  value = c(f_smap_epd$nmse[-14], 
            
            as.numeric(f_mdr_epd_tp_temp$NMSE),
            as.numeric(f_mdr_epd_tp_sal$NMSE),
            as.numeric(f_mdr_epd_tp_ph$NMSE),
            as.numeric(f_mdr_epd_tp_tur$NMSE),
            as.numeric(f_mdr_epd_tp_np$NMSE),
            as.numeric(f_mdr_epd_tp_sil$NMSE)),
  
  variable = rep(c("Chl-a",  "Chl-a + Temp","Chl-a + Sal",
                   "Chl-a + pH","Chl-a + Tur","Chl-a + N/P", "Chl-a + Sil"), 
                 times = c(length(f_smap_epd$nmse[-14]),  
                           
                           length(as.numeric(f_mdr_epd_tp_temp$NMSE)),
                           length(as.numeric(f_mdr_epd_tp_sal$NMSE)),
                           length(as.numeric(f_mdr_epd_tp_ph$NMSE)),
                           length(as.numeric(f_mdr_epd_tp_tur$NMSE)),
                           length(as.numeric(f_mdr_epd_tp_np$NMSE)),
                           length(as.numeric(f_mdr_epd_tp_sil$NMSE)))
                           ))


mdr_epd_sina$variable <- factor(mdr_epd_sina$variable, 
                                levels = c("Chl-a",  "Chl-a + Temp","Chl-a + Sal",
                                           "Chl-a + pH","Chl-a + Tur",
                                           "Chl-a + N/P", "Chl-a + Sil"))



# Create the sina plot
g24 <- ggplot(mdr_epd_sina, aes(x = variable, y = value)) +
  geom_boxplot(outlier.shape = NA) +  # Boxplot for the summary
  geom_jitter(aes(color = variable), size = 2, alpha = 0.5, width = 0.2) +  # Jittered points
  scale_color_manual(values = c( "brown", "blue4","olivedrab","gray",
                                 "brown", "blue4", "olivedrab"))  +
  labs(title = "", y = "NMSE", x = "") +
  ylim(0.35, 1) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),  # Remove grid lines
    axis.line = element_line(color = "black"),  # Keep axis lines
    legend.position = "none" 
  )


ggsave("~/2025hk_chl/uic_sur_output//word_mdrall.jpg", 
       plot = g24, width = 10, height = 5, dpi = 300) # Adjust dimensions as needed





data <- data.frame(
  SMAP = f_smap_epd_tp_tempsal$F_SMAP,
  MDR_temp = as.numeric(f_mdr_epd_tp_tempsal$NMSE)
)

mean(f_smap_epd_tp_tempsal$F_SMAP) #0.6668392
mean(as.numeric(f_mdr_epd_tp_tempsal$NMSE)) # 0.6563261

g23 <- ggplot(data, aes(x = SMAP, y = MDR_temp)) +
  geom_point(shape = 1, size = 3) +  # Add hollow circles
  geom_abline(slope = 1, intercept = 0) +  # Add the reference line
  xlim(0.4, 1) +  # Set x-axis limits
  ylim(0.4, 1) +  # Set y-axis limits
  labs(
    x = "NMSE of S-map",  # Label for x-axis
    y = "NMSE of MDR S-map",  # Label for y-axis
    title = expression("(a) " )  # Title of the plot
  ) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.background = element_blank()) +  # Optional: Remove panel background
  theme(axis.line = element_line())  # Keep axis lines

mdr_epd_sina <- data.frame(
  value = c(f_smap_epd_tp_temp$F_SMAP, 
           
            as.numeric(f_mdr_epd_tp_tempsal$NMSE)),
  
  variable = rep(c("Chl", 
                   "Chl+Temp+Sal"), 
                 times = c(length(f_smap_epd_tp_temp$F_SMAP),  
                          
                           length(as.numeric(f_mdr_epd_tp_tempsal$NMSE)))
  ))


mdr_epd_sina$variable <- factor(mdr_epd_sina$variable, 
                                levels = c("Chl", 
                                           "Chl+Temp+Sal"))



# Create the sina plot
g25 <- ggplot(mdr_epd_sina, aes(x = variable, y = value)) +
  geom_boxplot(outlier.shape = NA) +  # Boxplot for the summary
  geom_jitter(aes(color = variable), size = 2, alpha = 0.5, width = 0.2) +  # Jittered points
  scale_color_manual(values = c("blue4", "brown"))  +
  labs(title = "(b)", y = "NMSE", x = "") +
  ylim(0.4, 1) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),  # Remove grid lines
    axis.line = element_line(color = "black"),  # Keep axis lines
    legend.position = "none" 
  )


word_mdrepdtempsal <- grid.arrange(g23,g25, ncol = 2) # Change ncol to 1 for vertical arrangement



ETE_DO$site_code <- colnames(d_chl)[2:77]
ETE_PH$site_code <- colnames(d_chl)[2:77]
ETE_SAL$site_code <- colnames(d_chl)[2:77]
ETE_sil$site_code <- colnames(d_chl)[2:77]
ETE_temp$site_code <- colnames(d_chl)[2:77]
ETE_TN$site_code <- colnames(d_chl)[2:77]
ETE_TP$site_code <- colnames(d_chl)[2:77]
ETE_tur$site_code <- colnames(d_chl)[2:77]

ETE_DO_oo <- ETE_DO %>%
  filter(site_code %in% f_mdr_epd$site_code)

ETE_PH_oo <- ETE_PH %>%
  filter(site_code %in% f_mdr_epd$site_code)
ETE_SAL_oo <- ETE_SAL %>%
  filter(site_code %in% f_mdr_epd$site_code)
ETE_sil_oo <- ETE_sil %>%
  filter(site_code %in% f_mdr_epd$site_code)
ETE_temp_oo <- ETE_temp %>%
  filter(site_code %in% f_mdr_epd$site_code)
ETE_TN_oo <- ETE_TN %>%
  filter(site_code %in% f_mdr_epd$site_code)
ETE_TP_oo <- ETE_TP %>%
  filter(site_code %in% f_mdr_epd$site_code)
ETE_tur_oo <- ETE_tur %>%
  filter(site_code %in% f_mdr_epd$site_code)

ETE_oo <- data.frame(ETE_DO_oo$ete,
                     ETE_PH_oo$ete,
                     ETE_SAL_oo$ete,
                     ETE_sil_oo$ete,
                     ETE_temp_oo$ete,
                     ETE_TN_oo$ete,
                     ETE_TP_oo$ete,
                     ETE_tur_oo$ete,
                     comp_epd)
view(ETE_oo)
