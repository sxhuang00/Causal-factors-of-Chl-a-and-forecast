library(tidyverse); packageVersion("tidyverse") # 2.0.0
library(rEDM); packageVersion("rEDM") # 1.15.3
library(ggmap)
library(rUIC)
library(pbapply)
library(zoo)
library(ggsci) 
library(ggplot2)
library(gridExtra)

##------------------------------------------------------------------------
##first we do defalut uic using original chl ts
##----------------------------------------------------------------------------
d2 <- read_csv("2025hk_chl/data_modis/2002_2022_modis_wholeHK_tryagain_with_colname.csv")

view(d2)
##now see how many NA for each site
na_counts <- colSums(is.na(d2_o[2:102]))
na_counts
mean(na_counts)


d2$year_month <- ymd(paste0(d2$year_month, "01"))
d2[is.na(d2)] <- 0
head(d2) ###modis 
tail(d2)

d_SST <- read_csv("2025hk_chl/data_modis/MODIS_SST.csv")

head(d_SST)


columns_to_remove <- c(3:12, 17, 19, 20, 22:29, 35, 38:45, 53:61, 68:78, 87:94, 102,
                       106:108, 115:118, 121:124, 130:134)
SST_o <- d_SST[, -columns_to_remove]


d2_o <- d2[, -columns_to_remove]


# Loop through the specified range
UIC_SST <- data.frame()  # Initialize an empty data frame to store results

uic_SST <- list()

for (i in 2:102) {
  set.seed(123)
  
  # Step 3: Create block_MODIS data frame
  block_MODIS <- data.frame(d2_o_i = d2_o[[i]], SST_o_i = SST_o[[i]])
  block_MODIS <- apply(block_MODIS, 2, function(x) as.numeric(scale(x)))
  block_MODIS <- as.data.frame(block_MODIS)
  # Step 4: Conduct UIC analysis
  uic_SST[[i]] <- uic.optimal(block_MODIS, lib_var = "d2_o_i", tar_var = "SST_o_i",
                              E=0:16, tau=1, tp=-4:0)
  
  # Step 5: Select the row with the minimum pval
  min_pval_row <- uic_SST[[i]][which.min(uic_SST[[i]]$pval), ]
  
  # Step 7: Store results in UM_SST
  UIC_SST <- rbind(UIC_SST, min_pval_row)
}

# Output the results data frame

uic_SST
uic_SST[[2]]


saveRDS(uic_SST, file = "~/2025hk_chl/uic_sur_output/uic_SST")
uic_SST <- readRDS("~/2025hk_chl/uic_sur_output/uic_SST")
###-----------------------------------------------------------------------------
##here is surrogate
##----------------------------------------------------------------------------------
process_ts <- function(i) {
  set.seed(123)
  ts <- SST_o[[i]]  # Access the column directly
  ts <- as.numeric(ts)
  
  # Replace NA with the mean of the non-NA values
  ts_filled_mean <- na.fill(ts, fill = mean(ts, na.rm = TRUE))
  
  # Generate surrogate data
  surrogate_results <- make_surrogate_data(ts_filled_mean, method = "seasonal", 
                                           num_surr = 1000, T_period = 12)
  
  # Initialize a list to store te values for different tp
  uic_op <- list()
  count_values <- list()
  te_values <- list() 
  sum_counts <- numeric(5)  # Initialize a vector to store the sums for each k
  
  # Loop through each column of surrogate_results
  for (j in 1:1000) {
    
    block_MODIS <- data.frame(d2 = d2_o[[i]], SST_sur = surrogate_results[, j])
    block_MODIS <- apply(block_MODIS, 2, function(x) as.numeric(scale(x)))
    block_MODIS <- as.data.frame(block_MODIS)
    
    
    
    uic_op[[j]] <- uic.optimal(block_MODIS, lib_var = "d2", tar_var = "SST_sur",
                               E = 0:16, tau = 1, tp = -4:0, num_surr = 1)
    
    # Initialize a list to store te values for the current j
    te_values[[j]] <- numeric(5)  # Preallocate space for 5 values
    
    # Store te values for each tp
    for (k in 1:5) {
      te_values[[j]][k] <- uic_op[[j]]$te[k]  # Store the specific te value for each tp
    }
    
    # Initialize a vector to store counts for the current j
    count_values[[j]] <- numeric(5)  # Preallocate space for 5 counts
    
    # Count how many of the 2 values are less than uic_results[[i]]$te[k]
    for (k in 1:5) {
      count_values[[j]][k] <- sum(te_values[[j]][k] < uic_SST[[i]]$te[k])  # Count values less than
      # Update the sum_counts for the current k
      sum_counts[k] <- sum_counts[k] + count_values[[j]][k]  # Add to the sum for k
    }
  }
  
  # Calculate the ratios_tp based on sum_counts for the current i
  ratios_tp <- 1 - (sum_counts / 1000)  # Calculate the ratios for each k
  return(ratios_tp)  # Return the result for the current i
}

# Use pblapply to process each i in SSTallel
F_SUR_SST <- pblapply(2:102, process_ts, cl = 64)
F_SUR_SST
saveRDS(F_SUR_SST, file = "~/2025hk_chl/uic_sur_output/F_SUR_SST")
F_SUR_SST <- readRDS("~/2025hk_chl/uic_sur_output//F_SUR_SST")
F_MIN_SUR_SST <- sapply(F_SUR_SST, min)
F_MIN_SUR_SST


ff_MIN_SUR_SST <- data.frame(site_code = colnames(d2_o)[2:102], 
                             F_MIN_SUR_SST = F_MIN_SUR_SST)

d2_meta <- readRDS("2025hk_chl/data_modis/df_transposed.obj")
d2_meta
ff_MIN_SUR_SST <- data.frame(ff_MIN_SUR_SST, d2_meta[-c(2:11, 16, 18, 19, 21:28, 
                                                      34, 37:44, 52:60, 67:77,
                                                      86:93, 101, 105:107, 114:117,
                                                      120:123, 129:133), ])
################3------------------------------------------------------------
##-----------------------------NOW CHECK ETE-----------------------------------
########################################-------------------------------------
##since we're going to use ete based on pval based on SUR
##we're making an index based on SUR pval
F_SUR_SST <- readRDS("~/2025hk_chl/uic_sur_output/F_SUR_SST")
uic_SST <- readRDS("~/2025hk_chl/uic_sur_output//uic_SST")
pval_indices_sst <- sapply(F_SUR_SST, function(x) {
  min_value <- min(x)  # Find the minimum value
  last_index <- max(which(x == min_value))  # Find the last occurrence of the minimum value
  return(last_index)  # Return the index of the last occurrence
})
pval_indices_sst

##and then we PICK the ete row from uic_TN on above index
# Initialize an empty list to store the selected rows
selected_rows <- list()

# Loop through each index from 2 to 102
for (i in 1:101) {
  # Get the row index from pval_indices
  j <- pval_indices_sst[[i]]  # Adjusting for 1-based indexing in R
  
  # Extract the specified row from the dataframe
  selected_row <- uic_SST[[i + 1]][j, ]
  
  # Append the selected row to the list
  selected_rows[[i]] <- selected_row  # Store in the list, adjusting index for 1-based indexing
}


# Combine all selected rows into a new dataframe
ETE_SST <- do.call(rbind, selected_rows)
view(ETE_SST)
ETE_SST$tp
ETE_SST$ete
ff_MIN_SUR_SST

ff_MIN_SUR_SST$ete <- ETE_SST$ete
ff_MIN_SUR_SST$sig <- ifelse(ff_MIN_SUR_SST$F_MIN_SUR_SST > 0.05, 
                                       0.05, ff_MIN_SUR_SST$F_MIN_SUR_SST)

saveRDS(ff_MIN_SUR_SST, file = "~/2025hk_chl/uic_sur_output/ff_MIN_SUR_SST")



# Map type names
map_types <- c("stamen_terrain", "stamen_toner", "stamen_toner_lite",
               "stamen_watercolor", "alidade_smooth", "alidade_smooth_dark",
               "outdoors", "stamen_terrain_background", "stamen_toner_background",
               "stamen_terrain_labels", "stamen_terrain_lines",
               "stamen_toner_labels", "stamen_toner_lines")

# Draw Hong Kong Map
hk_lonlat <- c(left = 113.8, right = 114.5, top = 22.6, bottom = 22.1)
hk_detail <- get_stadiamap(hk_lonlat, zoom = 10, maptype = map_types[1])


g10 <-ggmap(hk_detail) + 
  geom_point(data = ff_MIN_SUR_SST, 
             aes(x = lon, y = lat, color = ete, size = sig, 
                 shape = ifelse(F_MIN_SUR_SST > 0.05, "Hollow", "Filled"))) +  # Use shape based on significance
  scale_color_gradient2(mid = "darkblue", high = "yellow3", name = "Effect Strength") +
  scale_size(range = c(5, 2), name = "Significance") +  # Reverse the size range
  scale_shape_manual(values = c("Hollow" = 1, "Filled" = 16), labels = NULL) +  # Remove shape labels
  guides(size = guide_legend(order = 1),
         shape = "none") +  # Remove the shape legend+  # Control legend order
  labs(title = expression("(b) Remote sensing data")) +
  labs(x = "Longitude", y = "Latitude") + # Control the order of legends
  theme(
    plot.title = element_text(size = 10),  # Adjust title size
    axis.title.x = element_text(size = 8),  # Adjust x-axis title size
    axis.title.y = element_text(size = 8),  # Adjust y-axis title size
    legend.text = element_text(size = 6),    # Adjust legend text size
    legend.title = element_text(size = 6)  # Adjust legend text size
  )


    # Control the order of legends

############################################






############################################################################
##------------------------here is surrogate FOR EPD-------------------------------
##################################################################################
d_DO <- read_csv("~/2025hk_chl/data_epd/EPD_DO.csv")
d_pH <- read_csv("~/2025hk_chl/data_epd/EPD_pH.csv")
d_sal <- read_csv("~/2025hk_chl/data_epd/EPD_salinity.csv")
d_temp <- read_csv("~/2025hk_chl/data_epd/EPD_temp.csv")
d_TN <- read_csv("~/2025hk_chl/data_epd/EPD_TN.csv")
d_TP <- read_csv("~/2025hk_chl/data_epd/EPD_TP.csv")

d_tur <- read_csv("~/2025hk_chl/data_epd/EPD_tur.csv")
d_sil <- read_csv("~/2025hk_chl/data_epd/EPD_sil.csv")



d_chl <- read_csv("~/2025hk_chl/data_epd/EPD_WholeHK.csv")

d_chl <- d_chl[-1, ]

# Delete the last five rows to keep time series is 200208-202207
d_chl <- d_chl[1:(nrow(d_chl) - 5), ]
head(d_chl)
tail(d_chl)
d_chl[is.na(d_chl)] <- 0
############################################################################
##------------------------TEMP-------------------------------
##################################################################################
###-----------------------------------------------------------------------------
##here is surrogate FOR EPD
##----------------------------------------------------------------------------------
UIC_TEMP <- data.frame()  # Initialize an empty data frame to store results
uic_TEMP <- list()
for (i in 2:77) {
  set.seed(123)
  
  # Step 3: Create block_MODIS data frame
  block_MODIS <- data.frame(d_chl_i = d_chl[[i]], d_temp_i = d_temp[[i]])
  block_MODIS <- apply(block_MODIS, 2, function(x) as.numeric(scale(x)))
  block_MODIS <- as.data.frame(block_MODIS)
  
  # Step 4: Conduct UIC analysis
  uic_TEMP[[i]] <- uic.optimal(block_MODIS, lib_var = "d_chl_i", tar_var = "d_temp_i",
                               E=0:16, tau=1, tp=-4:0)
  
  # Step 5: Select the row with the minimum pval
  min_pval_row <- uic_TEMP[[i]][which.min(uic_TEMP[[i]]$pval), ]
  
  # Step 7: Store results in UM_SST
  UIC_TEMP <- rbind(UIC_TEMP, min_pval_row)
}

uic_TEMP[[2]]

saveRDS(uic_TEMP, file = "~/2025hk_chl/uic_sur_output/uic_TEMP")


process_ts <- function(i) {
  set.seed(123)
  ts <- d_temp[[i]]  # Access the column directly
  ts <- as.numeric(ts)
  
  # Replace NA with the mean of the non-NA values
  ts_filled_mean <- na.fill(ts, fill = mean(ts, na.rm = TRUE))
  
  # Generate surrogate data
  surrogate_results <- make_surrogate_data(ts_filled_mean, method = "seasonal", 
                                     num_surr = 1000, T_period = 12)
  
  # Initialize a list to store te values for different temp
  count_values <- list()
  te_values <- list() 
  sum_counts <- numeric(5)  # Initialize a vector to store the sums for each k
  
  # Loop through each column of surrogate_results
  for (j in 1:1000) {
    
    block_MODIS <- data.frame(d_chl = d_chl[[i]], temp_sur = surrogate_results[, j])
   
    
    block_MODIS <- apply(block_MODIS, 2, function(x) as.numeric(scale(x)))
    block_MODIS <- as.data.frame(block_MODIS)
    
    
    
     uic_op[[j]] <- uic.optimal(block_MODIS, lib_var = "d_chl", tar_var = "temp_sur",
                               E = 0:16, tau = 1, tp = -4:0, num_surr = 1)
    
    # Initialize a list to store te values for the current j
    te_values[[j]] <- numeric(5)  # Preallocate space for 5 values
    
    # Store te values for each temp
    for (k in 1:5) {
      te_values[[j]][k] <- uic_op[[j]]$te[k]  # Store the specific te value for each temp
    }
    
    # Initialize a vector to store counts for the current j
    count_values[[j]] <- numeric(5)  # Preallocate space for 5 counts
    
    # Count how many of the 2 values are less than uic_results[[i]]$te[k]
    for (k in 1:5) {
      count_values[[j]][k] <- sum(te_values[[j]][k] < uic_TEMP[[i]]$te[k])  # Count values less than
      # Update the sum_counts for the current k
      sum_counts[k] <- sum_counts[k] + count_values[[j]][k]  # Add to the sum for k
    }
  }
  
  # Calculate the ratios_temp based on sum_counts for the current i
  ratios_temp <- 1 - (sum_counts / 1000)  # Calculate the ratios for each k
  return(ratios_temp)  # Return the result for the current i
}

# Use pblapply to process each i in parallel
F_SUR_temp <- pblapply(2:77, process_ts, cl = 64)
F_SUR_temp

saveRDS(F_SUR_temp, file = "~/2025hk_chl/uic_sur_output/F_SUR_temp")
F_SUR_temp <- readRDS("~/2025hk_chl/uic_sur_output/F_SUR_temp")

MIN_SUR_temp <- sapply(F_SUR_temp, min)
MIN_SUR_temp <- data.frame(site_code = colnames(d_chl)[2:77], 
                           pval = unlist(MIN_SUR_temp))

d_EPD_site <- readRDS("~/2025hk_chl/data_epd/d_EPD_site.obj")
MIN_SUR_temp <- merge(MIN_SUR_temp, d_EPD_site, by = "site_code")

################3------------------------------------------------------------
##-----------------------------NOW CHECK ETE-----------------------------------
########################################-------------------------------------
##since we're going to use ete based on pval based on SUR
##we're making an index based on SUR pval
F_SUR_temp

uic_TEMP <- readRDS("~/2025hk_chl/uic_sur_output/uic_TEMP")

pval_indices_temp <- sapply(F_SUR_temp, function(x) {
  min_value <- min(x)  # Find the minimum value
  last_index <- max(which(x == min_value))  # Find the last occurrence of the minimum value
  return(last_index)  # Return the index of the last occurrence
})
pval_indices_temp

##and then we tempK the ete row from uic_TN on above index
# Initialize an empty list to store the selected rows
selected_rows <- list()

# Loop through each index from 2 to 102
for (i in 1:76) {
  # Get the row index from pval_indices
  j <- pval_indices_temp[[i]]  # Adjusting for 1-based indexing in R
  
  # Extract the specified row from the dataframe
  selected_row <- uic_TEMP[[i + 1]][j, ]
  
  # Append the selected row to the list
  selected_rows[[i]] <- selected_row  # Store in the list, adjusting index for 1-based indexing
}


# Combine all selected rows into a new dataframe
ETE_temp <- do.call(rbind, selected_rows)
ETE_temp$tp
ETE_temp$ete
MIN_SUR_temp

MIN_SUR_temp$ete <- ETE_temp$ete

MIN_SUR_temp$sig <- ifelse(MIN_SUR_temp$pval > 0.05, 
                             0.05, MIN_SUR_temp$pval)

saveRDS(MIN_SUR_temp, file = "~/2025hk_chl/uic_sur_output/MIN_SUR_temp")


ggmap(hk_detail) + 
  geom_point(data = MIN_SUR_temp, 
             aes(x = lon, y = lat, color = ete, size = pval)) +  # Add size aesthetic
  scale_color_gradient2( mid = "darkblue", high = "yellow3",
                         name = "temp_ETE") +
  scale_size(range = c(2, 5), name = "Significance") 

g11 <- ggmap(hk_detail) + 
  geom_point(data = MIN_SUR_temp, 
             aes(x = lon, y = lat, color = ete, size = sig, 
                 shape = ifelse(pval > 0.05, "Hollow", "Filled"))) +  # Use shape based on significance
  scale_color_gradient2(mid = "darkblue", high = "yellow3", name = "Effect Strength") +
  scale_size(range = c(5, 2), name = "Significance") +  # Reverse the size range
  scale_shape_manual(values = c("Hollow" = 1, "Filled" = 16), labels = NULL) +  # Remove shape labels
  guides(size = guide_legend(order = 2),
       #  color = guide_legend(order = 1), 
         shape = "none") +  #
  labs(title = expression("(a) " * italic("In situ") * " measurement data")) +
  labs(x = "Longitude", y = "Latitude") + # Control the order of legends
  theme(
    plot.title = element_text(size = 10),  # Adjust title size
    axis.title.x = element_text(size = 8),  # Adjust x-axis title size
    axis.title.y = element_text(size = 8),  # Adjust y-axis title size
    legend.text = element_text(size = 6),    # Adjust legend text size
    legend.title = element_text(size = 6)  # Adjust legend text size
  )

# Control the order of legends

word_sst <- grid.arrange(g11, g10, ncol = 2) # Change ncol to 1 for vertical arrangement

# Save the combined figure
ggsave("~/2025hk_chl/uic_sur_output/word_sst.jpg", 
       plot = word_sst, width = 10, height = 5, dpi = 300) # Adjust dimensions as needed






############################################################################
##------------------------SAL-------------------------------
##################################################################################
###-----------------------------------------------------------------------------
##here is surrogate FOR EPD
##----------------------------------------------------------------------------------
UIC_SAL <- data.frame()  # Initialize an empty data frame to store results
uic_SAL <- list()
for (i in 2:77) {
  set.seed(123)
  
  # Step 3: Create block_MODIS data frame
  block_MODIS <- data.frame(d_chl_i = d_chl[[i]], d_sal_i = d_sal[[i]])
  block_MODIS <- apply(block_MODIS, 2, function(x) as.numeric(scale(x)))
  block_MODIS <- as.data.frame(block_MODIS)
  
  # Step 4: Conduct UIC analysis
  uic_SAL[[i]] <- uic.optimal(block_MODIS, lib_var = "d_chl_i", tar_var = "d_sal_i",
                               E=0:16, tau=1, tp=-4:0)
  
  # Step 5: Select the row with the minimum pval
  min_pval_row <- uic_SAL[[i]][which.min(uic_SAL[[i]]$pval), ]
  
  # Step 7: Store results in UM_SST
  UIC_SAL <- rbind(UIC_SAL, min_pval_row)
}

uic_SAL[[2]]

saveRDS(uic_SAL, file = "~/2025hk_chl/uic_sur_output/uic_SAL")
uic_SAL <- readRDS("~/2025hk_chl/uic_sur_output/uic_SAL")

process_ts <- function(i) {
  set.seed(123)
  ts <- d_sal[[i]]  # Access the column directly
  ts <- as.numeric(ts)
  
  # Replace NA with the mean of the non-NA values
  ts_filled_mean <- na.fill(ts, fill = mean(ts, na.rm = TRUE))
  
  # Generate surrogate data
  surrogate_results <- make_surrogate_data(ts_filled_mean, method = "seasonal", 
                                           num_surr = 1000, T_period = 12)
  
  # Initialize a list to store te values for different SAL
  count_values <- list()
  te_values <- list() 
  sum_counts <- numeric(5)  # Initialize a vector to store the sums for each k
  
  # Loop through each column of surrogate_results
  for (j in 1:1000) {
    
    block_MODIS <- data.frame(d_chl = d_chl[[i]], SAL_sur = surrogate_results[, j])
    
    
    block_MODIS <- apply(block_MODIS, 2, function(x) as.numeric(scale(x)))
    block_MODIS <- as.data.frame(block_MODIS)
    
    
    
    uic_op[[j]] <- uic.optimal(block_MODIS, lib_var = "d_chl", tar_var = "SAL_sur",
                               E = 0:16, tau = 1, tp = -4:0, num_surr = 1)
    
    # Initialize a list to store te values for the current j
    te_values[[j]] <- numeric(5)  # Preallocate space for 5 values
    
    # Store te values for each SAL
    for (k in 1:5) {
      te_values[[j]][k] <- uic_op[[j]]$te[k]  # Store the specific te value for each SAL
    }
    
    # Initialize a vector to store counts for the current j
    count_values[[j]] <- numeric(5)  # Preallocate space for 5 counts
    
    # Count how many of the 2 values are less than uic_results[[i]]$te[k]
    for (k in 1:5) {
      count_values[[j]][k] <- sum(te_values[[j]][k] < uic_SAL[[i]]$te[k])  # Count values less than
      # Update the sum_counts for the current k
      sum_counts[k] <- sum_counts[k] + count_values[[j]][k]  # Add to the sum for k
    }
  }
  
  # Calculate the ratios_SAL based on sum_counts for the current i
  ratios_SAL <- 1 - (sum_counts / 1000)  # Calculate the ratios for each k
  return(ratios_SAL)  # Return the result for the current i
}

# Use pblapply to process each i in parallel
F_SUR_SAL <- pblapply(2:77, process_ts, cl = 64)
F_SUR_SAL

saveRDS(F_SUR_SAL, file = "~/2025hk_chl/uic_sur_output/F_SUR_SAL")
F_SUR_SAL <- readRDS("~/2025hk_chl/uic_sur_output/F_SUR_SAL")

MIN_SUR_SAL <- sapply(F_SUR_SAL, min)
MIN_SUR_SAL <- data.frame(site_code = colnames(d_chl)[2:77], 
                           pval = unlist(MIN_SUR_SAL))


MIN_SUR_SAL <- merge(MIN_SUR_SAL, d_EPD_site, by = "site_code")

################3------------------------------------------------------------
##-----------------------------NOW CHECK ETE-----------------------------------
########################################-------------------------------------
##since we're going to use ete based on pval based on SUR
##we're making an index based on SUR pval
F_SUR_SAL

uic_SAL <- readRDS("~/2025hk_chl/uic_sur_output/uic_SAL")

pval_indices_SAL <- sapply(F_SUR_SAL, function(x) {
  min_value <- min(x)  # Find the minimum value
  last_index <- max(which(x == min_value))  # Find the last occurrence of the minimum value
  return(last_index)  # Return the index of the last occurrence
})
pval_indices_SAL

##and then we SALK the ete row from uic_TN on above index
# Initialize an empty list to store the selected rows
selected_rows <- list()

# Loop through each index from 2 to 102
for (i in 1:76) {
  # Get the row index from pval_indices
  j <- pval_indices_SAL[[i]]  # Adjusting for 1-based indexing in R
  
  # Extract the specified row from the dataframe
  selected_row <- uic_SAL[[i + 1]][j, ]
  
  # Append the selected row to the list
  selected_rows[[i]] <- selected_row  # Store in the list, adjusting index for 1-based indexing
}


# Combine all selected rows into a new dataframe
ETE_SAL <- do.call(rbind, selected_rows)
ETE_SAL$tp
ETE_SAL$ete
MIN_SUR_SAL

MIN_SUR_SAL$ete <- ETE_SAL$ete
saveRDS(MIN_SUR_SAL, file = "~/2025hk_chl/uic_sur_output/MIN_SUR_SAL")
MIN_SUR_SAL$sig <- ifelse(MIN_SUR_SAL$pval > 0.05, 
                           0.05, MIN_SUR_SAL$pval)





g12 <-ggmap(hk_detail) + 
  geom_point(data = MIN_SUR_SAL, 
             aes(x = lon, y = lat, color = ete, size = sig, 
                 shape = ifelse(pval > 0.05, "Hollow", "Filled"))) +  # Use shape based on significance
  scale_color_gradient2(mid = "darkblue", high = "yellow3", name = "Effect Strength") +
  scale_size(range = c(5, 2), name = "Significance") +  # Reverse the size range
  scale_shape_manual(values = c("Hollow" = 1, "Filled" = 16), labels = NULL) +  # Remove shape labels
  guides(size = guide_legend(order = 2),
         #  color = guide_legend(order = 1), 
         shape = "none") +  #
  labs(title = expression("(b) Salinity")) +
  labs(x = "Longitude", y = "Latitude") + # Control the order of legends
  theme(
    plot.title = element_text(size = 10),  # Adjust title size
    axis.title.x = element_text(size = 8),  # Adjust x-axis title size
    axis.title.y = element_text(size = 8),  # Adjust y-axis title size
    legend.text = element_text(size = 6),    # Adjust legend text size
    legend.title = element_text(size = 6)  # Adjust legend text size
  )



############################################################################
##------------------------PH-------------------------------
##################################################################################
###-----------------------------------------------------------------------------
##here is surrogate FOR EPD
##----------------------------------------------------------------------------------
UIC_PH <- data.frame()  # Initialize an empty data frame to store results
uic_PH <- list()
for (i in 2:77) {
  set.seed(123)
  
  # Step 3: Create block_MODIS data frame
  block_MODIS <- data.frame(d_chl_i = d_chl[[i]], d_pH_i = d_pH[[i]])
  block_MODIS <- apply(block_MODIS, 2, function(x) as.numeric(scale(x)))
  block_MODIS <- as.data.frame(block_MODIS)
  
  # Step 4: Conduct UIC analysis
  uic_PH[[i]] <- uic.optimal(block_MODIS, lib_var = "d_chl_i", tar_var = "d_pH_i",
                              E=0:16, tau=1, tp=-4:0)
  
  # Step 5: Select the row with the minimum pval
  min_pval_row <- uic_PH[[i]][which.min(uic_PH[[i]]$pval), ]
  
  # Step 7: Store results in UM_SST
  UIC_PH <- rbind(UIC_PH, min_pval_row)
}

uic_PH[[2]]

saveRDS(uic_PH, file = "~/2025hk_chl/uic_sur_output/uic_PH")


process_ts <- function(i) {
  set.seed(123)
  ts <- d_pH[[i]]  # Access the column directly
  ts <- as.numeric(ts)
  
  # Replace NA with the mean of the non-NA values
  ts_filled_mean <- na.fill(ts, fill = mean(ts, na.rm = TRUE))
  
  # Generate surrogate data
  surrogate_results <- make_surrogate_data(ts_filled_mean, method = "seasonal", 
                                           num_surr = 1000, T_period = 12)
  
  # Initialize a list to store te values for different PH
  count_values <- list()
  te_values <- list() 
  sum_counts <- numeric(5)  # Initialize a vector to store the sums for each k
  
  # Loop through each column of surrogate_results
  for (j in 1:1000) {
    
    block_MODIS <- data.frame(d_chl = d_chl[[i]], PH_sur = surrogate_results[, j])
    
    
    block_MODIS <- apply(block_MODIS, 2, function(x) as.numeric(scale(x)))
    block_MODIS <- as.data.frame(block_MODIS)
    
    
    
    uic_op[[j]] <- uic.optimal(block_MODIS, lib_var = "d_chl", tar_var = "PH_sur",
                               E = 0:16, tau = 1, tp = -4:0, num_surr = 1)
    
    # Initialize a list to store te values for the current j
    te_values[[j]] <- numeric(5)  # Preallocate space for 5 values
    
    # Store te values for each PH
    for (k in 1:5) {
      te_values[[j]][k] <- uic_op[[j]]$te[k]  # Store the specific te value for each PH
    }
    
    # Initialize a vector to store counts for the current j
    count_values[[j]] <- numeric(5)  # Preallocate space for 5 counts
    
    # Count how many of the 2 values are less than uic_results[[i]]$te[k]
    for (k in 1:5) {
      count_values[[j]][k] <- sum(te_values[[j]][k] < uic_PH[[i]]$te[k])  # Count values less than
      # Update the sum_counts for the current k
      sum_counts[k] <- sum_counts[k] + count_values[[j]][k]  # Add to the sum for k
    }
  }
  
  # Calculate the ratios_PH based on sum_counts for the current i
  ratios_PH <- 1 - (sum_counts / 1000)  # Calculate the ratios for each k
  return(ratios_PH)  # Return the result for the current i
}

# Use pblapply to process each i in parallel
F_SUR_PH <- pblapply(2:77, process_ts, cl = 64)
F_SUR_PH

saveRDS(F_SUR_PH, file = "~/2025hk_chl/uic_sur_output/F_SUR_PH")
F_SUR_PH <- readRDS("~/2025hk_chl/uic_sur_output/F_SUR_PH")

MIN_SUR_PH <- sapply(F_SUR_PH, min)
MIN_SUR_PH <- data.frame(site_code = colnames(d_chl)[2:77], 
                          pval = unlist(MIN_SUR_PH))


MIN_SUR_PH <- merge(MIN_SUR_PH, d_EPD_site, by = "site_code")

################3------------------------------------------------------------
##-----------------------------NOW CHECK ETE-----------------------------------
########################################-------------------------------------
##since we're going to use ete based on pval based on SUR
##we're making an index based on SUR pval
F_SUR_PH

uic_PH <- readRDS("~/2025hk_chl/uic_sur_output/uic_PH")

pval_indices_PH <- sapply(F_SUR_PH, function(x) {
  min_value <- min(x)  # Find the minimum value
  last_index <- max(which(x == min_value))  # Find the last occurrence of the minimum value
  return(last_index)  # Return the index of the last occurrence
})
pval_indices_PH

##and then we PHK the ete row from uic_TN on above index
# Initialize an empty list to store the selected rows
selected_rows <- list()

# Loop through each index from 2 to 102
for (i in 1:76) {
  # Get the row index from pval_indices
  j <- pval_indices_PH[[i]]  # Adjusting for 1-based indexing in R
  
  # Extract the specified row from the dataframe
  selected_row <- uic_PH[[i + 1]][j, ]
  
  # Append the selected row to the list
  selected_rows[[i]] <- selected_row  # Store in the list, adjusting index for 1-based indexing
}


# Combine all selected rows into a new dataframe
ETE_PH <- do.call(rbind, selected_rows)
ETE_PH$tp
ETE_PH$ete
MIN_SUR_PH

MIN_SUR_PH$ete <- ETE_PH$ete

MIN_SUR_PH$sig <- ifelse(MIN_SUR_PH$pval > 0.05, 
                          0.05, MIN_SUR_PH$pval)
saveRDS(MIN_SUR_PH, file="~/2025hk_chl/uic_sur_output/MIN_SUR_PH")

g13 <-ggmap(hk_detail) + 
  geom_point(data = MIN_SUR_PH, 
             aes(x = lon, y = lat, color = ete, size = sig, 
                 shape = ifelse(pval > 0.05, "Hollow", "Filled"))) +  # Use shape based on significance
  scale_color_gradient2(mid = "darkblue", high = "yellow3", name = "Effect Strength") +
  scale_size(range = c(5, 2), name = "Significance") +  # Reverse the size range
  scale_shape_manual(values = c("Hollow" = 1, "Filled" = 16), labels = NULL) +  # Remove shape labels
  guides(size = guide_legend(order = 2),
         #  color = guide_legend(order = 1), 
         shape = "none") +  #
  labs(title = expression("(c) pH")) +
  labs(x = "Longitude", y = "Latitude") + # Control the order of legends
  theme(
    plot.title = element_text(size = 10),  # Adjust title size
    axis.title.x = element_text(size = 8),  # Adjust x-axis title size
    axis.title.y = element_text(size = 8),  # Adjust y-axis title size
    legend.text = element_text(size = 6),    # Adjust legend text size
    legend.title = element_text(size = 6)  # Adjust legend text size
  )






###########################################################################
##############-------------------SINA----------------------------------###############################
##################################################################



MIN_SUR_temp_sig <- MIN_SUR_temp %>%
  filter(pval <= 0.05)

MIN_SUR_SAL_sig <- MIN_SUR_SAL %>%
  filter(pval <= 0.05)
MIN_SUR_PH_sig <- MIN_SUR_PH %>%
  filter(pval <= 0.05)
MIN_SUR_tur_sig <- MIN_SUR_tur %>%
  filter(pval <= 0.05)
MIN_SUR_NP_sig <- MIN_SUR_NP %>%
  filter(pval <= 0.05)
MIN_SUR_sil_sig <- MIN_SUR_sil %>%
  filter(pval <= 0.05)

mean(MIN_SUR_NP_sig[4:12,]$ete)






chap22_sina <- data.frame(
  value = c(MIN_SUR_temp_sig$ete,
           
            MIN_SUR_SAL_sig$ete,
            MIN_SUR_PH_sig$ete, 
            MIN_SUR_tur_sig$ete,
            
            MIN_SUR_NP_sig$ete,
            MIN_SUR_sil_sig$ete),
  
  variable = rep(c(  "Temperature", "Salinity", "pH", "Turbidity", "N/P ratio", "Silica"), 
                 times = c(
                           length(MIN_SUR_temp_sig$ete),  
                           
                           length(MIN_SUR_SAL_sig$ete),  
                           length(MIN_SUR_PH_sig$ete),
                           length(MIN_SUR_tur_sig$ete),  
                         
                           length(MIN_SUR_NP_sig$ete),  
                           length(MIN_SUR_sil_sig$ete)))
)

chap22_sina$variable <- factor(chap22_sina$variable,
                               levels = c( "Temperature", "Salinity", "pH", "Turbidity", 
                                          "N/P ratio", "Silica"))



# Create the sina plot
g18 <- ggplot(chap22_sina, aes(x = variable, y = value)) +
  geom_boxplot(outlier.shape = NA) +  # Boxplot for the summary
  geom_jitter(aes(color = variable), size = 2, alpha = 0.5, width = 0.2) +  # Jittered points
  scale_color_manual(values = c( "brown","blue4", "gray", 
                                 "brown","blue4", "gray")) +
  labs(title = "(a)", y = "Causal effect strength", x = "Environmental factors") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),  # Remove grid lines
    axis.line = element_line(color = "black"),  # Keep axis lines
    legend.position = "none"  # Remove legend
  )



############################################################################
##------------------------tur-------------------------------
##################################################################################
###-----------------------------------------------------------------------------
##here is surrogate FOR EPD
##----------------------------------------------------------------------------------
str(d_tur)
d_tur <- d_tur[1:240, ]
UIC_tur <- data.frame()  # Initialize an empty data frame to store results
uic_tur <- list()
for (i in 2:77) {
  set.seed(123)
  
  # Step 3: Create block_MODIS data frame
  block_MODIS <- data.frame(d_chl_i = d_chl[[i]], d_tur_i = d_tur[[i]])
  block_MODIS <- apply(block_MODIS, 2, function(x) as.numeric(scale(x)))
  block_MODIS <- as.data.frame(block_MODIS)
  
  # Step 4: Conduct UIC analysis
  uic_tur[[i]] <- uic.optimal(block_MODIS, lib_var = "d_chl_i", tar_var = "d_tur_i",
                             E=0:16, tau=1, tp=-4:0)
  
  # Step 5: Select the row with the minimum pval
  min_pval_row <- uic_tur[[i]][which.min(uic_tur[[i]]$pval), ]
  
  # Step 7: Store results in UM_SST
  UIC_tur <- rbind(UIC_tur, min_pval_row)
}

uic_tur[[77]]

saveRDS(uic_tur, file = "~/2025hk_chl/uic_sur_output/uic_tur")


process_ts <- function(i) {
  set.seed(123)
  ts <- d_tur[[i]]  # Access the column directly
  ts <- as.numeric(ts)
  
  # Replace NA with the mean of the non-NA values
  ts_filled_mean <- na.fill(ts, fill = mean(ts, na.rm = TRUE))
  
  # Generate surrogate data
  surrogate_results <- make_surrogate_data(ts_filled_mean, method = "seasonal", 
                                           num_surr = 1000, T_period = 12)
  
  # Initialize a list to store te values for different tur
  count_values <- list()
  te_values <- list() 
  sum_counts <- numeric(5)  # Initialize a vector to store the sums for each k
  
  # Loop through each column of surrogate_results
  for (j in 1:1000) {
    
    block_MODIS <- data.frame(d_chl = d_chl[[i]], tur_sur = surrogate_results[, j])
    
    
    block_MODIS <- apply(block_MODIS, 2, function(x) as.numeric(scale(x)))
    block_MODIS <- as.data.frame(block_MODIS)
    
    
    
    uic_op[[j]] <- uic.optimal(block_MODIS, lib_var = "d_chl", tar_var = "tur_sur",
                               E = 0:16, tau = 1, tp = -4:0, num_surr = 1)
    
    # Initialize a list to store te values for the current j
    te_values[[j]] <- numeric(5)  # Preallocate space for 5 values
    
    # Store te values for each tur
    for (k in 1:5) {
      te_values[[j]][k] <- uic_op[[j]]$te[k]  # Store the specific te value for each tur
    }
    
    # Initialize a vector to store counts for the current j
    count_values[[j]] <- numeric(5)  # Preallocate space for 5 counts
    
    # Count how many of the 2 values are less than uic_results[[i]]$te[k]
    for (k in 1:5) {
      count_values[[j]][k] <- sum(te_values[[j]][k] < uic_tur[[i]]$te[k])  # Count values less than
      # Update the sum_counts for the current k
      sum_counts[k] <- sum_counts[k] + count_values[[j]][k]  # Add to the sum for k
    }
  }
  
  # Calculate the ratios_tur based on sum_counts for the current i
  ratios_tur <- 1 - (sum_counts / 1000)  # Calculate the ratios for each k
  return(ratios_tur)  # Return the result for the current i
}

# Use pblapply to process each i in parallel
F_SUR_tur <- pblapply(2:77, process_ts, cl = 72)
F_SUR_tur

saveRDS(F_SUR_tur, file = "~/2025hk_chl/uic_sur_output/F_SUR_tur")
F_SUR_tur <- readRDS("~/2025hk_chl/uic_sur_output/F_SUR_tur")

MIN_SUR_tur <- sapply(F_SUR_tur, min)
MIN_SUR_tur <- data.frame(site_code = colnames(d_chl)[2:77], 
                         pval = unlist(MIN_SUR_tur))


MIN_SUR_tur <- merge(MIN_SUR_tur, d_EPD_site, by = "site_code")

################3------------------------------------------------------------
##-----------------------------NOW CHECK ETE-----------------------------------
########################################-------------------------------------
##since we're going to use ete based on pval based on SUR
##we're making an index based on SUR pval
F_SUR_tur

uic_tur <- readRDS("~/2025hk_chl/uic_sur_output/uic_tur")

pval_indices_tur <- sapply(F_SUR_tur, function(x) {
  min_value <- min(x)  # Find the minimum value
  last_index <- max(which(x == min_value))  # Find the last occurrence of the minimum value
  return(last_index)  # Return the index of the last occurrence
})
pval_indices_tur

##and then we turK the ete row from uic_tur on above index
# Initialize an empty list to store the selected rows
selected_rows <- list()

# Loop through each index from 2 to 102
for (i in 1:76) {
  # Get the row index from pval_indices
  j <- pval_indices_tur[[i]]  # Adjusting for 1-based indexing in R
  
  # Extract the specified row from the dataframe
  selected_row <- uic_tur[[i + 1]][j, ]
  
  # Append the selected row to the list
  selected_rows[[i]] <- selected_row  # Store in the list, adjusting index for 1-based indexing
}


# Combine all selected rows into a new dataframe
ETE_tur <- do.call(rbind, selected_rows)
ETE_tur$tp
ETE_tur$ete
MIN_SUR_tur

MIN_SUR_tur$ete <- ETE_tur$ete

MIN_SUR_tur$sig <- ifelse(MIN_SUR_tur$pval > 0.05, 
                         0.05, MIN_SUR_tur$pval)

saveRDS(MIN_SUR_tur, file = "~/2025hk_chl/uic_sur_output/MIN_SUR_tur")


g15 <- ggmap(hk_detail) + 
  geom_point(data = MIN_SUR_tur, 
             aes(x = lon, y = lat, color = ete, size = sig, 
                 shape = ifelse(pval > 0.05, "Hollow", "Filled"))) +  # Use shape based on significance
  scale_color_gradient2(mid = "darkblue", high = "yellow3", name = "Effect Strength") +
  scale_size(range = c(5, 2), name = "Significance") +  # Reverse the size range
  scale_shape_manual(values = c("Hollow" = 1, "Filled" = 16), labels = NULL) +  # Remove shape labels
  guides(size = guide_legend(order = 2),
         #  color = guide_legend(order = 1), 
         shape = "none") +  #
  labs(title = expression("(d) Turbidity")) +
  labs(x = "Longitude", y = "Latitude") + # Control the order of legends
  theme(
    plot.title = element_text(size = 10),  # Adjust title size
    axis.title.x = element_text(size = 8),  # Adjust x-axis title size
    axis.title.y = element_text(size = 8),  # Adjust y-axis title size
    legend.text = element_text(size = 6),    # Adjust legend text size
    legend.title = element_text(size = 6)  # Adjust legend text size
  )
############################################################################
##------------------------sil-------------------------------
##################################################################################
###-----------------------------------------------------------------------------
##here is surrogate FOR EPD
##----------------------------------------------------------------------------------
str(d_sil)

UIC_sil <- data.frame()  # Initialize an empty data frame to store results
uic_sil <- list()
for (i in 2:77) {
  set.seed(123)
  
  # Step 3: Create block_MODIS data frame
  block_MODIS <- data.frame(d_chl_i = d_chl[[i]], d_sil_i = d_sil[[i]])
  block_MODIS <- apply(block_MODIS, 2, function(x) as.numeric(scale(x)))
  block_MODIS <- as.data.frame(block_MODIS)
  
  # Step 4: Conduct UIC analysis
  uic_sil[[i]] <- uic.optimal(block_MODIS, lib_var = "d_chl_i", tar_var = "d_sil_i",
                              E=0:16, tau=1, tp=-4:0)
  
  # Step 5: Select the row with the minimum pval
  min_pval_row <- uic_sil[[i]][which.min(uic_sil[[i]]$pval), ]
  
  # Step 7: Store results in UM_SST
  UIC_sil <- rbind(UIC_sil, min_pval_row)
}

uic_sil[[77]]

saveRDS(uic_sil, file = "~/2025hk_chl/uic_sur_output/uic_sil")


process_ts <- function(i) {
  set.seed(123)
  ts <- d_sil[[i]]  # Access the column directly
  ts <- as.numeric(ts)
  
  # Replace NA with the mean of the non-NA values
  ts_filled_mean <- na.fill(ts, fill = mean(ts, na.rm = TRUE))
  
  # Generate surrogate data
  surrogate_results <- make_surrogate_data(ts_filled_mean, method = "seasonal", 
                                           num_surr = 1000, T_period = 12)
  
  # Initialize a list to store te values for different sil
  count_values <- list()
  te_values <- list() 
  sum_counts <- numeric(5)  # Initialize a vector to store the sums for each k
  
  # Loop through each column of surrogate_results
  for (j in 1:1000) {
    
    block_MODIS <- data.frame(d_chl = d_chl[[i]], sil_sur = surrogate_results[, j])
    
    
    block_MODIS <- apply(block_MODIS, 2, function(x) as.numeric(scale(x)))
    block_MODIS <- as.data.frame(block_MODIS)
    
    
    
    uic_op[[j]] <- uic.optimal(block_MODIS, lib_var = "d_chl", tar_var = "sil_sur",
                               E = 0:16, tau = 1, tp = -4:0, num_surr = 1)
    
    # Initialize a list to store te values for the current j
    te_values[[j]] <- numeric(5)  # Preallocate space for 5 values
    
    # Store te values for each sil
    for (k in 1:5) {
      te_values[[j]][k] <- uic_op[[j]]$te[k]  # Store the specific te value for each sil
    }
    
    # Initialize a vector to store counts for the current j
    count_values[[j]] <- numeric(5)  # Preallocate space for 5 counts
    
    # Count how many of the 2 values are less than uic_results[[i]]$te[k]
    for (k in 1:5) {
      count_values[[j]][k] <- sum(te_values[[j]][k] < uic_sil[[i]]$te[k])  # Count values less than
      # Update the sum_counts for the current k
      sum_counts[k] <- sum_counts[k] + count_values[[j]][k]  # Add to the sum for k
    }
  }
  
  # Calculate the ratios_sil based on sum_counts for the current i
  ratios_sil <- 1 - (sum_counts / 1000)  # Calculate the ratios for each k
  return(ratios_sil)  # Return the result for the current i
}

# Use pblapply to process each i in parallel
F_SUR_sil <- pblapply(2:77, process_ts, cl = 72)
F_SUR_sil

saveRDS(F_SUR_sil, file = "~/2025hk_chl/uic_sur_output/F_SUR_sil")
F_SUR_sil <- readRDS("~/2025hk_chl/uic_sur_output/F_SUR_sil")

MIN_SUR_sil <- sapply(F_SUR_sil, min)
MIN_SUR_sil <- data.frame(site_code = colnames(d_chl)[2:77], 
                          pval = unlist(MIN_SUR_sil))


MIN_SUR_sil <- merge(MIN_SUR_sil, d_EPD_site, by = "site_code")

################3------------------------------------------------------------
##-----------------------------NOW CHECK ETE-----------------------------------
########################################-------------------------------------
##since we're going to use ete based on pval based on SUR
##we're making an index based on SUR pval
F_SUR_sil

uic_sil <- readRDS("~/2025hk_chl/uic_sur_output/uic_sil")

pval_indices_sil <- sapply(F_SUR_sil, function(x) {
  min_value <- min(x)  # Find the minimum value
  last_index <- max(which(x == min_value))  # Find the last occurrence of the minimum value
  return(last_index)  # Resiln the index of the last occurrence
})
pval_indices_sil

##and then we silK the ete row from uic_sil on above index
# Initialize an empty list to store the selected rows
selected_rows <- list()

# Loop through each index from 2 to 102
for (i in 1:76) {
  # Get the row index from pval_indices
  j <- pval_indices_sil[[i]]  # Adjusting for 1-based indexing in R
  
  # Extract the specified row from the dataframe
  selected_row <- uic_sil[[i + 1]][j, ]
  
  # Append the selected row to the list
  selected_rows[[i]] <- selected_row  # Store in the list, adjusting index for 1-based indexing
}


# Combine all selected rows into a new dataframe
ETE_sil <- do.call(rbind, selected_rows)
ETE_sil$tp
ETE_sil$ete
MIN_SUR_sil

MIN_SUR_sil$ete <- ETE_sil$ete

MIN_SUR_sil$sig <- ifelse(MIN_SUR_sil$pval > 0.05, 
                         0.05, MIN_SUR_sil$pval)
saveRDS(MIN_SUR_sil, file= "~/2025hk_chl/uic_sur_output/MIN_SUR_sil")



g17 <-ggmap(hk_detail) + 
  geom_point(data = MIN_SUR_sil, 
             aes(x = lon, y = lat, color = ete, size = sig, 
                 shape = ifelse(pval > 0.05, "Hollow", "Filled"))) +  # Use shape based on significance
  scale_color_gradient2(mid = "darkblue", high = "yellow3", name = "Effect Strength") +
  scale_size(range = c(5, 2), name = "Significance") +  # Reverse the size range
  scale_shape_manual(values = c("Hollow" = 1, "Filled" = 16), labels = NULL) +  # Remove shape labels
  guides(size = guide_legend(order = 2),
         #  color = guide_legend(order = 1), 
         shape = "none") +  #
  
labs(title = expression("(f) Silica")) +
  labs(x = "Longitude", y = "Latitude") + # Control the order of legends
  theme(
    plot.title = element_text(size = 10),  # Adjust title size
    axis.title.x = element_text(size = 8),  # Adjust x-axis title size
    axis.title.y = element_text(size = 8),  # Adjust y-axis title size
    legend.text = element_text(size = 6),    # Adjust legend text size
    legend.title = element_text(size = 6)  # Adjust legend text size
  )









############################################################################
##------------------------np-------------------------------
##################################################################################
###-----------------------------------------------------------------------------
##here is surrogate FOR EPD
##----------------------------------------------------------------------------------
head(d_TN)


# Combine the first column with the remaining columns to create d_NP
d_NP <- data.frame(d_TN$year_month, d_TN[2:77] / d_TP[2:77])

for (i in 2:77) {
  d_NP[[i]][!is.finite(d_NP[[i]])] <- NA
}



# Optionally, you can set the column names if needed
colnames(d_NP)[1] <- "year_month"

str(d_NP)

UIC_NP <- data.frame()  # Initialize an empty data frame to store results
uic_NP <- list()
for (i in 2:77) {
  set.seed(123)
  
  # Step 3: Create block_MODIS data frame
  block_MODIS <- data.frame(d_chl_i = d_chl[[i]], d_NP_i = d_NP[[i]])
  block_MODIS <- apply(block_MODIS, 2, function(x) as.numeric(scale(x)))
  block_MODIS <- as.data.frame(block_MODIS)
  
  # Step 4: Conduct UIC analysis
  uic_NP[[i]] <- uic.optimal(block_MODIS, lib_var = "d_chl_i", tar_var = "d_NP_i",
                              E=0:16, tau=1, tp=-4:0)
  
  # Step 5: Select the row with the minimum pval
  min_pval_row <- uic_NP[[i]][which.min(uic_NP[[i]]$pval), ]
  
  # Step 7: Store results in UM_SST
  UIC_NP <- rbind(UIC_NP, min_pval_row)
}

uic_NP[[36]]

saveRDS(uic_NP, file = "~/2025hk_chl/uic_sur_output/uic_NP")
uic_NP <- readRDS("~/2025hk_chl/uic_sur_output/uic_NP")



process_ts <- function(i) {
  set.seed(123)
  ts <- d_NP[[i]]  # Access the column directly
  ts <- as.numeric(ts)
  
  # Replace NA with the mean of the non-NA values
  ts_filled_mean <- na.fill(ts, fill = mean(ts, na.rm = TRUE))
  
  # Generate surrogate data
  surrogate_results <- make_surrogate_data(ts_filled_mean, method = "seasonal", 
                                           num_surr = 1000, T_period = 12)
  
  # Initialize a list to store te values for different NP
  count_values <- list()
  te_values <- list() 
  sum_counts <- numeric(5)  # Initialize a vector to store the sums for each k
  
  # Loop through each column of surrogate_results
  for (j in 1:1000) {
    
    block_MODIS <- data.frame(d_chl = d_chl[[i]], NP_sur = surrogate_results[, j])
    
    
    block_MODIS <- apply(block_MODIS, 2, function(x) as.numeric(scale(x)))
    block_MODIS <- as.data.frame(block_MODIS)
    
    
    
    uic_op[[j]] <- uic.optimal(block_MODIS, lib_var = "d_chl", tar_var = "NP_sur",
                               E = 0:16, tau = 1, tp = -4:0, num_surr = 1)
    
    # Initialize a list to store te values for the current j
    te_values[[j]] <- numeric(5)  # Preallocate space for 5 values
    
    # Store te values for each NP
    for (k in 1:5) {
      te_values[[j]][k] <- uic_op[[j]]$te[k]  # Store the specific te value for each NP
    }
    
    # Initialize a vector to store counts for the current j
    count_values[[j]] <- numeric(5)  # Preallocate space for 5 counts
    
    # Count how many of the 2 values are less than uic_results[[i]]$te[k]
    for (k in 1:5) {
      count_values[[j]][k] <- sum(te_values[[j]][k] < uic_NP[[i]]$te[k])  # Count values less than
      # Update the sum_counts for the current k
      sum_counts[k] <- sum_counts[k] + count_values[[j]][k]  # Add to the sum for k
    }
  }
  
  # Calculate the ratios_NP based on sum_counts for the current i
  ratios_NP <- 1 - (sum_counts / 1000)  # Calculate the ratios for each k
  return(ratios_NP)  # Return the result for the current i
}

# Use pblapply to process each i in parallel
F_SUR_NP <- pblapply(2:77, process_ts, cl = 72)
F_SUR_NP

saveRDS(F_SUR_NP, file = "~/2025hk_chl/uic_sur_output/F_SUR_NP")
F_SUR_NP <- readRDS("~/2025hk_chl/uic_sur_output/F_SUR_NP")

MIN_SUR_NP <- sapply(F_SUR_NP, min)
MIN_SUR_NP <- data.frame(site_code = colnames(d_chl)[2:77], 
                          pval = unlist(MIN_SUR_NP))


MIN_SUR_NP <- merge(MIN_SUR_NP, d_EPD_site, by = "site_code")

################3------------------------------------------------------------
##-----------------------------NOW CHECK ETE-----------------------------------
########################################-------------------------------------
##since we're going to use ete based on pval based on SUR
##we're making an index based on SUR pval
F_SUR_NP

uic_NP <- readRDS("~/2025hk_chl/uic_sur_output/uic_NP")

pval_indices_NP <- sapply(F_SUR_NP, function(x) {
  min_value <- min(x)  # Find the minimum value
  last_index <- max(which(x == min_value))  # Find the last occurrence of the minimum value
  return(last_index)  # ReNPn the index of the last occurrence
})
pval_indices_NP

##and then we NPK the ete row from uic_NP on above index
# Initialize an empty list to store the selected rows
selected_rows <- list()

# Loop through each index from 2 to 102
for (i in 1:76) {
  # Get the row index from pval_indices
  j <- pval_indices_NP[[i]]  # Adjusting for 1-based indexing in R
  
  # Extract the specified row from the dataframe
  selected_row <- uic_NP[[i + 1]][j, ]
  
  # Append the selected row to the list
  selected_rows[[i]] <- selected_row  # Store in the list, adjusting index for 1-based indexing
}


# Combine all selected rows into a new dataframe
ETE_NP <- do.call(rbind, selected_rows)
ETE_NP$tp
ETE_NP$ete
MIN_SUR_NP

MIN_SUR_NP$ete <- ETE_NP$ete

MIN_SUR_NP$sig <- ifelse(MIN_SUR_NP$pval > 0.05, 
                          0.05, MIN_SUR_NP$pval)
view(MIN_SUR_NP)
saveRDS(MIN_SUR_NP, file="~/2025hk_chl/uic_sur_output/MIN_SUR_NP")


g16 <- ggmap(hk_detail) + 
  geom_point(data = MIN_SUR_NP, 
             aes(x = lon, y = lat, color = ete, size = sig, 
                 shape = ifelse(pval > 0.05, "Hollow", "Filled"))) +  # Use shape based on significance
  scale_color_gradient2(mid = "darkblue", high = "yellow3", name = "Effect Strength") +
  scale_size(range = c(5, 2), name = "Significance") +  # Reverse the size range
  scale_shape_manual(values = c("Hollow" = 1, "Filled" = 16), labels = NULL) +  # Remove shape labels
  guides(size = guide_legend(order = 2),
         #  color = guide_legend(order = 1), 
         shape = "none") +  #
  
  labs(title = expression("(e) N/P ratio")) +
  labs(x = "Longitude", y = "Latitude") + # Control the order of legends
  theme(
    plot.title = element_text(size = 10),  # Adjust title size
    axis.title.x = element_text(size = 8),  # Adjust x-axis title size
    axis.title.y = element_text(size = 8),  # Adjust y-axis title size
    legend.text = element_text(size = 6),    # Adjust legend text size
    legend.title = element_text(size = 6)  # Adjust legend text size
  )


word_env <- grid.arrange(g18, g12,g13,g15,  g16, g17, ncol = 3, nrow = 2)



# Save the combined figure
ggsave("~/2025hk_chl/uic_sur_output/word_env.jpg", 
       plot = word_env, width = 14, height = 8, dpi = 300) # Adjust dimensions as needed







