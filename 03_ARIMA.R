# Load packages
library(tidyverse); packageVersion("tidyverse") # 2.0.0
library(forecast); packageVersion("forecast") # 8.22.0
library(stringr)
library(lubridate)
library(pbapply)
library(gridExtra)
library(ggplot2)
library(ggmap)
library(ggpubr)
# magrittr, lubridate, and stringr are included in tidyverse, so you do not need to load them

# Load data

dineof <- read.csv("~/2024_HKredtide/data/2002_2022_modis_wholeHK_tryagain_with_colname_dineof.csv")
d_chl <- read.csv("~/2024_HKredtide/06_DataCompileOut/EPD_WholeHK.csv")
# Delete the first row
d_chl <- d_chl[-1, ]

# Delete the last five rows to keep time series is 200208-202207
d_chl <- d_chl[1:(nrow(d_chl) - 5), ]
head(d_chl)
tail(d_chl)
d_chl[is.na(d_chl)] <- 0



##sampling in last few sites(from ET1 to WT3) are new sites starting from 2012
##and they are not monthly, so i divide ts_d1 into two loops



###############below is for ts_d1##########################
###########################################################
# Apply auto.arima to columns 1 to 76
colnames(d_chl)

nmse_ARIMA_EPD_1 <- list()

ARIMA_EPD <- function(i, d_chl_diff) {
  fit <- auto.arima(d_chl_diff[, i],
                    d = NA,
                    max.p = 11, max.q = 11,
                    max.P = 2, max.Q = 2,
                    max.order = 22,
                    seasonal = TRUE,
                    start.p = 0, start.q = 0,
                    start.P = 0, start.Q = 0,
                    parallel = TRUE,
                    stepwise = FALSE, trace = FALSE)
  
  # Calculate the NMSE
  fitted <- fitted(fit)
  nmse <- mean((d_chl_diff[, i] - fitted)^2, na.rm = TRUE) / var(d_chl_diff[, i], na.rm = TRUE)
  
  return(nmse)
}

nmse_ARIMA_EPD_1 <- pblapply(2:77, ARIMA_EPD, d_chl = d_chl, cl = 36)
nmse_ARIMA_EPD_1


saveRDS(nmse_ARIMA_EPD_1, file = "~/2025hk_chl/uic_sur_output/nmse_ARIMA_EPD_1")
nmse_ARIMA_EPD_1 <- readRDS("~/2025hk_chl/uic_sur_output/nmse_ARIMA_EPD_1")

f_smap_epd <- readRDS("~/2025hk_chl/uic_sur_output/f_smap_epd")
write.csv(f_smap_epd, file="~/2025hk_chl/f_smap_epd.csv")
f_smap_epd$nmse
plot(f_smap_epd$nmse[-14],nmse_ARIMA_EPD_1[-14])

arima_epd <- data.frame(smap=f_smap_epd$nmse[-14],
                   arima=unlist(nmse_ARIMA_EPD_1[-14]))

fffs1 <- ggplot(arima_epd, aes(x = arima, y = smap)) +
  geom_point(shape = 1, size = 3) +  # Add hollow circles
  geom_abline(slope = 1, intercept = 0) +  # Add the reference line
  xlim(0.4, 1) +  # Set x-axis limits
  ylim(0.4, 1) +  # Set y-axis limits
  labs(
    x = "NMSE of ARIMA ",  # Label for x-axis
    y = "NMSE of S-map ",  # Label for y-axis
    title = expression("(a) " * italic("In situ") * " measurement data" )  # Title of the plot
  ) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.background = element_blank()) + 
  # Optional: Remove panel background
  theme(axis.line = element_line(),
        axis.text.x = element_text(size = 14),  # Adjust x-axis title size
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),  # Adjust x-axis title size
        axis.title.y = element_text(size = 14),
        plot.title = element_text(size = 16)) 
################################below is for ts_d2#################
#############################################################
head(dineof)
nmse_ARIMA_DIN <- list()
ARIMA_DIN <- function(i, dineof) {
  fit <- auto.arima(dineof[, i],
                    d = NA,
                    max.p = 11, max.q = 11,
                    max.P = 2, max.Q = 2,
                    max.order = 22,
                    seasonal = TRUE,
                    start.p = 0, start.q = 0,
                    start.P = 0, start.Q = 0,
                    parallel = TRUE,
                    stepwise = FALSE, trace = FALSE)
  
  # Calculate the NMSE
  fitted <- fitted(fit)
  nmse <- mean((dineof[, i] - fitted)^2, na.rm = TRUE) / var(dineof[, i], na.rm = TRUE)
  
  return(nmse)
}


nmse_ARIMA_DIN <- pblapply(1:101, ARIMA_DIN, dineof = dineof, cl = 36)
nmse_ARIMA_DIN

saveRDS(nmse_ARIMA_DIN, file = "~/2025hk_chl/uic_sur_output/nmse_ARIMA_DIN")
nmse_ARIMA_DIN <- readRDS("~/2025hk_chl/uic_sur_output/nmse_ARIMA_DIN")

f_smap_modis_din <- readRDS("~/2025hk_chl/uic_sur_output/f_smap_modis_din")


arima_modis <- data.frame(smap=f_smap_modis_din$nmse,
                   arima=unlist(nmse_ARIMA_DIN))

fffs2 <-ggplot(arima_modis, aes(x = arima, y = smap)) +
  geom_point(shape = 1, size = 3) +  # Add hollow circles
  geom_abline(slope = 1, intercept = 0) +  # Add the reference line
  xlim(0.4, 1) +  # Set x-axis limits
  ylim(0.4, 1) +  # Set y-axis limits
  labs(
    x = "NMSE of ARIMA ",  # Label for x-axis
    y = "NMSE of S-map ",  # Label for y-axis
    title = expression("(b) Remote sensing data " )  # Title of the plot
  ) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.background = element_blank()) + 
  # Optional: Remove panel background
  theme(axis.line = element_line(),
        axis.text.x = element_text(size = 14),  # Adjust x-axis title size
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),  # Adjust x-axis title size
        axis.title.y = element_text(size = 14),
        plot.title = element_text(size = 16)) 

word_arima <- grid.arrange(fffs1,fffs2, ncol = 2) 
ggsave("~/2025hk_chl/uic_sur_output/word_arima.tif", 
       plot = word_arima, width = 10, height = 5, dpi = 300) # Adjust dimensions as needed

d_EPD_site <- readRDS("~/2024_HKredtide/06_DataCompileOut/d_EPD_site.obj")


d_EPD_site$color <- NA
d_EPD_site$color[1:7] <- "orange"
d_EPD_site$color[8:23] <- "lightblue"
d_EPD_site$color[24:32] <- "pink"
d_EPD_site$color[33:34] <- "darkseagreen"
d_EPD_site$color[35:39] <- "brown"
d_EPD_site$color[40:45] <- "gray56"
d_EPD_site$color[46:59] <- "gold2"
d_EPD_site$color[60:63] <- "gray87"
d_EPD_site$color[64:66] <- "blue4"
d_EPD_site$color[67:76] <- "purple"

all_region_colors <- c(
  "brown" = "brown", "blue4" = "blue4", "darkseagreen" = "darkseagreen",
  "gold2" = "gold2", "gray56" = "gray56", "pink" = "pink",
  "lightblue" = "lightblue", "orange" = "orange",
  "purple" = "purple", "gray87" = "gray87"
)

all_region_labels <- c(
  "brown" = "Deep Bay", "blue4" = "Eastern", "darkseagreen" = "Junk Bay",
  "gold2" = "Mirs Bay", "gray56" = "North Western", "pink" = "Port Shelter",
  "lightblue" = "Southern", "orange" = "Tolo Harbour",
  "purple" = "Victoria Harbour", "gray87" = "Western"
)


all_levels <- c("brown", "blue4", "darkseagreen", "gold2", "gray56", 
                "pink", "lightblue", "orange", "purple", "gray87")


fffg1 <- ggmap(hk_detail) + 
  geom_point(data = d_EPD_site, aes(x = lon, y = lat, color = color),
             size = 2) +
  scale_color_identity() + 
  labs(title = expression("(a) " * italic("In situ") * " measurement sites")) +
  labs(x = "Longitude", y = "Latitude")+
  scale_color_manual(name = "Site Region",
                     values = all_region_colors,
                     labels = all_region_labels,
                     drop = FALSE) +
  theme(
    plot.title = element_text(size = 14),  
    axis.title.x = element_text(size = 12),  
    axis.title.y = element_text(size = 12),  
    legend.text = element_text(size = 6),    
    legend.title = element_text(size = 6)  ,
    legend.position = "left"
  )





d2_meta_oo$color <- NA
d2_meta_oo$color[1:91] <- "black"

fffg2 <-ggmap(hk_detail) + 
  geom_point(data = d2_meta_oo, aes(x = lon, y = lat, color = color), 
             size = 1.8) + 
  scale_color_identity() + 
  labs(title = expression("(b) Remote sensing sites " )) +
  labs(x = "Longitude", y = "Latitude")+

  theme(
    plot.title = element_text(size = 14),  
    axis.title.x = element_text(size = 12),  
    axis.title.y = element_text(size = 12),  
    axis.text.x = element_text(size = 12),  
    axis.text.y = element_text(size = 12),  
    legend.text = element_text(size = 6),    
    legend.title = element_text(size = 6)  ,
    legend.position = "left"
  )


dummy_df <- data.frame(color = factor(all_levels, levels = all_levels), x = NA, y = NA)

dummy_df$x <- -100
dummy_df$y <- -100
dummy_legend_plot <- ggplot(dummy_df, aes(x = x, y = y, color = color)) +
  geom_point(size = 3) +
  scale_color_manual(name = "Site Region",
                     values = all_region_colors,
                     labels = all_region_labels,
                     drop = FALSE) +
  coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +  
  theme_void() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 16),    
        legend.text = element_text(size = 16)   )

fffg1_nolegend <- fffg1 + theme(legend.position = "none")
fffg2_nolegend <- fffg2 + theme(legend.position = "none") 

fffg12 <- ggarrange(
  ggarrange(fffg1_nolegend,fffg2_nolegend,
            ncol = 2, nrow = 1),
  dummy_legend_plot,
  ncol = 1,
  heights = c(10, 1)
)



ggsave("~/2025hk_chl/uic_sur_output/map.tiff", 
       plot = fffg12, width = 12, height = 6.5, dpi = 300) 



g1c <- data.frame(
  time = d_chl$year_month,  
  epd = d_chl$SM6,      
  modis = nn_pred$S0151   
)

g1c$time <- as.Date(g1c$time, format = "%Y-%m-%d")
fffg1c <- ggplot(g1c, aes(x = time)) +
  geom_line(aes(y = epd, color = "In situ measurement"), 
            size = 1) +  
  geom_line(aes(y = modis, color = "Remote sensing"),
            size = 1) +  
  
  labs(x = "Year", y = expression("Chl-a Concentration (mg m"^-3*")"),
      title = "An example of Chl-a time series", 
       color = "Legend") +  # Add legend title
  scale_color_manual(
    values = c("In situ measurement" = "black", "Remote sensing" = "brown"),
    labels = c(expression(italic("In situ measurement")), "Remote sensing")  
  ) +  # Set colors
  theme_minimal() +
  theme(
        panel.grid = element_blank(),  # Remove grid lines
        axis.line = element_line(color = "black"),
        legend.position = c(0.1, 1),  # Position legend at top left (x = 0.1, y = 0.9)
        legend.justification = c("left", "top") ,
        plot.title = element_text(size = 14),  
        axis.title.x = element_text(size = 12),  
        axis.title.y = element_text(size = 12),  
        legend.text = element_text(size = 12),    
        legend.title = element_text(size = 12)  ) 



ggsave("~/2025hk_chl/uic_sur_output/fffg1c.tiff", 
       plot = fffg1c, width = 12, height = 6, dpi = 300) 

