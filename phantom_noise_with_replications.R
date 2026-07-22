# libraries ---------------------------------------------------------------

library(tidyverse)
library(paletteer)
library(patchwork)

# load data ---------------------------------------------------------------

omni_phantom_data_filename <- r"(data/harm_contour_stats_4replications_omni-2025-06-26.csv)" #file contains extracted ROI data from 4 replications per duration

# Read CSV into dataframe
omni_phantom_df <- read.csv(omni_phantom_data_filename) |>
  select(-seriestime) |> #redundant, all scans at same time
  filter(!grepl("00", contourname)) |> # remove groups
  filter(!(str_detect(seriesdesc,"^180") & afd == 60000)) # 2026-07-22: fourth replication of 3 minute scan had incorrect duration, filtered out

# Clean data --------------------------------------------------------------

# Extract BKG ROI data and re-structure
omni_bg_df <- omni_phantom_df |> 
  filter(grepl("^Bkg",contourname)) |> # background contours for noise analysis
  separate_wider_regex(
    contourname,
    c(".+",cyl_index = "\\d{2}","-",slice_index ="0\\d")  #Extract location identifiers in case we want to know individual locations of bkg rois
  ) |> 
  mutate(
    cyl_index = as.integer(cyl_index),
    slice_index = as.integer(slice_index),
    type = case_when( #Extract more filterable recon type distinction from series description
      str_detect(seriesdesc,"AC") | str_detect(seriesdesc,"osem") ~ "OSEM",
      str_detect(seriesdesc,"HPDL") ~ "HPDL",
      str_detect(seriesdesc,"MPDL") ~ "MPDL",
      .default = "QCLEAR"),
    replication = str_extract(seriesdesc,"( \\d)",1) #pulls replication number from description
  )  |> 
  mutate(
    replication = ifelse( # initial replication had no number, so replace NA with 1
      is.na(replication),
      yes = 1,
      no = parse_number(replication)
    )
  ) |> 
  select(-rows,-fov,-voxelsize,-suvpeak,-seriesdesc) |> # remove unnecisary fields
  relocate(type)

# Noise Metrics -----------------------------------------------------------

# summary metrics from Noise Paper
# As I understand it, background variance (aka NEMA noise) is the Coefficient of Variation (sd/mean) for the ROIs
# Similarly, Image Roughness is the mean of the CoV for each ROI
noise_metrics <- omni_bg_df |> 
  group_by(replication,type,afd) |> 
  reframe(
    count = n(),
    group_stdev = sd(suvmean),
    group_mean = mean(suvmean),
    image_roughness = mean(stdev/suvmean),
    bkg_var_noise = group_stdev / group_mean
  ) 

# with replications, we can get min/max/mean of previous metrics
summary_noise_metrics <- noise_metrics |> 
  group_by(type,afd) |> 
  reframe(
    group_stdev_min = min(group_stdev),
    group_stdev_max = max(group_stdev),
    group_stdev_mean= mean(group_stdev),
    group_mean_mean = mean(group_mean),
    image_roughness_min = min(image_roughness),
    image_roughness_max = max(image_roughness),
    image_roughness_mean= mean(image_roughness),
    bkg_var_noise_min = min(bkg_var_noise),
    bkg_var_noise_max = max(bkg_var_noise),
    bkg_var_noise_mean= mean(bkg_var_noise)
  )

# Plots -------------------------------------------------------------------

# Background variance (Noise) by scan type & duration
# we are using stdev in mean (and not dividing by mean of mean) because the mean is very consistent
# Mean of replications
mean_stdev_plot <- summary_noise_metrics |>
  ggplot(aes(x = afd/1000, y = group_stdev_mean, color = type, shape = type)) +
  geom_point() +
  geom_line() +
  #geom_segment(aes(x= afd/1000, xend = afd/1000, y = group_stdev_min,yend = group_stdev_max,color = type),alpha=0.5)+
  scale_x_continuous(breaks = c(60,90,120,150,180), minor_breaks = NULL) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.03)) +
  labs(
    title = "Noise by scan type & duration",
    x = "Scan Duration (seconds)",
    y = "Noise (stdev in ROI mean)",
    color = "Reconstruction",
    shape = "Reconstruction",
    subtitle = "Mean values of replications"
  ) +
  theme_minimal()
mean_stdev_plot

# plot above, but using individual replications to show spread, position dodge for legibility
noise_metrics |>
  ggplot(aes(x = afd/1000, y = group_stdev, color = type, shape = type)) +
  geom_point(position = position_dodge(width = 5)) +
  scale_x_continuous(breaks = c(60,90,120,150,180), minor_breaks = NULL) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.03)) +
  labs(
    title = "Noise by scan type & duration",
    x = "Scan Duration (seconds)",
    y = "Noise (stdev in ROI mean)",
    color = "Reconstruction",
    shape = "Reconstruction",
    subtitle = "All replications, position dodge for legibility"
  ) +
  theme_minimal() 

noise_metrics |>
  ggplot(aes(x = afd/1000, y = group_stdev, color = type, shape = type)) +
  geom_point() +
  geom_smooth(se = F)+
  scale_x_continuous(breaks = c(60,90,120,150,180), minor_breaks = NULL) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.03)) +
  labs(
    title = "Noise by scan type & duration",
    x = "Scan Duration (seconds)",
    y = "Noise (stdev in ROI mean)",
    color = "Reconstruction",
    shape = "Reconstruction",
    subtitle = "All replications, smoothed conditional mean trendline"
  ) +
  theme_minimal() 


# Image Roughness - Mean ROI CoV (stdev/mean)
img_rough_plot <- noise_metrics |> 
  ggplot(aes(x = afd/1000, y = image_roughness, color = type, shape = type)) +
  geom_point() +
  geom_smooth(se = F) +
  scale_x_continuous(breaks = c(60,90,120,150,180), minor_breaks = NULL) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.1)) +
  labs(
    title = "Image Roughness in ROI group by scan type/duration",
    x = "Scan Duration (seconds)",
    y = "Mean of ROI CoV",
    color = "Reconstruction",
    shape = "Reconstruction",
    subtitle = "Smoothed mean of replications"
  ) +
  theme_minimal()

# Comparison of Noise Metrics to see if there are any patterns
noise_chull <- noise_metrics |> group_by(type) |> slice(chull(image_roughness,group_stdev))
noise_metrics |> 
  ggplot(aes(x = image_roughness, y = group_stdev, color = type, fill=type, shape = type)) +
  geom_point() +
  geom_polygon(data = noise_chull, alpha = 0.4) +
  scale_fill_viridis_d()+
  scale_color_viridis_d()+
  labs(
    title = "Stdev of SUVmean vs mean CoV in SUVmean",
    x = "mean of ROI CoV",
    y = "Stdev of ROI means",
    color = "Reconstruction",
    shape = "Reconstruction",
    fill  = "Reconstruction"
  )

# Exploratory plots -------------------------------------------------------

# are the slices identical?
omni_bg_df |> 
  ggplot(aes(x = suvmean,y = as.factor(slice_index),fill=as.factor(slice_index))) +
  ggridges::geom_density_ridges(alpha=0.8,quantile_lines = T, quantiles = 2) +
  labs(
    x = "SUVmean",
    y = "Slice Index",
    fill = "Slice Index",
    title = "Distribution of SUVmean in background slices",
    subtitle = "Median indicated"
  ) + theme_minimal() +
  theme(legend.position = "none")

# any patterns in location (slice vs cyl)
omni_bg_df |> 
  ggplot(
    aes(x = as.factor(slice_index), y = as.factor(cyl_index), fill = suvmean)
  ) +
  geom_tile() +
  scale_fill_viridis_c()+
  coord_equal()

# Compare first acquisition to with average
first_noise_metrics <- omni_bg_df |> 
  filter(replication==1) |> 
  group_by(type,afd) |> 
  reframe(
    count = n(),
    group_stdev = sd(suvmean),
    group_mean = mean(suvmean),
    image_roughness = mean(stdev/suvmean),
    bkg_var_noise = group_stdev / group_mean
  )

first_stdev_plot <- first_noise_metrics |>
  ggplot(aes(x = afd/1000, y = group_stdev, color = type, shape = type)) +
  geom_point() +
  geom_line() +
  scale_x_continuous(breaks = c(60,90,120,150,180), minor_breaks = NULL) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.03)) +
  labs(
    title = "Noise by scan type & duration",
    x = "Scan Duration (seconds)",
    y = "Noise (stdev in ROI mean)",
    color = "Reconstruction",
    shape = "Reconstruction",
    subtitle = "Initial reconstruction"
  ) +
  theme_minimal()

first_img_rough_plot <- first_noise_metrics |> 
  ggplot(aes(x = afd/1000, y = image_roughness, color = type, shape = type)) +
  geom_point() +
  geom_line() +
  scale_x_continuous(breaks = c(60,90,120,150,180), minor_breaks = NULL) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.1)) +
  labs(
    title = "Image Roughness in ROI group by scan type/duration",
    x = "Scan Duration (seconds)",
    y = "Mean of ROI CoV",
    color = "Reconstruction",
    shape = "Reconstruction",
    subtitle = "Initial reconstruction"
  ) +
  theme_minimal()

comparison_plot <- (first_stdev_plot + (mean_stdev_plot+labs(title="")) + plot_layout(axes = "collect")) / (first_img_rough_plot + (img_rough_plot + labs(title = "") +theme(legend.position = "none")) + plot_layout(axes = "collect")) + plot_layout(guides = "collect")
comparison_plot

