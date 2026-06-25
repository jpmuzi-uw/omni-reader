
# libraries ---------------------------------------------------------------


library(tidyverse)


# import & clean data -----------------------------------------------------

# ROI data for phantom recons
phantom_data_filename <- r"(data/20250626_NEMA_IQ_Stats_2.csv)"
phantom_data <- read.csv(phantom_data_filename) |>
  select(-seriestime) 

# Focusing on the 50 background regions
bg_data <- phantom_data |>
  filter(grepl("^Bkg",contourname)) |>
  separate_wider_regex(
    contourname,
    c(".+",cyl_index = "\\d{2}","-",slice_index ="0\\d")  #Extract location identifiers in case we want to know individual locations of bkg rois
    ) |> 
  mutate(
    cyl_index = as.integer(cyl_index),
    slice_index = as.integer(slice_index),
    type = case_when( #Extract more filterable recon type distinction from series description
      str_detect(seriesdesc,"AC") ~ "AC",
      str_detect(seriesdesc,"HPDL") ~ "HPDL",
      str_detect(seriesdesc,"MPDL") ~ "MPDL",
      .default = "QCLEAR"
      )
  )  |> 
  select(-rows,-fov,-voxelsize,-suvpeak,-seriesdesc) |> 
  relocate(type)

# summary metrics from Noise Paper
# As I understand it, background variance is the Coefficient of Variation (sd/mean) for the ROIs
# Similarly, Image Roughness is the mean of the CoV for each ROI
noise_metrics <- bg_data |> 
  group_by(type,afd) |> 
  summarize(
    count = n(),
    group_stdev = sd(suvmean),
    group_mean = mean(suvmean),
    image_roughness = mean(stdev/suvmean),
    bkg_var = group_stdev / group_mean
  ) 
# Plots -------------------------------------------------------------------

# Background variance by scan type & duration

noise_metrics |> 
  ggplot(aes(x = afd/1000, y = bkg_var, color = type, shape = type)) +
  geom_point() +
  geom_line() +
  scale_x_continuous(breaks = c(90,120,150,180)) +
  #scale_y_continuous(expand = c(0, 0), limits = c(0, 0.03)) +
  labs(
    title = "Background variance by scan type & duration",
    x = "Scan Duration (seconds)",
    y = "Background variance: CoV (stdev / mean)",
    color = "Reconstruction",
    shape = "Reconstruction"
  ) +
  theme_minimal()

# Image Roughness
noise_metrics |> 
  ggplot(aes(x = afd/1000, y = image_roughness, color = type, shape = type)) +
  geom_point() +
  geom_line() +
  scale_x_continuous(breaks = c(90,120,150,180)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  labs(
    title = "Image Roughness in ROI group by scan type/duration",
    x = "Scan Duration (seconds)",
    y = "Image Roughness",
    color = "Reconstruction",
    shape = "Reconstruction"
  ) +
  theme_minimal()
