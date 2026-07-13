
# load data ---------------------------------------------------------------

source("load_data_for_dmi_omni_compare.R")

# sphere analysis ---------------------------------------------------------

omni_spheres <- omni_phantom_df |> filter(grepl("^sph",contourname,ignore.case = T))
dmi_spheres <- dmi_phantom_df |> filter(grepl("^sph",contourname,ignore.case = T))

sphere_data <- bind_rows(omni = omni_spheres,dmi = dmi_spheres, .id = "scanner") |> 
  mutate(seriesdesc = case_when(
    seriesdesc == "150s FDG AC" ~ "OSEM",
    seriesdesc == "150s FDG Qclear" ~ "QClear",
    seriesdesc == "150s FDG Qclear - HPDL" ~ "Qclear HPDL",
    seriesdesc == "150s FDG Qclear - MPDL" ~ "Qclear MPDL",
    seriesdesc == "150s TOF AC" ~ "ToF OSEM",
    seriesdesc == "150s QClear AC" ~ "QClear"
  )) |> 
  mutate(cov = stdev/suvmean) |> 
  pivot_longer(cols = starts_with("suv",ignore.case = T) | stdev | cov, names_to = "measure")

# Plots comparing SUV measures
sphere_suv_plot <- sphere_data |> 
  filter(grepl("suv",measure)) |> 
  ggplot(aes(x = volume, y = value,color = seriesdesc)) +
  geom_line(aes(linetype = scanner),linewidth = 1) +
  geom_point(size = 2) +
  scale_color_paletteer_d("lisa::GustavKlimt") +
  labs(
    x = "Sphere Volume",
    y = "SUV",
    title = "SUV measure comparison for Omni & DMI",
    subtitle = "NEMA IQ Phantom - Scan duration 150s",
    color = "Reconstruction Type",
    linetype = "Scanner"
  ) +
  facet_wrap(~measure)

ggsave(
  plot = sphere_suv_plot,
  filename = "plots/omni_dmi_sphere_suv.png",
  device = "png"
)
# stdev plot
stdev_sphere_plot <- sphere_data |> 
  filter(measure == "stdev") |>   
  ggplot(aes(x = volume, y = value,color = seriesdesc)) +
  geom_line(aes(linetype = scanner),linewidth = 1) +
  geom_point(size = 2) +
  scale_color_paletteer_d("lisa::GustavKlimt") +
  labs(
    x = "Sphere Volume",
    y = "Standard Deviation",
    title = "St. dev. in SUV for Omni & DMI",
    subtitle = "NEMA IQ Phantom - Scan duration 150s",
    color = "Reconstruction Type",
    linetype = "Scanner"
  )
ggsave(
  plot = stdev_sphere_plot,
  filename = "plots/omni_dmi_stdev_sphere.png",
  device = "png"
)

# Background (Noise) ------------------------------------------------------


bg_data <- bind_rows(
  omni = omni_phantom_df,
  dmi = dmi_phantom_df, 
  .id = "scanner"
  )|>
  filter(grepl("^Bkg",contourname)) |>
  separate_wider_regex(
    contourname,
    c(".+",cyl_index = "\\d{2}","-",slice_index ="0\\d")  #Extract location identifiers in case we want to know individual locations of bkg rois
  ) |> 
  mutate(
    cyl_index = as.integer(cyl_index),
    slice_index = as.integer(slice_index),
    type = case_when( #Extract more filterable recon type distinction from series description
      str_detect(seriesdesc,"FDG AC") ~ "OSEM",
      str_detect(seriesdesc,"HPDL")   ~ "QClear HPDL",
      str_detect(seriesdesc,"MPDL")   ~ "QClear MPDL",
      str_detect(seriesdesc,"TOF")    ~ "ToF OSEM",
      .default = "QCLEAR"
    )
  )  |> 
  select(-rows,-fov,-voxelsize,-suvpeak,-seriesdesc,-afd) |> 
  relocate(type)

# summary metrics from Noise Paper
# As I understand it, background variance is the Coefficient of Variation (sd/mean) for the ROIs
# Similarly, Image Roughness is the mean of the CoV for each ROI
noise_metrics <- bg_data |> 
  group_by(scanner,type) |> 
  summarize(
    count = n(),
    group_stdev = sd(suvmean),
    group_mean = mean(suvmean),
    image_roughness = mean(stdev/suvmean),
    bkg_var = group_stdev / group_mean
  ) 
noise_metrics

noise_plot <- noise_metrics |> 
  ggplot(aes(x = type, y = group_stdev, color = scanner)) +
  geom_point(size = 3) +
  scale_x_discrete() +
  scale_y_continuous(limits = c(0, 0.03)) +
  labs(
    title = "Noise in NEMA Phantom Background",
    subtitle = "Noise is defined as stdev in SUVmean between 50 Bkg ROIs",
    x = "Reconstruction Type",
    y = "Noise",
    color = "Scanner",
  )
ggsave(
  plot = noise_plot,
  filename = "plots/omni_dmi_bkg_noise.png",
  device = "png"
)
