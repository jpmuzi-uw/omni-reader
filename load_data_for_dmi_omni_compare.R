
# libraries ---------------------------------------------------------------

library(tidyverse)
library(paletteer)

# load data ---------------------------------------------------------------

omni_phantom_data_filename <- r"(data/20250626_NEMA_IQ_Stats_2.csv)"
omni_phantom_df <- read.csv(omni_phantom_data_filename) |>
  select(-seriestime) |> 
  filter(grepl("150",seriesdesc))

dmi_phantom_df <- read_csv("data/harm_contour_stats_DMI-2025-06-26.csv")|>
  select(-seriestime) |> 
  filter(grepl(" AC", seriesdesc)) |> 
  filter(!grepl("00", contourname))

