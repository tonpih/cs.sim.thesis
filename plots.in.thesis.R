library(dplyr)
library(tidyr)
library(ggplot2)
library(sf)
setwd("~/Library/CloudStorage/OneDrive-NTNU/NARM/MSNARM/Footprint/Data/R.master/leirelv")



merged_scenarios <- readRDS("saved/occupancy/merged_scenarios.rds")
head(merged_scenarios)

moose_baseline <- st_read("moose4.shp")

#View(moose_baseline)
#----- plot dominant habitat types in study area-----------------
#subset to studyarea
moose_studyarea <- moose_baseline[moose_baseline$studyarea != 0, ]
head(moose_studyarea)

library(sf)
library(dplyr)
library(tidyr)
library(ggplot2)

hab_cols <- c(
  "areatype10",
  "areatype20",
  "areatype30",
  "areatype50",
  "areatype60",
  "areatype81"
)

moose_dom <- moose_studyarea %>%
  st_drop_geometry() %>%
  mutate(cell_id = row_number()) %>%
  pivot_longer(
    cols = all_of(hab_cols),
    names_to = "habitat",
    values_to = "area"
  ) %>%
  group_by(cell_id) %>%
  slice_max(area, n = 1, with_ties = FALSE) %>%
  ungroup()

moose_plot <- moose_studyarea %>%
  mutate(cell_id = row_number()) %>%
  left_join(
    moose_dom %>% select(cell_id, habitat),
    by = "cell_id"
  )

moose_plot <- moose_plot %>%
  mutate(
    habitat = recode(
      habitat,
      areatype10 = "Developed areas",
      areatype20 = "Agriculture",
      areatype30 = "Forest",
      areatype50 = "Open terrain",
      areatype60 = "Peatland",
      areatype81 = "Freshwater"
    )
  )

moose_plot$habitat <- factor(
  moose_plot$habitat,
  levels = c(
    "Developed areas", 
    "Agriculture",
    "Forest",
    "Open terrain",
    "Peatland",
    "Freshwater"
  )
)

ggplot() +
  geom_sf(
    data = moose_plot,
    aes(fill = habitat),
    color = NA
  ) +
  geom_sf(
    data = moose_plot %>% filter(studyarea == 2),
    fill = NA,
    color = "black",
    linewidth = 0.3
  ) +
  scale_fill_manual(
    values = c(
      "Developed areas" = "#bdbdbd",  # grey
      "Agriculture"     = "#d9d976",  # yellow-green
      "Forest"          = "#1b7837",  # dark green
      "Open terrain"    = "#e69f00",  # strong ochre/orange (high contrast)
      "Peatland"        = "#8c6d31",  # brown
      "Freshwater"      = "#2b83ba"   # blue
    ),
    name = "Dominant habitat",
    drop = FALSE
  ) +
  labs(
    title = "Dominant habitat type in each grid cell",
    subtitle = "Corridor cells outlined in black"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "grey85"),
    legend.position = "right"
  )


hab_cols <- c("areatype10","areatype20","areatype30",
              "areatype50","areatype60","areatype81")

moose_plot <- moose_plot %>%
  mutate(
    total_area = rowSums(across(all_of(hab_cols))),
    max_area   = pmax(!!!syms(hab_cols)),
    dominance  = max_area / total_area
  )

geom_sf(
  data = moose_plot,
  aes(fill = habitat, alpha = dominance),
  color = NA
) +
  scale_alpha(range = c(0.4, 1), guide = "none")



#-------Sampling summaries ---------------

overview_table <- merged_scenarios %>%
  group_by(scenario_type, n_samples) %>%
  summarise(
    mean_inside = mean(mean_sampled_inside),
    mean_outside = mean(mean_sampled_outside),
    .groups = "drop"
  )

overview_table


# 1. Compute average proportions per scenario type
prop_by_scenario <- merged_scenarios %>%
  dplyr::select(
    scenario_type,
    mean_sampled_inside,
    mean_sampled_outside
  ) %>%
  group_by(scenario_type) %>%
  summarise(
    mean_inside = mean(mean_sampled_inside),
    mean_outside = mean(mean_sampled_outside),
    .groups = "drop"
  ) %>%
  mutate(
    total = mean_inside + mean_outside,
    pct_inside = mean_inside / total,
    pct_outside = mean_outside / total
  ) %>%
  dplyr::select(scenario_type, pct_inside, pct_outside)

# 2. Convert to long format for pie plotting
plot_data <- prop_by_scenario %>%
  pivot_longer(
    cols = c(pct_inside, pct_outside),
    names_to = "location",
    values_to = "proportion"
  ) %>%
  mutate(
    location = recode(location,
                      pct_inside = "Inside corridor",
                      pct_outside = "Outside corridor"
    ),
    label = paste0(round(proportion * 100, 1), "%")
  )

# 3. Pie chart with labels
ggplot(plot_data, aes(x = "", y = proportion, fill = location)) +
  geom_col(color = "white", width = 1) +
  geom_text(aes(label = label),
            position = position_stack(vjust = 0.5),
            size = 4, color = "white", fontface = "bold") +
  coord_polar(theta = "y") +
  facet_wrap(~ scenario_type, ncol = 3) +
  scale_fill_manual(values = c("#D55E00", "#0072B2")) +
  theme_void(base_size = 14) +
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "bottom"
  ) +
  labs(
    title = "Average Sampling Inside vs Outside Corridor",
    fill = "Location"
  )



#-------Detection summaries - pie chart:-----

# 1) Compute average detection proportions per scenario type
prop_by_scenario_det <- merged_scenarios %>%
  dplyr::select(
    scenario_type,
    mean_detected_inside,
    mean_detected_outside
  ) %>%
  group_by(scenario_type) %>%
  summarise(
    mean_inside = mean(mean_detected_inside, na.rm = TRUE),
    mean_outside = mean(mean_detected_outside, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    total = mean_inside + mean_outside,
    pct_inside  = ifelse(total > 0, mean_inside / total, NA_real_),
    pct_outside = ifelse(total > 0, mean_outside / total, NA_real_)
  ) %>%
  dplyr::select(scenario_type, pct_inside, pct_outside)

# 2) Convert to long format for pie plotting
plot_data_det <- prop_by_scenario_det %>%
  pivot_longer(
    cols = c(pct_inside, pct_outside),
    names_to = "location",
    values_to = "proportion"
  ) %>%
  mutate(
    location = recode(location,
                      pct_inside = "Inside corridor",
                      pct_outside = "Outside corridor"),
    label = ifelse(is.na(proportion), "NA", paste0(round(proportion * 100, 1), "%"))
  )

# 3) Pie chart with labels
ggplot(plot_data_det, aes(x = "", y = proportion, fill = location)) +
  geom_col(color = "white", width = 1) +
  geom_text(aes(label = label),
            position = position_stack(vjust = 0.5),
            size = 4, color = "white", fontface = "bold") +
  coord_polar(theta = "y") +
  facet_wrap(~ scenario_type, ncol = 3) +
  scale_fill_manual(values = c("#D55E00", "#0072B2")) +
  theme_void(base_size = 14) +
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "bottom"
  ) +
  labs(
    title = "Average Detections Inside vs Outside Corridor",
    fill = "Location"
  )


#----sampling and detection summaries combined --------
library(dplyr)
library(tidyr)
library(ggplot2)

scenario_order <- c("baseline", "strava", "core.sites", "targeted", "5050", "clustered")
plot_sampling <- plot_data %>%
  mutate(
    summary_type = "Sampling",
    scenario_type = factor(scenario_type, levels = scenario_order)
  )
plot_detection <- plot_data_det %>%
  mutate(
    summary_type = "Detection",
    scenario_type = factor(scenario_type, levels = scenario_order)
  )

plot_combined <- bind_rows(plot_sampling, plot_detection)

ggplot(plot_combined, aes(x = "", y = proportion, fill = location)) +
  geom_col(color = "white", width = 1) +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size = 4,
    color = "white",
    fontface = "bold",
    na.rm = TRUE
  ) +
  coord_polar(theta = "y") +
  facet_grid(
    summary_type ~ scenario_type
  ) +
  scale_fill_manual(
    values = c(
      "Inside corridor" = "#D55E00",
      "Outside corridor" = "#0072B2"
    )
  ) +
  theme_void(base_size = 14) +
  theme(
    strip.text.x = element_text(size = 13, face = "bold"),
    strip.text.y = element_text(size = 14, face = "bold"),
    legend.position = "bottom"
  ) +
  labs(
    title = "   Average Sampling and Detection Inside vs Outside Corridor  ",
    subtitle = "",
    fill = "Location"
  )



library(dplyr)
library(tidyr)

prop_by_scenario_both <- merged_scenarios %>%
  dplyr::select(
    scenario_type,
    mean_detected_inside, mean_detected_outside,
    mean_sampled_inside,  mean_sampled_outside
  ) %>%
  dplyr::group_by(scenario_type) %>%
  dplyr::summarise(
    det_in  = mean(mean_detected_inside, na.rm = TRUE),
    det_out = mean(mean_detected_outside, na.rm = TRUE),
    samp_in  = mean(mean_sampled_inside, na.rm = TRUE),
    samp_out = mean(mean_sampled_outside, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    det_total  = det_in + det_out,
    samp_total = samp_in + samp_out,
    det_in_pct  = dplyr::if_else(det_total  > 0, det_in  / det_total,  NA_real_),
    det_out_pct = dplyr::if_else(det_total  > 0, det_out / det_total, NA_real_),
    samp_in_pct  = dplyr::if_else(samp_total > 0, samp_in  / samp_total, NA_real_),
    samp_out_pct = dplyr::if_else(samp_total > 0, samp_out / samp_total, NA_real_)
  )

det_long <- prop_by_scenario_both %>%
  dplyr::select(scenario_type, inside = det_in_pct, outside = det_out_pct) %>%
  tidyr::pivot_longer(c(inside, outside), names_to = "where", values_to = "pct") %>%
  dplyr::mutate(ring = "Detections")

samp_long <- prop_by_scenario_both %>%
  dplyr::select(scenario_type, inside = samp_in_pct, outside = samp_out_pct) %>%
  tidyr::pivot_longer(c(inside, outside), names_to = "where", values_to = "pct") %>%
  dplyr::mutate(ring = "Sampled cells")

#donut plot

prop_det <- merged_scenarios %>%
  dplyr::select(
    scenario_type,
    mean_detected_inside,
    mean_detected_outside
  ) %>%
  group_by(scenario_type) %>%
  summarise(
    inside = mean(mean_detected_inside, na.rm = TRUE),
    outside = mean(mean_detected_outside, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    total = inside + outside,
    pct_inside  = inside / total,
    pct_outside = outside / total,
    metric = "Detections"
  ) %>%
  tidyr::pivot_longer(
    cols = c(pct_inside, pct_outside),
    names_to = "location",
    values_to = "pct"
  )

prop_samp <- merged_scenarios %>%
  dplyr::select(
    scenario_type,
    mean_sampled_inside,
    mean_sampled_outside
  ) %>%
  group_by(scenario_type) %>%
  summarise(
    inside = mean(mean_sampled_inside, na.rm = TRUE),
    outside = mean(mean_sampled_outside, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    total = inside + outside,
    pct_inside  = inside / total,
    pct_outside = outside / total,
    metric = "Sampled cells"
  ) %>%
  tidyr::pivot_longer(
    cols = c(pct_inside, pct_outside),
    names_to = "location",
    values_to = "pct"
  )

pie_df <- dplyr::bind_rows(prop_det, prop_samp) %>%
  mutate(
    location = recode(
      location,
      pct_inside  = "Inside corridor",
      pct_outside = "Outside corridor"
    ),
    metric = factor(metric, levels = c("Detections", "Sampled cells"))
  )

library(dplyr)

pie_df <- pie_df %>%
  group_by(scenario_type, metric) %>%
  arrange(metric, location) %>%
  mutate(
    label = scales::percent(pct, accuracy = 1),
    ypos  = cumsum(pct) - 0.5 * pct
  ) %>%
  ungroup()


library(ggplot2)
library(ggplot2)

ggplot(
  pie_df,
  aes(
    x = metric,
    y = pct,
    fill = location
  )
) +
  geom_col(
    width = 1,
    color = "white"
  ) +
  geom_text(
    aes(
      y = ypos,
      label = label
    ),
    color = "black",
    size = 3
  ) +
  coord_polar(theta = "y") +
  facet_wrap(~ scenario_type) +
  scale_x_discrete(
    limits = c("Detections", "Sampled cells")
  ) +
  scale_fill_manual(
    values = c(
      "Inside corridor"  = "#3182BD",
      "Outside corridor" = "#9ECAE1"
    )
  ) +
  theme_void() +
  labs(
    title = "Proportion inside and outside the corridor",
    subtitle = "Inner ring: detections · Outer ring: sampled cells",
    fill = NULL
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10)
  )



#########

pie_df <- pie_df %>%
  mutate(
    # this order controls stacking AND label placement
    location = factor(location, levels = c("Inside corridor", "Outside corridor")),
    metric = factor(metric, levels = c("Detections", "Sampled cells")),
    label = scales::percent(pct, accuracy = 1)
  )

ggplot(pie_df, aes(x = metric, y = pct, fill = location)) +
  geom_col(width = 1, color = "white") +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),  # <-- key
    color = "black",
    size = 3
  ) +
  coord_polar(theta = "y") +
  facet_wrap(~ scenario_type) +
  scale_fill_manual(
    values = c("Inside corridor" = "#2171B5", "Outside corridor" = "#BDD7E7")
  ) +
  theme_void() +
  labs(
    title = "Proportion of sampled cells and detections inside and outside the corridor",
    subtitle = "Inner ring: detections · Outer ring: sampled cells",
    fill = NULL
  )

pie_df <- pie_df %>%
  dplyr::mutate(
    location = factor(location, levels = c("Inside corridor", "Outside corridor")),
    metric   = factor(metric, levels = c("Detections", "Sampled cells")),
    label    = scales::percent(pct, accuracy = 1)
  )

ggplot(pie_df, aes(x = metric, y = pct, fill = location)) +
  geom_col(width = 1, color = "white") +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    color = "black",
    size = 4.2,            # larger labels (try 4–5)
    fontface = "bold"      # bold percentages
  ) +
  coord_polar(theta = "y") +
  facet_wrap(~ scenario_type) +
  scale_fill_manual(
    values = c("Inside corridor" = "#3B8BC4", "Outside corridor" = "#9ECAE1")
  ) +
  theme_void() +
  labs(
    title = "Proportion of sampled sites and detections inside and outside the corridor",
    subtitle = "Inner ring: detections | Outer ring: sampled sites \n",
    fill = NULL
  ) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 14),
    strip.text = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 11)
  )


##--##
scenario_order <- c("baseline", "strava", "core.sites", "targeted", "5050", "clustered")

pie_df <- pie_df %>%
  dplyr::mutate(
    scenario_type = factor(scenario_type, levels = scenario_order),
    location = factor(location, levels = c("Inside corridor", "Outside corridor")),
    metric   = factor(metric, levels = c("Detections", "Sampled cells")),
    label    = scales::percent(pct, accuracy = 1)
  )

ggplot(pie_df, aes(x = metric, y = pct, fill = location)) +
  geom_col(width = 1, color = "white") +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    color = "black",
    size = 4.2,
    fontface = "bold"
  ) +
  coord_polar(theta = "y") +
  facet_wrap(~ scenario_type) +   # will now use scenario_order
  scale_fill_manual(values = c("Inside corridor" = "#3B8BC4",
                               "Outside corridor" = "#9ECAE1")) +
  theme_void() +
  labs(
    title = "Proportion of sampled sites and detections. Inside and outside the corridor",
    subtitle = "Inner ring: detections | Outer ring: sampled sites\n",
    fill = NULL
  ) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 14),
    strip.text = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 11)
  )
# install.packages("ggtext")  # if needed
library(ggtext)

ggplot(pie_df, aes(x = metric, y = pct, fill = location)) +
  geom_col(width = 1, color = "white") +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    color = "black",
    size = 4.2,
    fontface = "bold"
  ) +
  coord_polar(theta = "y") +
  facet_wrap(~ scenario_type) +
  scale_fill_manual(values = c("Inside corridor" = "#3B8BC4", "Outside corridor" = "#9ECAE1")) +
  theme_void() +
  labs(
    title = "Proportion of sampled cells and detections inside and outside the corridor",
    subtitle = "inner ring = detections   |   outer ring = sampled cells",
    fill = NULL
  ) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = ggtext::element_textbox_simple(
      size = 13,
      face = "bold",
      color = "black",
      fill = "grey95",
      box.color = "grey95",
      linewidth = 0.6,
      padding = margin(4, 8, 4, 8),
      margin = margin(6, 0, 10, 0)
    ),
    strip.text = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 11)
  )

#-------Mean total detections per total visits -----
#scenario_colors <- c(
#  "baseline"   = "#66C2A5",
#  "strava"     = "#FC8D62",
#  "core.sites" = "#8DA0CB",
#  "targeted"   = "#E78AC3",
#  "5050"       = "#A6D854",
#  "clustered"  = "#FFD92F"
#)

scenario_order <- c("baseline", "strava", "core.sites", "targeted", "5050", "clustered")

# These are the colors you want (from your palette)
scenario_colors <- c(
  baseline   = "#FC8D62",  # salmon
  strava     = "#B3A200",  # olive/gold (adjust if your exact strava bar is different)
  core.sites = "#00BA38",  # bright green
  targeted   = "#00BFC4",  # cyan/teal
  `5050`     = "#619CFF",  # blue
  clustered  = "#F564E3"   # magenta
)


scenario_order <- c("baseline", "strava", "core.sites", "targeted", "5050", "clustered")

detections_visits_all <- merged_scenarios %>%
  dplyr::mutate(
    scenario_type = factor(scenario_type, levels = scenario_order)
  ) %>%
  dplyr::select(
    scenario_type,
    total_visits,
    mean_detected_inside,
    mean_detected_outside
  ) %>%
  dplyr::mutate(total_detected = mean_detected_inside + mean_detected_outside) %>%
  dplyr::group_by(scenario_type, total_visits) %>%
  dplyr::summarise(
    mean_total_detected = mean(total_detected),
    .groups = "drop"
  ) %>%
  dplyr::arrange(total_visits, scenario_type)


ggplot(detections_visits_all,
       aes(x = factor(total_visits), y = mean_total_detected, fill = scenario_type)) +
  geom_col(position = position_dodge(width = 0.8),
           width = 0.7, color = "black") +
  scale_fill_manual(values = scenario_colors,
                    breaks = scenario_order,   # forces legend order
                    name = "Scenario type") +
  theme_bw(base_size = 14) +
  labs(title = "Mean Total Detections per Total Visits",
       x = "Total visits", y = "Mean detections")




#-----RMSE across scenario types ------
scenario_order <- c("baseline", "strava", "core.sites", "targeted", "5050", "clustered")

merged_scenarios <- merged_scenarios %>%
  dplyr::mutate(
    scenario_type = factor(scenario_type, levels = scenario_order)
  )

ggplot(merged_scenarios,
       aes(x = scenario_type, y = rmse, fill = scenario_type)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.8) +
  theme_minimal() +
  scale_fill_manual(values = scenario_colors, breaks = scenario_order) +
  labs(
    title = "RMSE across scenario types",
    x = "Scenario type",
    y = "RMSE"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

#zoomed in 
ggplot(merged_scenarios,
       aes(x = scenario_type, y = rmse, fill = scenario_type)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.8) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_minimal() +
  scale_fill_manual(values = scenario_colors, breaks = scenario_order) +
  labs(
    title = "RMSE across scenario types (zoomed in)",
    x = "Scenario type",
    y = "RMSE"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )



#-----RMSE across sampling effort-----
#----visits
ggplot(merged_scenarios,
       aes(x = total_visits, y = rmse, color = scenario_type)) +
  stat_summary(fun = median, geom = "line", linewidth = 1) +
  stat_summary(fun = median, geom = "point", size = 2) +
  theme_minimal() +
  labs(
    title = "RMSE Across Sampling Effort (total visits)",
    x = "Number of visits",
    y = "RMSE",
    color = "Scenario type"
  )

#----sample size 
ggplot(merged_scenarios,
       aes(x = n_samples, y = rmse, color = scenario_type)) +
  stat_summary(fun = median, geom = "line", linewidth = 1) +
  stat_summary(fun = median, geom = "point", size = 2) +
  theme_minimal() +
  labs(
    title = "RMSE Across Number of Sampled Sites",
    x = "Number of sampled sites",
    y = "RMSE",
    color = "Scenario type"
  )

#------RMSE across detection probability --------
ggplot(merged_scenarios,
       aes(x = p_detect, y = rmse, color = scenario_type)) +
  stat_summary(fun = median, geom = "line", linewidth = 1) +
  stat_summary(fun = median, geom = "point", size = 2) +
  theme_minimal() +
  labs(
    title = "RMSE Across Detection Probability",
    x = "Detection probability",
    y = "RMSE",
    color = "Scenario type"
  )

ggplot(merged_scenarios,
       aes(x = as.numeric(as.character(p_detect)),
           y = rmse,
           color = scenario_type,
           group = scenario_type)) +
  stat_summary(fun = median, geom = "line", linewidth = 1) +
  stat_summary(fun = median, geom = "point", size = 2) +
  theme_minimal() +
  labs(
    title = "RMSE Across Detection Probability (p_detect)",
    x = "Detection probability",
    y = "RMSE",
    color = "Scenario type"
  )


#------BIAS across scenario types--------


#adding zero line
ggplot(merged_scenarios,
       aes(x = scenario_type, y = bias, fill = scenario_type)) +
  geom_hline(yintercept = 0, 
             linetype = "dashed", 
             color = "gray30", 
             linewidth = 0.8) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.8) +
  theme_minimal() +
  scale_fill_manual(values = scenario_colors) +
  labs(
    title = "Bias across Scenario Types",
    x = "Scenario type",
    y = "Bias"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggplot(merged_scenarios,
       aes(x = scenario_type, y = bias, fill = scenario_type)) +
  geom_hline(yintercept = 0, 
             linetype = "dashed", 
             color = "gray30", 
             linewidth = 0.8) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.8) +
  coord_cartesian(ylim = c(-0.1, 0.1)) +
  theme_minimal() +
  scale_fill_manual(values = scenario_colors) +
  labs(
    title = "Bias across Scenario Types (zoomed in)",
    x = "Scenario type",
    y = "Bias"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


#------BIAS across sampling effort ----------
#--- adding zero line:

#samples
ggplot(merged_scenarios,
       aes(x = n_samples, y = bias, color = scenario_type)) +
  geom_hline(yintercept = 0,
             linetype = "dashed",
             color = "gray30",
             linewidth = 0.8) +
  stat_summary(fun = median, geom = "line", linewidth = 1) +
  stat_summary(fun = median, geom = "point", size = 2) +
  theme_minimal() +
  labs(
    title = "Bias Across Number of Sampled Sites",
    x = "Number of sampled sites",
    y = "Bias",
    color = "Scenario type"
  )

#visits
ggplot(merged_scenarios,
       aes(x = total_visits, y = bias, color = scenario_type)) +
  geom_hline(yintercept = 0,
             linetype = "dashed",
             color = "gray30",
             linewidth = 0.8) +
  stat_summary(fun = median, geom = "line", linewidth = 1) +
  stat_summary(fun = median, geom = "point", size = 2) +
  theme_minimal() +
  labs(
    title = "Bias Across Sampling Effort (total visits)",
    x = "Total visits",
    y = "Bias",
    color = "Scenario type"
  )

# visits
ggplot(merged_scenarios,
       aes(x = total_visits, y = bias, color = scenario_type, group = scenario_type)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30", linewidth = 0.8) +
  stat_summary(fun = median, geom = "line", linewidth = 1) +
  stat_summary(fun = median, geom = "point", size = 2) +
  theme_minimal() +
  labs(
    title = "Bias Across Sampling Effort (total_visits)",
    x = "Total visits",
    y = "Bias",
    color = "Scenario type"
  )


#-------RMSE across different conditions of effort and detection ------
#------Comparing effort and detection -------------

rmse_drop_effort <- merged_scenarios %>%
  filter(total_visits %in% c(400, 2000)) %>%
  group_by(scenario_type, total_visits) %>%
  summarise(rmse = mean(rmse), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = total_visits, values_from = rmse,
                     names_prefix = "visits_") %>%
  mutate(rmse_drop = visits_400 - visits_2000)

rmse_drop_detect <- merged_scenarios %>%
  filter(p_detect %in% c(0.6, 0.9)) %>%
  group_by(scenario_type, p_detect) %>%
  summarise(rmse = mean(rmse), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = p_detect, values_from = rmse,
                     names_prefix = "p_") %>%
  mutate(rmse_drop = p_0.6 - p_0.9)



# Effort drop plot
p1 <- ggplot(rmse_drop_effort,
             aes(x = scenario_type, y = rmse_drop, fill = scenario_type)) +
  geom_col(alpha = 0.9) +
  scale_fill_manual(values = scenario_colors) +
  theme_minimal() +
  labs(
    title = "Reduction in RMSE from Low to High Sampling Effort",
    subtitle = "Difference between total_visits = 400 and 2000",
    x = "Scenario type",
    y = "RMSE drop"
  ) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

# Detection drop plot
p2 <- ggplot(rmse_drop_detect,
             aes(x = scenario_type, y = rmse_drop, fill = scenario_type)) +
  geom_col(alpha = 0.9) +
  scale_fill_manual(values = scenario_colors) +
  theme_minimal() +
  labs(
    title = "Reduction in RMSE from Low to High Detection Probability",
    subtitle = "Difference between p_detect = 0.6 and 0.9",
    x = "Scenario type",
    y = "RMSE drop"
  ) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

p1
p2


rmse_diff_p <- merged_scenarios %>%
  filter(p_detect %in% c(0.6, 0.9)) %>%
  group_by(scenario_type, p_detect) %>%
  summarise(median_rmse = median(rmse), .groups = "drop") %>%
  pivot_wider(
    names_from = p_detect,
    values_from = median_rmse,
    names_prefix = "p_"
  ) %>%
  mutate(rmse_drop = p_0.6 - p_0.9)

rmse_diff_p

rmse_diff <- merged_scenarios %>%
  filter(total_visits %in% c(400, 2000)) %>%
  group_by(scenario_type, total_visits) %>%
  summarise(median_rmse = median(rmse), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from = total_visits,
    values_from = median_rmse,
    names_prefix = "visits_"
  ) %>%
  mutate(rmse_drop = visits_400 - visits_2000)

rmse_diff

coord_cartesian(ylim = c(0, 1.6))
ggplot(rmse_diff,
       aes(x = scenario_type, y = rmse_drop, fill = scenario_type)) +
  geom_col(alpha = 0.8) +
  coord_cartesian(ylim = c(0, 1.7)) +   # <<< SAME SCALE
  theme_minimal() +
  labs(
    title = "Reduction in RMSE from Low to High Sampling Effort",
    subtitle = "Difference between total_visits = 400 and 2000 (aggregated across species and p_detect)",
    x = "Scenario type",
    y = "RMSE drop (400 − 2000)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  scale_fill_manual(values = scenario_colors)

ggplot(rmse_diff_p,
       aes(x = scenario_type, y = rmse_drop, fill = scenario_type)) +
  geom_col(alpha = 0.8) +
  coord_cartesian(ylim = c(0, 1.7)) +   # <<< SAME SCALE
  theme_minimal() +
  labs(
    title = "Reduction in RMSE from Low to High Detection Probability",
    subtitle = "Difference between p_detect = 0.6 and 0.9 (aggregated across species and sampling effort)",
    x = "Scenario type",
    y = "RMSE drop (0.6 − 0.9)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  scale_fill_manual(values = scenario_colors)


#statistical model

#treating p_detect and total_visits as continuous numerical variables, kind off wrong... 
mod <- lm(rmse ~ p_detect + total_visits + scenario_type, data = merged_scenarios)
summary(mod)

#as factor
merged_scenarios$p_detect <- factor(merged_scenarios$p_detect)
merged_scenarios$total_visits <- factor(merged_scenarios$total_visits)

mod2 <- lm(rmse ~ p_detect + total_visits + scenario_type, data = merged_scenarios)
summary(mod2)

mod3 <- lm(rmse ~ total_visits * scenario_type + p_detect, data = merged_scenarios) #with interaction visits and scenario type
summary(mod3)


table(merged_scenarios$p_detect)
table(merged_scenarios$total_visits)
xtabs(~ p_detect + total_visits + scenario_type, data = merged_scenarios)


#----RMSE across 3 levels of effort and detection conditions
# 1. Create numeric setup labels
plot_data <- merged_scenarios %>%
  mutate(
    setup_num = case_when(
      p_detect == 0.6 & total_visits == 400  ~ "1",   # low det, low effort
      p_detect == 0.9 & total_visits == 400  ~ "2",   # high det, low effort
      p_detect == 0.6 & total_visits == 2000 ~ "3",   # low det, high effort
      TRUE ~ NA_character_
    ),
    setup_desc = case_when(
      setup_num == "1" ~ "Low det, low effort",
      setup_num == "2" ~ "High det, low effort",
      setup_num == "3" ~ "Low det, high effort"
    )
  ) %>%
  filter(!is.na(setup_num))

# 2. Ensure numeric setup order
plot_data$setup_num <- factor(plot_data$setup_num, 
                              levels = c("1", "2", "3"))

# 3. ORDER scenario types:
plot_data$scenario_type <- factor(
  plot_data$scenario_type,
  levels = c("baseline", "strava", "core.sites", "targeted", "5050", "clustered")
)

# 4. Plot
ggplot(plot_data,
       aes(x = setup_num, y = rmse, fill = setup_num)) +
  geom_col() +
  facet_wrap(~ scenario_type, scales = "free_y") +
  theme_bw(base_size = 13) +
  labs(
    title = "RMSE Across Effort and Detection Conditions",
    x = "Sampling Setup",
    y = "RMSE",
    fill = "Setup"
  ) +
  scale_fill_manual(
    values = c(
      "1" = "#5AA1D6",   # darker pastel blue
      "2" = "#7BBCE8",   # medium pastel blue
      "3" = "#A6D1F0"    # light pastel blue
    ),
    labels = c(
      "1 = Low detection, low effort",
      "2 = High detection, low effort",
      "3 = Low detection, high effort"
    )
  ) +
  theme(
    strip.text = element_text(size = 12),
    axis.text.x = element_text(size = 11)
  )

ggplot(plot_data, aes(x = setup_num, y = rmse, fill = setup_num)) +
  geom_col() +
  facet_wrap(~ scenario_type) +   # <- remove scales="free_y"
  theme_bw(base_size = 13) +
  labs(
    title = "RMSE Across Effort and Detection Conditions",
    x = "Sampling Setup",
    y = "RMSE",
    fill = "Setup"
  ) +
  scale_fill_manual(
    values = c("1" = "#5AA1D6", "2" = "#7BBCE8", "3" = "#A6D1F0"),
    labels = c(
      "1 = Low detection, low effort",
      "2 = High detection, low effort",
      "3 = Low detection, high effort"
    )
  ) +
  theme(
    strip.text = element_text(size = 12),
    axis.text.x = element_text(size = 11)
  )


#-------Spatial bias --- scenario types ----
scenario_order <- c("baseline", "strava", "core.sites", "targeted", "5050", "clustered")

plot_df <- merged_scenarios %>%
  dplyr::filter(total_visits == 800, p_detect == 0.7) %>%
  dplyr::mutate(scenario_type = factor(scenario_type, levels = scenario_order)) %>%
  dplyr::select(species, scenario_type, bias, sd_beta_hat, rmse) %>%
  tidyr::pivot_longer(
    c(bias, sd_beta_hat, rmse),
    names_to = "metric",
    values_to = "value"
  )

ggplot(plot_df, aes(x = scenario_type, y = value, fill = scenario_type)) +
  geom_col(show.legend = FALSE) +
  facet_grid(metric ~ species, scales = "free_y") +
  theme_bw() +
  scale_fill_manual(values = scenario_colors, breaks = scenario_order) +
  labs(
    title = "Bias, RMSE, and SD(beta_hat) Across Scenario Types at total_visits = 800, p_detect = 0.7",
    x = "Scenario type",
    y = NULL
  )

head(merged_scenarios)

#----average-
scenario_order <- c("baseline", "strava", "core.sites", "targeted", "5050", "clustered")

plot_df <- merged_scenarios %>%
  dplyr::mutate(scenario_type = factor(scenario_type, levels = scenario_order)) %>%
  dplyr::group_by(species, scenario_type) %>%
  dplyr::summarise(
    bias = mean(bias, na.rm = TRUE),
    sd_beta_hat = mean(sd_beta_hat, na.rm = TRUE),
    rmse = mean(rmse, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  tidyr::pivot_longer(
    cols = c(bias, sd_beta_hat, rmse),
    names_to = "metric",
    values_to = "value"
  )
plot_df <- plot_df %>%
  dplyr::mutate(
    metric = factor(
      metric,
      levels = c("bias", "rmse", "sd_beta_hat"),
      labels = c("Bias", "RMSE", "SD")
    )
  )


ggplot(plot_df, aes(x = scenario_type, y = value, fill = scenario_type)) +
  geom_col(show.legend = FALSE) +
  facet_grid(metric ~ species, scales = "free_y") +
  theme_bw() +
  scale_fill_manual(values = scenario_colors, breaks = scenario_order) +
  labs(
    title = "Average Bias, RMSE, and SD(estimated corridor effect) Across Scenario Types",
    x = "Scenario type",
    y = NULL
  )

#bias per scenario type per species per detection and effort--------



plot_df <- merged_scenarios %>%
  dplyr::mutate(
    scenario_type = factor(scenario_type, levels = scenario_order)
  ) %>%
  dplyr::group_by(species, scenario_type, total_visits, p_detect) %>%
  dplyr::summarise(bias = mean(bias, na.rm = TRUE), .groups = "drop")

ggplot(plot_df, aes(x = total_visits, y = p_detect, fill = bias)) +
  geom_tile() +
  facet_grid(species ~ scenario_type) +
  theme_bw() +
  scale_fill_gradient2(midpoint = 0) +
  labs(
    title = "Bias across total_visits and p_detect (tiles)",
    x = "Total visits",
    y = "p_detect",
    fill = "Bias"
  )

bias_df <- merged_scenarios %>%
  dplyr::filter(!is.na(bias))

bias_df %>%
  dplyr::summarise(
    n = dplyr::n(),
    min = min(bias),
    q01 = quantile(bias, 0.01),
    q05 = quantile(bias, 0.05),
    q10 = quantile(bias, 0.10),
    q25 = quantile(bias, 0.25),
    median = median(bias),
    q75 = quantile(bias, 0.75),
    q90 = quantile(bias, 0.90),
    q95 = quantile(bias, 0.95),
    q99 = quantile(bias, 0.99),
    max = max(bias)
  )

qs <- c(0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99)
qv <- quantile(merged_scenarios$bias, probs = qs, na.rm = TRUE)

ggplot(merged_scenarios, aes(x = bias)) +
  geom_histogram(bins = 50, color = "white") +
  geom_vline(xintercept = 0, linewidth = 0.4) +
  geom_vline(xintercept = qv, linetype = "dashed", linewidth = 0.3) +
  theme_bw() +
  labs(title = "Bias distribution (all scenarios)", x = "Bias", y = "Count")


library(dplyr)
library(ggplot2)

scenario_order <- c("baseline", "strava", "core.sites", "targeted", "5050", "clustered")

bias_breaks <- c(-Inf, -0.10, -0.05, -0.02, -0.01, 0.01, 0.02, 0.05, 0.10, Inf)
bias_labels <- c("≤ -0.10", "≤ -0.05", "≤ -0.02", "≤ -0.01",
                 "|bias| ≤ 0.01",
                 "≤ 0.02", "≤ 0.05", "≤ 0.10", "> 0.10")

plot_df <- merged_scenarios %>%
  mutate(
    scenario_type = factor(scenario_type, levels = scenario_order),
    total_visits = factor(total_visits),
    p_detect = factor(p_detect)
  ) %>%
  group_by(species, scenario_type, total_visits, p_detect) %>%
  summarise(bias = mean(bias, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    bias_level = cut(bias, breaks = bias_breaks, labels = bias_labels, right = TRUE),
    bias_level = factor(bias_level, levels = bias_labels)  # keep legend ordered
  )

ggplot(plot_df, aes(x = p_detect, y = total_visits, fill = bias_level)) +
  geom_tile(color = "white", linewidth = 0.3) +
  facet_grid(species ~ scenario_type) +
  theme_bw() +
  labs(
    title = "Bias across effort (total_visits) and detection probability",
    x = "Detection probability",
    y = "Total visits (effort)",
    fill = "Bias"
  )


library(dplyr)
library(ggplot2)

scenario_order <- c("baseline", "strava", "core.sites", "targeted", "5050", "clustered")
eps <- 0.02  # threshold for "close to zero"

plot_df <- merged_scenarios %>%
  mutate(
    scenario_type = factor(scenario_type, levels = scenario_order),
    total_visits  = factor(total_visits),
    p_detect      = factor(p_detect)
  ) %>%
  group_by(species, scenario_type, total_visits, p_detect) %>%
  summarise(bias = mean(bias, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    bias_level = case_when(
      bias < -eps ~ "Underestimating",
      bias >  eps ~ "Overestimating",
      TRUE        ~ "Close to 0"
    ),
    bias_level = factor(bias_level, levels = c("Underestimating", "Close to 0", "Overestimating"))
  )


ggplot(plot_df, aes(x = p_detect, y = total_visits, fill = bias_level)) +
  geom_tile(color = "white", linewidth = 0.3) +
  facet_grid(species ~ scenario_type) +
  theme_bw() +
  scale_fill_manual(
    values = c(
      "Underestimating" = "#D7191C",  # red
      "Close to 0"      = "#F7F7F7",  # light/neutral
      "Overestimating"  = "#2C7BB6"   # blue
    )
  ) +
  labs(
    title = paste0("Bias direction across effort and detection (|bias| ≤ ", eps, " treated as ~0)"),
    x = "Detection probability",
    y = "Total visits (effort)",
    fill = "Bias"
  )

ggplot(plot_df, aes(x = p_detect, y = total_visits, fill = bias_level)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.3f", bias)), size = 2) +
  facet_grid(species ~ scenario_type) +
  theme_bw()

scenario_order <- c("baseline", "strava", "core.sites", "targeted", "5050", "clustered")
eps <- 0.02

plot_df <- merged_scenarios %>%
  dplyr::mutate(
    scenario_type = factor(scenario_type, levels = scenario_order),
    total_visits  = factor(total_visits),
    p_detect      = factor(p_detect)
  ) %>%
  dplyr::group_by(species, scenario_type, total_visits, p_detect) %>%
  dplyr::summarise(bias = mean(bias, na.rm = TRUE), .groups = "drop")

# clamp scale to avoid outliers washing out the colors
lims <- stats::quantile(plot_df$bias, probs = c(0.01, 0.99), na.rm = TRUE)
m <- max(abs(lims))
bmin <- -m
bmax <-  m

ggplot2::ggplot(plot_df, ggplot2::aes(x = p_detect, y = total_visits, fill = bias)) +
  ggplot2::geom_tile(color = "white", linewidth = 0.3) +
  ggplot2::facet_grid(species ~ scenario_type) +
  ggplot2::theme_bw() +
  ggplot2::scale_fill_gradientn(
    colours = c(
      "#67000D",  # dark red
      "#FCAE91",  # light red
      "#00A650",  # green near 0
      "#9ECAE1",  # light blue
      "#08306B"   # dark blue
    ),
    values = scales::rescale(c(bmin, -eps, 0, eps, bmax)),
    limits = c(bmin, bmax),
    oob = scales::squish
  ) +
  ggplot2::labs(
    title = paste0("Bias across effort and detection (green if bias ≈ 0)"),
    x = "Detection probability",
    y = "Total visits (effort)",
    fill = "Bias"
  )


## --- bias plot 2
scenario_order <- c("baseline", "strava", "core.sites", "targeted", "5050", "clustered")
eps <- 0.02

plot_df <- merged_scenarios %>%
  dplyr::mutate(
    scenario_type = factor(scenario_type, levels = scenario_order),
    total_visits  = factor(total_visits),
    p_detect      = factor(p_detect)
  ) %>%
  dplyr::group_by(species, scenario_type, total_visits, p_detect) %>%
  dplyr::summarise(bias = mean(bias, na.rm = TRUE), .groups = "drop")

# clamp scale to avoid outliers washing out the colors
lims <- stats::quantile(plot_df$bias, probs = c(0.01, 0.99), na.rm = TRUE)
m <- max(abs(lims))
bmin <- -m
bmax <-  m

ggplot2::ggplot(plot_df, ggplot2::aes(x = p_detect, y = total_visits, fill = bias)) +
  ggplot2::geom_tile(color = "white", linewidth = 0.3) +
  ggplot2::facet_grid(species ~ scenario_type) +
  ggplot2::theme_bw() +
  ggplot2::scale_fill_gradientn(
    colours = c(
      "#67000D",
      "#FCAE91",
      "#F7F7F7",  # zero = white/grey
      "#9ECAE1",
      "#08306B"
    ),
    values = scales::rescale(c(bmin, -eps, 0, eps, bmax)),
    limits = c(bmin, bmax),
    oob = scales::squish
  ) +
  ggplot2::theme(
    strip.background = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(color = "black")
  ) +
  ggplot2::labs(
    title = "Bias across effort and detection (white indicates bias ≈ 0)",
    x = "Detection probability",
    y = "Total visits (effort)",
    fill = "Bias"
  )




#-------RMSE performance across all combinations-----------


merged_scenarios <- merged_scenarios %>%
  mutate(RMSE_level = case_when(
    rmse <= 0.20 ~ "≤ 0.20",
    rmse <= 0.25 ~ "≤ 0.25",
    rmse <= 0.30 ~ "≤ 0.30",
    rmse <= 0.40 ~ "≤ 0.40",
    rmse <= 0.50 ~ "≤ 0.50",
    TRUE         ~ "> 0.50"
  ))

merged_scenarios <- merged_scenarios %>%
  mutate(
    RMSE_level = factor(
      RMSE_level,
      levels = c("≤ 0.20", "≤ 0.25", "≤ 0.30", "≤ 0.40", "≤ 0.50", "> 0.50")  # best → worst
    )
  )

merged_scenarios <- merged_scenarios %>%
  mutate(
    scenario_type = factor(
      scenario_type,
      levels = c("baseline", "strava", "core.sites", "targeted", "5050", "clustered")
    )
  )

merged_scenarios %>%
  ggplot(aes(x = p_detect,
             y = factor(total_visits),
             fill = RMSE_level)) +
  geom_tile(color = "white") +
  facet_grid(species ~ scenario_type) +
  scale_fill_manual(
    values = c(
      "≤ 0.20" = "#006837" , #(medium green)
      "≤ 0.25" = "#8BC34A", #(yellow-green)
      "≤ 0.30" = "#FFEB3B", #(pastel yellow)
      "≤ 0.40" = "#FB8C00", #(orange)
      "≤ 0.50" = "#E53935", #(orange-red)
      "> 0.50" =  "#7F0000" #(deep dark red)
    ),
    name = "RMSE"
  ) +
  labs(
    title = "RMSE performance across effort and detection probability",
    x = "Detection probability",
    y = "Total visits (effort)"
  ) +
  theme_minimal(base_size = 12)


#----Model performance, stability?------

###
# p_rmse summary table (no colors needed here, but keep scenario order)
p_rmse_summary <- merged_scenarios %>%
  group_by(scenario_type) %>%
  summarise(
    mean_p_rmse = mean(p_rmse, na.rm = TRUE),
    sd_p_rmse   = sd(p_rmse, na.rm = TRUE),
    n           = n(),
    .groups = "drop"
  )

print(p_rmse_summary)

# Detection probability RMSE (boxplot) — USE scenario_colors
ggplot(merged_scenarios,
       aes(x = scenario_type, y = p_rmse, fill = scenario_type)) +
  geom_boxplot(outlier.alpha = 0.3) +
  theme_minimal(base_size = 14) +
  scale_fill_manual(values = scenario_colors, breaks = scenario_order, drop = FALSE) +
  labs(
    title = "Distribution of detection probability RMSE across scenario types",
    x = "Scenario type",
    y = "RMSE of detection probability"
  ) +
  theme(legend.position = "none")


# Fail rate (mean) — USE scenario_colors
merged_scenarios %>%
  group_by(scenario_type) %>%
  summarise(mean_fail_rate = mean(fail_rate, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = scenario_type, y = mean_fail_rate, fill = scenario_type)) +
  geom_col(color = "black") +
  theme_minimal(base_size = 14) +
  scale_fill_manual(values = scenario_colors, breaks = scenario_order, drop = FALSE) +
  labs(
    title = "Mean model failure rate across scenario types",
    x = "Scenario type",
    y = "Mean failure rate"
  ) +
  theme(legend.position = "none")


# Error types (stacked) — scenario_colors CANNOT be used here as-is
# because fill is mapped to top_error, not scenario_type.
# What you *can* do is:
#   A) keep scenario order on x (already handled by factor levels), OR
#   B) facet by scenario_type and use a separate palette for top_error.

error_summary <- merged_scenarios %>%
  filter(!is.na(top_error)) %>%
  group_by(scenario_type, top_error) %>%
  summarise(n = n(), .groups = "drop")

# Option A (keep as stacked by error, ordered scenarios on x)
ggplot(error_summary,
       aes(x = scenario_type, y = n, fill = top_error)) +
  geom_col(position = "stack") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Types of fitting errors across scenario types",
    x = "Scenario type",
    y = "Number of errors",
    fill = "Error message"
  )




head(merged_scenarios)

library(dplyr)
library(tidyr)

error_overview_top <- merged_scenarios %>%
  filter(n_failed > 0, !is.na(top_error)) %>%
  count(scenario_type, top_error, wt = n_failed, name = "n_failed") %>%  # sum failures
  arrange(desc(n_failed))

error_overview_top


error_overview_top_wide <- error_overview_top %>%
  tidyr::pivot_wider(
    names_from = scenario_type,
    values_from = n_failed,
    values_fill = 0
  )

error_overview_top_wide

