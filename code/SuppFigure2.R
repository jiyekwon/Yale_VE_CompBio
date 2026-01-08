# ================================================
# INSTRUCTIONS: run this after Figure2.R file 
# ================================================


# Supp.Figure S2 -----------------------------# 
library(slider)

df_processed <- yale_drop1 %>%
  dplyr::arrange(variants_ft.x, week.date) %>%
  group_by(variants_ft.x) %>%
  
  mutate(rolling_sd = slide_index_dbl(
    .x = bnt_distance,
    .i = week.date,
    .f = ~sd(.x, na.rm = TRUE),
    .before = days(30))) %>% ungroup()

# Plot 
sd_run_plot<- ggplot(df_processed, 
                     aes(x = week.date, y = rolling_sd, color = variants_ft.x)) +
  
  geom_line(alpha = 0.8) + 
  # geom_smooth(method = "loess", span = 0.3, color = "black", size = 0.5) +
  scale_color_manual(values = color_map_variants, name = "Variants",
                     breaks = names(color_map_variants))+
  facet_wrap(~variants_ft.x, scales = "free_x") +
  theme_minimal(base_size = 14) +
  scale_x_date(date_breaks = "2 month", date_labels =  "%b %Y")+
  theme(axis.text.x=element_text(angle=60, hjust=1),
        legend.position="top")+
  guides(fill = guide_legend(override.aes = list(alpha = 1.0),
                             nrow = 1))+
  labs(
    y = "Rolling SD of AA substitutions",
    x = "Date")

ggsave(plot = sd_run_plot, filename = "./24draft_figures/RevisionFigure_sd_run.jpeg",
       device = "jpeg", width = 13, height = 8, dpi = 250)



# Supp.Figure S2 Alternative (coefficient of variance version)----------------# 
library(tidyverse)
library(zoo)

df_processed <- yale_drop1 %>%
  dplyr::arrange(variants_ft.x, week.date) %>%
  group_by(variants_ft.x) %>%
  mutate(
    # 30-day rolling SD
    rolling_sd = slide_index_dbl(
      .x = bnt_distance,
      .i = week.date,
      .f = ~sd(.x, na.rm = TRUE),
      .before = days(30)),
    
    # 30-day rolling Mean
    rolling_mean = slide_index_dbl(
      .x = bnt_distance,
      .i = week.date,
      .f = ~mean(.x, na.rm = TRUE),
      .before = days(30)),
    
    # Calculate CV (SD / Mean)
    rolling_cv = (rolling_sd / rolling_mean) * 100
  ) %>% 
  ungroup()

# Plotting the CV instead of just SD
cv_run_plot <- ggplot(df_processed, 
                      aes(x = week.date, y = rolling_cv, color = variants_ft.x)) +
  geom_line(alpha = 0.4, size = 1) + 
  scale_color_manual(values = color_map_variants, name = "Variants") +
  theme_minimal(base_size = 14) +
  scale_x_date(date_breaks = "2 month", date_labels = "%b %Y") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        legend.position = "top") +
  labs(
    y = "Rolling CV of AA substitutions (%)",
    x = "Date",
    subtitle = "Relative variation (SD/Mean) using a 30-day window")
