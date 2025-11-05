library(ggsci)
library(dplyr)
library(zoo)
library(ggplot2)
library(ggpubr)
# load data --------------
data_plot_2a<- readRDS("../Yale_VE_CompBio/Yale_VE_CompBio/data/Fig2A_data_plot_2a.RDS")


# global setting --------------
emergence_period <- 
  data.frame(date_start=as.Date(c("2021-06-03", "2021-11-30", "2022-02-23", "2022-06-09")),
             date_end = as.Date(c("2021-07-08","2022-01-04", "2022-03-30", "2022-07-14")))


emergence_period$date_start <- floor_date(emergence_period$date_start, unit='week')
emergence_period$date_end <- floor_date(emergence_period$date_end, unit='week')
reference_group <- "unvaccinated"

case_summary_timesince <- data_plot_2a %>%
  #filter(week.date <= "2022-09-30") %>%
  group_by(week.date, time_since_vax_factor) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(week.date) %>%
  mutate(percent = count / sum(count) * 100) %>%
  ungroup()


# Panel 2: 14-90days 
short_immunity<- case_summary_timesince %>% filter(time_since_vax_factor %in% c("14-90days"))
short_immunity<- short_immunity %>%
  tidyr::complete(week.date=seq.Date(min(week.date, na.rm=T), max(week.date, na.rm=T), 'week'),
                  fill = list(time_since_vax_factor = "14-90days", count = 0, percent = NA))

short_immunity <- short_immunity %>%
  arrange(week.date) %>%
  mutate(rolling_vax_prop = rollmean(percent, k = 3, fill = NA, align = "center")) %>%
  ungroup()


# Vaccinated groups  
plotting_data_for_facet <- data_plot_2a %>%
  filter(time_since_vax_factor != reference_group)

# Unvaccinated 
reference_data <- data_plot_2a %>%
  filter(time_since_vax_factor == reference_group) %>%
  rename(reference_bnt_distance = bnt_distance)  %>%
  select(-time_since_vax_factor) # Keep only necessary columns

# Subset to group of interest: 14-90days
select_data<- plotting_data_for_facet %>% 
  filter(time_since_vax_factor == "14-90days")


#---------------------------------#
# Figure 
#---------------------------------#
max_percent<- 100
vax_group_comparison_distance<- 
  ggplot(reference_data, aes(x = week.date, y =reference_bnt_distance)) +
  geom_rect(data = emergence_period, inherit.aes = FALSE,
            aes(xmin = date_start, xmax = date_end, ymin = -Inf, ymax = Inf),
            fill = "lightgrey", alpha= 0.3)+
  
  # Unvaccinated and vaccinated(14-90 days in Red)
  geom_point( col = "#F3AA4F14", fill = "#F3AA4F33", alpha = 0.2)+
  geom_point(data = select_data, 
             aes(y = bnt_distance, color = time_since_vax_factor), col = "#AD002AFF",alpha = 0.4) +
  
  # Vaccine update
  geom_vline(xintercept=as.Date("2022-09-15"), 
             col='#8B8682', lty=2,linewidth = 0.8)+
  
  
  # Rolling average of vaccinated controls as line and points
  geom_line(data = short_immunity,
            aes(x = week.date, y = rolling_vax_prop * max_percent / 100), 
            linetype = "solid",
            size = 0.8, alpha = 0.5, col ="#AD002AFF")+
  
  ggpubr::theme_pubr() +
  facet_wrap(~time_since_vax_factor, scales = "fixed", ncol = 1) + 
  ylab("Amino acid (AA) substitutions, \n original vaccine formulation") +
  xlab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_y_continuous(sec.axis = sec_axis(~ . * 100 / max_percent,
                                         name = "3-Week Rolling Avg.\n Vaccinated (%)"))+
  scale_x_date(limits = as.Date(c("2021-03-01", "2023-02-01")),
               date_breaks = "1 month", date_labels = "%b %Y") 



