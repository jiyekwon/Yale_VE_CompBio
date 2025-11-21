library(ggsci)
library(dplyr)
library(zoo)
library(ggplot2)
library(ggpubr)
# load data --------------
data_plot_2a<- readRDS("../Yale_VE_CompBio/Yale_VE_CompBio/data/Fig2A_data_plot_2a.RDS")
yale_drop1 <- readRDS("../Yale_VE_CompBio/Yale_VE_CompBio/data/yale_drop1.RDS")

# global setting --------------
color_map_variants <- c("Alpha" = "#17154FFF", "Delta" = "#1F6E9CFF", 
                        "Other" = "#818181", 
                        "BA.1" ="#BF3729FF", "BA.2" = "#E48171FF", 
                        "BA.4" = "#6C5D9EFF", "BA.5" = "#B57BA2", 
                        "XBB.1" = "#C2A43C", "XBB.2" = "#DEC895",
                        "JN.1" = "#24693D",
                        "KP.1" = "#ABCADF", "KP.2" = "#848FBE", "KP.3" = "#6D2D83")

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





# Figure 2B -----------------------------# 
Spike.version1<- yale_drop1 %>% 
  # filter(week.date<= "2022-09-30") %>% 
  ggplot(aes(x=week.date, y=bnt_distance, 
             col = variants_ft.x, fill = variants_ft.x)) +
  geom_rect(data = emergence_period, inherit.aes = FALSE,
            aes(xmin = date_start, xmax = date_end, ymin = -Inf, ymax = Inf),
            fill = "grey", alpha= 0.3)+
  geom_point(alpha = 0.4) +
  ylab("Amino acid (AA) differences, \n original vaccine formulation") +
  xlab("") +
  ylim(0, 45) + 
  scale_fill_manual(values = color_map_variants, name = "Variants",
                    breaks = names(color_map_variants))+
  scale_color_manual(values = color_map_variants, name = "Variants",
                     breaks = names(color_map_variants))+
  theme(panel.spacing= unit(1,'lines') , 
        axis.text.x=element_text(angle=90))+ 
  theme_classic()+
  geom_hline(yintercept=0, col='gray', lty=2)+
  geom_vline(xintercept=as.Date("2022-09-15"), 
             col='#8B8682', lty=2,linewidth = 0.8)+
  scale_x_date(limits = as.Date(c("2021-03-04", "2023-02-01")),
               date_breaks = "1 month", date_labels =  "%b %Y")+
  theme(axis.text.x=element_text(angle=60, hjust=1),
        legend.position="top")+ #top
  guides(fill = guide_legend(override.aes = list(alpha = 1.0),
                             nrow = 1))+
  # add annotation
  annotate("text", x = as.Date("2022-09-15"), 
           y = max(yale_drop1$bnt_distance, na.rm = TRUE) * 0.25, 
           label = "Vaccine update \n BA.4/BA.5", fontface = "bold", 
           color = "#6C7B8B", angle = 90, hjust = 0.8, vjust = 1.5)


