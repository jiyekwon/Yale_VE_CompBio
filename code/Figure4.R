library(ggplot2)
library(lubridate)
library(patchwork)


# load data --------------
drop_full_version2Vax<- readRDS(file = "../Yale_VE_CompBio/Yale_VE_CompBio/data/Figure4_drop_full_version2Vax.RDS")
weekly_means_ba4_rolling<- readRDS(file ="../Yale_VE_CompBio/Yale_VE_CompBio/data/bivalent_vaccine_eta_predictions.RDS")
weekly_means_xbb<- read.csv("../Yale_VE_CompBio/Yale_VE_CompBio/data/2023_24vax_CT_rolling_avg_3weekGD.csv", row.names = "X")
weekly_means_xbb$week.date<- as.Date(weekly_means_xbb$week.date, format = "%Y-%m-%d")



# Figure 4A -----------------------------# 
# CASES ONLY: Spike distance change against bivalent vaccine update
# the minimum distance to either ancestral or omicron ba4/5 vaccine
p_spike.version2<- drop_full_version2Vax %>%
  
  # Filter to period of interest 
  filter(Collection.date >= "2022-05-15" & Collection.date < "2023-09-30")%>%
  ggplot(aes(x=week.date, y=min_distance, col = variants, fill = variants)) +
  
  # Grey box  
  geom_rect(aes(xmin = as.Date("2022-03-01"), xmax = as.Date("2022-09-15"), 
                ymin = -Inf,  ymax = Inf), fill = "#F0F0F0", alpha = 0.3, colour = NA)+
  geom_point(alpha = 0.5) +
  ylab("Minimum amino acid (AA) substitutions, \n 2022/23 bivalent BA.4/BA.5 vaccine formulation") +
  xlab("Collection date (week)") +
  ylim(0,45)+
  
  # aesthetics 
  scale_fill_manual(values = color_map_variants, name = "Variants",
                    breaks = names(color_map_variants))+
  scale_color_manual(values = color_map_variants, name = "Variants",
                     breaks = names(color_map_variants))+
  theme(panel.spacing= unit(1,'lines'),axis.text.x=element_text(angle=90))+ 
  theme_classic()+
  geom_hline(yintercept=0, col='gray', lty=2)+
  geom_vline(xintercept=as.Date("2022-09-15"), col='#8B8682', lty=2,linewidth = 0.8)+
  geom_vline(xintercept=as.Date("2023-09-15"), col='#8B8682', lty=2,linewidth = 0.8)+
  scale_x_date(limits = as.Date(c("2022-03-01", "2023-10-01")),
               date_breaks = "1 month", date_labels =  "%b %Y")+
  theme(axis.text.x=element_text(angle=60, hjust=1), legend.position="top")+
  guides(fill = guide_legend(override.aes = list(alpha = 1.0), nrow = 1))+
  
  # add annotation
  annotate("text", x = as.Date("2022-09-15"), 
           y = max(gisaid_aug$bnt_distance_spike_vax2, na.rm = TRUE), 
          fontface = "bold", label = "Vaccine update \n BA.4/BA.5", 
           color = "#6C7B8B", angle = 90, hjust = 0.6, vjust = 1.5)+
  annotate("text", x = as.Date("2023-09-15"), 
           y = max(gisaid_aug$bnt_distance_spike_vax2, na.rm = TRUE), 
           label = "Vaccine update \n XBB 1.5", fontface = "bold",
           color = "#6C7B8B", angle = 90, hjust = 0.6, vjust = 1.5)

# Overlaying % VE reduction line 
p_spike.version2.2 <- p_spike.version2 +
  geom_line( data = weekly_means_ba4_rolling, inherit.aes = FALSE,
             aes(y=ba4_mean_eta_spike_100, x = week.date), size=1, 
             color = "indianred3", alpha = 0.5) +
  
  scale_y_continuous(
    # Add a second axis and specify its features
    sec.axis = sec_axis( ~.,name="Predicted % VE reduction",  
                         breaks = c( 0, 20, 40, 60)),
    limits = c(0, 63),
    breaks = c( 0, 20, 40, 60))
# ggsave(plot = p_spike.version2.2, filename = "./24draft_figures/Figure2.3_Spike_BA4_5_Formula_wReduction_lim65.jpeg",
#        device = "jpeg", width = 8, height = 6, dpi = 250)





# Figure 4B -----------------------------# 
# CASES ONLY: Spike distance change against XBB vaccine update
p_spike.version3<- gisaid_aug %>%
  filter(Collection.date >= "2023-05-15")%>%
  ggplot(aes(x=week.date, y=bnt_distance_spike_vax3, 
             col = variant_ft_short, fill = variant_ft_short)) +
  geom_rect(aes(xmin = as.Date("2023-03-01"), 
                xmax = as.Date("2023-09-15"), 
                ymin = -Inf, 
                ymax = Inf), 
            fill = "#F0F0F0", alpha = 0.3, colour = NA)+
  geom_point(alpha = 0.5) +
  ylab("Amino acid (AA) substitutions, \n 2023/24 XBB 1.5 vaccine formulation") +
  xlab("Collection date (week)") +
  ylim(0,45)+
  scale_fill_manual(values = color_map_variants, name = "Variants",
                    breaks = names(color_map_variants))+
  scale_color_manual(values = color_map_variants, name = "Variants",
                     breaks = names(color_map_variants))+
  theme(panel.spacing= unit(1,'lines') , 
        axis.text.x=element_text(angle=90))+ 
  theme_classic()+
  geom_hline(yintercept=0, col='gray', lty=2)+
  geom_vline(xintercept=as.Date("2023-09-15"), 
             col='#8B8682', lty=2,linewidth = 0.8)+
  scale_x_date(limits = as.Date(c("2023-03-01", "2024-08-01")),
               date_breaks = "1 month", date_labels =  "%b %Y")+
  theme(axis.text.x=element_text(angle=60, hjust=1),
        legend.position="top")+
  guides(fill = guide_legend(override.aes = list(alpha = 1.0),
                             nrow = 1))+
  # add annotation
  annotate("text", x = as.Date("2023-09-15"), 
           y = max(gisaid_aug$bnt_distance_spike_vax3, na.rm = TRUE), 
           label = "Vaccine update \n XBB 1.5", fontface = "bold",  
           color = "#6C7B8B", angle = 90, hjust = 0.6, vjust = 1.5)
p_spike.version3


# Overlaying % VE reduction line 
p_spike.version3.2 <- p_spike.version3 +
  geom_line( data = weekly_means_xbb, inherit.aes = FALSE,
             aes(y=xbb_mean_eta_spike_100, x = week.date), size=1, 
             color = "indianred3", alpha = 0.5) +
  
  scale_y_continuous(
    # Add a second axis and specify its features
    sec.axis = sec_axis( ~.,name="Predicted % VE reduction",  
                         breaks = c( 0, 20, 40, 60) ),
    limits = c(0, 63),
    breaks = c( 0, 20, 40, 60))

p_spike.version2.2/p_spike.version3.2
# ggsave(plot = p_spike.version2.2, filename = "./24draft_figures/Figure2.3_Spike_BA4_5_Formula_wReduction_lim65.jpeg",
#        device = "jpeg", width = 8, height = 6, dpi = 250)


