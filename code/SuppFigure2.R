library(seqinr)
library(stringr)
library(dplyr)
library(lubridate)
library(ggplot2)
library(ggsci)
library(stringr)
library(data.table)
library(tidyr)


# load data --------------
yale_drop1 <- readRDS("../Yale_VE_CompBio/Yale_VE_CompBio/data/yale_drop1.RDS")

# colors --------------
emergence_period <- 
  data.frame(date_start=as.Date(c("2021-06-03", "2021-11-30", "2022-02-23", "2022-06-09")),
             date_end = as.Date(c("2021-07-08","2022-01-04", "2022-03-30", "2022-07-14")))

color_map_variants <- c("Alpha" = "#17154FFF", "Delta" = "#1F6E9CFF", 
                        "Other" = "#818181", 
                        "BA.1" ="#BF3729FF", "BA.2" = "#E48171FF", 
                        "BA.4" = "#6C5D9EFF", "BA.5" = "#B57BA2", 
                        "XBB.1" = "#C2A43C", "XBB.2" = "#DEC895",
                        "JN.1" = "#24693D",
                        "KP.1" = "#ABCADF", "KP.2" = "#848FBE", "KP.3" = "#6D2D83")


#############################
# Figures 
# A. Correlation plots  ((1) spike AA vs  spike NT & ) 
# B. Temporal distance plots (spike NT, spike AA, RBD AA)
#############################

# Correlation plots --------------------
plot_nt_corr<- yale_drop1 %>% ggplot(., aes( x = ref_distance_spike_nt1, y = bnt_distance)) +
  geom_point(aes( color = variants_ft.x, fill= variants_ft.x), alpha = 0.4) +
  scale_fill_manual(values = color_map_variants, name = "Variants",
                    breaks = names(color_map_variants))+
  scale_color_manual(values = color_map_variants, name = "Variants",
                     breaks = names(color_map_variants))+
  ylim(0, 50)+
  # xlim(0,75)+
  theme_minimal()+ 
  labs( x = "Nucleotide differences to reference sequence (MN908947)*",
        y = "Amino acid (AA) differences, \n original vaccine formulation", 
        caption = "*'-' or nucleotide difference; spearman correlation") +
  theme_minimal()+
  stat_cor( aes(label = after_stat(rr.label)), 
            method = "spearman",  r.digits = 3, label.x.npc = 'left',  
            label.y.npc = 'top', geom = "text")


plot_rbd_corr<- yale_drop1 %>% ggplot(., aes( x = bnt_distance_rbd, y = bnt_distance)) +
  geom_point(aes( color = variants_ft.x, fill= variants_ft.x), alpha = 0.4) +
  scale_fill_manual(values = color_map_variants, name = "Variants",
                    breaks = names(color_map_variants))+
  scale_color_manual(values = color_map_variants, name = "Variants",
                     breaks = names(color_map_variants))+
  ylim(0, 50)+
  theme_minimal()+ 
  labs( x = "RBD differences, \n original vaccine formulation",
        y = "Amino acid (AA) differences, \n original vaccine formulation", 
        caption = "") +
  theme_minimal()+
  stat_cor( aes(label = after_stat(rr.label)), 
            method = "spearman",  r.digits = 3, label.x.npc = 'left',  
            label.y.npc = 'top', geom = "text")

# Distance plots --------------------
Spike_nt.version1<- yale_drop1 %>% 
  # filter(week.date<= "2022-09-30") %>% 
  ggplot(aes(x=week.date, y=ref_distance_spike_nt1, 
             col = variants_ft.x, fill = variants_ft.x)) +
  geom_rect(data = emergence_period, inherit.aes = FALSE,
            aes(xmin = date_start, xmax = date_end, ymin = -Inf, ymax = Inf),
            fill = "grey", alpha= 0.3)+
  geom_point(alpha = 0.4) +
  ylab("Nucleotide differences, \n reference sequence (MN908947)*") +
  xlab("") +
  # ylim(0, 45) + 
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


RBD.version1<- yale_drop1 %>% 
  # filter(week.date<= "2022-09-30") %>% 
  ggplot(aes(x=week.date, y=bnt_distance_rbd, 
             col = variants_ft.x, fill = variants_ft.x)) +
  geom_rect(data = emergence_period, inherit.aes = FALSE,
            aes(xmin = date_start, xmax = date_end, ymin = -Inf, ymax = Inf),
            fill = "grey", alpha= 0.3)+
  geom_point(alpha = 0.4) +
  ylab("Amino acid (AA) differences, \n original vaccine formulation") +
  xlab("Collection date (week)") +
  ylim(0,25)+
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
  scale_x_date(limits = as.Date(c("2021-03-01", "2023-02-01")),
               date_breaks = "1 month", date_labels =  "%b %Y")+
  theme(axis.text.x=element_text(angle=60, hjust=1),
        legend.position="top")+
  guides(fill = guide_legend(override.aes = list(alpha = 1.0),
                             nrow = 1))+
  # add annotation
  annotate("text", x = as.Date("2022-09-15"), 
           y = max(yale_drop1$bnt_distance_rbd, na.rm = TRUE) * 0.25, 
           label = "Vaccine update\nBA.4/BA.5",
           fontface = "bold", 
           color = "#6C7B8B", angle = 90, hjust = 0.8, vjust = 1.5)



# Plot collection -------------------------------------#
distance_plot<- Spike_nt.version1/ Spike.version1 / RBD.version1
A_nucleotide<- (Spike_nt.version1+
                  theme(legend.position = "none", plot.tag = element_text(face = "bold")) +
                  plot_nt_corr+
                  theme(legend.position = "none", plot.tag = element_text(face = "bold")))
B_RBD<- RBD.version1+
  theme(plot.tag = element_text(face = "bold"))+ 
  plot_rbd_corr+
  theme(legend.position = "none", plot.tag = element_text(face = "bold"))
distance_plot<- B_RBD/A_nucleotide+ plot_annotation(tag_levels = "A")
distance_plot
# ggsave(plot = distance_plot, filename = "./24draft_figures/RevisionFigure_RBD_nt_ogFormula.jpeg", device = "jpeg", width = 12, height = 9, dpi = 200)
