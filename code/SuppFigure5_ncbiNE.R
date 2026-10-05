library(seqinr)
library(stringr)
library(dplyr)
library(lubridate)
library(ggplot2)
library(ggsci)
library(stringr)
library(data.table)
library(tidyr)
library(patchwork)

# load data --------------
ne_AA <- seqinr::read.alignment(file ="../Yale_VE_CompBio/Yale_VE_CompBio/data/NE.aligned_S_aa.fasta", format="fasta")
ne_meta <- read.csv("../Yale_VE_CompBio/Yale_VE_CompBio/data/NE_metadata_ncbi251103.csv")
ne_clade <- read.table("../Yale_CDCBAA/NCBI_NE_test/output/nextclade.tsv", sep = "\t", header = 1)


# names of the sequences (49269 rows )
name.list<- as.data.frame(ne_AA$nam)
rownames(name.list)<- gsub(":.*$","",name.list$`ne_AA$nam`)
name.list$ID<- rownames(name.list)


# setting names for the sequences
seq.list<- as.list(ne_AA$seq)
names(seq.list)<- name.list$ID

# transform to one letter one cell format
seq.list2<- lapply(seq.list, function(x){ str_split_fixed(x, "", max(nchar(x)))})

# Convert the list to a dataframe
seq.list.df <- as.data.frame(do.call(rbind, seq.list2))
rownames(seq.list.df)<- names(seq.list2)


# Quality control: count Xs and -s across sites 
x_c<- colSums(seq.list.df =="x")
gaps_c<- colSums(seq.list.df =="-")
xs_c<- colSums(seq.list.df =="x" | seq.list.df =="-")
qc_data.column<- as.data.frame(cbind(x_c, gaps_c, xs_c))


# back to the data set 
spike_hamming<- seq.list.df
spike_hamming<- t(spike_hamming)
dim(spike_hamming) # 1274 * 49269 (until 49266 cases)
spike_hamming<- as.data.frame(spike_hamming)


# Northeast sequences vs third vaccine 
# if it is an "X" -- not considered as a difference 
# if it is a gap it is a difference 
hamming.sim = apply(spike_hamming [,c(1:49266)], 2, function(x) 
  case_when(x == "x" ~ 1,
            x != "x" ~ as.integer(str_equal(x, spike_hamming$`BNT162b2_S_XBB.1.5`, ignore_case = TRUE))))

# counting FALSEs (not equal to the vaccine sequence)
hamming.diff.NE<- as.data.frame(colSums(hamming.sim==0))
names(hamming.diff.NE)<- "bnt_distance_spike_vax3"
hamming.diff.NE$Accession <- rownames(hamming.diff.NE)


# merge with the metadata 
cases_NE= left_join(ne_meta, hamming.diff.NE)
cases_NE$Collection_Date<- as.Date(cases_NE$Collection_Date,format = "%m/%d/%y")
cases_NE$week.date<- floor_date(cases_NE$Collection_Date, "week")
  
clade<- ne_clade %>% select(seqName, clade_display) %>% rename(Accession = seqName)
cases_NE <- left_join(cases_NE, clade)



# Plot -- using nextclade clade assignment
colors<- c(colorRampPalette(c(c("#D0641EFF", "#E68E54FF", "#F9AB0EFF", "#9E4058FF", "#C2697FFF",  "#FBC559FF", "#589E40FF", "#7FC269FF", "#2C3778FF", "#4151B0FF", "#513965FF", "#785596FF")))( 26))

v4<- cases_NE %>%
  ggplot(aes(x=week.date, y=bnt_distance_spike_vax3,  fill = clade_display, col = clade_display)) +
  geom_rect(aes(xmin = as.Date("2023-03-01"), 
                xmax = as.Date("2023-09-15"), 
                ymin = -Inf, 
                ymax = Inf), 
            fill = "#F0F0F0", alpha = 0.3, colour = NA)+
  geom_point(alpha = 0.6) +
  ylab("Amino acid substitutions, \n 2023-2024 XBB 1.5 vaccine formulation") +
  xlab("Collection date (week)") +
  scale_fill_manual(values = colors,  "Nextclade \n assignment")+ 
  scale_color_manual(values = colors,  "Nextclade \n assignment")+ 
  theme(panel.spacing= unit(1,'lines') , axis.text.x=element_text(angle=90))+ 
  theme_classic()+
  geom_hline(yintercept=0, col='gray', lty=2)+
  geom_vline(xintercept=as.Date("2023-09-15"), 
             col='#8B8682', lty=2,linewidth = 0.8)+
  scale_x_date(limits = as.Date(c("2023-03-01", "2024-09-01")),
               date_breaks = "1 month", date_labels =  "%b %Y")+
  theme(axis.text.x=element_text(angle=60, hjust=1), legend.position="top",
        legend.spacing.x = unit(0.5, "cm"))+
  guides(col = guide_legend(override.aes = list(alpha = 1.0), ncol = 8, direction = "horizontal" ),
         fill = guide_legend(override.aes = list(alpha = 1.0), ncol = 8, direction = "horizontal" ))+
  ylim(0,60) 

v4


###################################################
# per week summary (3 week rolling average)
################################################### 
ct_mean<- read.csv("./2023_24vax_CT_rolling_avg_3weekGD.csv")

# 3-week rolling average
gd_mean_NE<- cases_NE%>%
  group_by(week.date) %>% 
  summarise( avg_distance = mean(bnt_distance_spike_vax3, na.rm = TRUE), .groups = 'drop') %>%
  as.data.table() %>%
  .[, `New England` := frollmean(avg_distance, n = 3, fill = NA)] 


gd_mean_CT<- ct_mean %>% select(week.date, avg_distance_rolling, avg_distance) %>% 
  rename(avg_CT = avg_distance,
         Yale = avg_distance_rolling)
# %>%mutate(location = "Yale")
gd_mean_CT$week.date<- as.Date(gd_mean_CT$week.date, "%Y-%m-%d")
gd_combined<- full_join(gd_mean_NE, gd_mean_CT)


df_long <- gd_combined %>%
  select(-c(avg_distance, avg_CT))%>%
  pivot_longer( cols = `New England`:Yale, names_to = "location",
    values_to = "avg_distance_rolling")


v_gd<- ggplot(df_long, aes(x = week.date, y = avg_distance_rolling, color = location)) +
  geom_rect(aes(xmin = as.Date("2023-03-01"), 
                xmax = as.Date("2023-09-15"), 
                ymin = -Inf, 
                ymax = Inf), 
            fill = "#F0F0F0", alpha = 0.3, colour = NA)+
  geom_line(size = 1) +
  geom_point(alpha = 0.6) + # Optional: add points for data clarity
  scale_color_manual(values = c("New England" = "#9E4058FF", "Yale" = "#4151B0FF")) +
  labs(title = "",
       y = "Amino acid substitutions, \n 2023-2024 XBB 1.5 vaccine formulation \n (3-week rolling average)",
       color = "Location") +
  ylim(0,60)+
  theme_classic()+
  geom_hline(yintercept=0, col='gray', lty=2)+
  geom_vline(xintercept=as.Date("2023-09-15"), 
             col='#8B8682', lty=2,linewidth = 0.8)+
  scale_x_date(limits = as.Date(c("2023-03-01", "2024-09-01")),
               date_breaks = "1 month", date_labels =  "%b %Y")+
  theme(axis.text.x=element_text(angle=60, hjust=1), legend.position="top",
        legend.spacing.x = unit(0.5, "cm"))


ne_plot<- v4 +v_gd + plot_annotation(tag_levels = "A")



###################################################
# State level 
################################################### 

# 3-week rolling average (state level)
gd_mean_states <- cases_NE %>%
  group_by(week.date, USA) %>% # Group by both week and state
  summarise( avg_distance = mean(bnt_distance_spike_vax3, na.rm = TRUE),
             .groups = 'drop') %>%
  ungroup()

gd_mean_states<- gd_mean_states %>% 
  group_by(USA)%>% 
  mutate( avg_distance_rolling = data.table::frollmean(avg_distance, n = 3, fill = NA) ) %>%
  ungroup() %>%
  rename(avg_state = avg_distance,
         avg_distance_rolling = avg_distance_rolling, 
         location = USA ) %>%
  select(week.date, location, avg_distance_rolling)

# combine with Yale info
gd_mean_CT2<- ct_mean %>% select(week.date, avg_distance_rolling) %>% 
  rename(
    avg_distance_rolling = avg_distance_rolling ) %>%
  mutate(location = "Yale")
gd_mean_CT2$week.date<- as.Date(gd_mean_CT2$week.date, "%Y-%m-%d")
gd_combined2<- full_join(gd_mean_states, gd_mean_CT2)

# making sure each week is represented, if missing then NA 
gd_combined2<- gd_combined2 %>% group_by(location) %>% 
  tidyr::complete(week.date=seq.Date(min(week.date, na.rm=T), 
                                     max(week.date, na.rm=T), 'week'), 
                fill=list(gd_combined2=NA))


# Separate the Yale data (the baseline line)
yale_data <- gd_combined2 %>%
  filter(location == "Yale") %>%
  rename(yale_distance = avg_distance_rolling) %>%
  ungroup() %>%
  dplyr::select(week.date, yale_distance)

# all non-Yale locations)
facet_data <- gd_combined2 %>%
  filter(location != "Yale")
df_facet_compare <- left_join(facet_data, yale_data, by = "week.date")


# facet label 
n_data <- cases_NE %>%
  group_by(USA) %>%
  summarise(n = n(), .groups = 'drop') %>%
  rename(location = USA) %>%
  mutate(facet_label = paste0(location, " (n=", n, ")"))


df_facet_labeled <- left_join(df_facet_compare, n_data, by = "location") %>%
  mutate(facet_label = factor(facet_label, levels = unique(facet_label[order(location)])))




# Plot C -----------------------# 

v_facet_gd_n <- ggplot(df_facet_labeled, aes(x = as.Date(week.date))) +
  geom_rect(aes(xmin = as.Date("2023-03-01"),
                xmax = as.Date("2023-09-15"),
                ymin = -Inf,
                ymax = Inf),
            fill = "#F0F0F0", alpha = 0.3, colour = NA) +
  
  geom_line(aes(y = yale_distance, color = "Yale"),
            size = 1, alpha = 1, linetype = "solid") +

  geom_line(aes(y = avg_distance_rolling, color = location),
            size = 1.2) +
  geom_point(aes(y = avg_distance_rolling, color = location),
             alpha = 0.8) +
  
  facet_wrap(~ facet_label, ncol = 3) +  

  scale_color_manual(
    name = "Location",
    values = c("Yale" = "#4151B0FF",
               "CT" = "#98768EFF","MA" =  "#B08BA5FF", "ME" = "#C7A2B6FF",
               "NH" =  "#DFBBC8FF", "NJ" = "#FFC680FF", "NY" =  "#FFB178FF",
               "PA"= "#DB8872FF", "RI" = "#A56457FF",  "VT" = "#C2CAE3FF") ) +
  labs( y = "Amino acid substitutions, \n 2023-2024 XBB 1.5 vaccine formulation \n (3-week rolling average)",
    x = "Collection Date") +
  ylim(0, 60) +
  theme_classic(base_size = 15)+
  geom_hline(yintercept=0, col='gray', lty=2)+
  geom_vline(xintercept=as.Date("2023-09-15"), 
             col='#8B8682', lty=2,linewidth = 0.8)+
  scale_x_date(limits = as.Date(c("2023-03-01", "2024-09-01")),
               date_breaks = "2 months", date_labels =  "%b %Y")+
  theme(axis.text.x=element_text(angle=60, hjust=1), legend.position="bottom",
        strip.background = element_rect(fill = "gray95"), 
        legend.spacing.x = unit(0.5, "cm"))

ne_final<- ne_plot/v_facet_gd_n+ plot_annotation(tag_levels = "A") +plot_layout(heights = c(2,3))

ggsave(plot = ne_final, filename = "./24draft_figures/RevisionFigure_NE_XBBFormula.jpeg",
       device = "jpeg", width = 16, height = 20, dpi = 250)
