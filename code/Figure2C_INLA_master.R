library(INLA)
library(fastDummies)
library(ggplot2)
library(dplyr)
library(tidyr)
library(lubridate)
ds1<- readRDS(file = "../Yale_VE_CompBio/Yale_VE_CompBio/data/EXPORT_inla_ds1.RDS") 
# trimmed_list<- readRDS(file = "../Yale_VE_CompBio/Yale_VE_CompBio/data/beta_export_240908.RDS") # AR1 version
# trimmed_list<- readRDS(file = "../Yale_VE_CompBio/Yale_VE_CompBio/data/beta_export_rw2_251121.RDS")# RW2 version

#######################################
# INLA master code
# NEVER TAKE THE FUNCTIONS OF MEANS; 
#######################################

######################################################
# Model Selection  ----------------------------------#
######################################################

# mod4: adding routine_flag term (best model)
# add constr = T to all of the ARs 
form4  <- as.formula(case_status ~ 1 + 
                       f(time_point_week, model = "ar1", hyper = list(prec = list(prior = "loggamma", param = c(3, 2))), constr = T) +
                       f(weekN1, vax_cat1, model = "ar1", hyper = list(prec = list(prior = "loggamma", param = c(3, 2))), constr = T) + 
                       f(weekN2, vax_cat2, model = "ar1", hyper = list(prec = list(prior = "loggamma", param = c(3, 2))), constr = T) +
                       f(weekN3, vax_cat3, model = "ar1", hyper = list(prec = list(prior = "loggamma", param = c(3, 2))), constr = T) +
                       f(weekN4, vax_cat4, model = "ar1", hyper = list(prec = list(prior = "loggamma", param = c(3, 2))), constr = T) +
                       
                       vax_cat1+
                       vax_cat2+
                       vax_cat3+
                       vax_cat4+
                       
                       AgeGroup+
                       income_std+
                       gender+
                       prev_vax_dose+
                       routine_flag)

# summary of vax_cat1 & 2 over time 
# fixed effect + random effects --> vaccine effectiveness time varying 

mod4 <- inla(form4,
             family = "binomial",
             data = ds1,
             control.fixed = list(
               mean = 0,
               prec = 1,
               mean.intercept=0,
               prec.intercept=1), 
             control.compute=list( waic=TRUE, config = TRUE),
             control.predictor=list(compute=FALSE)
)
summary(mod4) #  WAIC: 42986.77

# mod.r: revision, diff model rw2 (second-order random walk (RW2) model) instead of ar1 ----------------------------------# 
form.r  <- as.formula(case_status ~ 1 + 
                       f(time_point_week, model = "rw2",  constr = T) +
                       f(weekN1, vax_cat1, model = "rw2",  constr = T) + 
                       f(weekN2, vax_cat2, model = "rw2",  constr = T) +
                       f(weekN3, vax_cat3, model = "rw2",  constr = T) +
                       f(weekN4, vax_cat4, model = "rw2",  constr = T) +
                       
                       vax_cat1+
                       vax_cat2+
                       vax_cat3+
                       vax_cat4+
                       
                       AgeGroup+
                       income_std+
                       gender+
                       prev_vax_dose+
                       routine_flag)

# summary of vax_cat1 & 2 over time 
# fixed effect + random effects --> vaccine effectiveness time varying 

mod.r <- inla(form.r,
             family = "binomial",
             data = ds1,
             control.fixed = list(
               mean = 0,
               prec = 1,
               mean.intercept=0,
               prec.intercept=1), 
             control.compute=list( waic=TRUE, config = TRUE),
             control.predictor=list(compute=FALSE)
)
summary(mod.r) #  WAIC: 42882.98

######################################################
# INLA sampling 
######################################################
run_inla_sampling <- function(model, n_samples, num_cores, fixed_effects_names) {
  # Inner function to sample and process each sample in parallel
  
  sample_function <- function(n) {
    # Draw a single posterior sample
    posterior_sample <- inla.posterior.sample(1, model)
    
    # Extract the parameter sample and convert to data frame
    parameter_sample <- as.data.frame(posterior_sample[[1]]$latent)
    
    # Remove predictor rows
    rows_to_keep <- !grepl("^Predictor", rownames(parameter_sample))
    parameter_sample <- as.data.frame(parameter_sample[rows_to_keep, , drop = FALSE],  nm = "sim")
    
    # Add variable names
    parameter_sample$variable <- rownames(parameter_sample)
    
    # Separate fixed and random effects
    is_fixed_effect <- parameter_sample$variable %in% fixed_effects_names
    fixed_effects <- parameter_sample[is_fixed_effect, , drop = FALSE]
    random_effects <- parameter_sample[!is_fixed_effect, , drop = FALSE]
    
    # Return the list of fixed and random effects
    list(fixed_effects = fixed_effects, random_effects = random_effects)
  }
  
  # Draw and process posterior samples in parallel
  samples <- pbmclapply(1:n_samples, FUN = sample_function, mc.cores = num_cores)
  
  # Separate the samples into fixed and random effects lists
  fixed_effects_list <- lapply(samples, `[[`, "fixed_effects")
  random_effects_list <- lapply(samples, `[[`, "random_effects")
  
  return(list(fixed_effects = fixed_effects_list, random_effects = random_effects_list))
}


# Parameters for sampling
n_samples <- 5000
num_cores <- detectCores() - 1  # Adjust as needed'


fixed_effects_names <-  c("(Intercept):1", "vax_cat11:1", "vax_cat21:1", "vax_cat31:1",
                          "vax_cat41:1", "AgeGroup[40,55):1", "AgeGroup[55,65):1", "AgeGroup[65,75):1",
                          "AgeGroup[75,150]:1", "income_std:1", "genderMALE:1", "genderUNKNOWN:1",
                          "prev_vax_dose:1", "routine_flag1:1")


# Run the nested function
results <- run_inla_sampling(mod.r, n_samples, num_cores, fixed_effects_names) #mod4

# Accessing the fixed and random effects samples
fixed_effects_samples <- results$fixed_effects
random_effects_samples <- results$random_effects


# Fixed effects --------------------------------------------------------------#
# Convert each element of the list to a data frame using lapply
df_list <- lapply(seq_along(fixed_effects_samples), function(i) {
  data.frame(
    Simulation = i,
    FixedEffects = fixed_effects_samples[[i]]
  )
})

# Combine all data frames into one --> long format 
combined_df <- dplyr::bind_rows(df_list)

# Pivoting the data to wide format
fixed_wide_df <- combined_df %>%
  pivot_wider(
    names_from = FixedEffects.variable,
    values_from = FixedEffects.V1
  )


# Random effects --------------------------------------------------------------#
df_list_R <- lapply(seq_along(random_effects_samples), function(i) {
  data.frame(
    Simulation = i,
    RandomEffects = random_effects_samples[[i]]
  )
})

# Combine all data frames into one --> long format 
combined_df_R <- bind_rows(df_list_R)

# Pivoting the data to wide format
random_wide_df <- combined_df_R %>%
  pivot_wider(
    names_from = RandomEffects.variable,
    values_from = RandomEffects.V1
  )
# 
# saveRDS(fixed_wide_df, file = "./CovidBAA_INLA/fixed_wide_df_rw2_251121_batch2.RDS")
# saveRDS(random_wide_df, file = "./CovidBAA_INLA/random_wide_df_rw2_251121_batch2.RDS")
# saveRDS(fixed_wide_df, file = "./CovidBAA_INLA/fixed_wide_df_240908_batch2.RDS")
# saveRDS(random_wide_df, file = "./CovidBAA_INLA/random_wide_df_240908_batch2.RDS")

# RUN one more set.seed and proceed to INLA_viz.R

###############################################################################
# 1. combine two simulations to make a full 10,000 simulation dataframe (100%)
# 2. add fixed effect to each of the random effects 
# 3. Retain 95% crI for each week. (trimmed_list) represents betas 
# 4. Visualize the results in the form of exp(beta) = RR
###############################################################################


###############################################################################
# Visualization : INLA_viz.R
###############################################################################

# load in data
# random1<- readRDS("./CovidBAA_INLA/random_wide_df_rw2_251121_batch1.RDS")
# random2<- readRDS("./CovidBAA_INLA/random_wide_df_rw2_251121_batch2.RDS")
# random2$Simulation<- random2$Simulation+5000
# fixed1<- readRDS("./CovidBAA_INLA/fixed_wide_df_rw2_251121_batch1.RDS")
# fixed2<- readRDS("./CovidBAA_INLA/fixed_wide_df_rw2_251121_batch2.RDS")
# fixed2$Simulation<- fixed2$Simulation+5000
# dates_label<- readRDS("./cloud_INLA/week_dates.rds")
# names(dates_label)[2]<-"TimePoint"

# saveRDS(trimmed_list, file = "./CovidBAA_INLA/beta_export_240908.RDS") # AR1 version 
# saveRDS(trimmed_list, file = "./CovidBAA_INLA/beta_export_rw2_251121.RDS") # RW2 version


# relative risk : exp(beta) --------------------------------------------------#
# Function to exponentiate values in a data frame except for the "Simulation" column
exponentiate_values <- function(df) {
  df %>%
    mutate(across(-Simulation, exp))
}
# Apply the function to each data frame in the list
df_list_exp <- lapply(trimmed_list, exponentiate_values)

# Function to exponentiate values in a data frame except for the "Simulation" column
ve_values <- function(df) {
  df %>%
    mutate(across(-Simulation, ~ 1 - .))
}


# Apply the function to each data frame in the list
df_list_ve <- lapply(df_list_exp, ve_values)





# Plotting 95% CrI   --------------------------------------------------------------#
# flagging weeks of emergence
emergence_period <- data.frame(date_start= 
                                 as.Date(c("2021-06-03", "2021-11-30", "2022-02-23", "2022-06-09")),
                               date_end = as.Date(c("2021-07-08","2022-01-04", "2022-03-30", "2022-07-14")))

#"round" the date down to 
emergence_period$date_start <- floor_date(emergence_period$date_start, unit='week')
emergence_period$date_end <- floor_date(emergence_period$date_end, unit='week')


df <- df_list_ve[[2]] # main group of interest

# Convert to long format
df_long <- df %>%
  select(-Simulation) %>% 
  pivot_longer(cols = everything(), names_to = "TimePoint", values_to = "Value") %>% 
  mutate(TimePoint = as.numeric(sub(".*\\.", "", TimePoint))) %>% 
  #   mutate(TimePoint = as.numeric(sub(".*:", "", TimePoint)))%>% 
  left_join(dates_label, "TimePoint")

# Plotting
color_mean_ci<-  "#8E44AD"
timeVE_vaxGroup2_mean_ribbon <- df_long %>%
  group_by(week.date) %>%
  summarize(
    mean_value = mean(Value * 100, na.rm = TRUE),
    upper_ci = max(Value * 100, na.rm = TRUE),
    lower_ci = min(Value * 100, na.rm = TRUE)
  ) %>%
  ggplot(aes(x = week.date, y = mean_value)) +
  # geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci),
  #             fill = "#FAE09366") + # Ribbon for CI
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci),
                width = 5, color = color_mean_ci, alpha = 0.4) + # 95% CI lines
  geom_point(color = "#8E44AD" , alpha = 0.6) + # Mean points
  labs(title = "", x = "Collection date (week)",
       y = expression(atop("Vaccine effectiveness", "1 - " ~exp(beta)))) +
  theme_classic() +
  geom_hline(yintercept = 0, color = "blue", linetype = 3) +
  scale_x_date(limits = as.Date(c("2021-03-04", "2023-02-01")),
               date_breaks = "1 month", date_labels = "%b %Y") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        legend.position = "none") +
  scale_y_continuous(limits = c(-25, 100), breaks = c(-25, 0, 25, 50, 75, 100)) +
  geom_rect(data = emergence_period, inherit.aes = FALSE,
            aes(xmin = date_start, xmax = date_end, ymin = -Inf, ymax = Inf),
            fill = "grey", alpha = 0.3) +
  geom_vline(xintercept = as.Date("2022-09-15"), col = '#6C7B8B', lty = 2, linewidth = 1.2)+
  annotate("text", x = as.Date("2022-09-15"), 
           y = 80, 
           label = "Vaccine update \n BA.4/BA.5", fontface = "bold", 
           color = "#6C7B8B", angle = 90, hjust = 0.6, vjust = 1.5)+
  geom_vline(xintercept=as.Date("2022-09-15"), col='#6C7B8B', lty=2, linewidth = 1.2) # bivalent update 


timeVE_vaxGroup2_mean_ribbon
# ggsave(plot = timeVE_vaxGroup2_mean_ribbon, filename = "./24draft_figures/RevisionFigure_Figure2C_INLA_RW2.jpeg", device = "jpeg", width = 12, height = 6, dpi = 200)

