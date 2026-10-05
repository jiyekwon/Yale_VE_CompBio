######################################################################
# This is the baseline meta-regression JAGS code 
# GOAL: eta_1 or \gamma_1 parameter estimation 
# 
######################################################################
library(rjags)
library(dplyr)
library(ggplot2)
library(bayesplot)
library(reshape2)
library(tidyr)
library(scoringRules)
#==============================================================================
# 1. Loading in the right data to run `2. JAGS proper` yourself.
# NOTE: if you want to replicate the results skip to `3. post-processing` section
#==============================================================================
# Datasets -----------------------------------------------------------------#
test_Ve<- readRDS("./cloud_INLA/y_ve_df.rds") 
test_beta<- readRDS("./CovidBAA_INLA/beta_export_240908.RDS")
# test_dist<- readRDS("./cloud_INLA/cases_los_vax_visit_2024Aug.RDS")
test_dist<- readRDS("./24draft_data/cases_los_ct_first_sample.RDS") # 8897 *67
dist_all<- test_dist %>% select(bnt_distance, bnt_distance_rbd, week.date, week_order)
dist_all = dist_all %>% rename(TimePoint = week_order)


# Datasets: from INLA --------------------------------------------------------#
y_beta_df<- as.data.frame(test_beta[[2]]) #15-90 days, 9500 sims per week 

# Convert to long format
y_beta_long <- y_beta_df %>%
  pivot_longer(cols = -Simulation, names_to = "TimePoint", values_to = "Value") %>%
  mutate(TimePoint = as.numeric(sub(".*\\.", "", TimePoint)))

# calculate the mean beta (from INLA) for each week 
beta_mean = y_beta_long %>% group_by(TimePoint) %>% 
  summarize(mean_beta = mean(Value, na.rm = T))

# calculate the variance of beta for each week 
beta_var = y_beta_long %>% group_by(TimePoint) %>% 
  summarize(var_beta = var(Value, na.rm = T),
            inv_var_beta = 1/var_beta)



# Datasets: from AA substitutions ---------------------------------------------#
# calculate the mean distance (AA) for each week 
dist_spike_mean = dist_all %>% group_by(TimePoint) %>% 
  summarize(mean_dist_spike = mean(bnt_distance, na.rm = T),
            mean_dist_rbd = mean(bnt_distance_rbd, na.rm = T))
dist_spike_mean$TimePoint<- as.numeric(dist_spike_mean$TimePoint)

# standardize the distance 
dist_spike_mean$std_spike<- scale(dist_spike_mean$mean_dist_spike)
dist_spike_mean$std_rbd<- scale(dist_spike_mean$mean_dist_rbd)


#==============================================================================
# 2. JAGS proper 
#==============================================================================
# Initialize an empty list to store keep_beta results
results <- data.frame(
  TimePoint = integer(),
  eta_0 = numeric(),
  eta_1 = numeric(),
  keep_beta_mean = numeric(),
  keep_beta_2.5 = numeric(),
  keep_beta_97.5 = numeric()
)

model_string_ve <- "
    model{
    # Define initial condition for alpha[1]
    alpha[1] ~ dnorm(0, ((1-(rho)^2)*inv_sigma)) # initial condition for alpha_1
    
    for (t in 2:N_week){
    # AR(1) process for alpha_t (starting from t=2)
    alpha[t] ~ dnorm( (rho * alpha[t-1]), inv_sigma) # for t >1 
    }
    
    for (t in 1:N_week){
    
    RR_true[t] = exp(beta_true[t])

    # Observe beta_obs for all time points
    beta_obs[t] ~ dnorm(beta_true[t], inv_var[t])

    # Conditional model for beta_true based on availability of GD
    beta_true[t] ~ dnorm(eta_0 + eta_1 * GD[t] + alpha[t], tau_inv)
    
    error[t] ~ dnorm(0, tau_inv)
    }

    # Priors for beta coefficients and initial theta values
    eta_0 ~ dnorm(0, 0.0001)
    eta_1 ~ dnorm(0, 0.0001)
    rho ~ dunif(-1, 1)
    inv_sigma ~ dgamma(0.01, 0.01) 
    tau_inv ~ dgamma(0.01, 0.01)
    }
"

# Filter the dataset to include data up to the final week (assuming week 75 as the last week)
dist_spike_mean_sub <- dist_spike_mean %>% filter(TimePoint <= 91)
beta_mean_sub <- beta_mean %>% filter(TimePoint <= 91)
beta_var_sub <- beta_var %>% filter(TimePoint <= 91)

# Prepare the data list for JAGS
data_list <- list(
  N_week = length(beta_mean_sub$mean_beta),   # Number of time points (growing)
  GD = c(dist_spike_mean_sub$std_spike),      # GD data up to the current week
  beta_obs = beta_mean_sub$mean_beta,         # Beta observed values
  inv_var = beta_var_sub$inv_var_beta         # Inverse variance for each week
)

# Initialize the JAGS model
jags_model <- jags.model(textConnection(model_string_ve), data = data_list, 
                         n.chains = 4, n.adapt = 500000)

# Update (burn-in)
update(jags_model, n.iter = 500000)

# Sample from the posterior
params <- c("eta_0", "eta_1", "beta_true")  # Get all beta_true for each week
samples <- coda.samples(jags_model, variable.names = params, n.iter = 100000) # change to 10k 

# Extract the summary for beta_true
beta_true_stats <- summary(samples)$quantiles
beta_true_mean <- summary(samples)$statistics[, "Mean"]
beta_true_2.5 <- beta_true_stats[,"2.5%"]
beta_true_97.5 <- beta_true_stats[,"97.5%"]

# Extract the mean values for eta_0 and eta_1
eta_0_mean <- summary(samples)$statistics["eta_0", "Mean"]
eta_1_mean <- summary(samples)$statistics["eta_1", "Mean"]

# Append the results for each week to the results data frame
for (week in 1:length(beta_mean_sub$mean_beta)) {
  results <- rbind(results, data.frame(
    TimePoint = week,
    eta_0 = eta_0_mean,
    eta_1 = eta_1_mean,
    keep_beta_mean = beta_true_mean[week],
    keep_beta_2.5 = beta_true_2.5[week],
    keep_beta_97.5 = beta_true_97.5[week]
  ))
}

#==============================================================================
# 3. Post-processing
#==============================================================================
# escape hatch, read in RDS file below.
# samples<- readRDS(file ="~/OneDrive - Yale University/Research/Yale_CDCBAA/24draft_data/JAGS_sampleswf3_output_250527.RDS")
# results<- readRDS("~/OneDrive - Yale University/Research/Yale_CDCBAA/24draft_data/JAGS_betaAllwf3_output_91weeks.RDS")


#Calculate the mean of eta_1 from the samples
eta_1_all<- as.data.frame(as.matrix(samples)[, c("eta_1")])
names(eta_1_all)[1]<- "eta_1"


# divide everything by GD_obs sd
gd_spike_sd<- sd(dist_spike_mean_sub$mean_dist_spike)
gd_spike_mean<- mean(dist_spike_mean_sub$mean_dist_spike)

eta_1_all$eta_spike<- eta_1_all$eta_1/gd_spike_sd 

# apply multiplier by the predictor we want 
eta_1_all$eta_spike_10<- (eta_1_all$eta_spike) * 10

# Exponentiate all 
eta_1_all$eta_spike_10_RR <- exp(eta_1_all$eta_spike_10)

# get the median
eta_spike_10RR_median<- median(eta_1_all$eta_spike_10_RR) # (1.161) ~16% change in RR with 10 GD increase

# get the mean 
eta_spike_10RR_mean<- mean(eta_1_all$eta_spike_10_RR) # (1.161) ~16% change in RR with 10 GD increase

# get the CrIs 
# Calculate the 95% credible interval
eta_spike_10RR_lower <- quantile(eta_1_all$eta_spike_10_RR, probs = 0.025) # 0.29% increase 
eta_spike_10RR_upper <- quantile(eta_1_all$eta_spike_10_RR, probs = 0.975) # 33% increase 
print(paste0(round((eta_spike_10RR_median*100)-100,2),
             ",[", round(((eta_spike_10RR_lower-1)*100),2), ",",
             round(((eta_spike_10RR_upper-1)*100),2),"]")) # [-1.02,35.55]

# 90% CrI 
eta_spike_10RR_lower <- quantile(eta_1_all$eta_spike_10_RR, probs = 0.05)  
eta_spike_10RR_upper <- quantile(eta_1_all$eta_spike_10_RR, probs = 0.95) 
print(paste0(round((eta_spike_10RR_median*100)-100,2),
             ",[", round(((eta_spike_10RR_lower-1)*100),2), ",",
             round(((eta_spike_10RR_upper-1)*100),2),"]")) # [1.6, 32.35]


