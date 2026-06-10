# ================================================================
# BAYES FACTOR: MAIN EFFECT OF EXERCISE INTENSITY
# ================================================================

library(brms)
library(rstan)
library(rstudioapi)
library(dplyr)

# --------------------------------------------------
# LOAD DATA
# --------------------------------------------------

sub_online_ratings <- read.csv(
  "C:/Users/user/Desktop/projects/HotSpin/Data/LogExperiment/MAIN/sub_online_heat_and_pressure_clean.csv",
  sep = ",",
  header = TRUE
)

# --------------------------------------------------
# CREATE DATASET
# --------------------------------------------------

mean_ratings_sub <- sub_online_ratings %>%
  group_by(subject, time, exercise_intensity, trial, block) %>%
  summarise_at("online_rating", mean, na.rm = TRUE)

sub_max_online_ratings <- mean_ratings_sub %>%
  group_by(subject, exercise_intensity, trial, block) %>%
  summarise_at("online_rating", max, na.rm = TRUE)

sub_max_online_ratings$online_rating[
  is.infinite(sub_max_online_ratings$online_rating)
] <- NA

# --------------------------------------------------
# FACTORS
# --------------------------------------------------

sub_max_online_ratings$subject <- factor(sub_max_online_ratings$subject)
sub_max_online_ratings$exercise_intensity <- factor(sub_max_online_ratings$exercise_intensity)
sub_max_online_ratings$trial <- factor(sub_max_online_ratings$trial)
sub_max_online_ratings$block <- factor(sub_max_online_ratings$block)

# --------------------------------------------------
# REMOVE MISSING DATA
# --------------------------------------------------

sub_max_online_ratings <- sub_max_online_ratings %>%
  filter(
    !is.na(online_rating),
    !is.na(subject),
    !is.na(exercise_intensity),
    !is.na(trial),
    !is.na(block)
  )

# --------------------------------------------------
# PRIORS
# --------------------------------------------------

priors <- c(
  prior(normal(0, 8), class = "b"),
  prior(normal(75, 30), class = "Intercept"),
  prior(student_t(3, 0, 12), class = "sigma"),
  prior(student_t(3, 0, 12), class = "sd")
)

# --------------------------------------------------
# MODEL 1: BASELINE MODEL
# --------------------------------------------------

fit_null_all <- brm(
  online_rating ~
    trial +
    block +
    (1 | subject),
  
  data = sub_max_online_ratings,
  family = gaussian(),
  prior = priors,
  
  chains = 4,
  iter = 4000,
  warmup = 2000,
  cores = 4,
  
  save_pars = save_pars(all = TRUE)
)

# --------------------------------------------------
# MODEL 2: EXERCISE MODEL
# --------------------------------------------------

fit_main_all <- brm(
  online_rating ~
    exercise_intensity +
    trial +
    block +
    (1 | subject),
  
  data = sub_max_online_ratings,
  family = gaussian(),
  prior = priors,
  
  chains = 4,
  iter = 4000,
  warmup = 2000,
  cores = 4,
  
  control = list(
    adapt_delta = 0.95,
    max_treedepth = 12
  ),
  
  save_pars = save_pars(all = TRUE)
)

# --------------------------------------------------
# CHECK CONVERGENCE
# --------------------------------------------------

summary(fit_null_all)
summary(fit_main_all)

# --------------------------------------------------
# BAYES FACTOR
# --------------------------------------------------

bf_main_all <- bayes_factor(
  fit_main_all,
  fit_null_all
)

bf_main_all


# ================================================================
# HEAT
# ================================================================

# Filter heat trials
sub_online_ratings_heat <- sub_online_ratings %>%
  filter(modality == 2)

# --------------------------------------------------
# CREATE DATASET
# --------------------------------------------------

mean_ratings_sub_heat <- sub_online_ratings_heat %>%
  group_by(subject, time, exercise_intensity, trial, block) %>%
  summarise_at("online_rating", mean, na.rm = TRUE)

sub_max_online_ratings_heat <- mean_ratings_sub_heat %>%
  group_by(subject, exercise_intensity, trial, block) %>%
  summarise_at("online_rating", max, na.rm = TRUE)

sub_max_online_ratings_heat$online_rating[
  is.infinite(sub_max_online_ratings_heat$online_rating)
] <- NA

# --------------------------------------------------
# FACTORS
# --------------------------------------------------

sub_max_online_ratings_heat$subject <- factor(sub_max_online_ratings_heat$subject)
sub_max_online_ratings_heat$exercise_intensity <- factor(sub_max_online_ratings_heat$exercise_intensity)
sub_max_online_ratings_heat$trial <- factor(sub_max_online_ratings_heat$trial)
sub_max_online_ratings_heat$block <- factor(sub_max_online_ratings_heat$block)

# --------------------------------------------------
# REMOVE MISSING DATA
# --------------------------------------------------

sub_max_online_ratings_heat <- sub_max_online_ratings_heat %>%
  filter(
    !is.na(online_rating),
    !is.na(subject),
    !is.na(exercise_intensity),
    !is.na(trial),
    !is.na(block)
  )

# --------------------------------------------------
# MODEL 1: BASELINE MODEL
# --------------------------------------------------

fit_null_heat <- brm(
  online_rating ~
    trial +
    block +
    (1 | subject),
  
  data = sub_max_online_ratings_heat,
  family = gaussian(),
  prior = priors,
  
  chains = 4,
  iter = 4000,
  warmup = 2000,
  cores = 4,
  
  save_pars = save_pars(all = TRUE)
)

# --------------------------------------------------
# MODEL 2: EXERCISE MODEL
# --------------------------------------------------

fit_main_heat <- brm(
  online_rating ~
    exercise_intensity +
    trial +
    block +
    (1 | subject),
  
  data = sub_max_online_ratings_heat,
  family = gaussian(),
  prior = priors,
  
  chains = 4,
  iter = 4000,
  warmup = 2000,
  cores = 4,
  
  control = list(
    adapt_delta = 0.95,
    max_treedepth = 12
  ),
  
  save_pars = save_pars(all = TRUE)
)

# --------------------------------------------------
# BAYES FACTOR
# --------------------------------------------------

bf_heat <- bayes_factor(
  fit_main_heat,
  fit_null_heat
)

bf_heat


# ================================================================
# PRESSURE
# ================================================================

# Filter pressure trials
sub_online_ratings_pressure <- sub_online_ratings %>%
  filter(modality == 1)

# --------------------------------------------------
# CREATE DATASET
# --------------------------------------------------

mean_ratings_sub_pressure <- sub_online_ratings_pressure %>%
  group_by(subject, time, exercise_intensity, trial, block) %>%
  summarise_at("online_rating", mean, na.rm = TRUE)

sub_max_online_ratings_pressure <- mean_ratings_sub_pressure %>%
  group_by(subject, exercise_intensity, trial, block) %>%
  summarise_at("online_rating", max, na.rm = TRUE)

sub_max_online_ratings_pressure$online_rating[
  is.infinite(sub_max_online_ratings_pressure$online_rating)
] <- NA

# --------------------------------------------------
# FACTORS
# --------------------------------------------------

sub_max_online_ratings_pressure$subject <- factor(sub_max_online_ratings_pressure$subject)
sub_max_online_ratings_pressure$exercise_intensity <- factor(sub_max_online_ratings_pressure$exercise_intensity)
sub_max_online_ratings_pressure$trial <- factor(sub_max_online_ratings_pressure$trial)
sub_max_online_ratings_pressure$block <- factor(sub_max_online_ratings_pressure$block)

# --------------------------------------------------
# REMOVE MISSING DATA
# --------------------------------------------------

sub_max_online_ratings_pressure <- sub_max_online_ratings_pressure %>%
  filter(
    !is.na(online_rating),
    !is.na(subject),
    !is.na(exercise_intensity),
    !is.na(trial),
    !is.na(block)
  )

# --------------------------------------------------
# MODEL 1: BASELINE MODEL
# --------------------------------------------------

fit_null_pressure <- brm(
  online_rating ~
    trial +
    block +
    (1 | subject),
  
  data = sub_max_online_ratings_pressure,
  family = gaussian(),
  prior = priors,
  
  chains = 4,
  iter = 4000,
  warmup = 2000,
  cores = 4,
  
  save_pars = save_pars(all = TRUE)
)

# --------------------------------------------------
# MODEL 2: EXERCISE MODEL
# --------------------------------------------------

fit_main_pressure <- brm(
  online_rating ~
    exercise_intensity +
    trial +
    block +
    (1 | subject),
  
  data = sub_max_online_ratings_pressure,
  family = gaussian(),
  prior = priors,
  
  chains = 4,
  iter = 4000,
  warmup = 2000,
  cores = 4,
  
  control = list(
    adapt_delta = 0.95,
    max_treedepth = 12
  ),
  
  save_pars = save_pars(all = TRUE)
)

# --------------------------------------------------
# BAYES FACTOR
# --------------------------------------------------

bf_pressure <- bayes_factor(
  fit_main_pressure,
  fit_null_pressure
)

bf_pressure



# ================================================================
# BAYES FACTOR: EXERCISE INTENSITY × STIMULUS INTENSITY
# ================================================================

# --------------------------------------------------
# PRIORS
# --------------------------------------------------

priors <- c(
  prior(normal(0, 8), class = "b"),
  prior(normal(75, 30), class = "Intercept"),
  prior(student_t(3, 0, 12), class = "sigma"),
  prior(student_t(3, 0, 12), class = "sd")
)

# ================================================================
# ACROSS MODALITIES
# ================================================================

# --------------------------------------------------
# CREATE DATASET
# --------------------------------------------------

mean_ratings_sub_interaction <- sub_online_ratings %>%
  group_by(
    subject,
    time,
    exercise_intensity,
    VAS_intensity,
    trial,
    block,
    modality
  ) %>%
  summarise(
    online_rating = mean(online_rating, na.rm = TRUE),
    .groups = "drop"
  )

sub_max_online_ratings_interaction <- mean_ratings_sub_interaction %>%
  group_by(
    subject,
    exercise_intensity,
    VAS_intensity,
    trial,
    block,
    modality
  ) %>%
  summarise(
    online_rating = max(online_rating, na.rm = TRUE),
    .groups = "drop"
  )

# Replace Inf values with NA

sub_max_online_ratings_interaction$online_rating[
  is.infinite(sub_max_online_ratings_interaction$online_rating)
] <- NA

# --------------------------------------------------
# FACTORS
# --------------------------------------------------

sub_max_online_ratings_interaction$subject <-
  factor(sub_max_online_ratings_interaction$subject)

sub_max_online_ratings_interaction$exercise_intensity <-
  factor(sub_max_online_ratings_interaction$exercise_intensity)

sub_max_online_ratings_interaction$VAS_intensity <-
  factor(sub_max_online_ratings_interaction$VAS_intensity)

sub_max_online_ratings_interaction$trial <-
  factor(sub_max_online_ratings_interaction$trial)

sub_max_online_ratings_interaction$block <-
  factor(sub_max_online_ratings_interaction$block)

sub_max_online_ratings_interaction$modality <-
  factor(sub_max_online_ratings_interaction$modality)

# --------------------------------------------------
# REMOVE MISSING DATA
# --------------------------------------------------

sub_max_online_ratings_interaction <-
  sub_max_online_ratings_interaction %>%
  filter(
    !is.na(online_rating),
    !is.na(subject),
    !is.na(exercise_intensity),
    !is.na(VAS_intensity),
    !is.na(trial),
    !is.na(block)
  )

# --------------------------------------------------
# MODEL 1: BASELINE MODEL
# --------------------------------------------------

fit_null_interaction <- brm(
  online_rating ~
    trial +
    block +
    (1 | subject),
  
  data = sub_max_online_ratings_interaction,
  family = gaussian(),
  prior = priors,
  
  chains = 4,
  iter = 4000,
  warmup = 2000,
  cores = 4,
  
  save_pars = save_pars(all = TRUE)
)

# --------------------------------------------------
# MODEL 2: MAIN EFFECTS MODEL
# --------------------------------------------------

fit_main_interaction <- brm(
  online_rating ~
    exercise_intensity +
    VAS_intensity +
    trial +
    block +
    (1 | subject),
  
  data = sub_max_online_ratings_interaction,
  family = gaussian(),
  prior = priors,
  
  chains = 4,
  iter = 4000,
  warmup = 2000,
  cores = 4,
  
  control = list(
    adapt_delta = 0.95,
    max_treedepth = 12
  ),
  
  save_pars = save_pars(all = TRUE)
)

# --------------------------------------------------
# MODEL 3: INTERACTION MODEL
# --------------------------------------------------

fit_full_interaction <- brm(
  online_rating ~
    exercise_intensity * VAS_intensity +
    trial +
    block +
    (1 | subject),
  
  data = sub_max_online_ratings_interaction,
  family = gaussian(),
  prior = priors,
  
  chains = 4,
  iter = 4000,
  warmup = 2000,
  cores = 4,
  
  control = list(
    adapt_delta = 0.95,
    max_treedepth = 12
  ),
  
  save_pars = save_pars(all = TRUE)
)

# --------------------------------------------------
# CHECK CONVERGENCE
# --------------------------------------------------

summary(fit_null_interaction)
summary(fit_main_interaction)
summary(fit_full_interaction)

# --------------------------------------------------
# BAYES FACTORS
# --------------------------------------------------

# Main effects model vs baseline model

bf_main_vs_null <- bayes_factor(
  fit_main_interaction,
  fit_null_interaction
)

bf_main_vs_null

# Interaction model vs main effects model
# (primary test of the interaction)

bf_interaction_vs_main <- bayes_factor(
  fit_full_interaction,
  fit_main_interaction
)

bf_interaction_vs_main

# Full interaction model vs baseline model

bf_interaction_vs_null <- bayes_factor(
  fit_full_interaction,
  fit_null_interaction
)

bf_interaction_vs_null


# ================================================================
# BAYES FACTOR: HEAT
# Exercise Intensity × VAS Intensity
# ================================================================

# --------------------------------------------------
# CREATE HEAT DATASET
# --------------------------------------------------

sub_online_ratings_heat <- sub_online_ratings %>%
  filter(modality == 2)

# --------------------------------------------------
# MEAN TIME COURSES
# --------------------------------------------------

mean_ratings_sub_heat_int <- sub_online_ratings_heat %>%
  group_by(
    subject,
    time,
    exercise_intensity,
    VAS_intensity,
    trial,
    block
  ) %>%
  summarise(
    online_rating = mean(online_rating, na.rm = TRUE),
    .groups = "drop"
  )

# --------------------------------------------------
# MAXIMUM RATING
# --------------------------------------------------

sub_max_online_ratings_heat_int <- mean_ratings_sub_heat_int %>%
  group_by(
    subject,
    exercise_intensity,
    VAS_intensity,
    trial,
    block
  ) %>%
  summarise(
    online_rating = max(online_rating, na.rm = TRUE),
    .groups = "drop"
  )

# Replace Inf values with NA

sub_max_online_ratings_heat_int$online_rating[
  is.infinite(sub_max_online_ratings_heat_int$online_rating)
] <- NA

# --------------------------------------------------
# FACTORS
# --------------------------------------------------

sub_max_online_ratings_heat_int$subject <-
  factor(sub_max_online_ratings_heat_int$subject)

sub_max_online_ratings_heat_int$exercise_intensity <-
  factor(sub_max_online_ratings_heat_int$exercise_intensity)

sub_max_online_ratings_heat_int$VAS_intensity <-
  factor(sub_max_online_ratings_heat_int$VAS_intensity)

sub_max_online_ratings_heat_int$trial <-
  factor(sub_max_online_ratings_heat_int$trial)

sub_max_online_ratings_heat_int$block <-
  factor(sub_max_online_ratings_heat_int$block)

# --------------------------------------------------
# REMOVE MISSING DATA
# --------------------------------------------------

sub_max_online_ratings_heat_int <-
  sub_max_online_ratings_heat_int %>%
  filter(
    !is.na(online_rating),
    !is.na(subject),
    !is.na(exercise_intensity),
    !is.na(VAS_intensity),
    !is.na(trial),
    !is.na(block)
  )

# --------------------------------------------------
# MODEL 1: BASELINE MODEL
# --------------------------------------------------

fit_null_heat_int <- brm(
  online_rating ~
    trial +
    block +
    (1 | subject),
  
  data = sub_max_online_ratings_heat_int,
  family = gaussian(),
  prior = priors,
  
  chains = 4,
  iter = 4000,
  warmup = 2000,
  cores = 4,
  
  save_pars = save_pars(all = TRUE)
)

# --------------------------------------------------
# MODEL 2: MAIN EFFECTS MODEL
# --------------------------------------------------

fit_main_heat_int <- brm(
  online_rating ~
    exercise_intensity +
    VAS_intensity +
    trial +
    block +
    (1 | subject),
  
  data = sub_max_online_ratings_heat_int,
  family = gaussian(),
  prior = priors,
  
  chains = 4,
  iter = 4000,
  warmup = 2000,
  cores = 4,
  
  control = list(
    adapt_delta = 0.95,
    max_treedepth = 12
  ),
  
  save_pars = save_pars(all = TRUE)
)

# --------------------------------------------------
# MODEL 3: INTERACTION MODEL
# --------------------------------------------------

fit_full_heat_int <- brm(
  online_rating ~
    exercise_intensity * VAS_intensity +
    trial +
    block +
    (1 | subject),
  
  data = sub_max_online_ratings_heat_int,
  family = gaussian(),
  prior = priors,
  
  chains = 4,
  iter = 4000,
  warmup = 2000,
  cores = 4,
  
  control = list(
    adapt_delta = 0.95,
    max_treedepth = 12
  ),
  
  save_pars = save_pars(all = TRUE)
)

# --------------------------------------------------
# SUMMARIES
# --------------------------------------------------

summary(fit_null_heat_int)
summary(fit_main_heat_int)
summary(fit_full_heat_int)

# --------------------------------------------------
# BAYES FACTORS
# --------------------------------------------------

bf_main_vs_null_heat <- bayes_factor(
  fit_main_heat_int,
  fit_null_heat_int
)

bf_interaction_vs_main_heat <- bayes_factor(
  fit_full_heat_int,
  fit_main_heat_int
)

bf_interaction_vs_null_heat <- bayes_factor(
  fit_full_heat_int,
  fit_null_heat_int
)

bf_main_vs_null_heat
bf_interaction_vs_main_heat
bf_interaction_vs_null_heat

# ================================================================
# BAYES FACTOR: PRESSURE
# Exercise Intensity × VAS Intensity
# ================================================================

# --------------------------------------------------
# CREATE PRESSURE DATASET
# --------------------------------------------------

sub_online_ratings_pressure <- sub_online_ratings %>%
  filter(modality == 1)

# --------------------------------------------------
# MEAN TIME COURSES
# --------------------------------------------------

mean_ratings_sub_pressure_int <- sub_online_ratings_pressure %>%
  group_by(
    subject,
    time,
    exercise_intensity,
    VAS_intensity,
    trial,
    block
  ) %>%
  summarise(
    online_rating = mean(online_rating, na.rm = TRUE),
    .groups = "drop"
  )

# --------------------------------------------------
# MAXIMUM RATING
# --------------------------------------------------

sub_max_online_ratings_pressure_int <- mean_ratings_sub_pressure_int %>%
  group_by(
    subject,
    exercise_intensity,
    VAS_intensity,
    trial,
    block
  ) %>%
  summarise(
    online_rating = max(online_rating, na.rm = TRUE),
    .groups = "drop"
  )

# Replace Inf values with NA

sub_max_online_ratings_pressure_int$online_rating[
  is.infinite(sub_max_online_ratings_pressure_int$online_rating)
] <- NA

# --------------------------------------------------
# FACTORS
# --------------------------------------------------

sub_max_online_ratings_pressure_int$subject <-
  factor(sub_max_online_ratings_pressure_int$subject)

sub_max_online_ratings_pressure_int$exercise_intensity <-
  factor(sub_max_online_ratings_pressure_int$exercise_intensity)

sub_max_online_ratings_pressure_int$VAS_intensity <-
  factor(sub_max_online_ratings_pressure_int$VAS_intensity)

sub_max_online_ratings_pressure_int$trial <-
  factor(sub_max_online_ratings_pressure_int$trial)

sub_max_online_ratings_pressure_int$block <-
  factor(sub_max_online_ratings_pressure_int$block)

# --------------------------------------------------
# REMOVE MISSING DATA
# --------------------------------------------------

sub_max_online_ratings_pressure_int <-
  sub_max_online_ratings_pressure_int %>%
  filter(
    !is.na(online_rating),
    !is.na(subject),
    !is.na(exercise_intensity),
    !is.na(VAS_intensity),
    !is.na(trial),
    !is.na(block)
  )

# --------------------------------------------------
# MODEL 1: BASELINE MODEL
# --------------------------------------------------

fit_null_pressure_int <- brm(
  online_rating ~
    trial +
    block +
    (1 | subject),
  
  data = sub_max_online_ratings_pressure_int,
  family = gaussian(),
  prior = priors,
  
  chains = 4,
  iter = 4000,
  warmup = 2000,
  cores = 4,
  
  save_pars = save_pars(all = TRUE)
)

# --------------------------------------------------
# MODEL 2: MAIN EFFECTS MODEL
# --------------------------------------------------

fit_main_pressure_int <- brm(
  online_rating ~
    exercise_intensity +
    VAS_intensity +
    trial +
    block +
    (1 | subject),
  
  data = sub_max_online_ratings_pressure_int,
  family = gaussian(),
  prior = priors,
  
  chains = 4,
  iter = 4000,
  warmup = 2000,
  cores = 4,
  
  control = list(
    adapt_delta = 0.95,
    max_treedepth = 12
  ),
  
  save_pars = save_pars(all = TRUE)
)

# --------------------------------------------------
# MODEL 3: INTERACTION MODEL
# --------------------------------------------------

fit_full_pressure_int <- brm(
  online_rating ~
    exercise_intensity * VAS_intensity +
    trial +
    block +
    (1 | subject),
  
  data = sub_max_online_ratings_pressure_int,
  family = gaussian(),
  prior = priors,
  
  chains = 4,
  iter = 4000,
  warmup = 2000,
  cores = 4,
  
  control = list(
    adapt_delta = 0.95,
    max_treedepth = 12
  ),
  
  save_pars = save_pars(all = TRUE)
)

# --------------------------------------------------
# SUMMARIES
# --------------------------------------------------

summary(fit_null_pressure_int)
summary(fit_main_pressure_int)
summary(fit_full_pressure_int)

# --------------------------------------------------
# BAYES FACTORS
# --------------------------------------------------

bf_main_vs_null_pressure <- bayes_factor(
  fit_main_pressure_int,
  fit_null_pressure_int
)

bf_interaction_vs_main_pressure <- bayes_factor(
  fit_full_pressure_int,
  fit_main_pressure_int
)

bf_interaction_vs_null_pressure <- bayes_factor(
  fit_full_pressure_int,
  fit_null_pressure_int
)

bf_main_vs_null_pressure
bf_interaction_vs_main_pressure
bf_interaction_vs_null_pressure


#===============================================================
# Bayesian Analysis of Exercise Intensity x SubGroup across all
# Stimulus intensity and VAS 70
#==============================================================


# ================================================================
# DATA PREPARATION
# ================================================================

summary_df <- df_combined_heat %>%
  group_by(subject, exercise_intensity, treatment_order, gender, group, trial, block) %>%
  summarise_at(c('pain_rating'),mean,na.rm = T)

dat <- summary_df %>%
  filter(
    !is.na(pain_rating),
    !is.na(exercise_intensity),
    !is.na(gender),
    !is.na(group)
  )

#================================================================
# Group VARIABLE
#
# group 1 = males previous study
# group 2 = females previous study
# group 3 = females current study
# ================================================================


# ================================================================
# FACTORS
# ================================================================

dat$subject <- factor(dat$subject)

dat$exercise_intensity <- factor(dat$exercise_intensity)

dat$gender <- factor(
  dat$gender,
  levels = c(1,0),
  labels = c("female","male")
)

dat$group <- factor(dat$group)

dat$trial <- factor(dat$trial)

dat$block <- factor(dat$block)


# ================================================================
# PRIORS
# ================================================================

priors <- c(
  
  prior(normal(50,20), class = "Intercept"),
  
  prior(normal(0,5), class = "b"),
  
  prior(student_t(3,0,10), class = "sigma"),
  
  prior(student_t(3,0,10), class = "sd")
  
)

# ================================================================
# MODEL 0: NULL MODEL
# ================================================================

fit_null <- brm(
  
  pain_rating ~
    
    trial +
    block +
    
    (1  | subject),
  
  data = dat,
  
  family = gaussian(),
  
  prior = priors,
  
  chains = 4,
  iter = 6000,
  warmup = 3000,
  cores = 4,
  
  control = list(
    adapt_delta = 0.99,
    max_treedepth = 15
  ),
  
  save_pars = save_pars(all = TRUE)
  
)

# ================================================================
# MODEL 1: EIH ONLY
# ================================================================

fit_EIH <- brm(
  
  pain_rating ~
    
    exercise_intensity +
    
    trial +
    block +
    
    (1 | subject),
  
  data = dat,
  
  family = gaussian(),
  
  prior = priors,
  
  chains = 4,
  iter = 6000,
  warmup = 3000,
  cores = 4,
  
  control = list(
    adapt_delta = 0.99,
    max_treedepth = 15
  ),
  
  save_pars = save_pars(all = TRUE)
  
)

# ================================================================
# MODEL 2: Group MODERATION
# ================================================================

fit_group <- brm(
  
  pain_rating ~
    
    exercise_intensity * group +
    
    trial +
    block +
    
    (1  | subject),
  
  data = dat,
  
  family = gaussian(),
  
  prior = priors,
  
  chains = 4,
  iter = 6000,
  warmup = 3000,
  cores = 4,
  
  control = list(
    adapt_delta = 0.99,
    max_treedepth = 15
  ),
  
  save_pars = save_pars(all = TRUE)
  
)


#----bayes factor

# Evidence for EIH

bf_EIH <- bayes_factor(
  fit_EIH,
  fit_null
)

# Group moderation

bf_group <- bayes_factor(
  fit_group,
  fit_EIH
)

bf_group_null <- bayes_factor(
  fit_group,
  fit_null
)

bf_EIH
bf_group
bf_group_null

# ================================================================
# DATA PREPARATION
# ================================================================

df_combined_heat_70 <- df_combined_heat[df_combined_heat$VAS == 70,]

summary_df_70 <- df_combined_heat_70 %>%
  group_by(subject, exercise_intensity, treatment_order, gender, group, trial, block) %>%
  summarise_at(c('pain_rating'),mean,na.rm = T)

dat <- summary_df_70 %>%
  filter(
    !is.na(pain_rating),
    !is.na(exercise_intensity),
    !is.na(gender),
    !is.na(group)
  )

#================================================================
# Group VARIABLE
#
# group 1 = males previous study
# group 2 = females previous study
# group 3 = females current study
# ================================================================


# ================================================================
# FACTORS
# ================================================================

dat$subject <- factor(dat$subject)

dat$exercise_intensity <- factor(dat$exercise_intensity)

dat$gender <- factor(
  dat$gender,
  levels = c(1,0),
  labels = c("female","male")
)

dat$group <- factor(dat$group)

dat$trial <- factor(dat$trial)

dat$block <- factor(dat$block)


# ================================================================
# PRIORS
# ================================================================

priors <- c(
  
  prior(normal(50,20), class = "Intercept"),
  
  prior(normal(0,5), class = "b"),
  
  prior(student_t(3,0,10), class = "sigma"),
  
  prior(student_t(3,0,10), class = "sd")
  
)

# ================================================================
# MODEL 0: NULL MODEL
# ================================================================

fit_null_70 <- brm(
  
  pain_rating ~
    
    trial +
    block +
    
    (1 | subject),
  
  data = dat,
  
  family = gaussian(),
  
  prior = priors,
  
  chains = 4,
  iter = 6000,
  warmup = 3000,
  cores = 4,
  
  control = list(
    adapt_delta = 0.99,
    max_treedepth = 15
  ),
  
  save_pars = save_pars(all = TRUE)
  
)

# ================================================================
# MODEL 1: EIH ONLY
# ================================================================

fit_EIH_70 <- brm(
  
  pain_rating ~
    
    exercise_intensity +
    
    trial +
    block +
    
    (1 | subject),
  
  data = dat,
  
  family = gaussian(),
  
  prior = priors,
  
  chains = 4,
  iter = 6000,
  warmup = 3000,
  cores = 4,
  
  control = list(
    adapt_delta = 0.99,
    max_treedepth = 15
  ),
  
  save_pars = save_pars(all = TRUE)
  
)

# ================================================================
# MODEL 2: Group MODERATION
# ================================================================

fit_group_70 <- brm(
  
  pain_rating ~
    
    exercise_intensity * group +
    
    trial +
    block +
    
    (1  | subject),
  
  data = dat,
  
  family = gaussian(),
  
  prior = priors,
  
  chains = 4,
  iter = 6000,
  warmup = 3000,
  cores = 4,
  
  control = list(
    adapt_delta = 0.99,
    max_treedepth = 15
  ),
  
  save_pars = save_pars(all = TRUE)
  
)


#----bayes factor

# Evidence for EIH

bf_EIH_70 <- bayes_factor(
  fit_EIH_70,
  fit_null_70
)

# Group moderation

bf_group_70 <- bayes_factor(
  fit_group_70,
  fit_EIH_70
)

bf_group_null_70 <- bayes_factor(
  fit_group_70,
  fit_null_70
)

bf_EIH_70
bf_group_70
bf_group_null_70

# ================================================================
# BAYESIAN ANALYSIS OF EIH
#
# QUESTIONS:
# 1. Is there evidence for EIH?
# 2. Does fitness (PWC) moderate EIH?
# 3. Does sex moderate EIH?
# 4. Does study/sample moderate EIH?
# ================================================================


# ================================================================
# DATA PREPARATION
# ================================================================

summary_df <- df_combined_heat %>%
  group_by(subject, exercise_intensity, pwc, treatment_order, gender, group, trial, block) %>%
  summarise_at(c('pain_rating'),mean,na.rm = T)

dat <- summary_df %>%
  filter(
    !is.na(pain_rating),
    !is.na(pwc),
    !is.na(gender),
    !is.na(group)
  )

# ================================================================
# STUDY VARIABLE
#
# group 1 = males previous study
# group 2 = females previous study
# group 3 = females current study
# ================================================================

dat$study <- ifelse(
  dat$group == 3,
  "current",
  "previous"
)

# ================================================================
# FACTORS
# ================================================================

dat$subject <- factor(dat$subject)

dat$exercise_intensity <- factor(dat$exercise_intensity)

dat$gender <- factor(
  dat$gender,
  levels = c(1,0),
  labels = c("female","male")
)

dat$study <- factor(dat$study)

dat$trial <- factor(dat$trial)

dat$block <- factor(dat$block)


# ===================================================================
# LMER model
#==================================================================
lme_model_main <- lmer(pain_rating ~
  
  exercise_intensity * pwc * gender * study +

  
  trial +
  block +
  
  (1 | subject), data = dat)

summary(lme_model_main)
# ================================================================
# PRIORS
# ================================================================

priors <- c(
  
  prior(normal(50,20), class = "Intercept"),
  
  prior(normal(0,5), class = "b"),
  
  prior(student_t(3,0,10), class = "sigma"),
  
  prior(student_t(3,0,10), class = "sd")
  
)

# ================================================================
# MODEL 0: NULL MODEL
# ================================================================

fit_null <- brm(
  
  pain_rating ~
    
    trial +
    block +
    
    (1 | subject),
  
  data = dat,
  
  family = gaussian(),
  
  prior = priors,
  
  chains = 4,
  iter = 6000,
  warmup = 3000,
  cores = 4,
  
  control = list(
    adapt_delta = 0.99,
    max_treedepth = 15
  ),
  
  save_pars = save_pars(all = TRUE)
  
)

# ================================================================
# MODEL 1: EIH ONLY
# ================================================================

fit_EIH <- brm(
  
  pain_rating ~
    
    exercise_intensity +
    
    trial +
    block +
    
    (1 | subject),
  
  data = dat,
  
  family = gaussian(),
  
  prior = priors,
  
  chains = 4,
  iter = 6000,
  warmup = 3000,
  cores = 4,
  
  control = list(
    adapt_delta = 0.99,
    max_treedepth = 15
  ),
  
  save_pars = save_pars(all = TRUE)
  
)

# ================================================================
# MODEL 2: FITNESS MODERATION
# ================================================================

fit_fitness <- brm(
  
  pain_rating ~
    
    exercise_intensity * pwc +
    
    trial +
    block +
    
    (1 | subject),
  
  data = dat,
  
  family = gaussian(),
  
  prior = priors,
  
  chains = 4,
  iter = 6000,
  warmup = 3000,
  cores = 4,
  
  control = list(
    adapt_delta = 0.99,
    max_treedepth = 15
  ),
  
  save_pars = save_pars(all = TRUE)
  
)

# ================================================================
# MODEL 3: SEX MODERATION
# ================================================================

fit_sex <- brm(
  
  pain_rating ~
    
    exercise_intensity * gender +
    
    trial +
    block +
    
    (1 | subject),
  
  data = dat,
  
  family = gaussian(),
  
  prior = priors,
  
  chains = 4,
  iter = 6000,
  warmup = 3000,
  cores = 4,
  
  control = list(
    adapt_delta = 0.99,
    max_treedepth = 15
  ),
  
  save_pars = save_pars(all = TRUE)
  
)

# ================================================================
# MODEL 4: STUDY MODERATION
# ================================================================

fit_study <- brm(
  
  pain_rating ~
    
    exercise_intensity * study +
    
    trial +
    block +
    
    (1 | subject),
  
  data = dat,
  
  family = gaussian(),
  
  prior = priors,
  
  chains = 4,
  iter = 6000,
  warmup = 3000,
  cores = 4,
  
  control = list(
    adapt_delta = 0.99,
    max_treedepth = 15
  ),
  
  save_pars = save_pars(all = TRUE)
  
)

# ================================================================
# MODEL 5: FULL MODEL
# ================================================================

fit_full <- brm(
  
  pain_rating ~
    
    exercise_intensity * pwc +
    exercise_intensity * gender +
    exercise_intensity * study +

    trial +
    block +
    
    (1 | subject),
  
  data = dat,
  
  family = gaussian(),
  
  prior = priors,
  
  chains = 4,
  iter = 6000,
  warmup = 3000,
  cores = 4,
  
  control = list(
    adapt_delta = 0.995,
    max_treedepth = 15
  ),
  
  save_pars = save_pars(all = TRUE)
  
)

# ================================================================
# BAYES FACTORS
# ================================================================

# Evidence for EIH

bf_EIH <- bayes_factor(
  fit_EIH,
  fit_null
)

# Fitness moderation

bf_fitness <- bayes_factor(
  fit_fitness,
  fit_EIH
)

# Sex moderation

bf_sex <- bayes_factor(
  fit_sex,
  fit_EIH
)

# Study moderation

bf_study <- bayes_factor(
  fit_study,
  fit_EIH
)

# Full model

bf_full_vs_EIH <- bayes_factor(
  fit_full,
  fit_EIH
)

bf_full_vs_null <- bayes_factor(
  fit_full,
  fit_null
)

bf_EIH
bf_fitness
bf_sex
bf_study
bf_full_vs_EIH
bf_full_vs_null

# ================================================================
# POSTERIOR DISTRIBUTIONS
# ================================================================

posterior_draws <- as_draws_df(fit_full)

# inspect names

names(posterior_draws)

# ================================================================
# POSTERIOR SUMMARIES
# ================================================================

posterior_summary <- posterior_summary(fit_full)

posterior_summary

# ================================================================
# KEY EIH PARAMETERS
# ================================================================

mean(
  posterior_draws$b_exercise_intensity1 < 0
)

mean(
  posterior_draws$`b_exercise_intensity1:pwc` < 0
)

mean(
  posterior_draws$`b_exercise_intensity1:gendermale` < 0
)

mean(
  posterior_draws$`b_exercise_intensity1:studyprevious` < 0
)

# ================================================================
# CREDIBLE INTERVALS
# ================================================================

quantile(
  posterior_draws$b_exercise_intensity1,
  c(.025,.5,.975)
)

quantile(
  posterior_draws$`b_exercise_intensity1:pwc`,
  c(.025,.5,.975)
)

quantile(
  posterior_draws$`b_exercise_intensity1:gendermale`,
  c(.025,.5,.975)
)

quantile(
  posterior_draws$`b_exercise_intensity1:studyprevious`,
  c(.025,.5,.975)
)


# =======================================================
# Visualisation
#======================================================


# ---- define expected parameters ----
param_map <- tibble(
  parameter = c(
    "b_exercise_intensity1",
    "b_exercise_intensity1:pwc",
    "b_exercise_intensity1:gendermale",
    "b_exercise_intensity1:group3"
  ),
  label = c(
    "Exercise effect (EIH)",
    "Fitness moderation",
    "Sex moderation",
    "Previous-study moderation"
  )
)

# ---- keep only parameters that actually exist ----
existing_params <- param_map$parameter[param_map$parameter %in% names(posterior_draws)]

param_map <- param_map %>%
  filter(parameter %in% existing_params)

plot_df <- posterior_draws %>%
  select(all_of(param_map$parameter)) %>%
  pivot_longer(
    cols = everything(),
    names_to = "parameter",
    values_to = "beta"
  ) %>%
  left_join(param_map, by = "parameter")

p1 <- ggplot(plot_df, aes(x = beta)) +
  geom_density(fill = "grey80", color = "#222222", alpha = 0.8,size  =0.8) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.8) +
  facet_wrap(~label, scales = "free") +
  theme_classic() +
  theme(legend.key.size = unit(0.25, 'cm'),axis.title = element_text(size = axis_title_size,family="Helvetica"),
        axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),
        plot.title = element_text(size = plot_title_size,family="Helvetica"),
        legend.title = element_text(size = legend_title_size),
        legend.text = element_text(size = legend_title_size,family="Helvetica"),
        strip.text.x = element_text(size = legend_text_size,family="Helvetica")) +
  labs(x = "Posterior beta", y = "Density")

p1 

coef_summary <- map_dfr(param_map$parameter, function(p) {
  tibble(
    parameter = p,
    Mean = mean(posterior_draws[[p]]),
    Lower = quantile(posterior_draws[[p]], 0.025),
    Upper = quantile(posterior_draws[[p]], 0.975)
  )
}) %>%
  left_join(param_map, by = "parameter")


p2 <- ggplot(coef_summary, aes(y = label, x = Mean)) +
  geom_vline(xintercept = 0, linetype = 2, color = "black") +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper),
                 height = 0.15,
                 colour = "black",
                 linewidth = 0.5) +
  geom_point(size = 3, color = "grey20") +
  theme_classic() +
  theme(legend.key.size = unit(0.25, 'cm'),axis.title = element_text(size = axis_title_size,family="Helvetica"),
        axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),
        plot.title = element_text(size = plot_title_size,family="Helvetica"),
        legend.title = element_text(size = legend_title_size),
        legend.text = element_text(size = legend_title_size,family="Helvetica"),
        strip.text.x = element_text(size = legend_text_size,family="Helvetica")) +
  labs(x = "Posterior estimate (95% CrI)", y = NULL)

p2

prob_df <- map_dfr(param_map$parameter, function(p) {
  tibble(
    parameter = p,
    P_negative = mean(posterior_draws[[p]] < 0)
  )
}) %>%
  left_join(param_map, by = "parameter")

p3 <- ggplot(prob_df, aes(x = label, y = P_negative)) +
  geom_col(fill = "grey70", color = "black", alpha = 0.8) +
  coord_flip() +
  theme_classic() +
  theme(legend.key.size = unit(0.25, 'cm'),axis.title = element_text(size = axis_title_size,family="Helvetica"),
        axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),
        plot.title = element_text(size = plot_title_size,family="Helvetica"),
        legend.title = element_text(size = legend_title_size),
        legend.text = element_text(size = legend_title_size,family="Helvetica"),
        strip.text.x = element_text(size = legend_text_size,family="Helvetica")) +
  labs(y = "Posterior probability β < 0", x = NULL)

p3


combined_plot <- ggarrange(
  ggarrange(p1,ncol = 1,labels = 'D',font.label = list(size = 11), align = 'hv',heights = 0.8),
  ggarrange(
  p2,
  p3,
  ncol = 2,labels = c('E','F'),font.label = list(size = 11), align = 'hv',heights = c(0.8,0.8)),ncol = 1,
  widths = c(2.2, 1.3, 1)
)

combined_plot

ggsave(paste(save_path,'bayes_full_model.svg'), width = 15, height = 10, units = "cm")
ggsave(paste(save_path,'bayes_full_model.png'), width = 15, height = 10, units = "cm")


# ============================================================
#Female only dataset
#=============================================================
female_df <- dat %>%
  filter(gender == "female")

female_df$subject <- factor(female_df$subject)
female_df$group <- factor(female_df$group)
female_df$exercise_intensity <- factor(female_df$exercise_intensity)
female_df$trial <- factor(female_df$trial)
female_df$block <- factor(female_df$block)

female_df <- female_df %>%
  filter(group %in% c("2","3"))

fit_null_female <- brm(
  pain_rating ~
    trial +
    block +
    (1 + exercise_intensity | subject),
  
  data = female_df,
  family = gaussian(),
  prior = priors,
  
  chains = 4,
  iter = 8000,
  warmup = 4000,
  cores = 4,
  
  control = list(
    adapt_delta = 0.995,
    max_treedepth = 15
  ),
  
  save_pars = save_pars(all = TRUE)
)

# --------------------------------
# EIH
#--------------------------------

fit_EIH_female <- brm(
  pain_rating ~
    exercise_intensity +
    trial +
    block +
    (1 + exercise_intensity | subject),
  
  data = female_df,
  family = gaussian(),
  prior = priors,
  
  chains = 4,
  iter = 6000,
  warmup = 3000,
  cores = 4,
  
  control = list(
    adapt_delta = 0.995,
    max_treedepth = 15
  ),
  
  save_pars = save_pars(all = TRUE)
)


# --------------------------------------
# Study moderation
#---------------------------------------
fit_study_female <- brm(
  pain_rating ~
    exercise_intensity * study +
    trial +
    block +
    (1 + exercise_intensity | subject),
  
  data = female_df,
  family = gaussian(),
  prior = priors,
  
  chains = 4,
  iter = 6000,
  warmup = 3000,
  cores = 4,
  
  control = list(
    adapt_delta = 0.995,
    max_treedepth = 15
  ),
  
  save_pars = save_pars(all = TRUE)
)



# ----------------------------
# Baeys Factor
#------------------------------

# EIH in females
bf_EIH_female <- bayes_factor(
  fit_EIH_female,
  fit_null_female
)

bf_EIH_female

# # sutdy moderated EIH
bf_study_female <- bayes_factor(
  fit_study_female,
  fit_EIH_female
)

bf_study_female

# Full model
bf_study_vs_null_female <- bayes_factor(
  fit_study_female,
  fit_null_female
)


#---------------
# Posteriror Distributions
#-----------------
posterior_draws <- as_draws_df(fit_study_female)
grep(
  "exercise_intensity.*group",
  names(posterior_draws),
  value = TRUE
)


study_beta <- posterior_draws$`b_exercise_intensity1:group3`

mean(study_beta)

quantile(
  study_beta,
  probs = c(.025,.5,.975)
)

mean(study_beta > 0)
mean(study_beta < 0)
bf_study_vs_null_female



