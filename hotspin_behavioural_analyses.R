######################################################################
#
#                           HotSpin
#                     Behavioural Analyses
#
#####################################################################

#----------------------------------------------------------------
# Analysis script and visualisation for cortical results
# Before running this script, the script PEEP_LOAD_IN_DATA.R 
# should be run to load in data and settings
#
#
# Copyright Janne Nold 11-01-2024 (UKE, Hamburg)
#---------------------------------------------------------------

## -------------------------
## Load in Packages
## ------------------------
rm(list=ls())

library(lintr)
library(svglite)
library(reporttools)
library(ggridges)
library(ggdist)
library(rstan)
library(ggplot2)
library(dplyr)
library(MASS)
library(pwr)
library(tidyverse)
library(rstatix)
library('lme4')
library('lmerTest')
library(ggplot2)
library(dplyr)
library(ggpubr)
library(ez)
library(reshape2)
library(ggpubr)
library(data.table)
library(patchwork)
library(viridis)
library(lm.beta)
library(ggplot2)
library(dplyr)
library(MASS)
library(pwr)
library(tidyverse)
library(rstatix)
library(lme4)
library(lmerTest)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(ez)
library(reshape2)
library(ggpubr)
library(data.table)
library(car)
library(emmeans)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(viridis)
library(lm.beta)
library(ggExtra)
library(sjPlot)
library(sjmisc)
library(sjlabelled)
library(merTools)
library(ggeffects)
library(glmmTMB)
library(stargazer)
library(sjPlot)
library(sjmisc)
library(MuMIn)
library(plotrix)
library(extrafont)
library(multcomp)
library(emmeans)
library(simr)
library(effsize)
library(robustlmm)
library(brms)
library(tidyverse)
library(tidybayes)
library(ggpubr)

packages <- c("plyr", "lattice", "ggplot2", "dplyr", "readr", 
              "ggplot2","rmarkdown","Rmisc", "tidyr", "gghalves")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}
lapply(packages, require, character.only = TRUE)

## --------------------------------
## Set Paths and Figure Sizes 
## -------------------------------
save_path <- "C:/Users/user/Desktop/projects/HotSpin/Visualisation/graphs_V4.0/"
w_size = 11
h_size = 9
w_size_l = 13
h_size_l = 11
w_size_s = 5
h_size_s = 9
axis_text_size   = 9
axis_title_size  = 10
plot_title_size  = 9
legend_text_size = 9
legend_title_size= 9

## ---------------------------------
# Load in Calibration Data
# ---------------------------------
hotspin_calib  <- read.csv('C:/Users/user/Desktop/projects/HotSpin/data/LogCalibration/MAIN/hotspin_calib_data_cleaned.csv')

# Convert data to long format for heat and pressure VAS intensities 
data_heat_long <- hotspin_calib %>%
  select(subject, h_30, h_50, h_70) %>%
  pivot_longer(cols = starts_with("h_"), names_to = "VAS_int", values_to = "intensity") %>%
  mutate(modality = 2, VAS_int = as.numeric(sub("h_", "", VAS_int)))

data_heat_long$VAS_int <- as.factor(data_heat_long$VAS_int)

data_pressure_long <- hotspin_calib %>%
  select(subject, p_30, p_50, p_70) %>%
  pivot_longer(cols = starts_with("p_"), names_to = "VAS_int", values_to = "intensity") %>%
  mutate(modality = 1, VAS_int = as.numeric(sub("p_", "", VAS_int)))

data_pressure_long$VAS_int <- as.factor(data_pressure_long$VAS_int)

# summary statistics
summary_lmer_models_main_h<- data_heat_long %>%
  group_by(VAS_int)%>%
  summarise_at(c('intensity'),mean,na.rm = T)

# summary statistics
summary_lmer_models_main_p<- data_pressure_long %>%
  group_by(VAS_int)%>%
  summarise_at(c('intensity'),mean,na.rm = T)

  # Raincloud plot for heat data
  raincloud_heat <- ggplot(data_heat_long, aes(x = VAS_int, y = intensity, fill = VAS_int)) +
    PupillometryR::geom_flat_violin(position = position_nudge(x = .2, y = 0), adjust = 1.5, trim = T, alpha = .6, size = 0.25) +
    geom_point(position = position_jitter(width = .05), size = .5, alpha = .6, shape = 21) +
    geom_boxplot(width = .1, outlier.shape = NA, alpha = .5, size = 0.25) +
    #stat_summary(fun = mean, geom = "point", shape = 23, size = .7, fill = "white") +
    theme_classic() +
  scale_fill_manual(values = c("#ecd761", '#fa9c0e', "#af0404"))+
    theme(legend.key.size = unit(0.25, 'cm'),axis.title = element_text(size = axis_title_size,family="Helvetica"),axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),plot.title = element_text(size = plot_title_size,family="Helvetica"),legend.title = element_text(size = axis_title_size,family="Helvetica"),legend.text = element_text(size = legend_title_size,family="Helvetica"),strip.text.x = element_text(size = legend_text_size,family="Helvetica")) +
    geom_segment(data=summary_lmer_models_main_h,aes(x =c(0,0,0), xend = as.factor(VAS_int), y = c(intensity), yend = c(intensity)),size = 0.5,linetype = 3, colour = "#666666")+
    ylab('Pain Intensity [?C]') + xlab('Stimulus Intensity [VAS]')  + ggtitle('Heat')

  raincloud_heat

  # Raincloud plot for pressure data
  raincloud_pressure <- ggplot(data_pressure_long, aes(x = VAS_int, y = intensity, fill = VAS_int)) +
    PupillometryR::geom_flat_violin(position = position_nudge(x = .2, y = 0), adjust = 1.5, trim = T, alpha = .6, size = 0.25) +
    geom_point(position = position_jitter(width = .05), size = .5, alpha = .6, shape = 21) +
    geom_boxplot(width = .1, outlier.shape = NA, alpha = .5, size = 0.25) +
    #stat_summary(fun = mean, geom = "point", shape = 23, size = .7, fill = "white") +
    theme_classic() +
    scale_fill_manual(values = c("#969BF2", '#1F248C', "#020659")) +
    theme(legend.key.size = unit(0.25, 'cm'),axis.title = element_text(size = axis_title_size,family="Helvetica"),axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),plot.title = element_text(size = plot_title_size,family="Helvetica"),legend.title = element_text(size = axis_title_size,family="Helvetica"),legend.text = element_text(size = legend_title_size,family="Helvetica"),strip.text.x = element_text(size = legend_text_size,family="Helvetica")) +
    geom_segment(data=summary_lmer_models_main_p,aes(x =c(0,0,0), xend = as.factor(VAS_int), y = c(intensity), yend = c(intensity)),size = 0.5,linetype = 3, colour = "#666666")+
    ylab('Pain Intensity [kPa]') + xlab('Stimulus Intensity [VAS]') + ylim(0, 100) + ggtitle('Pressure')

  raincloud_pressure


# ------------- Plot FTP (weight corrected)
# Load in data
hotspin_calib_pwc <- hotspin_calib %>%
  select(subject, pwc) %>%
  mutate(pwc = as.numeric(pwc))

# Raincloud plot for pwc data
raincloud_pwc <- ggplot(hotspin_calib_pwc, aes(x = "", y = pwc)) +
  PupillometryR::geom_flat_violin(position = position_nudge(x = .2, y = 0), adjust = 1.5, trim = T, alpha = .6, fill = "#fcac00", size = 0.25) +
  geom_point(position = position_jitter(width = .1), size = .5, alpha = .6,shape = 21, fill = "#fcac00") +
  geom_boxplot(width = .1, outlier.shape = NA, alpha = .5, fill = "#fcac00", size = 0.25) +
  #stat_summary(fun = mean, geom = "point", shape = 23, size = 2, fill = "#000000") +
  theme_classic() +
  scale_fill_manual(values = c("#fcac00")) +
   theme(legend.key.size = unit(0.25, 'cm'),axis.title = element_text(size = axis_title_size,family="Helvetica"),axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),plot.title = element_text(size = plot_title_size,family="Helvetica"),legend.title = element_text(size = axis_title_size,family="Helvetica"),legend.text = element_text(size = legend_title_size,family="Helvetica"),strip.text.x = element_text(size = legend_text_size,family="Helvetica")) +
  ylab('FTP [W/kg]') + xlab('') + ggtitle('FTP (weight-corrected) Distribution') +
  coord_flip()

raincloud_pwc

#--------------------------------------------
# Distribution of mensutral caycle phase
#----------------------------------------------
hotspin_menstr_ccle  <- read.csv('C:/Users/user/Desktop/projects/HotSpin/analyses/01_data_cleaning/hotspin_cycle_phase_data.csv')


# exclude subject 4
hotspin_menstr_ccle <- hotspin_menstr_ccle %>% filter(subject != 4)



# ------------ statistics

# Calculate chi-square test for the distribution of cycle phase across days
cycle_phase_table <- table(hotspin_menstr_ccle$cycle_phase, hotspin_menstr_ccle$day)
chi_square_test <- chisq.test(cycle_phase_table)

# Print the results
print(chi_square_test)

# caöcuölate distribution of cycle phase within both days
cycle_phase_distribution_day1 <- table(hotspin_menstr_ccle$cycle_phase[hotspin_menstr_ccle$day == 1])
cycle_phase_distribution_day2 <- table(hotspin_menstr_ccle$cycle_phase[hotspin_menstr_ccle$day == 2])

chi_square_test1 <- chisq.test(cycle_phase_distribution_day1)
chi_square_test2 <- chisq.test(cycle_phase_distribution_day2)

print(chi_square_test1)
print(chi_square_test2)


# Plot the distribution of the cycle phase for both days
cycle_phase_distribution <- ggplot(hotspin_menstr_ccle, aes(x = cycle_phase, fill = as.factor(day))) +
  geom_bar(position = "dodge", alpha = 0.7, color = 'black') +
  scale_fill_manual(values = c("#1C02C7", "#C75302"), labels = c("Calibration (Day 1)", "Experiment (Day 2)")) +
  theme_classic() +
  theme(legend.position = "bottom", legend.key.size = unit(0.25, 'cm'), axis.title = element_text(size = axis_title_size, family = "Helvetica"), axis.text = element_text(size = axis_text_size, colour = 'black', family = "Helvetica"), plot.title = element_text(size = plot_title_size, family = "Helvetica"), legend.title = element_blank(), legend.text = element_text(size = legend_text_size, family = "Helvetica"), strip.text.x = element_text(size = legend_text_size, family = "Helvetica")) +
  xlab('Cycle Phase') +
  ylab('Count') +
  ggtitle('Distribution of Cycle Phase') +
  scale_x_continuous(breaks = 0:3, labels = c("Hormonal\ncontraceptives", "Follicular", "Ovulation", "Luteal")) 
  #annotate("text", x = 3.5, y = max(table(hotspin_menstr_ccle$cycle_phase, hotspin_menstr_ccle$day)) + 1, 
  #         label = paste("Chi-squared (", chi_square_test$parameter, ") =", round(chi_square_test$statistic, 2), "p =", round(chi_square_test$p.value, 2)), 
  #         size = 2, hjust = 1)

cycle_phase_distribution



##############################################
#      Group Comparisons
###############################################

# -----------------------
# Compare FTP between samples:
#----------------------
hotspin_peep_ftp  <- read.csv('C:/Users/user/Desktop/projects/HotSpin/analyses/01_data_cleaning/df_threshold_combined.csv',header = T, sep = ',')

hotspin_peep_ftp$group <- as.factor(hotspin_peep_ftp$group)
hotspin_peep_ftp$group <- ifelse(hotspin_peep_ftp$group == 1,2,1)
hotspin_peep_ftp$group[hotspin_peep_ftp$gender == 0] <- 0

hotspin_peep_ftp$group <- as.factor(hotspin_peep_ftp$group)


summary_ftp<- hotspin_peep_ftp %>%
  group_by(group)%>%
  summarise_at(c('pwc'),mean,na.rm = T)

  sd_group1 <- sd(hotspin_peep_ftp$pwc[hotspin_peep_ftp$group == 1], na.rm = TRUE)
  print(sd_group1)

# plot raincloud plots of the pwc with hue by group
raincloud_pwc_group <- ggplot(hotspin_peep_ftp, aes(x = group, y = pwc, fill = group)) +
  PupillometryR::geom_flat_violin(position = position_nudge(x = .2, y = 0), adjust = 1.5, trim = T, alpha = .8, size = 0.25,show.legend = FALSE) +
  geom_point(position = position_jitter(width = .1), size = .5, alpha = .6, shape = 21,show.legend = FALSE) +
  geom_boxplot(width = .1, outlier.shape = NA, alpha = .5, size = 0.25,show.legend = FALSE) +
  theme_classic() +
  scale_fill_manual(values = c("#C75302", "#8480F2", "#1C02C7")) +
  theme(legend.key.size = unit(0.25, 'cm'), axis.title = element_text(size = axis_title_size, family = "Helvetica"), axis.text = element_text(size = axis_text_size, colour = 'black', family = "Helvetica"), plot.title = element_text(size = plot_title_size, family = "Helvetica"), legend.title = element_text(size = axis_title_size, family = "Helvetica"), legend.text = element_text(size = legend_text_size, family = "Helvetica"), strip.text.x = element_text(size = legend_text_size, family = "Helvetica")) +
  ylab('FTP [W/kg]') + xlab('Group') + ggtitle('FTP (weight-corrected)') +
  coord_flip()+
  geom_segment(data=summary_ftp,aes(x =c(0,0,0), xend = as.factor(group), y = c(pwc), yend = c(pwc)),size = 0.5,linetype = 3,colour = '#666666') +
  scale_x_discrete(labels = c("0" = "Males (Previous Study)", "1" = "Females (Previous Study)", "2" = "Females (Current Study)"))+
  xlab('')

raincloud_pwc_group

# compare 
t.test(hotspin_peep_ftp$pwc[hotspin_peep_ftp$group==0], hotspin_peep_ftp$pwc[hotspin_peep_ftp$group==1], paired = F, alternative = "two.sided")
t.test(hotspin_peep_ftp$pwc[hotspin_peep_ftp$group==0], hotspin_peep_ftp$pwc[hotspin_peep_ftp$group==2], paired = F, alternative = "two.sided")
t.test(hotspin_peep_ftp$pwc[hotspin_peep_ftp$group==1], hotspin_peep_ftp$pwc[hotspin_peep_ftp$group==2], paired = F, alternative = "two.sided")

#---------------------------
# Training hours per week
#---------------------------
hotspin_training_hours  <- read.csv('C:/Users/user/Desktop/projects/HotSpin/analyses/01_data_cleaning/hotspin_training_hours.csv',header = T, sep = ',')

hotspin_training_hours <- hotspin_training_hours[hotspin_training_hours$SubID != 4, ]

# Calculate mean and standard deviation of training hours
mean_training_hours <- mean(hotspin_training_hours$sum_training_h_w, na.rm = TRUE)
sd_training_hours <- sd(hotspin_training_hours$sum_training_h_w, na.rm = TRUE)

# Print the results
print(paste("Mean training hours:", mean_training_hours))
print(paste("Standard deviation of training hours:", sd_training_hours))

# combine 
data_day1 <- read.csv('C:/Users/user/Desktop/PEEP/fMRI/Data/Questionnaires/data/results-survey_day1.csv',header = TRUE,sep = ',')
training_volume <- data_day1[,c('SubID','allgbefinden14','allgbefinden15')]
exclude_subjects <- c("testtet","sub999","sub005","sub008","sub010","sub012","sub016","sub020","sub033","sub037","sub042","test","")
training_volume<- as.data.frame(subset(training_volume, !(training_volume$SubID %in% exclude_subjects)))
training_volume <- training_volume[!duplicated(training_volume),]  
training_volume <- training_volume[-c(14),]  
training_volume$SubID[33]  <- "sub041"

travol <- c(0,3,3,2.25,2,0,2.5,2.5,4.5,4.5,4.5,3,0,6,6,5,5,17.5,11,5,2,5.5,3,11,0,4.5,3.75,7.5,5.25,3.25,6.75,3.75,3,6,5,4.5,3.5,3.5,4)
training_volume$hours_week <- travol

summary_lmer_models_2<- complete_data %>%
  group_by(subject,gender)%>%
  summarise_at(c('pain_rating'),mean,na.rm = T)

training_volume$gender <- summary_lmer_models_2$gender

training_volume$SubID <- as.numeric(str_extract(training_volume$SubID, "\\d+"))

peep_training_hours <- training_volume[, c("SubID", "hours_week","gender")]
hotspin_training_hours$SubID <- hotspin_training_hours$SubID + 100
hotspin_training_hours$gender <- 1
hotspin_training_hours$group <- 2
names(hotspin_training_hours) <- c("SubID","hours_week","gender","group")
peep_training_hours$group <- ifelse(peep_training_hours$gender==0,0,1)

# Combine data frames
training_hours_combined <- rbind(peep_training_hours, hotspin_training_hours[, names(peep_training_hours)])
names(training_hours_combined) <-  c("subject","hours_week","gender","group")
training_hours_combined$group <-as.factor(training_hours_combined$group)
training_hours_combined$subject <-as.factor(training_hours_combined$subject)
training_hours_combined$gender <-as.factor(training_hours_combined$gender)

summary_h_w<- training_hours_combined %>%
  group_by(group)%>%
  summarise_at(c('hours_week'),mean,na.rm = T)

sd_h_w <- hotspin_training_hours %>%
  summarise(
    sd_hours_week = sd(hours_week, na.rm = TRUE)
  )

# plot raincloud plots of the hours_week with hue by group
raincloud_training_h_w_group <- ggplot(training_hours_combined, aes(x = group, y = hours_week, fill = group)) +
  PupillometryR::geom_flat_violin(position = position_nudge(x = .2, y = 0), adjust = 1.5, trim = T, alpha = .8, size = 0.25) +
  geom_point(position = position_jitter(width = .1), size = .5, alpha = .6, shape = 21) +
  geom_boxplot(width = .1, outlier.shape = NA, alpha = .5, size = 0.25) +
  theme_classic() +
  scale_fill_manual(values = c("#C75302", "#8480F2", "#1C02C7")) +
  theme(legend.key.size = unit(0.25, 'cm'), axis.title = element_text(size = axis_title_size, family = "Helvetica"), axis.text = element_text(size = axis_text_size, colour = 'black', family = "Helvetica"), plot.title = element_text(size = plot_title_size, family = "Helvetica"), legend.title = element_text(size = axis_title_size, family = "Helvetica"), legend.text = element_text(size = legend_text_size, family = "Helvetica"), strip.text.x = element_text(size = legend_text_size, family = "Helvetica")) +
  ylab('Training Volume [h/week]') + xlab('Group') + ggtitle('Training Volume (h/week)') +
  coord_flip()+
  geom_segment(data=summary_h_w,aes(x =c(0,0,0), xend = as.factor(group), y = c(hours_week), yend = c(hours_week)),size = 0.5,linetype = 3,colour = '#666666') +
  scale_x_discrete(labels = c("0" = "Males (Previous Study)", "1" = "Females (Previous Study)", "2" = "Females (Current Study)"))+
  xlab('')

raincloud_training_h_w_group


# compare training volumne hotspin training hours with males and females of other study
t.test(training_hours_combined[training_hours_combined$group == 2,]$hours_week, training_hours_combined[training_hours_combined$group == 1,]$hours_week, paired = F, alternative = "two.sided")
t.test( training_hours_combined[training_hours_combined$group == 2,]$hours_week,  training_hours_combined[training_hours_combined$group == 0,]$hours_week, paired = F, alternative = "two.sided")
t.test( training_hours_combined[training_hours_combined$group == 0,]$hours_week,  training_hours_combined[training_hours_combined$group == 1,]$hours_week, paired = F, alternative = "two.sided")

#-------------------------
# Compare samples Height and wight 
#---------------------------------------------------

hotspin_peep_weight_height  <- read.csv('C:/Users/user/Desktop/projects/HotSpin/analyses/01_data_cleaning/hotspin_peep_weight_height.csv',header = T, sep = ';', dec = ".", stringsAsFactors = T)
hotspin_peep_weight_height$group <- as.factor(hotspin_peep_weight_height$group)



summary_weight<- hotspin_peep_weight_height %>%
  group_by(group)%>%
  summarise_at(c('weight'),mean,na.rm = T)

  sd_group1 <- sd(hotspin_peep_weight_height$weight[hotspin_peep_weight_height$group == 1], na.rm = TRUE)
  print(sd_group1)

summary_height <- hotspin_peep_weight_height %>%
  group_by(group)%>%
  summarise_at(c('height'),mean,na.rm = T)

  sd_group1 <- sd(hotspin_peep_weight_height$height[hotspin_peep_weight_height$group == 1], na.rm = TRUE)
  print(sd_group1)


  raincloud_weight_group <- ggplot(hotspin_peep_weight_height, aes(x = group, y = weight, fill = group)) +
    PupillometryR::geom_flat_violin(position = position_nudge(x = .2, y = 0), adjust = 1.5, trim = T, alpha = .8, size = 0.25) +
    geom_point(position = position_jitter(width = .1), size = .5, alpha = .6, shape = 21) +
    geom_boxplot(width = .1, outlier.shape = NA, alpha = .5, size = 0.25) +
    theme_classic() +
    scale_fill_manual(values = c("#C75302", "#8480F2", "#1C02C7")) +
    theme(legend.key.size = unit(0.25, 'cm'), axis.title = element_text(size = axis_title_size, family = "Helvetica"), axis.text = element_text(size = axis_text_size, colour = 'black', family = "Helvetica"), plot.title = element_text(size = plot_title_size, family = "Helvetica"), legend.title = element_text(size = axis_title_size, family = "Helvetica"), legend.text = element_text(size = legend_text_size, family = "Helvetica"), strip.text.x = element_text(size = legend_text_size, family = "Helvetica")) +
    ylab('Weight (kg)') + xlab('Group') + ggtitle('Weight (kg)') +
    coord_flip()+
    geom_segment(data=summary_weight,aes(x =c(0,0,0), xend = as.factor(group), y = c(weight), yend = c(weight)),size = 0.5,linetype = 3,colour = '#666666') +
    scale_x_discrete(labels = c("0" = "Males (Previous Study)", "1" = "Females (Previous Study)", "2" = "Females (Current Study)"))+
    xlab('')
  
  raincloud_weight_group
  
  


# calculate two sample t test for weight
t.test(hotspin_peep_weight_height$weight[hotspin_peep_weight_height$group == 2], hotspin_peep_weight_height$weight[hotspin_peep_weight_height$group == 1], paired = F, alternative = "two.sided")
cohen.d(hotspin_peep_weight_height$weight[hotspin_peep_weight_height$group == 1], hotspin_peep_weight_height$weight[hotspin_peep_weight_height$group == 2], paired=F, na.rm = T)



# plot raincloud plots of the pwc with hue by group
raincloud_height_group <- ggplot(hotspin_peep_weight_height, aes(x = group, y = height, fill = group)) +
  PupillometryR::geom_flat_violin(position = position_nudge(x = .2, y = 0), adjust = 1.5, trim = T, alpha = .8, size = 0.25) +
  geom_point(position = position_jitter(width = .1), size = .5, alpha = .6, shape = 21) +
  geom_boxplot(width = .1, outlier.shape = NA, alpha = .5, size = 0.25) +
  theme_classic() +
  theme(legend.key.size = unit(0.25, 'cm'), axis.title = element_text(size = axis_title_size, family = "Helvetica"), axis.text = element_text(size = axis_text_size, colour = 'black', family = "Helvetica"), plot.title = element_text(size = plot_title_size, family = "Helvetica"), legend.title = element_text(size = axis_title_size, family = "Helvetica"), legend.text = element_text(size = legend_text_size, family = "Helvetica"), strip.text.x = element_text(size = legend_text_size, family = "Helvetica")) +
   ylab('Height (cm)') + xlab('Group') + ggtitle('Height (cm)') +
  scale_fill_manual(values = c("#C75302", "#8480F2", "#1C02C7")) +
  coord_flip()+
  geom_segment(data=summary_height,aes(x =c(0,0,0), xend = as.factor(group), y = c(height), yend = c(height)),size = 0.5,linetype = 3,colour = '#666666') +
  scale_x_discrete(labels = c("0" = "Males (Previous Study)", "1" = "Females (Previous Study)", "2" = "Females (Current Study)"))+
  xlab('')

raincloud_height_group

# calculate two sample t test for height
t.test(hotspin_peep_weight_height$height[hotspin_peep_weight_height$group == 2], hotspin_peep_weight_height$height[hotspin_peep_weight_height$group == 1], paired = F, alternative = "two.sided")


# Calculate BMI for each subject
hotspin_peep_weight_height <- hotspin_peep_weight_height %>%
  mutate(bmi = weight / (height/100)^2)

# Calculate mean and standard deviation of BMI for each group
summary_bmi <- hotspin_peep_weight_height %>%
  group_by(group)%>%
  summarise_at(c('bmi'),mean,na.rm = T)

sd(hotspin_peep_weight_height$bmi[hotspin_peep_weight_height$group == 1])
sd(hotspin_peep_weight_height$bmi[hotspin_peep_weight_height$group == 2])

# calculate two sample t test for height
t.test(hotspin_peep_weight_height$bmi[hotspin_peep_weight_height$group == 1], hotspin_peep_weight_height$bmi[hotspin_peep_weight_height$group == 2], paired = F, alternative = "two.sided")
cohen.d(hotspin_peep_weight_height$bmi[hotspin_peep_weight_height$group == 1], hotspin_peep_weight_height$bmi[hotspin_peep_weight_height$group == 2], paired=F, na.rm = T)

# height

#-------------------------
# Compare samples age (females only)
#---------------------------------------------------


# Filter data for gender = 1
#hotspin_peep_age_female <- hotspin_peep_age %>% filter(gender == 1)

# Calculate summary statistics for age
summary_age <- hotspin_peep_weight_height %>%
  group_by(group) %>%
  summarise_at(c('age'), mean, na.rm = TRUE)

  # Calculate age range for group 1
  age_range_group1 <- range(hotspin_peep_weight_height$age[hotspin_peep_weight_height$group == 2], na.rm = TRUE)
  print(age_range_group1)

  
  # plot raincloud plots of the pwc with hue by group
  raincloud_age_group <- ggplot(hotspin_peep_weight_height, aes(x = group, y = age, fill = group)) +
    PupillometryR::geom_flat_violin(position = position_nudge(x = .2, y = 0), adjust = 1.5, trim = T, alpha = .8, size = 0.25) +
    geom_point(position = position_jitter(width = .1), size = .5, alpha = .6, shape = 21) +
    geom_boxplot(width = .1, outlier.shape = NA, alpha = .5, size = 0.25) +
    theme_classic() +
    scale_fill_manual(values = c("#C75302", "#8480F2", "#1C02C7")) +
    theme(legend.key.size = unit(0.25, 'cm'), axis.title = element_text(size = axis_title_size, family = "Helvetica"), axis.text = element_text(size = axis_text_size, colour = 'black', family = "Helvetica"), plot.title = element_text(size = plot_title_size, family = "Helvetica"), legend.title = element_text(size = axis_title_size, family = "Helvetica"), legend.text = element_text(size = legend_text_size, family = "Helvetica"), strip.text.x = element_text(size = legend_text_size, family = "Helvetica")) +
    ylab('Age (years)') + xlab('Group') + ggtitle('Age (years)') +
    coord_flip()+
    geom_segment(data=summary_age,aes(x =c(0,0,0), xend = as.factor(group), y = c(age), yend = c(age)),size = 0.5,linetype = 3,colour = '#666666') +
    scale_x_discrete(labels = c("0" = "Males (Previous Study)", "1" = "Females (Previous Study)", "2" = "Females (Current Study)"))+
    xlab('')
  
  raincloud_age_group

# two sample t test for age
t.test(hotspin_peep_weight_height$age[hotspin_peep_weight_height$group == 2], hotspin_peep_weight_height$age[hotspin_peep_weight_height$group == 1], paired = F, alternative = "two.sided")
t.test(hotspin_peep_weight_height$age[hotspin_peep_weight_height$group == 2], hotspin_peep_weight_height$age[hotspin_peep_weight_height$group == 0], paired = F, alternative = "two.sided")


##--------------------
## Load exercise data
##------------------------
hotspin_complete_data_hr_watt <- read.csv('C:/Users/user/Desktop/projects/HotSpin/analyses/01_data_cleaning/exercise_summary.csv')

##------------------------
## Exercise Power
##------------------------

#---------------- Statistics 
summary_data_watt <- hotspin_complete_data_hr_watt %>%
  group_by(subject,exercise_intensity)%>%
  summarise_at(c('watt_mean'),mean,na.rm = T) %>%
    # Convert exercise_intensity to integers
  mutate(exercise_intensity = ifelse(exercise_intensity == "Low", 0, 1))

names(summary_data_watt) <- c('subject','exercise_intensity','watt')

t.test(summary_data_watt$watt[summary_data_watt$exercise_intensity == 1],summary_data_watt$watt[summary_data_watt$exercise_intensity == 0],paired = T,alternative = "two.sided")
cohen.d(summary_data_watt$watt[summary_data_watt$exercise_intensity == 1],summary_data_watt$watt[summary_data_watt$exercise_intensity == 0],paired=TRUE,na.rm = T)


#------------------------Visualisation
plot_data<- summary_data_watt
plot_data_add <- spread(plot_data, exercise_intensity, watt)
plot_data_add$diff_rating <- plot_data_add$'1' -  plot_data_add$'0'  # Low Intensity - High Intensity Data
mean_diff = abs(mean(plot_data_add$diff_rating))
#mean_diff <- c(mean_diff,mean_diff)
se_diff = plotrix::std.error(plot_data_add$diff_rating)
se_diff <- c(se_diff,se_diff)

###################### CHANGE THINGS HERE FOR YOUR DATA #############################
# DATA: just an example data set...use your own here!
df <-  plot_data

# Data peparation ... change variable names here to match thos in your data frame 
dataLeft    <-  plot_data$watt[plot_data$exercise_intensity == 0]
dataRight   <-  plot_data$watt[plot_data$exercise_intensity== 1]

# Change y axis range
y_lim_min   <-  0
y_lim_max   <-  250

# Change colors for plot

leftColor         <-  "#005C53"
rightColor        <-  "#3d0092"
singleLineColor   <-  "gray"
meanLineColor     <-  "grey50"



####################### ACTUAL PLOT: no changes needed, except axis labels etc. #######
n <- length(dataLeft)
d <- data.frame(y = c(dataLeft, dataRight),
                x = rep(c(1,2), each=n),
                id = factor(rep(1:n,2)))

set.seed(321)
d$xj <- jitter(d$x, amount=.09)

score_mean_1 <- mean(d$y[d$x ==1])
score_mean_2 <- mean(d$y[d$x == 2])
score_median1 <- median(d$y[d$x ==1])
score_median2 <- median(d$y[d$x == 2])
score_sd_1 <- sd(d$y[d$x ==1])
score_sd_2 <- sd(d$y[d$x == 2])
score_se_1 <- score_sd_1/sqrt(n) 
score_se_2 <- score_sd_2/sqrt(n) 
score_ci_1 <- CI(d$y[d$x ==1], ci = 0.95)
score_ci_2 <- CI(d$y[d$x == 2], ci = 0.95)
#Create data frame with 2 rows and 7 columns containing the descriptives
group <- c("x", "z")
N <- c(n, n)
score_mean <- c(score_mean_1, score_mean_2)
score_median <- c(score_median1, score_median2)
sd <- c(score_sd_1, score_sd_2)
se <- c(score_se_1, score_se_2)
ci <- c((score_ci_1[1] - score_ci_1[3]), (score_ci_2[1] - score_ci_2[3]))
#Create the dataframe
summary_df <- data.frame(group, N, score_mean, score_median, sd, se, ci,se_diff)

x_tick_means <- c(.87, 2.13)
x_tick_diff_bar <- 1.5

power_cycling <- ggplot(data = d, aes(y = y)) +
  
  #Add geom_() objects
  geom_point(data = d %>% filter(x =="1"), aes(x = xj), color = leftColor, size = 0.5, 
             alpha = .6) +
  geom_point(data = d %>% filter(x =="2"), aes(x = xj), color = rightColor, size = 0.5, 
             alpha = .6) +
  geom_line(aes(x = xj, group = id), color = singleLineColor, alpha = .3) +
  
  geom_half_boxplot(
    data = d %>% filter(x=="1"), aes(x=x, y = y), position = position_nudge(x = -.28), 
    side = "r",outlier.shape = NA, center = TRUE, errorbar.draw = FALSE, width = .2, 
    fill = leftColor,alpha=0.9) +
  
  geom_half_boxplot(
    data = d %>% filter(x=="2"), aes(x=x, y = y), position = position_nudge(x = .18), 
    side = "r",outlier.shape = NA, center = TRUE, errorbar.draw = FALSE, width = .2, 
    fill = rightColor,alpha = 0.9) +
  
  geom_half_violin(
    data = d %>% filter(x=="1"),aes(x = x, y = y), position = position_nudge(x = -.3), 
    side = "l", fill = leftColor,alpha = 0.9) +
  
  geom_half_violin(
    data = d %>% filter(x=="2"),aes(x = x, y = y), position = position_nudge(x = .3), 
    side = "r", fill = rightColor,alpha = 0.9) +
  
  #Add a line connecting the two means
  geom_line(data = summary_df, aes(x = x_tick_means, y = score_mean), color = meanLineColor, 
            size = 0.5) +
  
  geom_point(data = summary_df, aes(x = c(1,2), y = score_mean), position = position_nudge(x = c(-.13,.13)), color = c(leftColor,rightColor), alpha = .6, size = 1) +
  
  geom_errorbar(data = summary_df, aes(x = c(1,2), y = score_mean, 
                                       ymin = score_mean-se, ymax = score_mean+se),position = position_nudge(x = c(-.13,.13)),
                color = c(leftColor,rightColor), width = 0.05, size = 0.4, alpha = .6)+
  
  
  scale_x_continuous(breaks=c(1,2), labels=c("Low", "High"), limits=c(0, 3)) +
  xlab("Exercise Intensity") + ylab("Power [Watt]") +
  theme_classic()+
  theme(axis.title = element_text(size = axis_title_size,family="Helvetica"),axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),plot.title = element_text(size = plot_title_size,family="Helvetica"),legend.title = element_text(size = axis_title_size,family="Helvetica"),legend.text = element_text(size = legend_title_size,family="Helvetica"),strip.text.x = element_text(size = legend_text_size,family="Helvetica")) +
  coord_cartesian(ylim=c(11, y_lim_max))+
  ggtitle("")

power_cycling

##------------------------
## Heart Rate 
##------------------------

#---------------- Statistics 
summary_data_hr_plot <- hotspin_complete_data_hr_watt %>%
  group_by(subject,exercise_intensity)%>%
  summarise_at(c('hr_mean'),mean,na.rm = T) %>%
    # Convert exercise_intensity to integers
  mutate(exercise_intensity = ifelse(exercise_intensity == "Low", 0, 1))

names(summary_data_hr_plot) <- c('subject','exercise_intensity','hr')

t.test(summary_data_hr_plot$hr[summary_data_hr_plot$exercise_intensity == 1],summary_data_hr_plot$hr[summary_data_hr_plot$exercise_intensity == 0],paired = T,alternative = "two.sided")
cohen.d(summary_data_hr_plot$hr[summary_data_hr_plot$exercise_intensity == 1],summary_data_hr_plot$hr[summary_data_hr_plot$exercise_intensity == 0],paired=TRUE,na.rm = T)


plot_data<- summary_data_hr_plot
plot_data_add <- spread(plot_data, exercise_intensity, hr)
plot_data_add$diff_rating <- plot_data_add$'1' -  plot_data_add$'0'  # Low Intensity - High Intensity Data
mean_diff = abs(mean(plot_data_add$diff_rating))
#mean_diff <- c(mean_diff,mean_diff)
se_diff = plotrix::std.error(plot_data_add$diff_rating)
se_diff <- c(se_diff,se_diff)

###################### CHANGE THINGS HERE FOR YOUR DATA #############################
# DATA: just an example data set...use your own here!
df <-  plot_data

# Data peparation ... change variable names here to match thos in your data frame 
dataLeft    <-  plot_data$hr[plot_data$exercise_intensity == 0]
dataRight   <-  plot_data$hr[plot_data$exercise_intensity== 1]

# Change y axis range
y_lim_min   <-  0
y_lim_max   <-  200

# Change colors for plot
leftColor         <-  "#005C53"
rightColor        <-  "#3C008E"
singleLineColor   <-  "gray"
meanLineColor     <-  "grey50"



####################### ACTUAL PLOT: no changes needed, except axis labels etc. #######
n <- length(dataLeft)
d <- data.frame(y = c(dataLeft, dataRight),
                x = rep(c(1,2), each=n),
                id = factor(rep(1:n,2)))

set.seed(321)
d$xj <- jitter(d$x, amount=.09)

score_mean_1 <- mean(d$y[d$x ==1])
score_mean_2 <- mean(d$y[d$x == 2])
score_median1 <- median(d$y[d$x ==1])
score_median2 <- median(d$y[d$x == 2])
score_sd_1 <- sd(d$y[d$x ==1])
score_sd_2 <- sd(d$y[d$x == 2])
score_se_1 <- score_sd_1/sqrt(n) 
score_se_2 <- score_sd_2/sqrt(n) 
score_ci_1 <- CI(d$y[d$x ==1], ci = 0.95)
score_ci_2 <- CI(d$y[d$x == 2], ci = 0.95)
#Create data frame with 2 rows and 7 columns containing the descriptives
group <- c("x", "z")
N <- c(n, n)
score_mean <- c(score_mean_1, score_mean_2)
score_median <- c(score_median1, score_median2)
sd <- c(score_sd_1, score_sd_2)
se <- c(score_se_1, score_se_2)
ci <- c((score_ci_1[1] - score_ci_1[3]), (score_ci_2[1] - score_ci_2[3]))
#Create the dataframe
summary_df <- data.frame(group, N, score_mean, score_median, sd, se, ci,se_diff)

x_tick_means <- c(.87, 2.13)
x_tick_diff_bar <- 1.5

hr_cycling <- ggplot(data = d, aes(y = y)) +
  
  #Add geom_() objects
  geom_point(data = d %>% filter(x =="1"), aes(x = xj), color = leftColor, size = 0.5, 
             alpha = .6) +
  geom_point(data = d %>% filter(x =="2"), aes(x = xj), color = rightColor, size = 0.5, 
             alpha = .6) +
  geom_line(aes(x = xj, group = id), color = singleLineColor, alpha = .3) +
  
  geom_half_boxplot(
    data = d %>% filter(x=="1"), aes(x=x, y = y), position = position_nudge(x = -.28), 
    side = "r",outlier.shape = NA, center = TRUE, errorbar.draw = FALSE, width = .2, 
    fill = leftColor,alpha=0.9) +
  
  geom_half_boxplot(
    data = d %>% filter(x=="2"), aes(x=x, y = y), position = position_nudge(x = .18), 
    side = "r",outlier.shape = NA, center = TRUE, errorbar.draw = FALSE, width = .2, 
    fill = rightColor,alpha = 0.9) +
  
  geom_half_violin(
    data = d %>% filter(x=="1"),aes(x = x, y = y), position = position_nudge(x = -.3), 
    side = "l", fill = leftColor,alpha = 0.9) +
  
  geom_half_violin(
    data = d %>% filter(x=="2"),aes(x = x, y = y), position = position_nudge(x = .3), 
    side = "r", fill = rightColor,alpha = 0.9) +
  
  #Add a line connecting the two means
  geom_line(data = summary_df, aes(x = x_tick_means, y = score_mean), color = meanLineColor, 
            size = 0.5) +
  
  geom_point(data = summary_df, aes(x = c(1,2), y = score_mean), position = position_nudge(x = c(-.13,.13)), color = c(leftColor,rightColor), alpha = .6, size = 1) +
  
  geom_errorbar(data = summary_df, aes(x = c(1,2), y = score_mean, 
                                       ymin = score_mean-se, ymax = score_mean+se),position = position_nudge(x = c(-.13,.13)),
                color = c(leftColor,rightColor), width = 0.05, size = 0.4, alpha = .6)+
  
  
  scale_x_continuous(breaks=c(1,2), labels=c("Low", "High"), limits=c(0, 3)) +
  xlab("Exercise Intensity") + ylab("Heart Rate [bpm]") +
  theme_classic()+
  theme(axis.title = element_text(size = axis_title_size,family="Helvetica"),axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),plot.title = element_text(size = plot_title_size,family="Helvetica"),legend.title = element_text(size = axis_title_size,family="Helvetica"),legend.text = element_text(size = legend_title_size,family="Helvetica"),strip.text.x = element_text(size = legend_text_size,family="Helvetica")) +
  coord_cartesian(ylim=c(11, y_lim_max))+
  ggtitle("")

hr_cycling


##----------------------
## Exercise Ratings
## ------------------------

# --------- Statistics 

exercise_data_plot <- hotspin_complete_data_hr_watt %>%
  group_by(subject,exercise_intensity)%>%
  summarise_at(c('exercise_rating_mean'),mean,na.rm = T) %>%
    # Convert exercise_intensity to integers
  mutate(exercise_intensity = ifelse(exercise_intensity == "Low", 0, 1))

names(exercise_data_plot) <- c('subject','exercise_intensity','exercise_rating')


# By Intensity
low_data <- exercise_data_plot[exercise_data_plot$exercise_intensity==0,]
high_data <- exercise_data_plot[exercise_data_plot$exercise_intensity==1,]

# paired t test
t.test(exercise_data_plot$exercise_rating[exercise_data_plot$exercise_intensity == 1],exercise_data_plot$exercise_rating[exercise_data_plot$exercise_intensity == 0],paired = T,alternative = "two.sided")
cohen.d(exercise_data_plot$exercise_rating[exercise_data_plot$exercise_intensity == 1],exercise_data_plot$exercise_rating[exercise_data_plot$exercise_intensity == 0],paired=TRUE)


# ---------- RAINCLOUD PLOT ------------

plot_data<- exercise_data_plot
plot_data_add <- spread(plot_data, exercise_intensity, exercise_rating)
plot_data_add$diff_rating <- plot_data_add$'1' -  plot_data_add$'0'  # Low Intensity - High Intensity Data
mean_diff = abs(mean(plot_data_add$diff_rating))
se_diff = plotrix::std.error(plot_data_add$diff_rating)
se_diff <- c(se_diff,se_diff)

###################### CHANGE THINGS HERE FOR YOUR DATA #############################
# DATA: just an example data set...use your own here!
df <-  plot_data

# Data peparation ... change variable names here to match thos in your data frame 
dataLeft    <-  plot_data$exercise_rating[plot_data$exercise_intensity == 0]
dataRight   <-  plot_data$exercise_rating[plot_data$exercise_intensity== 1]

# Change y axis range
y_lim_min   <-  6
y_lim_max   <-  20

# Change colors for plot
leftColor         <-  "#005C53"
rightColor        <-  "#3C008E"
singleLineColor   <-  "gray"
meanLineColor     <-  "grey50"



####################### ACTUAL PLOT: no changes needed, except axis labels etc. #######
n <- length(dataLeft)
d <- data.frame(y = c(dataLeft, dataRight),
                x = rep(c(1,2), each=n),
                id = factor(rep(1:n,2)))

set.seed(321)
d$xj <- jitter(d$x, amount=.09)

score_mean_1 <- mean(d$y[d$x ==1])
score_mean_2 <- mean(d$y[d$x == 2])
score_median1 <- median(d$y[d$x ==1])
score_median2 <- median(d$y[d$x == 2])
score_sd_1 <- sd(d$y[d$x ==1])
score_sd_2 <- sd(d$y[d$x == 2])
score_se_1 <- score_sd_1/sqrt(n) 
score_se_2 <- score_sd_2/sqrt(n) 
score_ci_1 <- CI(d$y[d$x ==1], ci = 0.95)
score_ci_2 <- CI(d$y[d$x == 2], ci = 0.95)
#Create data frame with 2 rows and 7 columns containing the descriptives
group <- c("x", "z")
N <- c(n, n)
score_mean <- c(score_mean_1, score_mean_2)
score_median <- c(score_median1, score_median2)
sd <- c(score_sd_1, score_sd_2)
se <- c(score_se_1, score_se_2)
ci <- c((score_ci_1[1] - score_ci_1[3]), (score_ci_2[1] - score_ci_2[3]))
#Create the dataframe
summary_df <- data.frame(group, N, score_mean, score_median, sd, se, ci,se_diff)

x_tick_means <- c(.87, 2.13)
x_tick_diff_bar <- 1.5

borg_rating <- ggplot(data = d, aes(y = y)) +
  
  #Add geom_() objects
  geom_point(data = d %>% filter(x =="1"), aes(x = xj), color = leftColor, size = 0.5, 
             alpha = .6) +
  geom_point(data = d %>% filter(x =="2"), aes(x = xj), color = rightColor, size = 0.5, 
             alpha = .6) +
  geom_line(aes(x = xj, group = id), color = singleLineColor, alpha = .3) +
  
  geom_half_boxplot(
    data = d %>% filter(x=="1"), aes(x=x, y = y), position = position_nudge(x = -.28), 
    side = "r",outlier.shape = NA, center = TRUE, errorbar.draw = FALSE, width = .2, 
    fill = leftColor,alpha=0.9) +
  
  geom_half_boxplot(
    data = d %>% filter(x=="2"), aes(x=x, y = y), position = position_nudge(x = .18), 
    side = "r",outlier.shape = NA, center = TRUE, errorbar.draw = FALSE, width = .2, 
    fill = rightColor,alpha = 0.9) +
  
  geom_half_violin(
    data = d %>% filter(x=="1"),aes(x = x, y = y), position = position_nudge(x = -.3), 
    side = "l", fill = leftColor,alpha = 0.9) +
  
  geom_half_violin(
    data = d %>% filter(x=="2"),aes(x = x, y = y), position = position_nudge(x = .3), 
    side = "r", fill = rightColor,alpha = 0.9) +
  
  #Add a line connecting the two means
  geom_line(data = summary_df, aes(x = x_tick_means, y = score_mean), color = meanLineColor, 
            size = 0.5) +
  
  geom_point(data = summary_df, aes(x = c(1,2), y = score_mean), position = position_nudge(x = c(-.13,.13)), color = c(leftColor,rightColor), alpha = .6, size = 1) +
  
  geom_errorbar(data = summary_df, aes(x = c(1,2), y = score_mean, 
                                       ymin = score_mean-se, ymax = score_mean+se),position = position_nudge(x = c(-.13,.13)),
                color = c(leftColor,rightColor), width = 0.05, size = 0.4, alpha = .6)+
  
  
  scale_x_continuous(breaks=c(1,2), labels=c("Low", "High"), limits=c(0, 3)) +
  xlab("Exercise Intensity") + ylab("Rating of Perceived Exertion\n[BORG 6-20]") +
  theme_classic()+
  theme(axis.title = element_text(size = axis_title_size,family="Helvetica"),axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),plot.title = element_text(size = plot_title_size,family="Helvetica"),legend.title = element_text(size = axis_title_size,family="Helvetica"),legend.text = element_text(size = legend_title_size,family="Helvetica"),strip.text.x = element_text(size = legend_text_size,family="Helvetica")) +
  coord_cartesian(ylim=c(6, y_lim_max))+
  ggtitle("")

borg_rating


##----------------------------------------
## Calculate relative maintained power during high intensitiy relative to FTP
##----------------------------------------

summary_data_watt <- hotspin_complete_data_hr_watt %>%
  group_by(subject,exercise_intensity,adj_ftp)%>%
  summarise_at(c('watt_mean'),mean,na.rm = T)

# calculate relative maintained power during low and high intensity relative to FTP
summary_data_watt$rel_power <- summary_data_watt$watt_mean/summary_data_watt$adj_ftp

summary_data_watt <- summary_data_watt %>%
  group_by(subject,exercise_intensity)%>%
  summarise_at(c('rel_power'),mean,na.rm = T)

  summary_data_watt <- summary_data_watt %>%
    mutate(exercise_intensity = ifelse(exercise_intensity == "High", 1, 0))
#---------------- Statistics 

t.test(summary_data_watt$rel_power[summary_data_watt$exercise_intensity == 1],summary_data_watt$rel_power[summary_data_watt$exercise_intensity == 0],paired = T,alternative = "two.sided")
cohen.d(summary_data_watt$rel_power[summary_data_watt$exercise_intensity == 1],summary_data_watt$rel_power[summary_data_watt$exercise_intensity == 0],paired=TRUE)


#------------------------Visualisation
plot_data<- summary_data_watt
plot_data_add <- spread(plot_data, exercise_intensity, rel_power)
plot_data_add$diff_rating <- plot_data_add$'1' -  plot_data_add$'0'  # Low Intensity - High Intensity Data
mean_diff = abs(mean(plot_data_add$diff_rating))
#mean_diff <- c(mean_diff,mean_diff)
se_diff = plotrix::std.error(plot_data_add$diff_rating)
se_diff <- c(se_diff,se_diff)

###################### CHANGE THINGS HERE FOR YOUR DATA #############################
# DATA: just an example data set...use your own here!
df <-  plot_data

# Data peparation ... change variable names here to match thos in your data frame 
dataLeft    <-  plot_data$rel_power[plot_data$exercise_intensity == 0]*100
dataRight   <-  plot_data$rel_power[plot_data$exercise_intensity== 1]*100

# Change y axis range
y_lim_min   <-  min(dataLeft)
y_lim_max   <-  max(dataRight)

# Change colors for plot

leftColor         <-  "#005C53"
rightColor        <-  "#3C008E"
singleLineColor   <-  "gray"
meanLineColor     <-  "grey50"



####################### ACTUAL PLOT: no changes needed, except axis labels etc. #######
n <- length(dataLeft)
d <- data.frame(y = c(dataLeft, dataRight),
                x = rep(c(1,2), each=n),
                id = factor(rep(1:n,2)))

set.seed(321)
d$xj <- jitter(d$x, amount=.09)

score_mean_1 <- mean(d$y[d$x ==1])
score_mean_2 <- mean(d$y[d$x == 2])
score_median1 <- median(d$y[d$x ==1])
score_median2 <- median(d$y[d$x == 2])
score_sd_1 <- sd(d$y[d$x ==1])
score_sd_2 <- sd(d$y[d$x == 2])
score_se_1 <- score_sd_1/sqrt(n) 
score_se_2 <- score_sd_2/sqrt(n) 
score_ci_1 <- CI(d$y[d$x ==1], ci = 0.95)
score_ci_2 <- CI(d$y[d$x == 2], ci = 0.95)
#Create data frame with 2 rows and 7 columns containing the descriptives
group <- c("x", "z")
N <- c(n, n)
score_mean <- c(score_mean_1, score_mean_2)
score_median <- c(score_median1, score_median2)
sd <- c(score_sd_1, score_sd_2)
se <- c(score_se_1, score_se_2)
ci <- c((score_ci_1[1] - score_ci_1[3]), (score_ci_2[1] - score_ci_2[3]))
#Create the dataframe
summary_df <- data.frame(group, N, score_mean, score_median, sd, se, ci,se_diff)

x_tick_means <- c(.87, 2.13)
x_tick_diff_bar <- 1.5

rel_power_cycling <- ggplot(data = d, aes(y = y)) +
  
  #Add geom_() objects
  geom_point(data = d %>% filter(x =="1"), aes(x = xj), color = leftColor, size = 0.5, 
             alpha = .6) +
  geom_point(data = d %>% filter(x =="2"), aes(x = xj), color = rightColor, size = 0.5, 
             alpha = .6) +
  geom_line(aes(x = xj, group = id), color = singleLineColor, alpha = .3) +
  
  geom_half_boxplot(
    data = d %>% filter(x=="1"), aes(x=x, y = y), position = position_nudge(x = -.28), 
    side = "r",outlier.shape = NA, center = TRUE, errorbar.draw = FALSE, width = .2, 
    fill = leftColor,alpha=0.9) +
  
  geom_half_boxplot(
    data = d %>% filter(x=="2"), aes(x=x, y = y), position = position_nudge(x = .18), 
    side = "r",outlier.shape = NA, center = TRUE, errorbar.draw = FALSE, width = .2, 
    fill = rightColor,alpha = 0.9) +
  
  geom_half_violin(
    data = d %>% filter(x=="1"),aes(x = x, y = y), position = position_nudge(x = -.3), 
    side = "l", fill = leftColor,alpha = 0.9) +
  
  geom_half_violin(
    data = d %>% filter(x=="2"),aes(x = x, y = y), position = position_nudge(x = .3), 
    side = "r", fill = rightColor,alpha = 0.9) +
  
  #Add a line connecting the two means
  geom_line(data = summary_df, aes(x = x_tick_means, y = score_mean), color = meanLineColor, 
            size = 0.5) +
  
  geom_point(data = summary_df, aes(x = c(1,2), y = score_mean), position = position_nudge(x = c(-.13,.13)), color = c(leftColor,rightColor), alpha = .6, size = 1) +
  
  geom_errorbar(data = summary_df, aes(x = c(1,2), y = score_mean, 
                                       ymin = score_mean-se, ymax = score_mean+se),position = position_nudge(x = c(-.13,.13)),
                color = c(leftColor,rightColor), width = 0.05, size = 0.4, alpha = .6)+
  
  
  scale_x_continuous(breaks=c(1,2), labels=c("Low", "High"), limits=c(0, 3)) +
  xlab("Exercise Intensity") + ylab("relative Power [%FTP]") +
  theme_classic()+
  theme(axis.title = element_text(size = axis_title_size,family="Helvetica"),axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),plot.title = element_text(size = plot_title_size,family="Helvetica"),legend.title = element_text(size = axis_title_size,family="Helvetica"),legend.text = element_text(size = legend_title_size,family="Helvetica"),strip.text.x = element_text(size = legend_text_size,family="Helvetica")) +
  coord_cartesian(ylim=c(y_lim_min, y_lim_max))+
  ggtitle("")

rel_power_cycling


# ======================================================
#   Prepare Dataframes PEEP and HotSpin
# =====================================================

# ---- PEEP (original Study)
peep_df<- complete_data_sal %>%
  group_by(subject,VAS,exercise_intensity,pwc,gender,modality,treatment_order,nr_pain_rating,exercise_block)%>%
  summarise_at(c('pain_rating'),mean,na.rm = T)

names(peep_df) <- c("subject","VAS","exercise_intensity","pwc","gender","modality","treatment_order","trial","block","pain_rating")
peep_df$group <- ifelse(peep_df$gender==0,1,2)
peep_df_heat <- peep_df[peep_df$modality == 2,]
peep_df$sample <- 1
peep_df$sample <- as.factor(peep_df$sample)
peep_df$trial <- as.factor(peep_df$trial)

#----- Hotspin Max Pain ratings and FTP

# Load in data
sub_online_ratings <- read.csv('C:/Users/user/Desktop/projects/HotSpin/Data/LogExperiment/MAIN/sub_online_heat_and_pressure_clean.csv',sep = ',',header = T)

# calculate the mean time courses for each subject
sub_mean_ratings <- sub_online_ratings %>%
  group_by(subject,VAS_intensity,exercise_intensity,pwc,modality,time,trial,block)%>%
  summarise_at(c('online_rating'),mean,na.rm = T)

# Calculate the max rating for each exercise intensity, VAS intensity, modality, and subject
sub_max_online_ratings <- sub_mean_ratings %>%
  group_by(subject,VAS_intensity,exercise_intensity,pwc,modality,trial,block) %>%
  summarise_at(c('online_rating'),max,na.rm = T)

# Replace Inf values with NA
sub_max_online_ratings$online_rating[is.infinite(sub_max_online_ratings$online_rating)] <- NA

# Rows with missing pwc
na_pwc_rows <- sub_max_online_ratings %>%
  filter(is.na(pwc))

# Unique subjects with missing pwc
na_pwc_subjects <- unique(na_pwc_rows$subject)

na_pwc_subjects

sub_max_online_ratings$pwc[
  sub_max_online_ratings$subject == 24 &
    is.na(sub_max_online_ratings$pwc)
] <- 1.47


# Reassign
hotspin_df <- sub_max_online_ratings

hotspin_df$group <- 3
hotspin_df$gender <- 1
hotspin_df$sample <- 2

# add artifical treatment order as 0 for stats
hotspin_df$treatment_order <- 0

# rename VAS and pain ratings column
names(hotspin_df) <-  c("subject", "VAS", "exercise_intensity", "pwc", "modality","trial","block","pain_rating", "group", "gender","sample","treatment_order")
hotspin_df$VAS <- as.numeric(as.character(hotspin_df$VAS))
hotspin_df$treatment_order <- as.factor(hotspin_df$treatment_order)
hotspin_df$gender <- as.factor(hotspin_df$gender)
hotspin_df$subject <- as.factor(hotspin_df$subject)
hotspin_df$sample <- as.factor(hotspin_df$sample)

# change subject number
hotspin_df$subject <- as.numeric(as.character(hotspin_df$subject))
hotspin_df$subject <-  hotspin_df$subject + 100
hotspin_df$subject <- as.factor(hotspin_df$subject)

hotspin_df$exercise_intensity <- as.factor(hotspin_df$exercise_intensity)
hotspin_df$gender <- as.factor(hotspin_df$gender)
hotspin_df$treatment_order <- as.factor(hotspin_df$treatment_order)
hotspin_df$modality <- as.factor(hotspin_df$modality)
hotspin_df$trial <- as.factor(hotspin_df$trial)
hotspin_df$block <- as.factor(hotspin_df$block)

# ----------- rescale pain rating Hotspin for comparability to 0 - 100 (instead of ) ----
hotspin_df$pain_rating_original <- hotspin_df$pain_rating
hotspin_df$pain_rating <- (hotspin_df$pain_rating - 50) / 100 * 100

# set pain rescaled below 0 to 0
# Set negative values to 0
hotspin_df$pain_rating <- pmax(hotspin_df$pain_rating, 0)



# -------------------------------
# Combine dataframes
# -------------------------------

# Combine data frames
df_combined <- rbind(peep_df, hotspin_df[, names(peep_df)])
df_combined$group <-as.factor(df_combined$group)
df_combined$gender <-as.factor(df_combined$gender)

# add training volume
df_combined <- df_combined %>%
  left_join(
    training_hours_combined %>% select(subject, hours_week),
    by = "subject"
  )

df_combined_heat <- df_combined[df_combined$modality == 2,]
df_combined_pressure <- df_combined[df_combined$modality == 1,]


#-----------------------------------------------------
#   Online Ratings across exercise intensities
#----------------------------------------------------

# Load in data
sub_online_ratings <- read.csv('C:/Users/user/Desktop/projects/HotSpin/Data/LogExperiment/MAIN/sub_online_heat_and_pressure_clean.csv',sep = ',',header = T)
#sub_online_ratings$rating <- as.numeric(sub_online_ratings$rating)


# Calculate summary statistic for better visualisation
summary_sub_online_ratings<- sub_online_ratings %>%
  group_by(time,modality,VAS_intensity)%>%
  summarise_at(c('online_rating'),mean,na.rm = T)


# Calculating SE for within design 
# (https://www.niklasjohannes.com/post/calculating-and-visualizing-error-bars-for-within-subjects-designs/)
sum_se <- summarySEwithin(sub_online_ratings, 
                          measurevar = "online_rating", 
                          withinvars  = c("time","modality","VAS_intensity"),
                          idvar = 'subject',
                          na.rm = T)

summary_sub_online_ratings$se <- sum_se$se
summary_sub_online_ratings$VAS_intensity <- as.factor(summary_sub_online_ratings$VAS_intensity)
summary_sub_online_ratings$modality <- as.factor(summary_sub_online_ratings$modality)
summary_sub_online_ratings$time <- as.numeric(summary_sub_online_ratings$time)


# Split summary data frames for heat and pressure
summary_sub_online_ratings_heat <- summary_sub_online_ratings %>% filter(modality == 2)
summary_sub_online_ratings_pressure <- summary_sub_online_ratings %>% filter(modality == 1)


# GGplot visualise online ratings
ratings_pressure <- ggplot(summary_sub_online_ratings_pressure,aes(time,online_rating,fill = VAS_intensity,colour = VAS_intensity))+
  #  geom_rect(data = data.frame(xmin = 0, xmax = 17.5, ymin = 0, ymax = 102 + 2),
  #          aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
  #          fill = "#cccccc", alpha = 0.5, inherit.aes = FALSE) +
  geom_ribbon(aes(ymin = online_rating-se,ymax = online_rating+se),alpha = 0.2,show.legend = F,colour = NA,legend = F)+
  #geom_point(size = 1,show.legend = F)+
  geom_line(size = 0.6,show.legend = F)+
  scale_color_manual(labels = c("30", "50","70"),
                     values = c("#969BF2", '#1F248C',"#020659"))+
  scale_fill_manual(labels = c("30", "50","70"),
                    values = c("#969BF2", '#1F248C',"#020659"))+
  theme_classic()+
  theme(axis.title = element_text(size = axis_title_size,family="Helvetica"),axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),plot.title = element_text(size = plot_title_size,family="Helvetica"),legend.title = element_text(size = axis_title_size,family="Helvetica"),
        legend.text = element_text(size = legend_title_size,family="Helvetica"),strip.text.x = element_text(size = legend_text_size,family="Helvetica"),
        plot.margin = margin(t = 0,  # Top margin
                             r = 0,  # Right margin
                             b = 0,  # Bottom margin
                             l = 0,  # Left margin
                             unit = "cm")) +
  # guides(fill=guide_legend(title="")) +
  scale_y_continuous(breaks = seq(0,150,10))+
  scale_x_continuous(breaks=seq(0,20,5))+
  xlab('Time [sec]')+
  ylab('Pain Rating [VAS]')+
  #geom_segment(data = summary_sub_online_ratings, aes(x = 0, xend =  17.5, y = 102 + 2, yend = 102 + 2), na.rm = TRUE, colour = "black", linewidth = 1) +
  # Adding tick marks at segment ends only in the correct facet
  #geom_segment(data = summary_sub_online_ratings, aes(x = 0, xend = 0, y = 99.5 + 2 , yend = 104.5 + 2 ), na.rm = TRUE, colour = "black", linewidth = 1) +
  #geom_segment(data = summary_sub_online_ratings, aes(x = 17.5, xend = 17.5, y = 99.5 + 2 , yend = 104.5 + 2 ), na.rm = TRUE, colour = "black", linewidth = 1) +
  #geom_text(data = summary_sub_online_ratings, aes(x = 10, y = 110, label = "Stimulus Duration"), size = 1, colour = "black")+ 
  geom_hline(yintercept = 0)+
  ylim(0,150)+
  ggtitle("")

 ratings_pressure 


ratings_heat <- ggplot(summary_sub_online_ratings_heat,aes(time,online_rating,fill = VAS_intensity,colour = VAS_intensity))+
  #  geom_rect(data = data.frame(xmin = 0, xmax = 17.5, ymin = 0, ymax = 102 + 2),
  #          aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
  #          fill = "#cccccc", alpha = 0.5, inherit.aes = FALSE) +
  geom_ribbon(aes(ymin = online_rating-se,ymax = online_rating+se),alpha = 0.2,show.legend = F,colour = NA,legend = F)+
  #geom_point(size = 1,show.legend = F)+
  geom_line(size = 0.6,show.legend = F)+
  scale_color_manual(labels = c("30", "50","70"),
                     values = c("#ecd761", '#fa9c0e', "#af0404"))+
  scale_fill_manual(labels = c("30", "50","70"),
                    values = c("#ecd761", '#fa9c0e', "#af0404"))+
  theme_classic()+
  theme(axis.title = element_text(size = axis_title_size,family="Helvetica"),axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),plot.title = element_text(size = plot_title_size,family="Helvetica"),legend.title = element_text(size = axis_title_size,family="Helvetica"),
        legend.text = element_text(size = legend_title_size,family="Helvetica"),strip.text.x = element_text(size = legend_text_size,family="Helvetica"),
        plot.margin = margin(t = 0,  # Top margin
                             r = 0,  # Right margin
                             b = 0,  # Bottom margin
                             l = 0,  # Left margin
                             unit = "cm")) +
  # guides(fill=guide_legend(title="")) +
  scale_y_continuous(breaks = seq(0,150,10))+
  scale_x_continuous(breaks=seq(0,20,5))+
  xlab('Time [sec]')+
  ylab('Pain Rating [VAS]')+
  #geom_segment(data = summary_sub_online_ratings, aes(x = 0, xend =  17.5, y = 102 + 2, yend = 102 + 2), na.rm = TRUE, colour = "black", linewidth = 1) +
  # Adding tick marks at segment ends only in the correct facet
  #geom_segment(data = summary_sub_online_ratings, aes(x = 0, xend = 0, y = 99.5 + 2 , yend = 104.5 + 2 ), na.rm = TRUE, colour = "black", linewidth = 1) +
  #geom_segment(data = summary_sub_online_ratings, aes(x = 17.5, xend = 17.5, y = 99.5 + 2 , yend = 104.5 + 2 ), na.rm = TRUE, colour = "black", linewidth = 1) +
  #geom_text(data = summary_sub_online_ratings, aes(x = 10, y = 110, label = "Stimulus Duration"), size = 1, colour = "black")+ 
  geom_hline(yintercept = 0)+
  ylim(0,150)+
  ggtitle("")

 ratings_heat 

# --------------------------------------------------------
#  Online Ratings and Exercise Intenstiy (across VAS)
#---------------------------------------------------------
 
 #-------- Statistics ------

 #load in data
 sub_online_ratings <- read.csv('C:/Users/user/Desktop/projects/HotSpin/Data/LogExperiment/MAIN/sub_online_heat_and_pressure_clean.csv',sep = ',',header = T)
 
 # Calculate summary statistic for better visualisation
 summary_sub_online_ratings<- sub_online_ratings %>%
   group_by(subject,time,exercise_intensity,modality)%>%
   summarise_at(c('online_rating'),mean,na.rm = T)
 
summary_sub_online_ratings$modality <- as.factor(summary_sub_online_ratings$modality)
 summary_sub_online_ratings$exercise_intensity <- as.factor(summary_sub_online_ratings$exercise_intensity)

 
 # Filter dataframes for heat and pressure
 summary_sub_online_ratings_heat <- summary_sub_online_ratings %>% filter(modality == 2)
 summary_sub_online_ratings_pressure <- summary_sub_online_ratings %>% filter(modality == 1)

 # Heat
 lme_model_heat_online_rating <- lmer(online_rating ~  exercise_intensity*time + (1|subject), data = summary_sub_online_ratings_heat)
 summary(lme_model_heat_online_rating)
 
 # Pressure
 lme_model_pressure_online_rating <- lmer(online_rating ~  exercise_intensity*time + (1|subject), data = summary_sub_online_ratings_pressure)
 summary(lme_model_pressure_online_rating)
 
 
 #-------- Visualisation ------

 
#load ad in data
sub_online_ratings <- read.csv('C:/Users/user/Desktop/projects/HotSpin/Data/LogExperiment/MAIN/sub_online_heat_and_pressure_clean.csv',sep = ',',header = T)

# Calculate summary statistic for better visualisation
summary_sub_online_ratings<- sub_online_ratings %>%
  group_by(time,exercise_intensity,modality)%>%
  summarise_at(c('online_rating'),mean,na.rm = T)


# Calculating SE for within design 
# (https://www.niklasjohannes.com/post/calculating-and-visualizing-error-bars-for-within-subjects-designs/)
sum_se <- summarySEwithin(sub_online_ratings, 
                          measurevar = "online_rating", 
                          withinvars  = c("time","exercise_intensity","modality"),
                          idvar = 'subject',
                          na.rm = T)

summary_sub_online_ratings$se <- sum_se$se
summary_sub_online_ratings$modality <- as.factor(summary_sub_online_ratings$modality)
summary_sub_online_ratings$exercise_intensity <- as.factor(summary_sub_online_ratings$exercise_intensity)
summary_sub_online_ratings$time <- as.numeric(summary_sub_online_ratings$time)

# Filter dataframes for heat and pressure
summary_sub_online_ratings_heat <- summary_sub_online_ratings %>% filter(modality == 2)
summary_sub_online_ratings_pressure <- summary_sub_online_ratings %>% filter(modality == 1)


heat_online_rating_across_VAS <- ggplot(summary_sub_online_ratings_heat, aes(time, online_rating, fill = exercise_intensity, colour = exercise_intensity)) +
  # Adding rect for heat pain plateau (modality '2')
  geom_rect(data = data.frame(xmin = 0.95, xmax = 16.1, ymin = 0, ymax = 90, modality = '2'),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "#7d7d7e", alpha = 0.5, inherit.aes = FALSE) +
  
  # Adding rect for stimulus duration
  geom_rect(data = data.frame(xmin = 0, xmax = 17.1, ymin = 0, ymax = 102 + 2),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "#cccccc", alpha = 0.5, inherit.aes = FALSE) +


  geom_ribbon(aes(ymin = online_rating - se, ymax = online_rating + se), alpha = 0.2, show.legend = FALSE, colour = NA) +
  geom_line(size = 1, show.legend = TRUE) +
  
  scale_color_manual(labels = c("Low", "High"),
                     values = c("#005C53", "#3C008E")) +
  scale_fill_manual(labels = c("Low", "High"),
                    values = c("#005C53", "#3C008E")) +
  theme_classic() +
  
  theme(axis.title = element_text(size = axis_title_size,family="Helvetica"),axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),plot.title = element_text(size = plot_title_size,family="Helvetica"),legend.title = element_text(size = axis_title_size,family="Helvetica"),
        legend.text = element_text(size = legend_title_size,family="Helvetica"),strip.text.x = element_text(size = legend_text_size,family="Helvetica"),
        plot.margin = margin(t = 0,  # Top margin
                             r = 0,  # Right margin
                             b = 0,  # Bottom margin
                             l = 0,  # Left margin
                             unit = "cm")) +
  
  ylim(0, 150) +
  xlab('Time [sec]') +
  ylab('Pain Rating [VAS]') +
  

   geom_segment(data = summary_sub_online_ratings, aes(x = ifelse(modality == '2', 0.95, NA), 
                                                      xend = ifelse(modality == '2', 16.1, NA), 
                                                      y = 90, yend = 90), na.rm = TRUE, colour = "black", linewidth = 1) +
  geom_segment(data = summary_sub_online_ratings, aes(x = 0, xend =  17.1, y = 102 + 2, yend = 102 + 2), na.rm = TRUE, colour = "black", linewidth = 1) +
  
  # Adding tick marks at segment ends only in the correct facet
    geom_segment(data = subset(summary_sub_online_ratings, modality == '2'), aes(x = 0.95, xend = 0.95, y = 88, yend = 92), na.rm = TRUE, colour = "black", linewidth = 1) +
  geom_segment(data = subset(summary_sub_online_ratings, modality == '2'), aes(x = 16, xend = 16, y = 88, yend = 92), na.rm = TRUE, colour = "black", linewidth = 1) +
  geom_segment(data = summary_sub_online_ratings, aes(x = 0, xend = 0, y = 99.5 + 2 , yend = 104.5 + 2 ), na.rm = TRUE, colour = "black", linewidth = 1) +
  geom_segment(data = summary_sub_online_ratings, aes(x = 17.1, xend = 17.1, y = 99.5 + 2 , yend = 104.5 + 2 ), na.rm = TRUE, colour = "black", linewidth = 1) +
  
  #geom_text(aes(x = 10, y = 104 + 3, label = "Stimulus Duration"), size = 1.25, colour = "black") +
  #geom_text(data = subset(summary_sub_online_ratings, modality == '1'), aes(x = 10, y = 93, label = "Plateau Duration"), size = 1, colour = "black") +

  geom_hline(yintercept = 0)+
  ggtitle('Heat')

heat_online_rating_across_VAS

# ---- pressure 
pressure_online_rating_across_VAS <- ggplot(summary_sub_online_ratings_pressure, aes(time, online_rating, fill = exercise_intensity, colour = exercise_intensity)) +
  # Adding rect for heat pain plateau (modality '2')
  geom_rect(data = data.frame(xmin = 1.5, xmax = 16.5, ymin = 0, ymax = 90, modality = '2'),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "#7d7d7e", alpha = 0.5, inherit.aes = FALSE) +
  
  # Adding rect for stimulus duration
  geom_rect(data = data.frame(xmin = 0, xmax = 17.1, ymin = 0, ymax = 102 + 2),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "#cccccc", alpha = 0.5, inherit.aes = FALSE) +


  geom_ribbon(aes(ymin = online_rating - se, ymax = online_rating + se), alpha = 0.2, show.legend = FALSE, colour = NA) +
  geom_line(size = 1, show.legend = TRUE) +
  
  scale_color_manual(labels = c("Low", "High"),
                     values = c("#005C53", "#3C008E")) +
  scale_fill_manual(labels = c("Low", "High"),
                    values = c("#005C53", "#3C008E")) +
  theme_classic() +
  
  theme(axis.title = element_text(size = axis_title_size,family="Helvetica"),axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),plot.title = element_text(size = plot_title_size,family="Helvetica"),legend.title = element_text(size = axis_title_size,family="Helvetica"),
        legend.text = element_text(size = legend_title_size,family="Helvetica"),strip.text.x = element_text(size = legend_text_size,family="Helvetica"),
        plot.margin = margin(t = 0,  # Top margin
                             r = 0,  # Right margin
                             b = 0,  # Bottom margin
                             l = 0,  # Left margin
                             unit = "cm")) +
  
  ylim(0, 150) +
  xlab('Time [sec]') +
  ylab('Pain Rating [VAS]') +
  

   geom_segment(data = summary_sub_online_ratings, aes(x = ifelse(modality == '1', 1.5, NA), 
                                                      xend = ifelse(modality == '1', 16.5, NA), 
                                                      y = 90, yend = 90), na.rm = TRUE, colour = "black", linewidth = 1) +
  geom_segment(data = summary_sub_online_ratings, aes(x = 0, xend =  17.1, y = 102 + 2, yend = 102 + 2), na.rm = TRUE, colour = "black", linewidth = 1) +
  
  # Adding tick marks at segment ends only in the correct facet
    geom_segment(data = subset(summary_sub_online_ratings, modality == '2'), aes(x = 1.5, xend = 1.5, y = 88, yend = 92), na.rm = TRUE, colour = "black", linewidth = 1) +
  geom_segment(data = subset(summary_sub_online_ratings, modality == '2'), aes(x = 16.5, xend = 16.5, y = 88, yend = 92), na.rm = TRUE, colour = "black", linewidth = 1) +
  geom_segment(data = summary_sub_online_ratings, aes(x = 0, xend = 0, y = 99.5 + 2 , yend = 104.5 + 2 ), na.rm = TRUE, colour = "black", linewidth = 1) +
  geom_segment(data = summary_sub_online_ratings, aes(x = 17.1, xend = 17.1, y = 99.5 + 2 , yend = 104.5 + 2 ), na.rm = TRUE, colour = "black", linewidth = 1) +
  
  #geom_text(aes(x = 10, y = 104 + 3, label = "Stimulus Duration"), size = 1.25, colour = "black") +
  #geom_text(data = subset(summary_sub_online_ratings, modality == '1'), aes(x = 10, y = 93, label = "Plateau Duration"), size = 1, colour = "black") +
  
  geom_hline(yintercept = 0)+
  ggtitle('Pressure')

pressure_online_rating_across_VAS


#----------------------------------------------------
# Online Ratings and Exercise Intenstiy
#-----------------------------------------------------


#-------- Statistics ------

#load in data
sub_online_ratings <- read.csv('C:/Users/user/Desktop/projects/HotSpin/Data/LogExperiment/MAIN/sub_online_heat_and_pressure_clean.csv',sep = ',',header = T)

# Calculate summary statistic for better visualisation
summary_sub_online_ratings<- sub_online_ratings %>%
  group_by(subject,time,exercise_intensity,modality,VAS_intensity)%>%
  summarise_at(c('online_rating'),mean,na.rm = T)

summary_sub_online_ratings$modality <- as.factor(summary_sub_online_ratings$modality)
summary_sub_online_ratings$exercise_intensity <- as.factor(summary_sub_online_ratings$exercise_intensity)


# Filter dataframes for heat and pressure
summary_sub_online_ratings_heat <- summary_sub_online_ratings %>% filter(modality == 2)
summary_sub_online_ratings_pressure <- summary_sub_online_ratings %>% filter(modality == 1)

# Heat
lme_model_heat_online_rating <- lmer(online_rating ~  exercise_intensity*time*VAS_intensity + (1|subject), data = summary_sub_online_ratings_heat)
summary(lme_model_heat_online_rating)

# Pressure
lme_model_pressure_online_rating <- lmer(online_rating ~  exercise_intensity*time*VAS_intensity + (1|subject), data = summary_sub_online_ratings_pressure)
summary(lme_model_pressure_online_rating)



# ------ Visualisation --------

# Load in data
sub_online_ratings <- read.csv('C:/Users/user/Desktop/projects/HotSpin/Data/LogExperiment/MAIN/sub_online_heat_and_pressure_clean.csv',sep = ',',header = T)

# Divide dataframes for both modalities 
sub_online_ratings_heat <- sub_online_ratings[sub_online_ratings$modality == 2,]
sub_online_ratings_pressure <- sub_online_ratings[sub_online_ratings$modality == 1,]


# Calculate summary statistic for better visualisation
summary_sub_online_ratings<- sub_online_ratings %>%
  group_by(time,exercise_intensity,modality,VAS_intensity)%>%
  summarise_at(c('online_rating'),mean,na.rm = T)


# Calculating SE for within design 
# (https://www.niklasjohannes.com/post/calculating-and-visualizing-error-bars-for-within-subjects-designs/)
sum_se <- summarySEwithin(sub_online_ratings, 
                          measurevar = "online_rating", 
                          withinvars  = c("time","exercise_intensity","modality","VAS_intensity"),
                          idvar = 'subject',
                          na.rm = T)

summary_sub_online_ratings$se <- sum_se$se
summary_sub_online_ratings$VAS_intensity <- as.factor(summary_sub_online_ratings$VAS_intensity)
summary_sub_online_ratings$modality <- as.factor(summary_sub_online_ratings$modality)
summary_sub_online_ratings$exercise_intensity <- as.factor(summary_sub_online_ratings$exercise_intensity)
summary_sub_online_ratings$time <- as.numeric(summary_sub_online_ratings$time)

(online_rating_VAS <-ggplot(summary_sub_online_ratings,aes(time,online_rating,fill = exercise_intensity,colour = exercise_intensity))+
   
    geom_ribbon(aes(ymin = online_rating-se,ymax = online_rating+se),alpha = 0.2,show.legend = F,colour = NA)+
    #geom_point(size = 1,show.legend = F)+
    geom_line(size = 0.6,show.legend = T)+
    
    scale_color_manual(labels = c("Low", "High"),
                       values = c("#005C53", "#3C008E"))+
    scale_fill_manual(labels = c("Low", "High"),
                      values = c("#005C53", "#3C008E"))+
    theme_classic()+
    
    theme(axis.title = element_text(size = axis_title_size,family="Helvetica"),axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),plot.title = element_text(size = plot_title_size,family="Helvetica"),legend.title = element_text(size = axis_title_size,family="Helvetica"),
          legend.text = element_text(size = legend_title_size,family="Helvetica"),        strip.text = element_blank(),
          strip.background = element_blank(),
          plot.margin = margin(t = 0,  # Top margin
                               r = 0,  # Right margin
                               b = 0,  # Bottom margin
                               l = 0,  # Left margin
                               unit = "cm")) +
    # guides(fill=guide_legend(title="")) +
    #scale_y_continuous(breaks = seq(0,150,10))+
    #scale_x_continuous(breaks=seq(0,20,5))+
    
    ylim(0, 150) +
    xlab('Time [sec]')+
    ylab('Pain Rating [VAS]')+
    facet_wrap(modality~VAS_intensity,labeller = as_labeller(c('1'='Pressure','2'= 'Heat','30'='30','50' = '50','70'='70')),ncol = 3,nrow = 2)+
    #geom_hline(yintercept = 50,size = 1)+
    
    geom_hline(yintercept = 0)+
    
    
    ggtitle("")
)

#-----------------------------------------------------------------------------------------------------------------
# Main effect exercise (peak pain ratings from subject ratings across VAS across modality)
#-----------------------------------------------------------------------------------------------------------------

# Load in data
sub_online_ratings <- read.csv('C:/Users/user/Desktop/projects/HotSpin/Data/LogExperiment/MAIN/sub_online_heat_and_pressure_clean.csv',sep = ',',header = T)

# calculate the mean time courses for each subject
mean_ratings_sub <- sub_online_ratings %>%
  group_by(subject,time,exercise_intensity)%>%
  summarise_at(c('online_rating'),mean,na.rm = T)

# Calculate the max rating for each exercise intensity, VAS intensity, modality, and subject
sub_max_online_ratings <- mean_ratings_sub %>%
  group_by(subject, exercise_intensity) %>%
  summarise_at(c('online_rating'),max,na.rm = T)

# Replace Inf values with NA
sub_max_online_ratings$online_rating[is.infinite(sub_max_online_ratings$online_rating)] <- NA

# calculat ethe mean of the max ratings across subjects
max_online_ratings <- sub_max_online_ratings %>%
  group_by(exercise_intensity) %>%
  summarise_at(c('online_rating'),mean,na.rm = T)


# Calculating SE for within design 
# (https://www.niklasjohannes.com/post/calculating-and-visualizing-error-bars-for-within-subjects-designs/)
sum_se <- summarySEwithin(sub_max_online_ratings, 
                          measurevar = "online_rating", 
                          withinvars  = c("exercise_intensity"),
                          idvar = 'subject',
                          na.rm = T)

max_online_ratings$se <- sum_se$se


max_online_ratings$exercise_intensity <- as.factor(max_online_ratings$exercise_intensity)
sub_max_online_ratings$exercise_intensity <- as.factor(sub_max_online_ratings$exercise_intensity)

# Plot as barplots with standard error
max_ratings_whole_across_VAS_across_modality <- ggplot(max_online_ratings, aes(x = exercise_intensity, y = online_rating, fill = exercise_intensity,colour = exercise_intensity)) +
  geom_jitter(data = sub_max_online_ratings,aes(x = exercise_intensity, y = online_rating, fill = exercise_intensity,color = exercise_intensity), shape = 21,alpha = 0.4,size = 0.5,position = position_jitterdodge(jitter.width = 0.1,dodge.width = 0.9),show.legend = F) +
  geom_bar(stat = "identity", position = position_dodge(), alpha = 0.7) +
  geom_errorbar(aes(ymin = online_rating - se, ymax = online_rating + se), 
                colour = "black",width  =0.2,size = 0.5,position = position_dodge(0.9))+
  scale_fill_manual(labels = c("Low", "High"),
                    values = c("#005c23", "#3C008E")) +
  scale_color_manual(labels = c("Low", "High"),
                     values = c("#005c23", "#3C008E")) +
  theme_classic() +
  theme(axis.title = element_text(size = axis_title_size, family = "Helvetica"),
        axis.text = element_text(size = axis_text_size, colour = 'black', family = "Helvetica"),
        plot.title = element_text(size = plot_title_size, family = "Helvetica"),
        legend.title = element_text(size = axis_title_size, family = "Helvetica"),
        legend.text = element_text(size = legend_text_size, family = "Helvetica"),
        strip.text = element_blank(),
        strip.background = element_blank(),
        legend.key.size = unit(0.25, 'cm')) +
  xlab('Stimulus Intensity [VAS]') +
  ylab('Max Pain Rating [VAS]') +
  ggtitle('') +
  ylim(0, 150)

max_ratings_whole_across_VAS_across_modality


#---------------- Statistics---------------------


# calculate the mean time courses for each subject
mean_ratings_sub <- sub_online_ratings %>%
  group_by(subject,time,exercise_intensity,trial,block)%>%
  summarise_at(c('online_rating'),mean,na.rm = T)


# Calculate the max rating for each exercise intensity, VAS intensity, modality, and subject
sub_max_online_ratings <- mean_ratings_sub %>%
  group_by(subject, exercise_intensity,trial,block) %>%
  summarise_at(c('online_rating'),max,na.rm = T)

# Replace Inf values with NA
sub_max_online_ratings$online_rating[is.infinite(sub_max_online_ratings$online_rating)] <- NA

sub_max_online_ratings$exercise_intensity <- as.factor(sub_max_online_ratings$exercise_intensity)

#  Exercise intensity  mean pain ratings
lme_model_main <- lmer(online_rating ~   exercise_intensity + (1|subject) + trial + block, data = sub_max_online_ratings)
summary(lme_model_main)
confint(lme_model_main)


# Post hoc t-tests for heat max ratings at all stimulus intensities based on the lmer model
# Perform pairwise comparisons for each VAS intensity

lme_model_main <- lmer(online_rating ~   exercise_intensity + (1|subject) + trial + block, data = sub_max_online_ratings)

emm_results_h <- emmeans(lme_model_main, ~ exercise_intensity)

# View the results for heat model
summary(emm_results_h)

# Perform pairwise comparisons for the heat model across all VAS intensity levels
pairwise_comparisons_h <- contrast(emm_results_h, interaction = "pairwise", adjust = "tukey")

# View pairwise results for heat model
summary(pairwise_comparisons_h)



#-----------------------------------------------------------------------------------------------------------------
# Main effect exercise (peak pain ratings from subject ratings across VAS for heat and pressure seperateyl)
#-----------------------------------------------------------------------------------------------------------------

# Load in data
sub_online_ratings <- read.csv('C:/Users/user/Desktop/projects/HotSpin/Data/LogExperiment/MAIN/sub_online_heat_and_pressure_clean.csv',sep = ',',header = T)

# calculate the mean time courses for each subject
mean_ratings_sub <- sub_online_ratings %>%
  group_by(subject,time,exercise_intensity,modality)%>%
  summarise_at(c('online_rating'),mean,na.rm = T)

# Calculate the max rating for each exercise intensity, VAS intensity, modality, and subject
sub_max_online_ratings <- mean_ratings_sub %>%
  group_by(subject, exercise_intensity, modality) %>%
  summarise_at(c('online_rating'),max,na.rm = T)

# Replace Inf values with NA
sub_max_online_ratings$online_rating[is.infinite(sub_max_online_ratings$online_rating)] <- NA

# calculat ethe mean of the max ratings across subjects
max_online_ratings <- sub_max_online_ratings %>%
  group_by(exercise_intensity, modality) %>%
  summarise_at(c('online_rating'),mean,na.rm = T)


# Calculating SE for within design 
# (https://www.niklasjohannes.com/post/calculating-and-visualizing-error-bars-for-within-subjects-designs/)
sum_se <- summarySEwithin(sub_max_online_ratings, 
                          measurevar = "online_rating", 
                          withinvars  = c("exercise_intensity","modality"),
                          idvar = 'subject',
                          na.rm = T)

max_online_ratings$se <- sum_se$se

max_online_ratings$modality <- as.factor(max_online_ratings$modality)
max_online_ratings$exercise_intensity <- as.factor(max_online_ratings$exercise_intensity)
sub_max_online_ratings$modality <- as.factor(sub_max_online_ratings$modality)
sub_max_online_ratings$exercise_intensity <- as.factor(sub_max_online_ratings$exercise_intensity)

# Convert to factor with labels
max_online_ratings$exercise_intensity <- factor(max_online_ratings$exercise_intensity, 
                                                levels = c(0, 1), 
                                                labels = c("Low", "High"))

sub_max_online_ratings$exercise_intensity <- factor(sub_max_online_ratings$exercise_intensity, 
                                                    levels = c(0, 1), 
                                                    labels = c("Low", "High"))

# Plot as barplots with standard error
max_ratings_whole_across_VAS <- ggplot(max_online_ratings, aes(x = exercise_intensity, y = online_rating, fill = exercise_intensity,colour = exercise_intensity)) +
  geom_jitter(data = sub_max_online_ratings,aes(x = exercise_intensity, y = online_rating, fill = exercise_intensity,color = exercise_intensity), shape = 21,alpha = 0.4,size = 0.5,position = position_jitterdodge(jitter.width = 0.1,dodge.width = 0.9),show.legend = F) +
  geom_bar(stat = "identity", position = position_dodge(), alpha = 0.7) +
  geom_errorbar(aes(ymin = online_rating - se, ymax = online_rating + se), 
                colour = "black",width  =0.2,size = 0.5,position = position_dodge(0.9))+
  scale_fill_manual(labels = c("Low", "High"),
                    values = c("#005c23", "#3C008E")) +
  scale_color_manual(labels = c("Low", "High"),
                     values = c("#005c23", "#3C008E")) +
  theme_classic() +
  theme(axis.title = element_text(size = axis_title_size, family = "Helvetica"),
        axis.text = element_text(size = axis_text_size, colour = 'black', family = "Helvetica"),
        plot.title = element_text(size = plot_title_size, family = "Helvetica"),
        legend.title = element_text(size = axis_title_size, family = "Helvetica"),
        legend.text = element_text(size = legend_text_size, family = "Helvetica"),
        strip.text = element_blank(),
        strip.background = element_blank(),
        legend.key.size = unit(0.25, 'cm')) +
  xlab('Exercise Intensity') +
  ylab('Max Pain Rating [VAS]') +
  ggtitle('') +
  ylim(0, 150)+
  facet_wrap(~modality, labeller = as_labeller(c('1' = 'Pressure', '2' = 'Heat')),ncol = 1)

max_ratings_whole_across_VAS



#---------------- Statistics---------------------

# Filter sub_online_ratings_plateau into new dataframes for heat and pressure pain separately
sub_online_ratings_max_heat <- sub_online_ratings[sub_online_ratings$modality == 2,]
sub_online_ratings_max_pressure <- sub_online_ratings[sub_online_ratings$modality == 1,]

# calculate the mean time courses for each subject
mean_ratings_sub_heat <- sub_online_ratings_max_heat %>%
  group_by(subject,time,exercise_intensity,trial,block)%>%
  summarise_at(c('online_rating'),mean,na.rm = T)

mean_ratings_sub_pressure <- sub_online_ratings_max_pressure %>%
  group_by(subject,time,exercise_intensity,trial,block)%>%
  summarise_at(c('online_rating'),mean,na.rm = T)

# Calculate the max rating for each exercise intensity, VAS intensity, modality, and subject
sub_max_online_ratings_heat <- mean_ratings_sub_heat %>%
  group_by(subject, exercise_intensity,trial,block) %>%
  summarise_at(c('online_rating'),max,na.rm = T)

sub_max_online_ratings_pressure <- mean_ratings_sub_pressure %>%
  group_by(subject, exercise_intensity,trial,block) %>%
  summarise_at(c('online_rating'),max,na.rm = T)

# Replace Inf values with NA
sub_max_online_ratings_heat$online_rating[is.infinite(sub_max_online_ratings_heat$online_rating)] <- NA
sub_max_online_ratings_pressure$online_rating[is.infinite(sub_max_online_ratings_pressure$online_rating)] <- NA

# Heat Exercise intensity and VAS intenstiy on mean pain ratings
lme_model_main_h <- lmer(online_rating ~   exercise_intensity + (1|subject)+ trial + block, data = sub_max_online_ratings_heat)
summary(lme_model_main_h)
confint(lme_model_main_h)

# Pressure Exercise intensity and VAS intenstiy on mean pain ratings
lme_model_main_p <- lmer(online_rating ~   exercise_intensity + (1|subject)+ trial + block, data = sub_max_online_ratings_pressure)
summary(lme_model_main_p)
confint(lme_model_main_p)


# Post hoc t-tests for heat max ratings at all stimulus intensities based on the lmer model
# Perform pairwise comparisons for each VAS intensity
sub_max_online_ratings_heat$exercise_intensity <- as.factor(sub_max_online_ratings_heat$exercise_intensity)
lme_model_main_h <- lmer(online_rating ~   exercise_intensity + (1|subject)+ trial + block, data = sub_max_online_ratings_heat)

emm_results_h <- emmeans(lme_model_main_h, ~ exercise_intensity)

# View the results for heat model
summary(emm_results_h)

# Perform pairwise comparisons for the heat model across all VAS intensity levels
pairwise_comparisons_h <- contrast(emm_results_h, interaction = "pairwise", adjust = "tukey")

# View pairwise results for heat model
summary(pairwise_comparisons_h)

#-------------- pressure 
# Post hoc t-tests for heat max ratings at all stimulus intensities based on the lmer model
# Perform pairwise comparisons for each VAS intensity
sub_max_online_ratings_pressure$exercise_intensity <- as.factor(sub_max_online_ratings_pressure$exercise_intensity)

lme_model_main_p <- lmer(online_rating ~   exercise_intensity + (1|subject)+ trial + block, data = sub_max_online_ratings_pressure)

emm_results_p <- emmeans(lme_model_main_p, ~ exercise_intensity)

# View the results for heat model
summary(emm_results_p)

# Perform pairwise comparisons for the heat model across all VAS intensity levels
pairwise_comparisons_p <- contrast(emm_results_p, interaction = "pairwise", adjust = "tukey")

# View pairwise results for heat model
summary(pairwise_comparisons_p)



#-----------------------------------------------------
# Calculate peak pain ratings from subject ratings at each VAS intensity
#----------------------------------------------------

# Load in data
sub_online_ratings <- read.csv('C:/Users/user/Desktop/projects/HotSpin/Data/LogExperiment/MAIN/sub_online_heat_and_pressure_clean.csv',sep = ',',header = T)

# calculate the mean time courses for each subject
mean_ratings_sub <- sub_online_ratings %>%
  group_by(subject,time,exercise_intensity,modality,VAS_intensity)%>%
  summarise_at(c('online_rating'),mean,na.rm = T)

# Calculate the max rating for each exercise intensity, VAS intensity, modality, and subject
sub_max_online_ratings <- mean_ratings_sub %>%
    group_by(subject, exercise_intensity, modality, VAS_intensity) %>%
    summarise_at(c('online_rating'),max,na.rm = T)

# Replace Inf values with NA
sub_max_online_ratings$online_rating[is.infinite(sub_max_online_ratings$online_rating)] <- NA

# calculat ethe mean of the max ratings across subjects
 max_online_ratings <- sub_max_online_ratings %>%
    group_by(exercise_intensity, modality, VAS_intensity) %>%
    summarise_at(c('online_rating'),mean,na.rm = T)


# Calculating SE for within design 
# (https://www.niklasjohannes.com/post/calculating-and-visualizing-error-bars-for-within-subjects-designs/)
sum_se <- summarySEwithin(sub_max_online_ratings, 
                          measurevar = "online_rating", 
                          withinvars  = c("exercise_intensity","modality","VAS_intensity"),
                          idvar = 'subject',
                          na.rm = T)

max_online_ratings$se <- sum_se$se
max_online_ratings$VAS_intensity <- as.factor(max_online_ratings$VAS_intensity)
max_online_ratings$modality <- as.factor(max_online_ratings$modality)
max_online_ratings$exercise_intensity <- as.factor(max_online_ratings$exercise_intensity)
sub_max_online_ratings$VAS_intensity <- as.factor(sub_max_online_ratings$VAS_intensity)
sub_max_online_ratings$modality <- as.factor(sub_max_online_ratings$modality)
sub_max_online_ratings$exercise_intensity <- as.factor(sub_max_online_ratings$exercise_intensity)


# Plot as barplots with standard error
max_ratings_whole <- ggplot(max_online_ratings, aes(x = VAS_intensity, y = online_rating, fill = exercise_intensity,colour = exercise_intensity)) +
  geom_jitter(data = sub_max_online_ratings,aes(x = VAS_intensity, y = online_rating, fill = exercise_intensity,color = exercise_intensity), shape = 21,alpha = 0.4,size = 0.5,position = position_jitterdodge(jitter.width = 0.1,dodge.width = 0.9),show.legend = F) +
  geom_bar(stat = "identity", position = position_dodge(), alpha = 0.7) +
  geom_errorbar(aes(ymin = online_rating - se, ymax = online_rating + se), 
                colour = "black",width  =0.2,size = 0.5,position = position_dodge(0.9))+
  scale_fill_manual(labels = c("Low", "High"),
                    values = c("#005c23", "#3C008E")) +
  scale_color_manual(labels = c("Low", "High"),
                     values = c("#005c23", "#3C008E")) +
  theme_classic() +
  theme(axis.title = element_text(size = axis_title_size, family = "Helvetica"),
        axis.text = element_text(size = axis_text_size, colour = 'black', family = "Helvetica"),
        plot.title = element_text(size = plot_title_size, family = "Helvetica"),
        legend.title = element_text(size = axis_title_size, family = "Helvetica"),
        legend.text = element_text(size = legend_text_size, family = "Helvetica"),
        strip.text = element_blank(),
        strip.background = element_blank(),
        legend.key.size = unit(0.25, 'cm')) +
  xlab('Stimulus Intensity [VAS]') +
  ylab('Max Pain Rating [VAS]') +
  ggtitle('') +
        ylim(0, 150)+
  facet_wrap(~modality, labeller = as_labeller(c('1' = 'Pressure', '2' = 'Heat')),ncol = 1)

  max_ratings_whole



#---------------- Statistics---------------------

# Filter sub_online_ratings_plateau into new dataframes for heat and pressure pain separately
sub_online_ratings_max_heat <- sub_online_ratings[sub_online_ratings$modality == 2,]
sub_online_ratings_max_pressure <- sub_online_ratings[sub_online_ratings$modality == 1,]

# calculate the mean time courses for each subject
mean_ratings_sub_heat <- sub_online_ratings_max_heat %>%
  group_by(subject,time,exercise_intensity,VAS_intensity,trial,block)%>%
  summarise_at(c('online_rating'),mean,na.rm = T)

mean_ratings_sub_pressure <- sub_online_ratings_max_pressure %>%
  group_by(subject,time,exercise_intensity,VAS_intensity,trial,block)%>%
  summarise_at(c('online_rating'),mean,na.rm = T)

# Calculate the max rating for each exercise intensity, VAS intensity, modality, and subject
sub_max_online_ratings_heat <- mean_ratings_sub_heat %>%
    group_by(subject, exercise_intensity, VAS_intensity,trial,block) %>%
    summarise_at(c('online_rating'),max,na.rm = T)

sub_max_online_ratings_pressure <- mean_ratings_sub_pressure %>%
    group_by(subject, exercise_intensity, VAS_intensity,trial,block) %>%
    summarise_at(c('online_rating'),max,na.rm = T)

# Replace Inf values with NA
sub_max_online_ratings_heat$online_rating[is.infinite(sub_max_online_ratings_heat$online_rating)] <- NA
sub_max_online_ratings_pressure$online_rating[is.infinite(sub_max_online_ratings_pressure$online_rating)] <- NA

# Heat Exercise intensity and VAS intenstiy on mean pain ratings
lme_model_main_h <- lmer(online_rating ~   exercise_intensity * VAS_intensity + (1|subject)+ trial + block, data = sub_max_online_ratings_heat)
summary(lme_model_main_h)
confint(lme_model_main_h)

# Pressure Exercise intensity and VAS intenstiy on mean pain ratings
lme_model_main_p <- lmer(online_rating ~   exercise_intensity * VAS_intensity + (1|subject)+ trial + block, data = sub_max_online_ratings_pressure)
summary(lme_model_main_p)
confint(lme_model_main_p)

# Post hoc t-tests for heat max ratings at all stimulus intensities based on the lmer model
# Perform pairwise comparisons for each VAS intensity
sub_max_online_ratings_heat$exercise_intensity <- as.factor(sub_max_online_ratings_heat$exercise_intensity)
sub_max_online_ratings_heat$VAS_intensity <- as.factor(sub_max_online_ratings_heat$VAS_intensity)
lme_model_main_h <- lmer(online_rating ~   exercise_intensity * VAS_intensity + (1|subject)+ trial + block, data = sub_max_online_ratings_heat)

emm_results_h <- emmeans(lme_model_main_h, ~ exercise_intensity | VAS_intensity)

# View the results for heat model
summary(emm_results_h)

# Perform pairwise comparisons for the heat model across all VAS intensity levels
pairwise_comparisons_h <- contrast(emm_results_h, interaction = "pairwise", adjust = "tukey")

# View pairwise results for heat model
summary(pairwise_comparisons_h)

#-------------- pressure 
# Post hoc t-tests for heat max ratings at all stimulus intensities based on the lmer model
# Perform pairwise comparisons for each VAS intensity
sub_max_online_ratings_pressure$exercise_intensity <- as.factor(sub_max_online_ratings_pressure$exercise_intensity)
sub_max_online_ratings_pressure$VAS_intensity <- as.factor(sub_max_online_ratings_pressure$VAS_intensity)
lme_model_main_p <- lmer(online_rating ~   exercise_intensity * VAS_intensity + (1|subject)+ trial + block, data = sub_max_online_ratings_pressure)

emm_results_p <- emmeans(lme_model_main_p, ~ exercise_intensity | VAS_intensity)

# View the results for heat model
summary(emm_results_p)

# Perform pairwise comparisons for the heat model across all VAS intensity levels
pairwise_comparisons_p <- contrast(emm_results_p, interaction = "pairwise", adjust = "tukey")

# View pairwise results for heat model
summary(pairwise_comparisons_p)



#---------------- Statistics across modality---------------------

# calculate the mean time courses for each subject
mean_ratings_sub <- sub_online_ratings %>%
  group_by(subject,time,exercise_intensity,VAS_intensity,trial,block,modality)%>%
  summarise_at(c('online_rating'),mean,na.rm = T)


# Calculate the max rating for each exercise intensity, VAS intensity, modality, and subject
sub_max_online_ratings<- mean_ratings_sub %>%
  group_by(subject, exercise_intensity, VAS_intensity,trial,block,modality) %>%
  summarise_at(c('online_rating'),max,na.rm = T)


# Replace Inf values with NA
sub_max_online_ratings$online_rating[is.infinite(sub_max_online_ratings$online_rating)] <- NA


# Heat Exercise intensity and VAS intenstiy on mean pain ratings
lme_model_main_h <- lmer(online_rating ~   exercise_intensity * VAS_intensity + (1|subject)+ trial + block, data = sub_max_online_ratings)
summary(lme_model_main_h)
confint(lme_model_main_h)


# ------------------
# Deltas 
# --------------

# Load in data
sub_online_ratings <- read.csv('C:/Users/user/Desktop/projects/HotSpin/Data/LogExperiment/MAIN/sub_online_heat_and_pressure_clean.csv',sep = ',',header = T)

# calculate the mean time courses for each subject
mean_ratings_sub <- sub_online_ratings %>%
  group_by(subject,time,exercise_intensity,modality,VAS_intensity)%>%
  summarise_at(c('online_rating'),mean,na.rm = T)

# Calculate the max rating for each exercise intensity, VAS intensity, modality, and subject
sub_max_online_ratings <- mean_ratings_sub %>%
  group_by(subject, exercise_intensity, modality, VAS_intensity) %>%
  summarise_at(c('online_rating'),max,na.rm = T)

# Replace Inf values with NA
sub_max_online_ratings$online_rating[is.infinite(sub_max_online_ratings$online_rating)] <- NA

# calculat ethe mean of the max ratings across subjects
max_online_ratings <- sub_max_online_ratings %>%
  group_by(exercise_intensity, modality, VAS_intensity) %>%
  summarise_at(c('online_rating'),mean,na.rm = T)

max_ratings_diff <- spread(max_online_ratings, exercise_intensity, online_rating)
max_ratings_diff$diff_ints <- max_ratings_diff$'0' -  max_ratings_diff$'1'  # Low Intensity - High Intensity Data


max_online_ratings <- sub_max_online_ratings %>%
  group_by(subject,exercise_intensity, modality, VAS_intensity) %>%
  summarise_at(c('online_rating'),mean,na.rm = T)

sub_max_online_ratings_diff <- spread(max_online_ratings, exercise_intensity, online_rating)
sub_max_online_ratings_diff$diff_ints <- sub_max_online_ratings_diff$'0' -  sub_max_online_ratings_diff$'1'  # Low Intensity - High Intensity Data

sub_max_online_ratings_diff$VAS_intensity <-
  factor(sub_max_online_ratings_diff$VAS_intensity)

# Pressure
lme_model_main_p <- lmer(diff_ints ~   VAS_intensity + (1|subject), data = sub_max_online_ratings_diff[sub_max_online_ratings_diff$modality==1,])
summary(lme_model_main_p)
emm_results_p <- emmeans(lme_model_main_p, ~ VAS_intensity)

## Heat 
lme_model_main_h <- lmer(diff_ints ~   VAS_intensity + (1|subject), data = sub_max_online_ratings_diff[sub_max_online_ratings_diff$modality==2,])
summary(lme_model_main_h)
emm_results_h <- emmeans(lme_model_main_h, ~ VAS_intensity)

p <- summary(emmeans(lme_model_main_p, ~ VAS_intensity))
h <- summary(emmeans(lme_model_main_h, ~ VAS_intensity))

max_ratings_diff$se <- c(p$SE, h$SE)
max_ratings_diff$VAS_intensity <- as.factor(max_ratings_diff$VAS_intensity)
sub_max_online_ratings_diff$VAS_intensity <- as.factor(sub_max_online_ratings_diff$VAS_intensity)

(deltas_param_max <-ggplot(max_ratings_diff,aes(VAS_intensity,diff_ints))+
    #geom_jitter(data = sub_max_online_ratings_diff,aes(x = VAS_intensity, y = diff_ints), shape = 21,alpha = 0.4,size = 0.5,position = position_jitter(width = 0.1),show.legend = F) +
    geom_bar(stat = 'identity',alpha = 0.8,width = 0.6,colour = 'black',position = position_dodge(0.7))+
    geom_errorbar(aes(VAS_intensity,ymin = diff_ints-se,ymax=diff_ints+se),colour = "black",width  =0.2,size = 0.5)+
    theme_classic()+
    theme(legend.key.size = unit(0.25, 'cm'),axis.title = element_text(size = axis_title_size,family="Helvetica"),axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),plot.title = element_text(size = plot_title_size,family="Helvetica"),
          legend.title = element_text(size = axis_title_size,family="Helvetica"),legend.text = element_text(size = legend_title_size,family="Helvetica"),        strip.text = element_blank(),
          strip.background = element_blank()) +
    guides(fill=guide_legend(title=""))+
    ylab('\u0394 Max Pain Rating [LI - HI exercise]')+xlab('Stimulus Intensity [VAS]')+
    #ggtitle('Difference between Max Pain Rating [LI - HI exercise] Whole Stimulus Duration') +
    geom_hline(yintercept = 0)+
    ylim(-5, 15)+
    facet_wrap(~modality, labeller = as_labeller(c('1' = 'Pressure', '2' = 'Heat')),ncol = 1)
)


# 
# #----------------------------------------
# # Peak Pain ratings and FTP
# #---------------------------------------
# # Calculate mean difference ratings by subject
# mean_diff_ratings_sub <- mean_ratings_diff_sub %>%
#   group_by(subject, modality) %>%
#   summarise_at(c('diff_ints'),mean,na.rm = T)
# 
# # Calculate mean difference ratings plateau by subject
# mean_diff_ratings_plateau_sub <- mean_ratings_diff_plateau_sub %>%
#   group_by(subject, modality) %>%
#   summarise_at(c('diff_ints'),mean,na.rm = T)
# 
# # Calculate max difference ratings by subject
# max_diff_ratings_sub <- sub_max_online_ratings_diff %>%
#   group_by(subject, modality) %>%
#   summarise_at(c('diff_ints'),mean,na.rm = T)
# 
#   # Merge data frames by subject and modality
# diff_ratings_combined <- mean_diff_ratings_sub
# diff_ratings_combined$diff_ints_plateau <- mean_diff_ratings_plateau_sub$diff_ints
# diff_ratings_combined$diff_ints_max <- max_diff_ratings_sub$diff_ints
#   
# # Add the hotspin_calib$pwc by subject
# diff_ratings_combined <- diff_ratings_combined %>%
#   left_join(hotspin_calib %>% select(subject, pwc), by = "subject")
# 
# diff_ratings_combined$modality <- as.factor(diff_ratings_combined$modality)
#   
#   # Plot correlations
#   correlation_plot <- ggplot(diff_ratings_combined, aes(x = pwc, y = diff_ints,  color = modality, fill = modality)) +
#   geom_point(shape = 21, size = 1,alpha = 0.5) +
#     geom_smooth(method = "lm", se = T,alpha = 0.2) +
#     facet_wrap(~modality, labeller = as_labeller(c('1' = 'Pressure', '2' = 'Heat'))) +
#     theme_classic() +
#      theme(legend.key.size = unit(0.25, 'cm'),axis.title = element_text(size = axis_title_size,family="Helvetica"),axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),plot.title = element_text(size = plot_title_size,family="Helvetica"),legend.title = element_text(size = axis_title_size,family="Helvetica"),legend.text = element_text(size = legend_title_size,family="Helvetica"),strip.text.x = element_text(size = legend_text_size,family="Helvetica")) +
#     
#     scale_color_manual(labels = c("Pressure", "Heat"),
#                        values = c("#1C02C7", "#C75302"))+
#     scale_fill_manual(labels = c("Pressure", "Heat"),
#                       values = c("#1C02C7", "#C75302"))+
#     ylab("Mean Difference Ratings [LI-HI]") +
#     xlab("FTP [W/kg]") +
#     ggtitle("Correlation between  FTP (weight-corrected) and Mean Difference Ratings (Stimulus Duration)")
# 
#   correlation_plot_plateau <- ggplot(diff_ratings_combined, aes(x = pwc, y = diff_ints_plateau,  color = modality, fill = modality)) +
#   geom_point(shape = 21, size = 1,alpha = 0.5) +
#     geom_smooth(method = "lm", se = T,alpha = 0.2) +
#     facet_wrap(~modality, labeller = as_labeller(c('1' = 'Pressure', '2' = 'Heat'))) +
#     theme_classic() +
#      theme(legend.key.size = unit(0.25, 'cm'),axis.title = element_text(size = axis_title_size,family="Helvetica"),axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),plot.title = element_text(size = plot_title_size,family="Helvetica"),legend.title = element_text(size = axis_title_size,family="Helvetica"),legend.text = element_text(size = legend_title_size,family="Helvetica"),strip.text.x = element_text(size = legend_text_size,family="Helvetica")) +
#     
#         scale_color_manual(labels = c("Pressure", "Heat"),
#                        values = c("#1C02C7", "#C75302"))+
#     scale_fill_manual(labels = c("Pressure", "Heat"),
#                       values = c("#1C02C7", "#C75302"))+
# 
#     xlab("FTP [W/kg]") +
#     ylab("Mean Difference Ratings [LI-HI]") +
#     ggtitle("Correlation between  FTP (weight-corrected) and Mean Difference Ratings (Plateau Duration)")
# 
# 
# correlation_plot_max <- ggplot(diff_ratings_combined, aes(x = pwc, y = diff_ints_max, color = modality, fill = modality)) +
#   geom_point(shape = 21, size = 1,alpha = 0.5) +
#   geom_smooth(method = "lm", se = T, alpha = 0.2) +
#   facet_wrap(~modality, labeller = as_labeller(c('1' = 'Pressure', '2' = 'Heat'))) +
#   theme_classic() +
#    theme(legend.key.size = unit(0.25, 'cm'),axis.title = element_text(size = axis_title_size,family="Helvetica"),axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),plot.title = element_text(size = plot_title_size,family="Helvetica"),legend.title = element_text(size = axis_title_size,family="Helvetica"),legend.text = element_text(size = legend_title_size,family="Helvetica"),strip.text.x = element_text(size = legend_text_size,family="Helvetica")) +
#     
#     scale_color_manual(labels = c("Pressure", "Heat"),
#                        values = c("#1C02C7", "#C75302"))+
#     scale_fill_manual(labels = c("Pressure", "Heat"),
#                       values = c("#1C02C7", "#C75302"))+
#   xlab("FTP [W/kg]") +
#   ylab("Max Difference Ratings [LI-HI]") +
#   ggtitle("Correlation between FTP (weight-corrected) and Max Difference Ratings")
# 
# 
#   # Calculate correlations and p-values for heat and pressure separately
# 
#   # Function to calculate correlation and p-value
#   correlation_test <- function(x, y) {
#     cor_test <- cor.test(x, y, use = "complete.obs")
#     list(correlation = cor_test$estimate, p_value = cor_test$p.value)
#   }
# 
#   # Correlation between FTP and mean difference ratings (whole stimulus duration)
#   correlation_mean_diff_heat <- correlation_test(diff_ratings_combined$pwc[diff_ratings_combined$modality == 2], diff_ratings_combined$diff_ints[diff_ratings_combined$modality == 2])
#   correlation_mean_diff_pressure <- correlation_test(diff_ratings_combined$pwc[diff_ratings_combined$modality == 1], diff_ratings_combined$diff_ints[diff_ratings_combined$modality == 1])
# 
#   # Correlation between FTP and mean difference ratings (plateau duration)
#   correlation_mean_diff_plateau_heat <- correlation_test(diff_ratings_combined$pwc[diff_ratings_combined$modality == 2], diff_ratings_combined$diff_ints_plateau[diff_ratings_combined$modality == 2])
#   correlation_mean_diff_plateau_pressure <- correlation_test(diff_ratings_combined$pwc[diff_ratings_combined$modality == 1], diff_ratings_combined$diff_ints_plateau[diff_ratings_combined$modality == 1])
# 
#   # Correlation between FTP and max difference ratings
#   correlation_max_diff_heat <- correlation_test(diff_ratings_combined$pwc[diff_ratings_combined$modality == 2], diff_ratings_combined$diff_ints_max[diff_ratings_combined$modality == 2])
#   correlation_max_diff_pressure <- correlation_test(diff_ratings_combined$pwc[diff_ratings_combined$modality == 1], diff_ratings_combined$diff_ints_max[diff_ratings_combined$modality == 1])
# 
#   # Print correlations and p-values
#   print(correlation_mean_diff_heat)
#   print(correlation_mean_diff_pressure)
#   print(correlation_mean_diff_plateau_heat)
#   print(correlation_mean_diff_plateau_pressure)
#   print(correlation_max_diff_heat)
#   print(correlation_max_diff_pressure)
#   
# 
# 
# #----------------------------------------
# # Peak Pain ratings and FTP at VAS 70
# #---------------------------------------
# 
# # Calculate mean difference ratings by subject
# mean_diff_ratings_sub <- mean_ratings_diff_sub %>%
#   group_by(subject, modality,VAS_intensity) %>%
#   summarise_at(c('diff_ints'),mean,na.rm = T)
# 
# # Calculate mean difference ratings plateau by subject
# mean_diff_ratings_plateau_sub <- mean_ratings_diff_plateau_sub %>%
#   group_by(subject, modality,VAS_intensity) %>%
#   summarise_at(c('diff_ints'),mean,na.rm = T)
# 
# # Calculate max difference ratings by subject
# max_diff_ratings_sub <- sub_max_online_ratings_diff %>%
#   group_by(subject, modality,VAS_intensity) %>%
#   summarise_at(c('diff_ints'),mean,na.rm = T)
# 
#   # Merge data frames by subject and modality
# diff_ratings_combined <- mean_diff_ratings_sub
# diff_ratings_combined$diff_ints_plateau <- mean_diff_ratings_plateau_sub$diff_ints
# diff_ratings_combined$diff_ints_max <- max_diff_ratings_sub$diff_ints
#   
# # Add the hotspin_calib$pwc by subject
# diff_ratings_combined <- diff_ratings_combined %>%
#   left_join(hotspin_calib %>% select(subject, pwc), by = "subject")
# 
# diff_ratings_combined$modality <- as.factor(diff_ratings_combined$modality)
# 
# diff_ratings_combined_VAS70 <- diff_ratings_combined %>% filter(VAS_intensity == "70")
# 
#   # Plot correlations
#   correlation_plot_VAS70 <- ggplot(diff_ratings_combined_VAS70, aes(x = pwc, y = diff_ints,  color = modality, fill = modality)) +
#   geom_point(shape = 21, size = 1,alpha = 0.5) +
#     geom_smooth(method = "lm", se = T,alpha = 0.2) +
#     facet_wrap(~modality, labeller = as_labeller(c('1' = 'Pressure', '2' = 'Heat'))) +
#     theme_classic() +
#      theme(legend.key.size = unit(0.25, 'cm'),axis.title = element_text(size = axis_title_size,family="Helvetica"),axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),plot.title = element_text(size = plot_title_size,family="Helvetica"),legend.title = element_text(size = axis_title_size,family="Helvetica"),legend.text = element_text(size = legend_title_size,family="Helvetica"),strip.text.x = element_text(size = legend_text_size,family="Helvetica")) +
#     
#     scale_color_manual(labels = c("Pressure", "Heat"),
#                        values = c("#1C02C7", "#C75302"))+
#     scale_fill_manual(labels = c("Pressure", "Heat"),
#                       values = c("#1C02C7", "#C75302"))+
#     ylab("Mean Difference Ratings [LI-HI]") +
#     xlab("FTP [W/kg]") +
#     ggtitle("Correlation between  FTP (weight-corrected) and Mean Difference Ratings (Stimulus Duration) at VAS 70")
# 
#   correlation_plot_plateau_VAS70 <- ggplot(diff_ratings_combined_VAS70, aes(x = pwc, y = diff_ints_plateau,  color = modality, fill = modality)) +
#   geom_point(shape = 21, size = 1,alpha = 0.5) +
#     geom_smooth(method = "lm", se = T,alpha = 0.2) +
#     facet_wrap(~modality, labeller = as_labeller(c('1' = 'Pressure', '2' = 'Heat'))) +
#     theme_classic() +
#      theme(legend.key.size = unit(0.25, 'cm'),axis.title = element_text(size = axis_title_size,family="Helvetica"),axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),plot.title = element_text(size = plot_title_size,family="Helvetica"),legend.title = element_text(size = axis_title_size,family="Helvetica"),legend.text = element_text(size = legend_title_size,family="Helvetica"),strip.text.x = element_text(size = legend_text_size,family="Helvetica")) +
#     
#         scale_color_manual(labels = c("Pressure", "Heat"),
#                        values = c("#1C02C7", "#C75302"))+
#     scale_fill_manual(labels = c("Pressure", "Heat"),
#                       values = c("#1C02C7", "#C75302"))+
# 
#     xlab("FTP [W/kg]") +
#     ylab("Mean Difference Ratings [LI-HI]") +
#     ggtitle("Correlation between  FTP (weight-corrected) and Mean Difference Ratings (Plateau Duration) at VAS 70")
# 
# 
# correlation_plot_max_VAS70 <- ggplot(diff_ratings_combined_VAS70, aes(x = pwc, y = diff_ints_max, color = modality, fill = modality)) +
#   geom_point(shape = 21, size = 1,alpha = 0.5) +
#   geom_smooth(method = "lm", se = T, alpha = 0.2) +
#   facet_wrap(~modality, labeller = as_labeller(c('1' = 'Pressure', '2' = 'Heat'))) +
#   theme_classic() +
#    theme(legend.key.size = unit(0.25, 'cm'),axis.title = element_text(size = axis_title_size,family="Helvetica"),axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),plot.title = element_text(size = plot_title_size,family="Helvetica"),legend.title = element_text(size = axis_title_size,family="Helvetica"),legend.text = element_text(size = legend_title_size,family="Helvetica"),strip.text.x = element_text(size = legend_text_size,family="Helvetica")) +
#     
#     scale_color_manual(labels = c("Pressure", "Heat"),
#                        values = c("#1C02C7", "#C75302"))+
#     scale_fill_manual(labels = c("Pressure", "Heat"),
#                       values = c("#1C02C7", "#C75302"))+
#   xlab("FTP [W/kg]") +
#   ylab("Max Difference Ratings [LI-HI]") +
#   ggtitle("Correlation between FTP (weight-corrected) and Max Difference Ratings at VAS 70")
# 
#   # Calculate correlations for each plot
#   # Calculate correlations and p-values for heat and pressure separately
# 
#   # Function to calculate correlation and p-value
#   correlation_test <- function(x, y) {
#     cor_test <- cor.test(x, y, use = "complete.obs")
#     list(correlation = cor_test$estimate, p_value = cor_test$p.value)
#   }
# 
#   # Correlation between FTP and mean difference ratings (whole stimulus duration)
#   correlation_mean_diff_heat <- correlation_test(diff_ratings_combined_VAS70$pwc[diff_ratings_combined_VAS70$modality == 2], diff_ratings_combined_VAS70$diff_ints[diff_ratings_combined_VAS70$modality == 2])
#   correlation_mean_diff_pressure <- correlation_test(diff_ratings_combined_VAS70$pwc[diff_ratings_combined_VAS70$modality == 1], diff_ratings_combined_VAS70$diff_ints[diff_ratings_combined_VAS70$modality == 1])
# 
#   # Correlation between FTP and mean difference ratings (plateau duration)
#   correlation_mean_diff_plateau_heat <- correlation_test(diff_ratings_combined_VAS70$pwc[diff_ratings_combined_VAS70$modality == 2], diff_ratings_combined_VAS70$diff_ints_plateau[diff_ratings_combined_VAS70$modality == 2])
#   correlation_mean_diff_plateau_pressure <- correlation_test(diff_ratings_combined_VAS70$pwc[diff_ratings_combined_VAS70$modality == 1], diff_ratings_combined_VAS70$diff_ints_plateau[diff_ratings_combined_VAS70$modality == 1])
# 
#   # Correlation between FTP and max difference ratings
#   correlation_max_diff_heat <- correlation_test(diff_ratings_combined_VAS70$pwc[diff_ratings_combined_VAS70$modality == 2], diff_ratings_combined_VAS70$diff_ints_max[diff_ratings_combined_VAS70$modality == 2])
#   correlation_max_diff_pressure <- correlation_test(diff_ratings_combined_VAS70$pwc[diff_ratings_combined_VAS70$modality == 1], diff_ratings_combined_VAS70$diff_ints_max[diff_ratings_combined_VAS70$modality == 1])
# 
#   # Print correlations and p-values
#   print(correlation_mean_diff_heat)
#   print(correlation_mean_diff_pressure)
#   print(correlation_mean_diff_plateau_heat)
#   print(correlation_mean_diff_plateau_pressure)
#   print(correlation_max_diff_heat)
#   print(correlation_max_diff_pressure)
#   



  # --------------------------------------------
  # Correlation Fitness on EIH HEAT
  #-------------------------------------------
  
  
  # ---- Fitness level On difference pain ratings:  
  summary_lmer_models_main<- df_combined_heat %>%
    group_by(subject,exercise_intensity,pwc,treatment_order,gender,group)%>%
    summarise_at(c('pain_rating'),mean,na.rm = T)
  
  lmer_model_df <- spread(summary_lmer_models_main, exercise_intensity, pain_rating)
  lmer_model_df$diff_hi_low_rating <- lmer_model_df$'0' -  lmer_model_df$'1'  # Low Intensity - High Intensity Data
  
  
  # fit regular model
  lme_model_main <- lm(diff_hi_low_rating ~   pwc + treatment_order, data = lmer_model_df)
  summary(lme_model_main)
  confint(lme_model_main)
  
  # Test Correlation
  cor.test(lmer_model_df$diff_hi_low_rating,lmer_model_df$pwc)
  
  #------- Fitness level and Sex On difference pain ratings:  

  # fit regular model
  lme_model_main <- lm(diff_hi_low_rating ~   pwc*gender + treatment_order, data = lmer_model_df)
  summary(lme_model_main)
  confint(lme_model_main)
  
  # Test Correlation
  cor.test(lmer_model_df$diff_hi_low_rating[lmer_model_df$gender==0],lmer_model_df$pwc[lmer_model_df$gender==0])
  cor.test(lmer_model_df$diff_hi_low_rating[lmer_model_df$gender==1],lmer_model_df$pwc[lmer_model_df$gender==1])
  
  
  #------- Fitness level and Group On difference pain ratings:  
  lmer_model_df$group <- as.integer(lmer_model_df$group)
  # fit regular model
  lme_model_main <- lm(diff_hi_low_rating ~   pwc*group + treatment_order, data = lmer_model_df)
  summary(lme_model_main)
  confint(lme_model_main)
  
  # Test Correlation
  cor.test(lmer_model_df$diff_hi_low_rating[lmer_model_df$group==1],lmer_model_df$pwc[lmer_model_df$group==1])
  cor.test(lmer_model_df$diff_hi_low_rating[lmer_model_df$group==2],lmer_model_df$pwc[lmer_model_df$group==2])
  cor.test(lmer_model_df$diff_hi_low_rating[lmer_model_df$group==3],lmer_model_df$pwc[lmer_model_df$group==3])
  
  
  
  #---------------- Visualisation
  
  (peep_hotspin_heat_saline_fitness <-  ggplot(lmer_model_df,aes(x = pwc,y=diff_hi_low_rating,colour = "#024873",fill = "#024873"))+
      geom_point(size = 0.5,aes(colour = exercise_intensity),shape = 21,colour = 'black',alpha = 0.6,show.legend = F)+
      geom_smooth(method = 'lm',alpha = 0.1,size = 1,se = T,show.legend = F)+
      theme_classic(base_size = 11) +
      scale_color_manual(labels = c("SAL"),
                         values = c("#024873")) +
      scale_fill_manual(labels = c("SAL"),
                        values = c("#024873"))+
      geom_hline(yintercept = 0,colour = 'black',size = 0.5)+
      theme(axis.title = element_text(size = axis_title_size,family="Helvetica"),
            axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),
            plot.title = element_text(size = plot_title_size,family="Helvetica"),legend.title = element_text(size = axis_title_size,family="Helvetica"),
            legend.text = element_text(size = legend_title_size,family="Helvetica"),
            strip.text.x = element_text(size = legend_text_size,family="Helvetica")) +
      #guides(fill=guide_legend(title=""))+
      ylab('\u0394 Pain Rating\n[LI - HI Exercise Pain Rating]')+xlab(expression(paste('FTP [Watt * kg'^-1,']')))+
      #stat_cor(method = "pearson",alternative = 'two.sided',label.sep = "\n", size = 2,show.legend = F,colour = 'black')+
      ggtitle('')
      )

  
  #------------- Visualisation by gender:
  
  (peep_hotspin_heat_saline_fitness_gender <-  ggplot(lmer_model_df,aes(x = pwc,y=diff_hi_low_rating,colour = gender,fill = gender))+
      geom_point(size = 0.5,aes(colour = exercise_intensity),shape = 21,colour = 'black',alpha = 0.6,show.legend = F)+
      geom_smooth(method = 'lm',alpha = 0.1,size = 1,se = T,show.legend = F)+
      theme_classic(base_size = 11) +
      scale_color_manual(labels = c("Female", "Male"),
                         values = c("#AF00C7", "#C75302"))+
      scale_fill_manual(labels = c("Female", "Male"),
                        values = c("#AF00C7", "#C75302"))+
      geom_hline(yintercept = 0,colour = 'black',size = 0.5)+
      theme(axis.title = element_text(size = axis_title_size,family="Helvetica"),
            axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),
            plot.title = element_text(size = plot_title_size,family="Helvetica"),legend.title = element_text(size = axis_title_size,family="Helvetica"),
            legend.text = element_text(size = legend_title_size,family="Helvetica"),
            strip.text.x = element_text(size = legend_text_size,family="Helvetica")) +
      #guides(fill=guide_legend(title=""))+
      ylab('\u0394 Pain Rating\n[LI - HI Exercise Pain Rating]')+xlab(expression(paste('FTP [Watt * kg'^-1,']')))+
      #stat_cor(method = "pearson",alternative = 'two.sided',label.sep = "\n", size = 2,show.legend = F,colour = 'black')+
      ggtitle('')
  )
  
  
  #------------- Visualisation by group:
  lmer_model_df$group <- as.factor(lmer_model_df$group)
  (peep_hotspin_heat_saline_fitness_group <-  ggplot(lmer_model_df,aes(x = pwc,y=diff_hi_low_rating,colour = group,fill = group))+
      geom_point(size = 0.5,aes(colour = exercise_intensity),shape = 21,colour = 'black',alpha = 0.6,show.legend = F)+
      geom_smooth(method = 'lm',alpha = 0.1,size = 1,se = T,show.legend = F)+
      theme_classic(base_size = 11) +
      scale_color_manual(labels = c("Females (current study), Females (previous study), Males (previous study)"),
                         values = c("#C75302", "#8480F2", "#1C02C7")) +
      scale_fill_manual(labels = c("Females (current study), Females (previous study), Males (previous study)"),
                        values = c("#C75302", "#8480F2", "#1C02C7")) +
      geom_hline(yintercept = 0,colour = 'black',size = 0.5)+
      theme(axis.title = element_text(size = axis_title_size,family="Helvetica"),
            axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),
            plot.title = element_text(size = plot_title_size,family="Helvetica"),legend.title = element_text(size = axis_title_size,family="Helvetica"),
            legend.text = element_text(size = legend_title_size,family="Helvetica"),
            strip.text.x = element_text(size = legend_text_size,family="Helvetica")) +
      #guides(fill=guide_legend(title=""))+
      ylab('\u0394 Pain Rating\n[LI - HI Exercise Pain Rating]')+xlab(expression(paste('FTP [Watt * kg'^-1,']')))+
      #stat_cor(method = "pearson",alternative = 'two.sided',label.sep = "\n", size = 2,show.legend = F,colour = 'black')+
      ggtitle('')
      )
  
  
  #----------------------------------------------
  # Calculate with training hours / week instead of FTP for heat
  #----------------------------------------------

  # ---- Training status On difference pain ratings:  
  summary_lmer_models_main<- df_combined_heat %>%
    group_by(subject,exercise_intensity,hours_week,treatment_order,gender,group)%>%
    summarise_at(c('pain_rating'),mean,na.rm = T)
  
  lmer_model_df <- spread(summary_lmer_models_main, exercise_intensity, pain_rating)
  lmer_model_df$diff_hi_low_rating <- lmer_model_df$'0' -  lmer_model_df$'1'  # Low Intensity - High Intensity Data
  
  
  # fit regular model
  lme_model_main <- lm(diff_hi_low_rating ~   hours_week + treatment_order, data = lmer_model_df)
  summary(lme_model_main)
  confint(lme_model_main)
  
  # Test Correlation
  cor.test(lmer_model_df$diff_hi_low_rating,lmer_model_df$hours_week)
  
  #------- TRAVOL and Sex On difference pain ratings:  
  
  # fit regular model
  lme_model_main <- lm(diff_hi_low_rating ~   hours_week*gender + treatment_order, data = lmer_model_df)
  summary(lme_model_main)
  confint(lme_model_main)
  
  # Test Correlation
  cor.test(lmer_model_df$diff_hi_low_rating[lmer_model_df$gender==0],lmer_model_df$hours_week[lmer_model_df$gender==0])
  cor.test(lmer_model_df$diff_hi_low_rating[lmer_model_df$gender==1],lmer_model_df$hours_week[lmer_model_df$gender==1])
  
  
  #------- Fitness level and Group On difference pain ratings:  
  lmer_model_df$group <- as.integer(lmer_model_df$group)
  # fit regular model
  lme_model_main <- lm(diff_hi_low_rating ~   hours_week*group + treatment_order, data = lmer_model_df)
  summary(lme_model_main)
  confint(lme_model_main)
  
  # Test Correlation
  cor.test(lmer_model_df$diff_hi_low_rating[lmer_model_df$group==1],lmer_model_df$hours_week[lmer_model_df$group==1])
  cor.test(lmer_model_df$diff_hi_low_rating[lmer_model_df$group==2],lmer_model_df$hours_week[lmer_model_df$group==2])
  cor.test(lmer_model_df$diff_hi_low_rating[lmer_model_df$group==3],lmer_model_df$hours_week[lmer_model_df$group==3])
  
  
  
  #---------------- Visualisation
  lmer_model_df$group <- as.factor(lmer_model_df$group)
  
  (peep_hotspin_heat_saline_fitness_travol <-  ggplot(lmer_model_df,aes(x = hours_week,y=diff_hi_low_rating,colour = "#024873",fill = "#024873"))+
      geom_point(size = 0.5,aes(colour = exercise_intensity),shape = 21,colour = 'black',alpha = 0.6,show.legend = F)+
      geom_smooth(method = 'lm',alpha = 0.1,size = 1,se = T,show.legend = F)+
      theme_classic()+
      scale_color_manual(labels = c("SAL"),
                         values = c("#024873")) +
      scale_fill_manual(labels = c("SAL"),
                        values = c("#024873"))+
      geom_hline(yintercept = 0,colour = 'black',size = 0.5)+
      theme(axis.title = element_text(size = axis_title_size,family="Helvetica"),
            axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),
            plot.title = element_text(size = plot_title_size,family="Helvetica"),legend.title = element_text(size = axis_title_size,family="Helvetica"),
            legend.text = element_text(size = legend_title_size,family="Helvetica"),
            strip.text.x = element_text(size = legend_text_size,family="Helvetica")) +
      #guides(fill=guide_legend(title=""))+
      ylab('\u0394 Pain Rating\n[LI - HI Exercise Pain Rating]')+xlab(expression(paste('Training Volume (h/week)')))+
      #stat_cor(method = "pearson",alternative = 'two.sided',label.sep = "\n", size = 2,show.legend = F,colour = 'black')+
      ggtitle('Heat')
  )
  
  
  #------------- Visualisation by gender:
  
  (peep_hotspin_heat_saline_fitness_gender_travol <-  ggplot(lmer_model_df,aes(x = hours_week,y=diff_hi_low_rating,colour = gender,fill = gender))+
      geom_point(size = 0.5,aes(colour = exercise_intensity),shape = 21,colour = 'black',alpha = 0.6,show.legend = F)+
      geom_smooth(method = 'lm',alpha = 0.1,size = 1,se = T,show.legend = F)+
      theme_classic()+
      scale_color_manual(labels = c("Female", "Male"),
                         values = c("#AF00C7", "#C75302"))+
      scale_fill_manual(labels = c("Female", "Male"),
                        values = c("#AF00C7", "#C75302"))+
      geom_hline(yintercept = 0,colour = 'black',size = 0.5)+
      theme(axis.title = element_text(size = axis_title_size,family="Helvetica"),
            axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),
            plot.title = element_text(size = plot_title_size,family="Helvetica"),legend.title = element_text(size = axis_title_size,family="Helvetica"),
            legend.text = element_text(size = legend_title_size,family="Helvetica"),
            strip.text.x = element_text(size = legend_text_size,family="Helvetica")) +
      #guides(fill=guide_legend(title=""))+
      ylab('\u0394 Pain Rating\n[LI - HI Exercise Pain Rating]')+xlab(expression(paste('Training Volume (h/week)')))+
      #stat_cor(method = "pearson",alternative = 'two.sided',label.sep = "\n", size = 2,show.legend = F,colour = 'black')+
      ggtitle('')
  )
  
  
  #------------- Visualisation by group:
  
  (peep_hotspin_heat_saline_fitness_group_travol <-  ggplot(lmer_model_df,aes(x = hours_week,y=diff_hi_low_rating,colour = group,fill = group))+
      geom_point(size = 0.5,aes(colour = exercise_intensity),shape = 21,colour = 'black',alpha = 0.6,show.legend = F)+
      geom_smooth(method = 'lm',alpha = 0.1,size = 1,se = T,show.legend = F)+
      theme_classic()+
      scale_color_manual(labels = c("Females (current study), Females (previous study), Males (previous study)"),
                         values = c("#C75302", "#8480F2", "#1C02C7")) +
      scale_fill_manual(labels = c("Females (current study), Females (previous study), Males (previous study)"),
                        values = c("#C75302", "#8480F2", "#1C02C7")) +
      geom_hline(yintercept = 0,colour = 'black',size = 0.5)+
      theme(axis.title = element_text(size = axis_title_size,family="Helvetica"),
            axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),
            plot.title = element_text(size = plot_title_size,family="Helvetica"),legend.title = element_text(size = axis_title_size,family="Helvetica"),
            legend.text = element_text(size = legend_title_size,family="Helvetica"),
            strip.text.x = element_text(size = legend_text_size,family="Helvetica")) +
      #guides(fill=guide_legend(title=""))+
      ylab('\u0394 Pain Rating\n[LI - HI Exercise Pain Rating]')+xlab(expression(paste('Training Volume (h/week)')))+
      #stat_cor(method = "pearson",alternative = 'two.sided',label.sep = "\n", size = 2,show.legend = F,colour = 'black')+
      ggtitle('')
  )
  

  
  # --------------------------------------------
  # Correlation Fitness on EIH Pressure
  #-------------------------------------------
  
  
  # ---- Fitness level On difference pain ratings:  
  summary_lmer_models_main<- df_combined_pressure %>%
    group_by(subject,exercise_intensity,pwc,treatment_order,gender,group)%>%
    summarise_at(c('pain_rating'),mean,na.rm = T)
  
  lmer_model_df <- spread(summary_lmer_models_main, exercise_intensity, pain_rating)
  lmer_model_df$diff_hi_low_rating <- lmer_model_df$'0' -  lmer_model_df$'1'  # Low Intensity - High Intensity Data
  
  
  # fit regular model
  lme_model_main <- lm(diff_hi_low_rating ~   pwc + treatment_order, data = lmer_model_df)
  summary(lme_model_main)
  confint(lme_model_main)
  
  # Test Correlation
  cor.test(lmer_model_df$diff_hi_low_rating,lmer_model_df$pwc)
  
  #------- Fitness level and Sex On difference pain ratings:  
  
  # fit regular model
  lme_model_main <- lm(diff_hi_low_rating ~   pwc*gender + treatment_order, data = lmer_model_df)
  summary(lme_model_main)
  confint(lme_model_main)
  
  # Test Correlation
  cor.test(lmer_model_df$diff_hi_low_rating[lmer_model_df$gender==0],lmer_model_df$pwc[lmer_model_df$gender==0])
  cor.test(lmer_model_df$diff_hi_low_rating[lmer_model_df$gender==1],lmer_model_df$pwc[lmer_model_df$gender==1])
  
  
  #------- Fitness level and Group On difference pain ratings:  
  lmer_model_df$group <- as.integer(lmer_model_df$group)
  # fit regular model
  lme_model_main <- lm(diff_hi_low_rating ~   pwc*group + treatment_order, data = lmer_model_df)
  summary(lme_model_main)
  confint(lme_model_main)
  
  # Test Correlation
  cor.test(lmer_model_df$diff_hi_low_rating[lmer_model_df$group==1],lmer_model_df$pwc[lmer_model_df$group==1])
  cor.test(lmer_model_df$diff_hi_low_rating[lmer_model_df$group==2],lmer_model_df$pwc[lmer_model_df$group==2])
  cor.test(lmer_model_df$diff_hi_low_rating[lmer_model_df$group==3],lmer_model_df$pwc[lmer_model_df$group==3])
  
  
  
  #---------------- Visualisation
  
  (peep_hotspin_pressure_saline_fitness <-  ggplot(lmer_model_df,aes(x = pwc,y=diff_hi_low_rating,colour = "#024873",fill = "#024873"))+
      geom_point(size = 0.5,aes(colour = exercise_intensity),shape = 21,colour = 'black',alpha = 0.6,show.legend = F)+
      geom_smooth(method = 'lm',alpha = 0.1,size = 1,se = T,show.legend = F)+
      theme_classic()+
      scale_color_manual(labels = c("SAL"),
                         values = c("#024873")) +
      scale_fill_manual(labels = c("SAL"),
                        values = c("#024873"))+
      geom_hline(yintercept = 0,colour = 'black',size = 0.5)+
      theme(axis.title = element_text(size = axis_title_size,family="Helvetica"),
            axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),
            plot.title = element_text(size = plot_title_size,family="Helvetica"),legend.title = element_text(size = axis_title_size,family="Helvetica"),
            legend.text = element_text(size = legend_title_size,family="Helvetica"),
            strip.text.x = element_text(size = legend_text_size,family="Helvetica")) +
      #guides(fill=guide_legend(title=""))+
      ylab('\u0394 Pain Rating\n[LI - HI Exercise Pain Rating]')+xlab(expression(paste('FTP [Watt * kg'^-1,']')))+
      #stat_cor(method = "pearson",alternative = 'two.sided',label.sep = "\n", size = 2,show.legend = F,colour = 'black')+
      ggtitle('')
  )
  
  
  #------------- Visualisation by gender:
  
  (peep_hotspin_pressure_saline_fitness_gender <-  ggplot(lmer_model_df,aes(x = pwc,y=diff_hi_low_rating,colour = gender,fill = gender))+
      geom_point(size = 0.5,aes(colour = exercise_intensity),shape = 21,colour = 'black',alpha = 0.6,show.legend = F)+
      geom_smooth(method = 'lm',alpha = 0.1,size = 1,se = T,show.legend = F)+
      theme_classic()+
      scale_color_manual(labels = c("Female", "Male"),
                         values = c("#AF00C7", "#C75302"))+
      scale_fill_manual(labels = c("Female", "Male"),
                        values = c("#AF00C7", "#C75302"))+
      geom_hline(yintercept = 0,colour = 'black',size = 0.5)+
      theme(axis.title = element_text(size = axis_title_size,family="Helvetica"),
            axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),
            plot.title = element_text(size = plot_title_size,family="Helvetica"),legend.title = element_text(size = axis_title_size,family="Helvetica"),
            legend.text = element_text(size = legend_title_size,family="Helvetica"),
            strip.text.x = element_text(size = legend_text_size,family="Helvetica")) +
      #guides(fill=guide_legend(title=""))+
      ylab('\u0394 Pain Rating\n[LI - HI Exercise Pain Rating]')+xlab(expression(paste('FTP [Watt * kg'^-1,']')))+
      #stat_cor(method = "pearson",alternative = 'two.sided',label.sep = "\n", size = 2,show.legend = F,colour = 'black')+
      ggtitle('')
  )
  
  
  #------------- Visualisation by group:
  lmer_model_df$group <- as.factor(lmer_model_df$group)
  (peep_hotspin_pressure_saline_fitness_group <-  ggplot(lmer_model_df,aes(x = pwc,y=diff_hi_low_rating,colour = group,fill = group))+
      geom_point(size = 0.5,aes(colour = exercise_intensity),shape = 21,colour = 'black',alpha = 0.6,show.legend = F)+
      geom_smooth(method = 'lm',alpha = 0.1,size = 1,se = T,show.legend = F)+
      theme_classic()+
      scale_color_manual(labels = c("Females (current study), Females (previous study), Males (previous study)"),
                         values = c("#C75302", "#8480F2", "#1C02C7")) +
      scale_fill_manual(labels = c("Females (current study), Females (previous study), Males (previous study)"),
                        values = c("#C75302", "#8480F2", "#1C02C7")) +
      geom_hline(yintercept = 0,colour = 'black',size = 0.5)+
      theme(axis.title = element_text(size = axis_title_size,family="Helvetica"),
            axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),
            plot.title = element_text(size = plot_title_size,family="Helvetica"),legend.title = element_text(size = axis_title_size,family="Helvetica"),
            legend.text = element_text(size = legend_title_size,family="Helvetica"),
            strip.text.x = element_text(size = legend_text_size,family="Helvetica")) +
      #guides(fill=guide_legend(title=""))+
      ylab('\u0394 Pain Rating\n[LI - HI Exercise Pain Rating]')+xlab(expression(paste('FTP [Watt * kg'^-1,']')))+
      #stat_cor(method = "pearson",alternative = 'two.sided',label.sep = "\n", size = 2,show.legend = F,colour = 'black')+
      ggtitle('')
  )
  
  
  #----------------------------------------------
  # Calculate with training hours / week instead of FTP for pressure
  #----------------------------------------------
  
  # ---- Training status On difference pain ratings:  
  summary_lmer_models_main<- df_combined_pressure %>%
    group_by(subject,exercise_intensity,hours_week,treatment_order,gender,group)%>%
    summarise_at(c('pain_rating'),mean,na.rm = T)
  
  lmer_model_df <- spread(summary_lmer_models_main, exercise_intensity, pain_rating)
  lmer_model_df$diff_hi_low_rating <- lmer_model_df$'0' -  lmer_model_df$'1'  # Low Intensity - High Intensity Data
  
  
  # fit regular model
  lme_model_main <- lm(diff_hi_low_rating ~   hours_week + treatment_order, data = lmer_model_df)
  summary(lme_model_main)
  confint(lme_model_main)
  
  # Test Correlation
  cor.test(lmer_model_df$diff_hi_low_rating,lmer_model_df$hours_week)
  
  #------- TRAVOL and Sex On difference pain ratings:  
  
  # fit regular model
  lme_model_main <- lm(diff_hi_low_rating ~   hours_week*gender + treatment_order, data = lmer_model_df)
  summary(lme_model_main)
  confint(lme_model_main)
  
  # Test Correlation
  cor.test(lmer_model_df$diff_hi_low_rating[lmer_model_df$gender==0],lmer_model_df$hours_week[lmer_model_df$gender==0])
  cor.test(lmer_model_df$diff_hi_low_rating[lmer_model_df$gender==1],lmer_model_df$hours_week[lmer_model_df$gender==1])
  
  
  #------- Fitness level and Group On difference pain ratings:  
  lmer_model_df$group <- as.integer(lmer_model_df$group)
  # fit regular model
  lme_model_main <- lm(diff_hi_low_rating ~   hours_week*group + treatment_order, data = lmer_model_df)
  summary(lme_model_main)
  confint(lme_model_main)
  
  # Test Correlation
  cor.test(lmer_model_df$diff_hi_low_rating[lmer_model_df$group==1],lmer_model_df$hours_week[lmer_model_df$group==1])
  cor.test(lmer_model_df$diff_hi_low_rating[lmer_model_df$group==2],lmer_model_df$hours_week[lmer_model_df$group==2])
  cor.test(lmer_model_df$diff_hi_low_rating[lmer_model_df$group==3],lmer_model_df$hours_week[lmer_model_df$group==3])
  
  
  
  #---------------- Visualisation
  lmer_model_df$group <- as.factor(lmer_model_df$group)
  
  (peep_hotspin_pressure_saline_fitness_travol <-  ggplot(lmer_model_df,aes(x = hours_week,y=diff_hi_low_rating,colour = "#024873",fill = "#024873"))+
      geom_point(size = 0.5,aes(colour = exercise_intensity),shape = 21,colour = 'black',alpha = 0.6,show.legend = F)+
      geom_smooth(method = 'lm',alpha = 0.1,size = 1,se = T,show.legend = F)+
      theme_classic()+
      scale_color_manual(labels = c("SAL"),
                         values = c("#024873")) +
      scale_fill_manual(labels = c("SAL"),
                        values = c("#024873"))+
      geom_hline(yintercept = 0,colour = 'black',size = 0.5)+
      theme(axis.title = element_text(size = axis_title_size,family="Helvetica"),
            axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),
            plot.title = element_text(size = plot_title_size,family="Helvetica"),legend.title = element_text(size = axis_title_size,family="Helvetica"),
            legend.text = element_text(size = legend_title_size,family="Helvetica"),
            strip.text.x = element_text(size = legend_text_size,family="Helvetica")) +
      #guides(fill=guide_legend(title=""))+
      ylab('\u0394 Pain Rating\n[LI - HI Exercise Pain Rating]')+xlab(expression(paste('Training Volume (h/week)')))+
      #stat_cor(method = "pearson",alternative = 'two.sided',label.sep = "\n", size = 2,show.legend = F,colour = 'black')+
      ggtitle('Pressure')
  )
  
  
  #------------- Visualisation by gender:
  
  (peep_hotspin_pressure_saline_fitness_gender_travol <-  ggplot(lmer_model_df,aes(x = hours_week,y=diff_hi_low_rating,colour = gender,fill = gender))+
      geom_point(size = 0.5,aes(colour = exercise_intensity),shape = 21,colour = 'black',alpha = 0.6,show.legend = F)+
      geom_smooth(method = 'lm',alpha = 0.1,size = 1,se = T,show.legend = F)+
      theme_classic()+
      scale_color_manual(labels = c("Female", "Male"),
                         values = c("#AF00C7", "#C75302"))+
      scale_fill_manual(labels = c("Female", "Male"),
                        values = c("#AF00C7", "#C75302"))+
      geom_hline(yintercept = 0,colour = 'black',size = 0.5)+
      theme(axis.title = element_text(size = axis_title_size,family="Helvetica"),
            axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),
            plot.title = element_text(size = plot_title_size,family="Helvetica"),legend.title = element_text(size = axis_title_size,family="Helvetica"),
            legend.text = element_text(size = legend_title_size,family="Helvetica"),
            strip.text.x = element_text(size = legend_text_size,family="Helvetica")) +
      #guides(fill=guide_legend(title=""))+
      ylab('\u0394 Pain Rating\n[LI - HI Exercise Pain Rating]')+xlab(expression(paste('Training Volume (h/week)')))+
      #stat_cor(method = "pearson",alternative = 'two.sided',label.sep = "\n", size = 2,show.legend = F,colour = 'black')+
      ggtitle('')
  )
  
  
  #------------- Visualisation by group:
  
  (peep_hotspin_pressure_saline_fitness_group_travol <-  ggplot(lmer_model_df,aes(x = hours_week,y=diff_hi_low_rating,colour = group,fill = group))+
      geom_point(size = 0.5,aes(colour = exercise_intensity),shape = 21,colour = 'black',alpha = 0.6,show.legend = F)+
      geom_smooth(method = 'lm',alpha = 0.1,size = 1,se = T,show.legend = F)+
      theme_classic()+
      scale_color_manual(labels = c("Females (current study), Females (previous study), Males (previous study)"),
                         values = c("#C75302", "#8480F2", "#1C02C7")) +
      scale_fill_manual(labels = c("Females (current study), Females (previous study), Males (previous study)"),
                        values = c("#C75302", "#8480F2", "#1C02C7")) +
      geom_hline(yintercept = 0,colour = 'black',size = 0.5)+
      theme(axis.title = element_text(size = axis_title_size,family="Helvetica"),
            axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),
            plot.title = element_text(size = plot_title_size,family="Helvetica"),legend.title = element_text(size = axis_title_size,family="Helvetica"),
            legend.text = element_text(size = legend_title_size,family="Helvetica"),
            strip.text.x = element_text(size = legend_text_size,family="Helvetica")) +
      #guides(fill=guide_legend(title=""))+
      ylab('\u0394 Pain Rating\n[LI - HI Exercise Pain Rating]')+xlab(expression(paste('Training Volume (h/week)')))+
      #stat_cor(method = "pearson",alternative = 'two.sided',label.sep = "\n", size = 2,show.legend = F,colour = 'black')+
      ggtitle('')
  )
  
  
  
  
  # ---------------------------------------------
  # Compare fit females and unfit females
  #----------------------------------------------
  
  # PEEP females vs. Hotspin females (SAL)
  
  df_combined_females <- df_combined_heat[df_combined_heat$gender==1,]
  
  f_vs_f <- df_combined_females %>%
    group_by(subject,exercise_intensity,treatment_order,group)%>%
    summarise_at(c('pain_rating'),mean,na.rm = T)
  
  # On Pain Rating:
  lme_model_main_f_vs_f <- lmer(pain_rating ~   exercise_intensity*group+ treatment_order+ (1|subject), data = f_vs_f)
  summary(lme_model_main_f_vs_f)
  confint(lme_model_main)
  
  # Post hoc pairwise comparisons
  emm <- emmeans(lme_model_main_f_vs_f, ~ exercise_intensity * group)
  
  # Pairwise comparisons WITHIN each group
  contrast_within <- contrast(emm, method = "pairwise", by = "group")
  summary(contrast_within, infer = TRUE)
  
  # Pairwise comparisons BETWEEN groups at each level of exercise_intensity
  contrast_between <- contrast(emm, method = "pairwise", by = "exercise_intensity")
  summary(contrast_between, infer = TRUE)
  
  
  # ---- Plotting
  summary_lmer_models_main_sub<- f_vs_f %>%
    group_by(subject,exercise_intensity,group)%>%
    summarise_at(c('pain_rating'),mean,na.rm = T)
  
  summary_lmer_models_main<- f_vs_f %>%
    group_by(exercise_intensity,group)%>%
    summarise_at(c('pain_rating'),mean,na.rm = T)
  
  sum_se <- summarySEwithin(summary_lmer_models_main_sub, 
                            measurevar = "pain_rating", 
                            withinvars  = c("group","exercise_intensity"),
                            idvar = 'subject',
                            na.rm = T)
  
  summary_lmer_models_main$se <- sum_se$se
  summary_lmer_models_main$group<- as.factor(summary_lmer_models_main$group)
  summary_lmer_models_main_sub$group<- as.factor(summary_lmer_models_main_sub$group)
  summary_lmer_models_main$exercise_intensity<- as.factor(summary_lmer_models_main$exercise_intensity)
  summary_lmer_models_main_sub$exercise_intensity<- as.factor(summary_lmer_models_main_sub$exercise_intensity)
  
  (f_vs_f_plot <-ggplot(summary_lmer_models_main,aes(group,pain_rating, fill = exercise_intensity, colour = exercise_intensity))+
      geom_jitter(data = summary_lmer_models_main_sub, aes(x = group,y = pain_rating),shape = 21,alpha = 0.5,size = 1,position = position_jitterdodge(jitter.width= 0.1,dodge.width = 0.7),colour = 'black')+
      geom_bar(stat = 'identity',position = position_dodge(0.7),alpha = 0.8,width = 0.6,colour = 'black')+
      geom_errorbar(aes(group,ymin = pain_rating-se,ymax=pain_rating+se),colour = "black",width  =0.2,size = 0.5,position = position_dodge(0.7))+
      theme_classic()+
      scale_fill_manual(labels = c("Low", "High"),
                        values = c("#005c23", "#3C008E")) +
      scale_color_manual(labels = c("Low", "High"),
                         values = c("#005c23", "#3C008E")) +
      guides(fill=guide_legend(title=""))+
      theme(legend.key.size = unit(0.25, 'cm'),axis.title = element_text(size = axis_title_size,family="Helvetica"),
            axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),
            plot.title = element_text(size = plot_title_size,family="Helvetica"),
            legend.title = element_text(size = legend_title_size),
            legend.text = element_text(size = legend_title_size,family="Helvetica"),
            strip.text.x = element_text(size = legend_text_size,family="Helvetica")) +
      ylab('Pain Rating [VAS]')+xlab('Exercise Intensity')+ylim(0,150)+ggtitle('')+
      scale_x_discrete(limits = c("2","3"),labels = c("Females (Previous Study)","Females (Current Study)"))
  )
  
  # ---- Plotting VAS 70
  f_vs_f_vas70 <- df_combined_females %>%
    group_by(subject,exercise_intensity,treatment_order,group,VAS)%>%
    summarise_at(c('pain_rating'),mean,na.rm = T)
  
  f_vs_f_vas70 <- f_vs_f_vas70[f_vs_f_vas70$VAS==70,]
  
  summary_lmer_models_main_sub<- f_vs_f_vas70 %>%
    group_by(subject,exercise_intensity,group)%>%
    summarise_at(c('pain_rating'),mean,na.rm = T)
  
  summary_lmer_models_main<- f_vs_f_vas70 %>%
    group_by(exercise_intensity,group)%>%
    summarise_at(c('pain_rating'),mean,na.rm = T)
  
  sum_se <- summarySEwithin(f_vs_f_vas70, 
                            measurevar = "pain_rating", 
                            withinvars  = c("group","exercise_intensity"),
                            idvar = 'subject',
                            na.rm = T)
  
  summary_lmer_models_main$se <- sum_se$se
  summary_lmer_models_main$group<- as.factor(summary_lmer_models_main$group)
  summary_lmer_models_main_sub$group<- as.factor(summary_lmer_models_main_sub$group)
  summary_lmer_models_main$exercise_intensity<- as.factor(summary_lmer_models_main$exercise_intensity)
  summary_lmer_models_main_sub$exercise_intensity<- as.factor(summary_lmer_models_main_sub$exercise_intensity)
  
  (f_vs_f_plot <-ggplot(summary_lmer_models_main,aes(group,pain_rating, fill = exercise_intensity, colour = exercise_intensity))+
      geom_jitter(data = summary_lmer_models_main_sub, aes(x = group,y = pain_rating),shape = 21,alpha = 0.5,size = 1,position = position_jitterdodge(jitter.width= 0.1,dodge.width = 0.7),colour = 'black')+
      geom_bar(stat = 'identity',position = position_dodge(0.7),alpha = 0.8,width = 0.6,colour = 'black')+
      geom_errorbar(aes(group,ymin = pain_rating-se,ymax=pain_rating+se),colour = "black",width  =0.2,size = 0.5,position = position_dodge(0.7))+
      theme_classic()+
      scale_fill_manual(labels = c("Low", "High"),
                        values = c("#005c23", "#3C008E")) +
      scale_color_manual(labels = c("Low", "High"),
                         values = c("#005c23", "#3C008E")) +
      guides(fill=guide_legend(title=""))+
      theme(legend.key.size = unit(0.25, 'cm'),axis.title = element_text(size = axis_title_size,family="Helvetica"),
            axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),
            plot.title = element_text(size = plot_title_size,family="Helvetica"),
            legend.title = element_text(size = legend_title_size),
            legend.text = element_text(size = legend_title_size,family="Helvetica"),
            strip.text.x = element_text(size = legend_text_size,family="Helvetica")) +
      ylab('Pain Rating [VAS]')+xlab('Exercise Intensity')+ylim(0,150)+ggtitle('')+
      scale_x_discrete(limits = c("2","3"),labels = c("Females (Previous Study)","Females (Current Study)"))
  )
  
  # On Pain Rating:
  lme_model_main_f_vs_f <- lmer(pain_rating ~   exercise_intensity*group+ treatment_order+ (1|subject), data = f_vs_f_vas70)
  summary(lme_model_main_f_vs_f)
  confint(lme_model_main)
  # Post hoc pairwise comparisons
  emm <- emmeans(lme_model_main_f_vs_f, ~ exercise_intensity * group)
  
  # Pairwise comparisons WITHIN each group
  contrast_within <- contrast(emm, method = "pairwise", by = "group")
  summary(contrast_within, infer = TRUE)
  
  # Pairwise comparisons BETWEEN groups at each level of exercise_intensity
  contrast_between <- contrast(emm, method = "pairwise", by = "exercise_intensity")
  summary(contrast_between, infer = TRUE)
  
  
  
  #-------- Media Split (pwc) fit vs. unfit females
  
  df_combined_females <- df_combined_heat[df_combined_heat$gender==1,]
  
  f_vs_f_med_split <- df_combined_females %>%
    group_by(subject,exercise_intensity,treatment_order,group,pwc)%>%
    summarise_at(c('pain_rating'),mean,na.rm = T)
  
  # Compute median
  median_val <- median(f_vs_f_med_split$pwc, na.rm = TRUE)
  
  # Create new column with median split
  f_vs_f_med_split$group_median_split_pwc <- ifelse(f_vs_f_med_split$pwc <= median_val, "0", "1")
  
  # On Pain Rating:
  lme_model_main_f_vs_f_median_split <- lmer(pain_rating ~   exercise_intensity*group_median_split_pwc+ treatment_order+ (1|subject), data = f_vs_f_med_split)
  summary(lme_model_main_f_vs_f_median_split)
  confint(lme_model_main_f_vs_f_median_split)
  
  # Post hoc pairwise comparisons
  emm <- emmeans(lme_model_main_f_vs_f_median_split, ~ exercise_intensity * group_median_split_pwc)
  
  # Pairwise comparisons WITHIN each group
  contrast_within <- contrast(emm, method = "pairwise", by = "group_median_split_pwc")
  summary(contrast_within, infer = TRUE)
  
  # Pairwise comparisons BETWEEN groups at each level of exercise_intensity
  contrast_between <- contrast(emm, method = "pairwise", by = "exercise_intensity")
  summary(contrast_between, infer = TRUE)
  
  
  # ---- Plotting
  summary_lmer_models_main_sub<- f_vs_f_med_split %>%
    group_by(subject,exercise_intensity,group_median_split_pwc)%>%
    summarise_at(c('pain_rating'),mean,na.rm = T)
  
  summary_lmer_models_main<- f_vs_f_med_split %>%
    group_by(exercise_intensity,group_median_split_pwc)%>%
    summarise_at(c('pain_rating'),mean,na.rm = T)
  
  sum_se <- summarySEwithin(f_vs_f_med_split, 
                            measurevar = "pain_rating", 
                            withinvars  = c("group_median_split_pwc","exercise_intensity"),
                            idvar = 'subject',
                            na.rm = T)
  
  summary_lmer_models_main$se <- sum_se$se
  summary_lmer_models_main$group_median_split_pwc<- as.factor(summary_lmer_models_main$group_median_split_pwc)
  summary_lmer_models_main_sub$group_median_split_pwc<- as.factor(summary_lmer_models_main_sub$group_median_split_pwc)
  summary_lmer_models_main$exercise_intensity<- as.factor(summary_lmer_models_main$exercise_intensity)
  summary_lmer_models_main_sub$exercise_intensity<- as.factor(summary_lmer_models_main_sub$exercise_intensity)
  
  (f_vs_f_plot_median_split <-ggplot(summary_lmer_models_main,aes(group_median_split_pwc,pain_rating, fill = exercise_intensity, colour = exercise_intensity))+
      geom_jitter(data = summary_lmer_models_main_sub, aes(x = group_median_split_pwc,y = pain_rating),shape = 21,alpha = 0.4,size = 0.5,position = position_jitterdodge(jitter.width= 0.1,dodge.width = 0.7),colour = 'black')+
      geom_bar(stat = 'identity',position = position_dodge(0.7),alpha = 0.8,width = 0.6,colour = 'black')+
      geom_errorbar(aes(group_median_split_pwc,ymin = pain_rating-se,ymax=pain_rating+se),colour = "black",width  =0.2,size = 0.5,position = position_dodge(0.7))+
      theme_classic()+
      scale_fill_manual(labels = c("Low", "High"),
                        values = c("#005c23", "#3C008E")) +
      scale_color_manual(labels = c("Low", "High"),
                         values = c("#005c23", "#3C008E")) +
      guides(fill=guide_legend(title=""))+
      theme(legend.key.size = unit(0.25, 'cm'),axis.title = element_text(size = axis_title_size,family="Helvetica"),
            axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),
            plot.title = element_text(size = plot_title_size,family="Helvetica"),
            legend.title = element_text(size = legend_title_size),
            legend.text = element_text(size = legend_title_size,family="Helvetica"),
            strip.text.x = element_text(size = legend_text_size,family="Helvetica")) +
      ylab('Pain Rating [VAS]')+xlab('Exercise Intensity')+ylim(0,100)+ggtitle('')+
      scale_x_discrete(limits = c("0","1"),labels = c(" unfit females","fit Females"))
  )
  
  #--------------------------------------------------------------
  # PEEP females vs. Hotspin females (SAL) vs. PEEP males
  #----------------------------------------------------------------
  
  
  f_vs_f_m <- df_combined_heat %>%
    group_by(subject,exercise_intensity,treatment_order,group,gender,trial,block)%>%
    summarise_at(c('pain_rating'),mean,na.rm = T)

  f_vs_f_m$group <- as.numeric(as.character(f_vs_f_m$group))
  f_vs_f_m$exercise_intensity <- as.factor(f_vs_f_m$exercise_intensity)
  f_vs_f_m$trial <- as.factor(f_vs_f_m$trial)
  f_vs_f_m$block <- as.factor(f_vs_f_m$block)
  
  # On Pain Rating:
  lme_model_main_f_vs_f <- lmer(pain_rating ~   exercise_intensity*group+ treatment_order+(1|subject)+trial + block, data = f_vs_f_m)
  summary(lme_model_main_f_vs_f)
  confint(lme_model_main)
  
  f_vs_f_m$group <- as.factor(f_vs_f_m$group)
  f_vs_f_m$exercise_intensity <- as.factor(f_vs_f_m$exercise_intensity)
  lme_model_main_f_vs_f <- lmer(pain_rating ~   exercise_intensity*group+ treatment_order+ (1|subject) + trial +  block, data = f_vs_f_m)

  # Post hoc pairwise comparisons
  f_vs_f_m$group <-as.factor(f_vs_f_m$group)
  emm <- emmeans(lme_model_main_f_vs_f, ~ exercise_intensity * group)
  
  # Pairwise comparisons WITHIN each group
  contrast_within <- contrast(emm, method = "pairwise", by = "group")
  summary(contrast_within, infer = TRUE)
  
  # Pairwise comparisons BETWEEN groups at each level of exercise_intensity
  contrast_between <- contrast(emm, method = "pairwise", by = "group")
  summary(contrast_between, infer = TRUE)
  
  
  # ---- Plotting
  f_vs_f_m <- df_combined_heat %>%
    group_by(subject,exercise_intensity,treatment_order,group,gender)%>%
    summarise_at(c('pain_rating'),mean,na.rm = T)
  
  f_vs_f_m$group <- as.numeric(as.character(f_vs_f_m$group))
  f_vs_f_m$exercise_intensity <- as.factor(f_vs_f_m$exercise_intensity)

  summary_lmer_models_main_sub<- f_vs_f_m %>%
    group_by(subject,exercise_intensity,group)%>%
    summarise_at(c('pain_rating'),mean,na.rm = T)
  
  summary_lmer_models_main<- f_vs_f_m %>%
    group_by(exercise_intensity,group)%>%
    summarise_at(c('pain_rating'),mean,na.rm = T)
  
  sum_se <- summarySEwithin(summary_lmer_models_main_sub, 
                            measurevar = "pain_rating", 
                            withinvars  = c("group","exercise_intensity"),
                            idvar = 'subject',
                            na.rm = T)
  
  summary_lmer_models_main$se <- sum_se$se
  summary_lmer_models_main$group<- as.factor(summary_lmer_models_main$group)
  summary_lmer_models_main_sub$group<- as.factor(summary_lmer_models_main_sub$group)
  summary_lmer_models_main$exercise_intensity<- as.factor(summary_lmer_models_main$exercise_intensity)
  summary_lmer_models_main_sub$exercise_intensity<- as.factor(summary_lmer_models_main_sub$exercise_intensity)
  
  (f_vs_f_vs_m_plot <-ggplot(summary_lmer_models_main,aes(group,pain_rating, fill = exercise_intensity, colour = exercise_intensity))+
      geom_jitter(data = summary_lmer_models_main_sub, aes(x = group,y = pain_rating),shape = 21,alpha = 0.5,size = 1,position = position_jitterdodge(jitter.width= 0.1,dodge.width = 0.7),colour = 'black')+
      geom_bar(stat = 'identity',position = position_dodge(0.7),alpha = 0.8,width = 0.6,colour = 'black')+
      geom_errorbar(aes(group,ymin = pain_rating-se,ymax=pain_rating+se),colour = "black",width  =0.2,size = 0.5,position = position_dodge(0.7))+
      theme_classic(base_size = 11) +
      scale_fill_manual(labels = c("Low", "High"),
                        values = c("#005c23", "#3C008E")) +
      scale_color_manual(labels = c("Low", "High"),
                         values = c("#005c23", "#3C008E")) +
      guides(fill=guide_legend(title=""))+
      theme(legend.key.size = unit(0.25, 'cm'),axis.title = element_text(size = axis_title_size,family="Helvetica"),
            axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),
            plot.title = element_text(size = plot_title_size,family="Helvetica"),
            legend.title = element_text(size = legend_title_size),
            legend.text = element_text(size = legend_title_size,family="Helvetica"),
            strip.text.x = element_text(size = legend_text_size,family="Helvetica")) +
      ylab('Pain Rating [VAS]')+xlab('')+ylim(0,100)+ggtitle('')+
      scale_x_discrete(limits = c("1","2","3"),labels = c("Males\n(Previous)", "Females\n(Previous)","Females\n(Current)"))+
      ggtitle('Heat (Saline, all stimulus intensities)')
  )
  
  
  # calculat ethe mean of the max ratings across subjects
  deltas_df_f_f_m <- f_vs_f_m %>%
    group_by(exercise_intensity,gender,group) %>%
    summarise_at(c('pain_rating'),mean,na.rm = T)
  
  sub_deltas_df_f_f_m <- f_vs_f_m %>%
    group_by(subject,exercise_intensity,gender,group) %>%
    summarise_at(c('pain_rating'),mean,na.rm = T)
  

  deltas_df_f_f_m_diff <- spread(deltas_df_f_f_m, exercise_intensity, pain_rating)
  deltas_df_f_f_m_diff$diff_ints <- deltas_df_f_f_m_diff$'0' -  deltas_df_f_f_m_diff$'1'  # Low Intensity - High Intensity Data
  
  sub_deltas_df_f_f_m_diff <- spread(sub_deltas_df_f_f_m, exercise_intensity, pain_rating)
  sub_deltas_df_f_f_m_diff$diff_ints <- sub_deltas_df_f_f_m_diff$'0' -  sub_deltas_df_f_f_m_diff$'1'  # Low Intensity - High Intensity Data
  
  # On Pain Rating:
  #lme_model_main_f_vs_f <- lmer(diff_ints ~   group + (1|subject), data = sub_deltas_df_f_f_m_diff)
  #summary(lme_model_main_f_vs_f)
  
  
  sum_se <- summarySE(sub_deltas_df_f_f_m_diff, 
                      measurevar = "diff_ints", 
                      groupvars = c("group",'gender'),
                      na.rm = T)
  
  deltas_df_f_f_m_diff$se <- sum_se$se

  (deltas_param_samples <-ggplot(deltas_df_f_f_m_diff,aes(group,diff_ints))+
      #geom_jitter(data = sub_max_online_ratings_diff,aes(x = VAS_intensity, y = diff_ints), shape = 21,alpha = 0.4,size = 0.5,position = position_jitter(width = 0.1),show.legend = F) +
      geom_bar(stat = 'identity',alpha = 0.8,width = 0.6,colour = 'black',position = position_dodge(0.7))+
      geom_errorbar(aes(group,ymin = diff_ints-se,ymax=diff_ints+se),colour = "black",width  =0.2,size = 0.5)+
      theme_classic(base_size = 11) +
      theme(legend.key.size = unit(0.25, 'cm'),axis.title = element_text(size = axis_title_size,family="Helvetica"),axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),plot.title = element_text(size = plot_title_size,family="Helvetica"),legend.title = element_text(size = axis_title_size,family="Helvetica"),legend.text = element_text(size = legend_title_size,family="Helvetica"),strip.text.x = element_text(size = legend_text_size,family="Helvetica")) +
      guides(fill=guide_legend(title=""))+
      ylab('\u0394 Pain Rating [LI - HI exercise]')+xlab('')+
      scale_x_discrete(limits = c("1","2","3"),labels = c("Males\n(Previous)", "Females\n(Previous)","Females\n(Current)"))+
      #ggtitle('Difference between Max Pain Rating [LI - HI exercise] Whole Stimulus Duration') +
      geom_hline(yintercept = 0)+
      ylim(-10,10)
  )
  
  #================ VAS 70 =======================
  f_vs_f_m <- df_combined_heat %>%
    group_by(subject,exercise_intensity,treatment_order,group,gender,VAS,trial,block)%>%
    summarise_at(c('pain_rating'),mean,na.rm = T)
  
  f_vs_f_m_70 <- f_vs_f_m[f_vs_f_m$VAS==70,]
  f_vs_f_m_70$group <- as.numeric(as.character(f_vs_f_m_70$group))
  f_vs_f_m_70$exercise_intensity <- as.factor(f_vs_f_m_70$exercise_intensity)
  f_vs_f_m_70$VAS <- as.factor(f_vs_f_m_70$VAS)
  f_vs_f_m_70$trial <- as.factor(f_vs_f_m_70$trial)
  f_vs_f_m_70$block <- as.factor(f_vs_f_m_70$block)
  
  # On Pain Rating:
  lme_model_main_f_vs_f <- lmer(pain_rating ~   exercise_intensity*group+ treatment_order+ trial + block +(1|subject), data = f_vs_f_m_70)
  summary(lme_model_main_f_vs_f)
  confint(lme_model_main)
  
  # Post hoc pairwise comparisons
  f_vs_f_m_70$group <- as.factor(f_vs_f_m_70$group)
  f_vs_f_m_70$exercise_intensity <- as.factor(f_vs_f_m_70$exercise_intensity)
  lme_model_main_f_vs_f_70 <- lmer(pain_rating ~   exercise_intensity*group+ treatment_order+ trial + block + (1|subject), data = f_vs_f_m_70)
  
  emm <- emmeans(lme_model_main_f_vs_f_70, ~ exercise_intensity * group)
  
  # Pairwise comparisons WITHIN each group
  contrast_within <- contrast(emm, method = "pairwise", by = "group")
  summary(contrast_within, infer = TRUE)
  
  # Pairwise comparisons BETWEEN groups at each level of exercise_intensity
  contrast_between <- contrast(emm, method = "pairwise", by = "exercise_intensity")
  summary(contrast_between, infer = TRUE)
  
  
  # ---- Plotting
  summary_lmer_models_main_sub<- f_vs_f_m_70 %>%
    group_by(subject,exercise_intensity,group)%>%
    summarise_at(c('pain_rating'),mean,na.rm = T)
  
  summary_lmer_models_main<- f_vs_f_m_70 %>%
    group_by(exercise_intensity,group)%>%
    summarise_at(c('pain_rating'),mean,na.rm = T)
  
  sum_se <- summarySEwithin(summary_lmer_models_main_sub, 
                            measurevar = "pain_rating", 
                            withinvars  = c("group","exercise_intensity"),
                            idvar = 'subject',
                            na.rm = T)
  
  summary_lmer_models_main$se <- sum_se$se
  summary_lmer_models_main$group<- as.factor(summary_lmer_models_main$group)
  summary_lmer_models_main_sub$group<- as.factor(summary_lmer_models_main_sub$group)
  summary_lmer_models_main$exercise_intensity<- as.factor(summary_lmer_models_main$exercise_intensity)
  summary_lmer_models_main_sub$exercise_intensity<- as.factor(summary_lmer_models_main_sub$exercise_intensity)
  
  (f_vs_f_vs_m_plot_vas70 <-ggplot(summary_lmer_models_main,aes(group,pain_rating, fill = exercise_intensity, colour = exercise_intensity))+
      geom_jitter(data = summary_lmer_models_main_sub, aes(x = group,y = pain_rating),shape = 21,alpha = 0.5,size = 1,position = position_jitterdodge(jitter.width= 0.1,dodge.width = 0.7),colour = 'black')+
      geom_bar(stat = 'identity',position = position_dodge(0.7),alpha = 0.8,width = 0.6,colour = 'black')+
      geom_errorbar(aes(group,ymin = pain_rating-se,ymax=pain_rating+se),colour = "black",width  =0.2,size = 0.5,position = position_dodge(0.7))+
      theme_classic(base_size = 11) +
      scale_fill_manual(labels = c("Low", "High"),
                        values = c("#005c23", "#3C008E")) +
      scale_color_manual(labels = c("Low", "High"),
                         values = c("#005c23", "#3C008E")) +
      guides(fill=guide_legend(title=""))+
      theme(legend.key.size = unit(0.25, 'cm'),axis.title = element_text(size = axis_title_size,family="Helvetica"),
            axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),
            plot.title = element_text(size = plot_title_size,family="Helvetica"),
            legend.title = element_text(size = legend_title_size),
            legend.text = element_text(size = legend_title_size,family="Helvetica"),
            strip.text.x = element_text(size = legend_text_size,family="Helvetica")) +
      ylab('Pain Rating [VAS]')+xlab('')+ylim(0,100)+ggtitle('')+
      scale_x_discrete(limits = c("1","2","3"),labels = c("Males\n(Previous)", "Females\n(Previous)","Females\n(Current)"))+
      ggtitle('Heat (Saline, VAS 70)')
  )
  
  # calculate the mean of the max ratings across subjects
  f_vs_f_m <- df_combined_heat %>%
    group_by(subject,exercise_intensity,treatment_order,group,gender,VAS)%>%
    summarise_at(c('pain_rating'),mean,na.rm = T)
  
  f_vs_f_m_70 <- f_vs_f_m[f_vs_f_m$VAS==70,]
  f_vs_f_m_70$group <- as.numeric(as.character(f_vs_f_m_70$group))
  f_vs_f_m_70$exercise_intensity <- as.factor(f_vs_f_m_70$exercise_intensity)
  
  deltas_df_f_f_m_70 <- f_vs_f_m_70 %>%
    group_by(exercise_intensity,gender,group) %>%
    summarise_at(c('pain_rating'),mean,na.rm = T)
  
  sub_deltas_df_f_f_m_70 <- f_vs_f_m_70 %>%
    group_by(subject,exercise_intensity,gender,group) %>%
    summarise_at(c('pain_rating'),mean,na.rm = T)
  
  
  deltas_df_f_f_m_diff <- spread(deltas_df_f_f_m_70, exercise_intensity, pain_rating)
  deltas_df_f_f_m_diff$diff_ints <- deltas_df_f_f_m_diff$'0' -  deltas_df_f_f_m_diff$'1'  # Low Intensity - High Intensity Data
  
  sub_deltas_df_f_f_m_diff <- spread(sub_deltas_df_f_f_m_70, exercise_intensity, pain_rating)
  sub_deltas_df_f_f_m_diff$diff_ints <- sub_deltas_df_f_f_m_diff$'0' -  sub_deltas_df_f_f_m_diff$'1'  # Low Intensity - High Intensity Data
  
  # On Pain Rating:
  #lme_model_main_f_vs_f <- lmer(diff_ints ~   group+ (1|subject), data = sub_deltas_df_f_f_m_diff)
  #summary(lme_model_main_f_vs_f)
  
  
  sum_se <- summarySE(sub_deltas_df_f_f_m_diff, 
                      measurevar = "diff_ints", 
                      groupvars = c("group",'gender'),
                      na.rm = T)
  
  deltas_df_f_f_m_diff$se <- sum_se$se
  
  
  (deltas_param_samples_vas70 <-ggplot(deltas_df_f_f_m_diff,aes(group,diff_ints))+
      #geom_jitter(data = sub_max_online_ratings_diff,aes(x = VAS_intensity, y = diff_ints), shape = 21,alpha = 0.4,size = 0.5,position = position_jitter(width = 0.1),show.legend = F) +
      geom_bar(stat = 'identity',alpha = 0.8,width = 0.6,colour = 'black',position = position_dodge(0.7))+
      geom_errorbar(aes(group,ymin = diff_ints-se,ymax=diff_ints+se),colour = "black",width  =0.2,size = 0.5)+
      theme_classic(base_size = 11) +
      theme(legend.key.size = unit(0.25, 'cm'),axis.title = element_text(size = axis_title_size,family="Helvetica"),axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),plot.title = element_text(size = plot_title_size,family="Helvetica"),legend.title = element_text(size = axis_title_size,family="Helvetica"),legend.text = element_text(size = legend_title_size,family="Helvetica"),strip.text.x = element_text(size = legend_text_size,family="Helvetica")) +
      guides(fill=guide_legend(title=""))+
      ylab('\u0394 Pain Rating [LI - HI exercise]')+xlab('')+
      scale_x_discrete(limits = c("1","2","3"),labels = c("Males\n(Previous)", "Females\n(Previous)","Females\n(Current)"))+
      #ggtitle('Difference between Max Pain Rating [LI - HI exercise] Whole Stimulus Duration') +
      geom_hline(yintercept = 0)+
      ylim(-10,10)
  )
  
  #------------------------- Pressure --------------------------
  
  f_vs_f_m <- df_combined_pressure %>%
    group_by(subject,exercise_intensity,treatment_order,group,gender,trial, block)%>%
    summarise_at(c('pain_rating'),mean,na.rm = T)
  
  f_vs_f_m$group <- as.numeric(as.character(f_vs_f_m$group))
  f_vs_f_m$exercise_intensity <- as.factor(f_vs_f_m$exercise_intensity)
  f_vs_f_m$trial <- as.factor(f_vs_f_m$trial)
  f_vs_f_m$block <- as.factor(f_vs_f_m$block)
  
  # On Pain Rating:
  lme_model_main_f_vs_f <- lmer(pain_rating ~   exercise_intensity*group+ treatment_order+(1|subject)+trial + block, data = f_vs_f_m)
  summary(lme_model_main_f_vs_f)
  confint(lme_model_main)
  
  f_vs_f_m$group <- as.factor(f_vs_f_m$group)
  f_vs_f_m$exercise_intensity <- as.factor(f_vs_f_m$exercise_intensity)
  lme_model_main_f_vs_f <- lmer(pain_rating ~   exercise_intensity*group+ treatment_order+ (1|subject) + trial +  block, data = f_vs_f_m)
  
  # Post hoc pairwise comparisons
  f_vs_f_m$group <-as.factor(f_vs_f_m$group)
  emm <- emmeans(lme_model_main_f_vs_f, ~ exercise_intensity * group)
  
  # Pairwise comparisons WITHIN each group
  contrast_within <- contrast(emm, method = "pairwise", by = "group")
  summary(contrast_within, infer = TRUE)
  
  # Pairwise comparisons BETWEEN groups at each level of exercise_intensity
  contrast_between <- contrast(emm, method = "pairwise", by = "group")
  summary(contrast_between, infer = TRUE)
  
  
  # ---- Plotting
  summary_lmer_models_main_sub<- f_vs_f_m %>%
    group_by(subject,exercise_intensity,group)%>%
    summarise_at(c('pain_rating'),mean,na.rm = T)
  
  summary_lmer_models_main<- f_vs_f_m %>%
    group_by(exercise_intensity,group)%>%
    summarise_at(c('pain_rating'),mean,na.rm = T)
  
  sum_se <- summarySEwithin(f_vs_f_m, 
                            measurevar = "pain_rating", 
                            withinvars  = c("group","exercise_intensity"),
                            idvar = 'subject',
                            na.rm = T)
  
  summary_lmer_models_main$se <- sum_se$se
  summary_lmer_models_main$group<- as.factor(summary_lmer_models_main$group)
  summary_lmer_models_main_sub$group<- as.factor(summary_lmer_models_main_sub$group)
  summary_lmer_models_main$exercise_intensity<- as.factor(summary_lmer_models_main$exercise_intensity)
  summary_lmer_models_main_sub$exercise_intensity<- as.factor(summary_lmer_models_main_sub$exercise_intensity)
  
  (f_vs_f_vs_m_plot_pressure <-ggplot(summary_lmer_models_main,aes(group,pain_rating, fill = exercise_intensity, colour = exercise_intensity))+
      geom_jitter(data = summary_lmer_models_main_sub, aes(x = group,y = pain_rating),shape = 21,alpha = 0.5,size = 1,position = position_jitterdodge(jitter.width= 0.1,dodge.width = 0.7),colour = 'black')+
      geom_bar(stat = 'identity',position = position_dodge(0.7),alpha = 0.8,width = 0.6,colour = 'black')+
      geom_errorbar(aes(group,ymin = pain_rating-se,ymax=pain_rating+se),colour = "black",width  =0.2,size = 0.5,position = position_dodge(0.7))+
      theme_classic()+
      scale_fill_manual(labels = c("Low", "High"),
                        values = c("#005c23", "#3C008E")) +
      scale_color_manual(labels = c("Low", "High"),
                         values = c("#005c23", "#3C008E")) +
      guides(fill=guide_legend(title=""))+
      theme(legend.key.size = unit(0.25, 'cm'),axis.title = element_text(size = axis_title_size,family="Helvetica"),
            axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),
            plot.title = element_text(size = plot_title_size,family="Helvetica"),
            legend.title = element_text(size = legend_title_size),
            legend.text = element_text(size = legend_title_size,family="Helvetica"),
            strip.text.x = element_text(size = legend_text_size,family="Helvetica")) +
      ylab('Pain Rating [VAS]')+xlab('')+ylim(0,100)+ggtitle('')+
      scale_x_discrete(limits = c("1","2","3"),labels = c("Males\n(Previous)", "Females\n(Previous)","Females\n(Current)"))+
      ggtitle('Pressure (Saline, all stimulus intensities)')
  )
  
  
  # calculat ethe mean of the max ratings across subjects
  deltas_df_f_f_m <- f_vs_f_m %>%
    group_by(exercise_intensity,gender,group) %>%
    summarise_at(c('pain_rating'),mean,na.rm = T)
  
  sub_deltas_df_f_f_m <- f_vs_f_m %>%
    group_by(subject,exercise_intensity,gender,group) %>%
    summarise_at(c('pain_rating'),mean,na.rm = T)
  
  
  deltas_df_f_f_m_diff <- spread(deltas_df_f_f_m, exercise_intensity, pain_rating)
  deltas_df_f_f_m_diff$diff_ints <- deltas_df_f_f_m_diff$'0' -  deltas_df_f_f_m_diff$'1'  # Low Intensity - High Intensity Data
  
  sub_deltas_df_f_f_m_diff <- spread(sub_deltas_df_f_f_m, exercise_intensity, pain_rating)
  sub_deltas_df_f_f_m_diff$diff_ints <- sub_deltas_df_f_f_m_diff$'0' -  sub_deltas_df_f_f_m_diff$'1'  # Low Intensity - High Intensity Data
  
  # On Pain Rating:
  #lme_model_main_f_vs_f <- lmer(diff_ints ~   group+ (1|subject), data = sub_deltas_df_f_f_m_diff)
  #summary(lme_model_main_f_vs_f)
  
  
  sum_se <- summarySE(sub_deltas_df_f_f_m_diff, 
                      measurevar = "diff_ints", 
                      groupvars = c("group",'gender'),
                      na.rm = T)
  
  deltas_df_f_f_m_diff$se <- sum_se$se
  
  (deltas_param_samples_pressure <-ggplot(deltas_df_f_f_m_diff,aes(group,diff_ints))+
      #geom_jitter(data = sub_max_online_ratings_diff,aes(x = VAS_intensity, y = diff_ints), shape = 21,alpha = 0.4,size = 0.5,position = position_jitter(width = 0.1),show.legend = F) +
      geom_bar(stat = 'identity',alpha = 0.8,width = 0.6,colour = 'black',position = position_dodge(0.7))+
      geom_errorbar(aes(group,ymin = diff_ints-se,ymax=diff_ints+se),colour = "black",width  =0.2,size = 0.5)+
      theme_classic()+
      theme(legend.key.size = unit(0.25, 'cm'),axis.title = element_text(size = axis_title_size,family="Helvetica"),axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),plot.title = element_text(size = plot_title_size,family="Helvetica"),legend.title = element_text(size = axis_title_size,family="Helvetica"),legend.text = element_text(size = legend_title_size,family="Helvetica"),strip.text.x = element_text(size = legend_text_size,family="Helvetica")) +
      guides(fill=guide_legend(title=""))+
      ylab('\u0394 Pain Rating [LI - HI exercise]')+xlab('')+
      scale_x_discrete(limits = c("1","2","3"),labels = c("Males\n(Previous)", "Females\n(Previous)","Females\n(Current)"))+
      #ggtitle('Difference between Max Pain Rating [LI - HI exercise] Whole Stimulus Duration') +
      geom_hline(yintercept = 0)+
      ylim(-10,10)
  )
  
  #================ VAS 70 =======================
  f_vs_f_m <- df_combined_pressure %>%
    group_by(subject,exercise_intensity,treatment_order,group,gender,VAS,trial,block)%>%
    summarise_at(c('pain_rating'),mean,na.rm = T)
  
  f_vs_f_m <- f_vs_f_m[f_vs_f_m$VAS==70,]
  f_vs_f_m$group <- as.numeric(as.character(f_vs_f_m$group))
  f_vs_f_m$exercise_intensity <- as.factor(f_vs_f_m$exercise_intensity)
  f_vs_f_m$VAS <- as.factor(f_vs_f_m$VAS)
  f_vs_f_m$trial <- as.factor(f_vs_f_m$trial)
  f_vs_f_m$block <- as.factor(f_vs_f_m$block)
  
  # On Pain Rating:
  lme_model_main_f_vs_f <- lmer(pain_rating ~   exercise_intensity*group+ treatment_order+ trial + block +(1|subject), data = f_vs_f_m)
  summary(lme_model_main_f_vs_f)
  confint(lme_model_main)
  
  # Post hoc pairwise comparisons
  f_vs_f_m$group <- as.factor(f_vs_f_m$group)
  f_vs_f_m$exercise_intensity <- as.factor(f_vs_f_m$exercise_intensity)
  lme_model_main_f_vs_f <- lmer(pain_rating ~   exercise_intensity*group+ treatment_order+ trial + block +(1|subject), data = f_vs_f_m)
  
  emm <- emmeans(lme_model_main_f_vs_f, ~ exercise_intensity * group)
  
  # Pairwise comparisons WITHIN each group
  contrast_within <- contrast(emm, method = "pairwise", by = "group")
  summary(contrast_within, infer = TRUE)
  
  # Pairwise comparisons BETWEEN groups at each level of exercise_intensity
  contrast_between <- contrast(emm, method = "pairwise", by = "exercise_intensity")
  summary(contrast_between, infer = TRUE)
  
  
  # ---- Plotting
  summary_lmer_models_main_sub<- f_vs_f_m %>%
    group_by(subject,exercise_intensity,group)%>%
    summarise_at(c('pain_rating'),mean,na.rm = T)
  
  summary_lmer_models_main<- f_vs_f_m %>%
    group_by(exercise_intensity,group)%>%
    summarise_at(c('pain_rating'),mean,na.rm = T)
  
  sum_se <- summarySEwithin(f_vs_f_m, 
                            measurevar = "pain_rating", 
                            withinvars  = c("group","exercise_intensity"),
                            idvar = 'subject',
                            na.rm = T)
  
  summary_lmer_models_main$se <- sum_se$se
  summary_lmer_models_main$group<- as.factor(summary_lmer_models_main$group)
  summary_lmer_models_main_sub$group<- as.factor(summary_lmer_models_main_sub$group)
  summary_lmer_models_main$exercise_intensity<- as.factor(summary_lmer_models_main$exercise_intensity)
  summary_lmer_models_main_sub$exercise_intensity<- as.factor(summary_lmer_models_main_sub$exercise_intensity)
  
  (f_vs_f_vs_m_plot_vas70_pressure <-ggplot(summary_lmer_models_main,aes(group,pain_rating, fill = exercise_intensity, colour = exercise_intensity))+
      geom_jitter(data = summary_lmer_models_main_sub, aes(x = group,y = pain_rating),shape = 21,alpha = 0.5,size = 1,position = position_jitterdodge(jitter.width= 0.1,dodge.width = 0.7),colour = 'black')+
      geom_bar(stat = 'identity',position = position_dodge(0.7),alpha = 0.8,width = 0.6,colour = 'black')+
      geom_errorbar(aes(group,ymin = pain_rating-se,ymax=pain_rating+se),colour = "black",width  =0.2,size = 0.5,position = position_dodge(0.7))+
      theme_classic()+
      scale_fill_manual(labels = c("Low", "High"),
                        values = c("#005c23", "#3C008E")) +
      scale_color_manual(labels = c("Low", "High"),
                         values = c("#005c23", "#3C008E")) +
      guides(fill=guide_legend(title=""))+
      theme(legend.key.size = unit(0.25, 'cm'),axis.title = element_text(size = axis_title_size,family="Helvetica"),
            axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),
            plot.title = element_text(size = plot_title_size,family="Helvetica"),
            legend.title = element_text(size = legend_title_size),
            legend.text = element_text(size = legend_title_size,family="Helvetica"),
            strip.text.x = element_text(size = legend_text_size,family="Helvetica")) +
      ylab('Pain Rating [VAS]')+xlab('')+ylim(0,100)+ggtitle('')+
      scale_x_discrete(limits = c("1","2","3"),labels = c("Males\n(Previous)", "Females\n(Previous)","Females\n(Current)"))+
      ggtitle('Pressure (Saline, VAS 70)')
  )
  
  # calculat ethe mean of the max ratings across subjects
  deltas_df_f_f_m <- f_vs_f_m %>%
    group_by(exercise_intensity, VAS,gender,group) %>%
    summarise_at(c('pain_rating'),mean,na.rm = T)
  
  sub_deltas_df_f_f_m <- f_vs_f_m %>%
    group_by(subject,exercise_intensity, VAS,gender,group) %>%
    summarise_at(c('pain_rating'),mean,na.rm = T)
  
  deltas_df_f_f_m_diff <- spread(deltas_df_f_f_m, exercise_intensity, pain_rating)
  deltas_df_f_f_m_diff$diff_ints <- deltas_df_f_f_m_diff$'0' -  deltas_df_f_f_m_diff$'1'  # Low Intensity - High Intensity Data
  
  sub_deltas_df_f_f_m_diff <- spread(sub_deltas_df_f_f_m, exercise_intensity, pain_rating)
  sub_deltas_df_f_f_m_diff$diff_ints <- sub_deltas_df_f_f_m_diff$'0' -  sub_deltas_df_f_f_m_diff$'1'  # Low Intensity - High Intensity Data
  
  
  sum_se <- summarySE(sub_deltas_df_f_f_m_diff, 
                      measurevar = "diff_ints", 
                      groupvars = c("group",'gender'),
                      na.rm = T)
  
  deltas_df_f_f_m_diff$se <- sum_se$se
  
  (deltas_param_samples_vas70_pressure <-ggplot(deltas_df_f_f_m_diff,aes(group,diff_ints))+
      #geom_jitter(data = sub_max_online_ratings_diff,aes(x = VAS_intensity, y = diff_ints), shape = 21,alpha = 0.4,size = 0.5,position = position_jitter(width = 0.1),show.legend = F) +
      geom_bar(stat = 'identity',alpha = 0.8,width = 0.6,colour = 'black',position = position_dodge(0.7))+
      geom_errorbar(aes(group,ymin = diff_ints-se,ymax=diff_ints+se),colour = "black",width  =0.2,size = 0.5)+
      theme_classic()+
      theme(legend.key.size = unit(0.25, 'cm'),axis.title = element_text(size = axis_title_size,family="Helvetica"),axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),plot.title = element_text(size = plot_title_size,family="Helvetica"),legend.title = element_text(size = axis_title_size,family="Helvetica"),legend.text = element_text(size = legend_title_size,family="Helvetica"),strip.text.x = element_text(size = legend_text_size,family="Helvetica")) +
      guides(fill=guide_legend(title=""))+
      ylab('\u0394 Pain Rating [LI - HI exercise]')+xlab('')+
      scale_x_discrete(limits = c("1","2","3"),labels = c("Males\n(Previous)", "Females\n(Previous)","Females\n(Current)"))+
      #ggtitle('Difference between Max Pain Rating [LI - HI exercise] Whole Stimulus Duration') +
      geom_hline(yintercept = 0)+
      ylim(-10,10)
  )
  
  
  # 
  # # --------------------------------------------------------------------
  # # Mixed Model with PEEP and Hotspin: Fitness Level on Diff Pain Ratings
  # #----------------------------------------------------------------------
  # 
  # df_combined_diff_summary<- df_combined_heat %>%
  #   group_by(subject,exercise_intensity,pwc,gender,treatment_order,group)%>%
  #   summarise_at(c('pain_rating'),mean,na.rm = T)
  # 
  # df_combined_diff <- spread(df_combined_diff_summary, exercise_intensity, pain_rating)
  # df_combined_diff$diff_ints <- df_combined_diff$'0' -  df_combined_diff$'1'  # Low Intensity - High Intensity Data
  # 
  # # statistics
  # lme_model_main_df_combined_diff <- lmer(diff_ints ~   pwc + treatment_order + (1|subject), data = df_combined_diff)
  # summary(lme_model_main_df_combined_diff)
  # confint(lme_model_main_df_combined_diff)
  # 
  # # Test Correlation
  # cor.test(df_combined_diff$diff_ints,df_combined_diff$pwc)
  # 
  # #- Plot
  # (corr_fitness_overall <- 
  #     ggplot(df_combined_diff,aes(x = pwc,y=diff_ints,colour = '#024873',fill = '#024873'))+
  #     geom_point(size = 1,shape = 21,alpha = 0.9,show.legend = F, colour = 'black')+
  #     geom_smooth(method = 'lm',alpha = 0.1,size = 1,se = T,show.legend = F, fullrange = F)+
  #     theme_classic()+
  #     geom_hline(yintercept = 0,colour = 'black',size = 0.5)+
  #     guides(fill=guide_legend(title=""))+
  #     scale_color_manual(labels = c("SAL"),
  #                        values = c("#024873"))+
  #     scale_fill_manual(labels = c("SAL"),
  #                       values = c("#024873"))+
  #     theme(legend.key.size = unit(0.25, 'cm'),axis.title = element_text(size = axis_title_size,family="Helvetica"),
  #           axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),
  #           plot.title = element_text(size = plot_title_size,family="Helvetica"),
  #           legend.title = element_text(size = legend_title_size),
  #           legend.text = element_text(size = legend_title_size,family="Helvetica"),
  #           strip.text.x = element_text(size = legend_text_size,family="Helvetica")) +
  #     ylab('\u0394 Pain Ratings\n[LI - HI Exercise Pain Rating]')+
  #     xlab(expression(paste('FTP [Watt * kg'^-1,']')))+
  #     ggtitle('')
  # )
  # 
  # 
  # # --------------------------------------------------------------------
  # # Mixed Model with PEEP and Hotspin: Fitness Level x Sex on Diff Pain Ratings
  # #----------------------------------------------------------------------
  # 
  #  # statistics
  # lme_model_main_df_combined_diff <- lmer(diff_ints ~   pwc*gender + treatment_order + (1|subject), data = df_combined_diff)
  # summary(lme_model_main_df_combined_diff)
  # confint(lme_model_main_df_combined_diff)
  # 
  # # Post hoc pairwise comparisons
  # emm <- emmeans(lme_model_main_df_combined_diff, ~  pwc*gender)
  # 
  # # Calculate correlations and p-values
  # cor1 <- cor.test(df_combined_diff$diff_ints[df_combined_diff$gender==0], df_combined_diff$pwc[df_combined_diff$gender==0])
  # cor2 <- cor.test(df_combined_diff$diff_ints[df_combined_diff$gender==1], df_combined_diff$pwc[df_combined_diff$gender==1])
  # 
  # # Create text labels for annotations
  # label1 <- sprintf("Males (PEEP): r = %.2f, p = %.3f", cor1$estimate, cor1$p.value)
  # label2 <- sprintf("Females (PEEP + HotSpin): r = %.2f, p = %.3f", cor2$estimate, cor2$p.value)
  # 
  # (corr_fitness_gender <- 
  #     ggplot(df_combined_diff,aes(x = pwc,y=diff_ints,fill = gender, colour = gender,group = gender))+
  #     geom_point(size = 1,shape = 21,alpha = 0.9,show.legend = F, colour = 'black')+
  #     geom_smooth(method = 'lm',alpha = 0.1,size = 1,se = T,show.legend = F, fullrange = F)+
  #     theme_classic()+
  #     scale_color_manual(labels = c("Female (PEEP + HotSpin)", "Male (PEEP)"),
  #                        values = c("#DC8166", "#00c4b4")) +
  #     scale_fill_manual(labels = c("Female (PEEP + HotSpin)", "Male (PEEP)"),
  #                       values = c("#DC8166", "#00c4b4")) +
  #     geom_hline(yintercept = 0,colour = 'black',size = 0.5)+
  #     guides(fill=guide_legend(title=""))+
  #     theme(legend.key.size = unit(0.25, 'cm'),axis.title = element_text(size = axis_title_size,family="Helvetica"),
  #           axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),
  #           plot.title = element_text(size = plot_title_size,family="Helvetica"),
  #           legend.title = element_text(size = legend_title_size),
  #           legend.text = element_text(size = legend_title_size,family="Helvetica"),
  #           strip.text.x = element_text(size = legend_text_size,family="Helvetica")) +
  #     ylab('\u0394 Pain Ratings\n[LI - HI Exercise Pain Rating]')+
  #     xlab(expression(paste('FTP [Watt * kg'^-1,']')))+
  #     annotate("text", x = Inf, y = Inf, hjust = 6.45, vjust = 4.0,
  #              label = label1, colour = "#00c4b4", size = 3.5, parse = FALSE) +
  #     annotate("text", x = Inf, y = Inf, hjust = 4.6, vjust = 6.0,
  #              label = label2, colour = "#DC8166", size = 3.5, parse = FALSE)+
  #     ggtitle('')
  # )
  # 
  # # --------------------------------------------------------------------
  # # Mixed Model with PEEP and Hotspin: Fitness Level x Group on Diff Pain Ratings
  # #----------------------------------------------------------------------
  # 
  # df_combined_diff$group <- as.numeric(as.character(df_combined_diff$group))
  # lme_model_main <- lmer(diff_ints ~   group*pwc+treatment_order+(1|subject), data = df_combined_diff)
  # summary(lme_model_main)
  # 
  # 
  # # Calculate correlations and p-values
  # cor1 <- cor.test(df_combined_diff$diff_ints[df_combined_diff$group==1], df_combined_diff$pwc[df_combined_diff$group==1])
  # cor2 <- cor.test(df_combined_diff$diff_ints[df_combined_diff$group==2], df_combined_diff$pwc[df_combined_diff$group==2])
  # cor3 <- cor.test(df_combined_diff$diff_ints[df_combined_diff$group==3], df_combined_diff$pwc[df_combined_diff$group==3])
  # 
  # # Create text labels for annotations
  # label1 <- sprintf("Females (PEEP): r = %.2f, p = %.3f", cor1$estimate, cor1$p.value)
  # label2 <- sprintf("Males (PEEP): r = %.2f, p = %.3f", cor2$estimate, cor2$p.value)
  # label3 <- sprintf("Females (HotSpin): r = %.2f, p = %.3f", cor3$estimate, cor3$p.value)
  # 
  # # Add annotations in the respective group colors
  # (corr_fitness_group <- 
  #     ggplot(df_combined_diff, aes(x = pwc, y = diff_ints, fill = group, colour = group, group = group)) +
  #     geom_point(size = 1, alpha = .9, shape = 21, show.legend = FALSE, colour = 'black') +
  #     geom_smooth(method = 'lm', alpha = 0.1, size = 1, se = TRUE, show.legend = FALSE, fullrange = FALSE) +
  #     theme_classic() +
  #     scale_color_manual(labels = c("Female (PEEP)", "Male (PEEP)", "Female (HotSpin)"),
  #                        values = c("#00c4b4", "#be00c4", "#fcac00")) +
  #     scale_fill_manual(labels = c("Female (PEEP)", "Male (PEEP)", "Female (HotSpin)"),
  #                       values = c("#00c4b4", "#be00c4", "#fcac00")) +
  #     geom_hline(yintercept = 0, colour = 'black', size = 0.5) +
  #     guides(fill = guide_legend(title = "")) +
  #     theme(
  #       legend.key.size = unit(0.25, 'cm'),
  #       axis.title = element_text(size = axis_title_size, family = "Helvetica"),
  #       axis.text = element_text(size = axis_text_size, colour = 'black', family = "Helvetica"),
  #       plot.title = element_text(size = plot_title_size, family = "Helvetica"),
  #       legend.title = element_text(size = legend_title_size),
  #       legend.text = element_text(size = legend_title_size, family = "Helvetica"),
  #       strip.text.x = element_text(size = legend_text_size, family = "Helvetica")
  #     ) +
  #     ylab('\u0394 Pain Ratings\n[LI - HI Exercise Pain Rating]') +
  #     xlab(expression(paste('FTP [Watt * kg'^-1,']'))) +
  #     ggtitle('') +
  #     annotate("text", x = Inf, y = Inf, hjust = 6.0, vjust = 2.0,
  #              label = label1, colour = "#be00c4", size = 3.5, parse = FALSE) +
  #     annotate("text", x = Inf, y = Inf, hjust = 6.5, vjust = 4.0,
  #              label = label2, colour = "#00c4b4", size = 3.5, parse = FALSE) +
  #   annotate("text", x = Inf, y = Inf, hjust = 5.6, vjust = 6.0,
  #            label = label3, colour = "#fcac00", size = 3.5, parse = FALSE)
  # )
  # 

#----------------------------------------------
# Expectation Exercise no Pain
#----------------------------------------------

# Load in data
data_exercise <- read.csv('C:/Users/user/Desktop/projects/HotSpin/analyses/01_data_cleaning/hotspin_exercise_expectation.csv', sep = ',', header = T)
names(data_exercise) <- c('SubID',"sport.SQ001.","sport.SQ002.","sport.SQ003.","sport.SQ004.", "sport.SQ005.",
                         "sport.SQ006.","sport.SQ007.","sport.SQ008.","sport.SQ009.","sport.SQ010.","sport.SQ011.", 
                         "sport.SQ012.","sport.SQ013.","sport.SQ014." )

# Extract the number from the SubID column and replace the original values
data_exercise$SubID <- as.numeric(gsub("[^0-9]", "", data_exercise$SubID))

# Only take items that regard pain expectation and exercise ()
data_expect_exercise_pain <- data_exercise[,c('SubID',"sport.SQ007.","sport.SQ009.","sport.SQ014." )]

# exclude subject 4
data_expect_exercise_pain <- data_expect_exercise_pain[data_expect_exercise_pain$SubID != 4,]

# Convert wide to long
data_exercise_pain_long <-  gather(data_expect_exercise_pain, question, rating, sport.SQ007.:sport.SQ014.)

# Center the ratings around 0
data_exercise_pain_long$rating <- data_exercise_pain_long$rating - 4

# Summarise over everything
data_exercise_pain_long_sum_sub<- data_exercise_pain_long %>%
  group_by(SubID,question)%>%
  summarise_at(c('rating'),mean,na.rm = T)

data_exercise_pain_long_sum<- data_exercise_pain_long %>%
  group_by(question)%>%
  summarise_at(c('rating'),mean,na.rm = T)


sum_se <- summarySEwithin(data_exercise_pain_long, 
                          measurevar = "rating", 
                          withinvars  = c("question"),
                          idvar = 'SubID',
                          na.rm = T)


data_exercise_pain_long_sum$se <- sum_se$se 
data_exercise_pain_long_sum$question <- as.factor(data_exercise_pain_long_sum$question)
data_exercise_pain_long_sum_sub$question <- as.factor(data_exercise_pain_long_sum_sub$question)

#---------- Statistics
# Perform t-tests 
t_test_sport_SQ007 <- t.test(data_exercise_pain_long$rating[data_exercise_pain_long$question == "sport.SQ007."], mu = 0)
print(t_test_sport_SQ007)

t_test_sport_SQ009 <- t.test(data_exercise_pain_long$rating[data_exercise_pain_long$question == "sport.SQ009."], mu = 0)
print(t_test_sport_SQ009)

t_test_sport_SQ014 <- t.test(data_exercise_pain_long$rating[data_exercise_pain_long$question == "sport.SQ014."], mu = 0)
print(t_test_sport_SQ014)

#---------- Plotting
pain_exp <- ggplot(data_exercise_pain_long_sum_sub, aes(question, rating)) +
  geom_jitter(width = 0.1, shape = 21, alpha = .5, size = 0.75, color = 'darkblue', fill = 'darkblue') +
  geom_boxplot(alpha = 0.5, width = 0.3, fill = '#1C02C7') +
  geom_errorbar(data = data_exercise_pain_long_sum, aes(x = question, ymin = rating - se, ymax = rating + se), position = position_dodge(5), width = .1, size = 0.75) +
  geom_point(data = data_exercise_pain_long_sum, aes(question, rating), shape = 21, colour = 'black', fill = 'black', size = 1) +
  ylim(-3, 3) +
  scale_x_discrete(limits = c("sport.SQ007.", "sport.SQ009.", "sport.SQ014."),
                   labels = c("Joint Pain", "Muscle Pain", "Whole-body Pain")) +
  scale_y_discrete(limits = c(-3, -2, -1, 0, 1, 2, 3), labels = c('greatly reduce (-3)', "-2", "-1", "0", "1", "2", "greatly increase (3)")) +
  xlab('Pain Type') + ylab('Expectation') +
  ggtitle('Expectation of acute exercise on...') +
  theme_classic() +
  geom_hline(yintercept = 0, size = 0.5) +
  theme(axis.title = element_text(size = axis_title_size, family = "Helvetica"),
        axis.text = element_text(size = axis_text_size, colour = 'black', family = "Helvetica"),
        plot.title = element_text(size = plot_title_size, family = "Helvetica"),
        legend.title = element_text(size = axis_title_size, family = "Helvetica"),
        legend.text = element_text(size = legend_text_size, family = "Helvetica"),
        strip.text.x = element_text(size = legend_text_size, family = "Helvetica")) +
  annotate("text", x = 1, y = 3, label = ifelse(t_test_sport_SQ007$p.value < 0.05, "*", "n.s."), size = 3, hjust = 0.5) +
  annotate("text", x = 2, y = 3, label = ifelse(t_test_sport_SQ009$p.value < 0.05, "*", "n.s."), size = 3, hjust = 0.5) +
  annotate("text", x = 3, y = 3, label = ifelse(t_test_sport_SQ014$p.value < 0.05, "*", "n.s."), size = 3, hjust = 0.5)

pain_exp

#----------------------------------------
#- Changes in mood (pre post) based on POMS
#----------------------------------

# Load in mood data
hotspin_mood <- read.csv('C:/Users/user/Desktop/projects/HotSpin/analyses/01_data_cleaning/hotspin_mood_data.csv', sep = ',', header = T)

# Remove POMS1 from question
hotspin_mood$question <- gsub("POMS1", "", hotspin_mood$question)

# Extract the number of the question from the 'question' column and convert it to numeric
hotspin_mood$question <- as.numeric(gsub("[^0-9]", "", hotspin_mood$question))

# Create a new column 'question_group' based on the question number
hotspin_mood <- hotspin_mood %>%
  mutate(dimension = case_when(
    question >= 1 & question <= 14 ~ 1,
    question >= 15 & question <= 21 ~ 2,
    question >= 22 & question <= 28 ~ 3,
    question >= 29 & question <= 35 ~ 4
  ))

# exclude subject 4
hotspin_mood <- hotspin_mood[hotspin_mood$subid != 4,]

# Calculate summary statistics for mood data
summary_mood <- hotspin_mood %>%
  group_by(pre_post,dimension) %>%
  summarise_at(c('score'), mean, na.rm = TRUE)

summary_mood_sub <- hotspin_mood %>%
  group_by(subid,pre_post,dimension) %>%
  summarise_at(c('score'), mean, na.rm = TRUE)


# Calculate standard error for mood data
sum_se_mood <- summarySEwithin(hotspin_mood, 
                               measurevar = "score", 
                               withinvars = c("pre_post", "dimension"),
                               idvar = 'subid',
                               na.rm = TRUE)

summary_mood$se <- sum_se_mood$se

summary_mood$pre_post <- as.factor(summary_mood$pre_post)
summary_mood_sub$pre_post <- as.factor(summary_mood_sub$pre_post)
summary_mood$dimension <- as.factor(summary_mood$dimension)
summary_mood_sub$dimension <- as.factor(summary_mood_sub$dimension)

#------------ Statistics
# Perform paired samples t-tests for each dimension between pre and post ratings
t_tests <- summary_mood_sub %>%
  pivot_wider(
    id_cols = c(subid, dimension),
    names_from = pre_post,
    values_from = score
  ) %>%
  split(.$dimension) %>%
  map(~ t.test(.x$`1`, .x$`2`, paired = TRUE))

# Print the results
t_tests[["1"]]  # Dejection
t_tests[["2"]]  # Fatigue
t_tests[["3"]]  # Discontent
t_tests[["4"]]  # Drive


# Plot bar plots comparing pre and post ratings for each dimension
mood_plot <- ggplot(summary_mood, aes(x = dimension, y = score, fill = pre_post)) +
  geom_bar(stat = "identity", position = position_dodge(), alpha = 0.7, color = 'black') +
  geom_errorbar(aes(ymin = score - se, ymax = score + se), 
                colour = "black", width = 0.2, size = 0.5, position = position_dodge(0.9)) +
  geom_jitter(data = summary_mood_sub, aes(x = dimension, y = score, fill = pre_post), shape = 21, alpha = 0.4, size = 0.5, position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.9), show.legend = F) +
  scale_fill_manual(values = c("#1C02C7", "#C75302"), labels = c("Pre", "Post")) +
  theme_classic() +
  theme(legend.key.size = unit(0.25, 'cm'), axis.title = element_text(size = axis_title_size, family = "Helvetica"), axis.text = element_text(size = axis_text_size, colour = 'black', family = "Helvetica"), plot.title = element_text(size = plot_title_size, family = "Helvetica"), legend.title = element_blank(), legend.text = element_text(size = legend_text_size, family = "Helvetica"), strip.text.x = element_text(size = legend_text_size, family = "Helvetica")) +
  xlab('Dimension (POMS)') +
  ylab('Rating') +
  ggtitle('Pre and Post Mood Ratings (by dimension)') +
  scale_x_discrete(labels = c('1' = 'Dejection', '2' = 'Fatigue', '3' = 'Discontent', '4' = 'Drive')) +
  scale_y_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6), labels = c('not at all (0)', '1', '2', '3', '4', '5', 'very strongly (6)')) +
  annotate("text", x = 1, y = max(summary_mood_sub$score) + 0.5, label = ifelse(t_tests[["1"]]$p.value < 0.05, "**", "n.s."), size = 3, hjust = 0.5) +
  annotate("text", x = 2, y = max(summary_mood_sub$score) + 0.5, label = ifelse(t_tests[["2"]]$p.value < 0.05, "*", "n.s."), size = 3, hjust = 0.5) +
  annotate("text", x = 3, y = max(summary_mood_sub$score) + 0.5, label = ifelse(t_tests[["3"]]$p.value < 0.05, "*", "n.s."), size = 3, hjust = 0.5) +
  annotate("text", x = 4, y = max(summary_mood_sub$score) + 0.5, label = ifelse(t_tests[["4"]]$p.value < 0.05, "*", "n.s."), size = 3, hjust = 0.5) +
  geom_segment(aes(x = 0.8, xend = 1.2, y = max(summary_mood_sub$score) + 0.3, yend = max(summary_mood_sub$score) + 0.3), size = 0.5) +
  geom_segment(aes(x = 1.8, xend = 2.2, y = max(summary_mood_sub$score) + 0.3, yend = max(summary_mood_sub$score) + 0.3), size = 0.5) +
  geom_segment(aes(x = 2.8, xend = 3.2, y = max(summary_mood_sub$score) + 0.3, yend = max(summary_mood_sub$score) + 0.3), size = 0.5) +
  geom_segment(aes(x = 3.8, xend = 4.2, y = max(summary_mood_sub$score) + 0.3, yend = max(summary_mood_sub$score) + 0.3), size = 0.5)

mood_plot

# # --------------------------------------------------------------------
# # Mixed Model with PEEP and Hotspin: Fitness Level on Diff Pain Ratings
# #----------------------------------------------------------------------
# training_hours_combined$subject <- training_hours_combined$SubID
# training_hours_combined$subject <- as.factor(training_hours_combined$subject)
# 
# df_combined_heat <- left_join(df_combined_heat, training_hours_combined %>% select(subject,hours_week), by = "subject")
# 
# df_combined_diff_summary<- df_combined_heat %>%
#   group_by(subject,exercise_intensity,hours_week,gender,treatment_order,group,VAS)%>%
#   summarise_at(c('pain_rating'),mean,na.rm = T)
# 
# df_combined_diff_summary <- df_combined_diff_summary[df_combined_diff_summary$VAS == 70,]
# 
# df_combined_diff <- spread(df_combined_diff_summary, exercise_intensity, pain_rating)
# df_combined_diff$diff_ints <- df_combined_diff$'0' -  df_combined_diff$'1'  # Low Intensity - High Intensity Data
# 
# df_combined_diff$subject <- as.numeric(as.character(df_combined_diff$subject))
# 
# # statistics
# lme_model_main_df_combined_diff <- lmrob(diff_ints ~   hours_week + treatment_order, data = df_combined_diff)
# summary(lme_model_main_df_combined_diff)
# confint(lme_model_main_df_combined_diff)
# 
# # Test Correlation
# cor.test(df_combined_diff$diff_ints,df_combined_diff$hours_week)
# 
# #- Plot
# (corr_fitness_overall <- 
#     ggplot(df_combined_diff,aes(x = hours_week,y=diff_ints,colour = '#024873',fill = '#024873'))+
#     geom_point(size = 1,shape = 21,alpha = 0.9,show.legend = F, colour = 'black')+
#     geom_smooth(method = 'lm',alpha = 0.1,size = 1,se = T,show.legend = F, fullrange = F)+
#     theme_classic()+
#     geom_hline(yintercept = 0,colour = 'black',size = 0.5)+
#     guides(fill=guide_legend(title=""))+
#     scale_color_manual(labels = c("SAL"),
#                        values = c("#024873"))+
#     scale_fill_manual(labels = c("SAL"),
#                       values = c("#024873"))+
#     theme(legend.key.size = unit(0.25, 'cm'),axis.title = element_text(size = axis_title_size,family="Helvetica"),
#           axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),
#           plot.title = element_text(size = plot_title_size,family="Helvetica"),
#           legend.title = element_text(size = legend_title_size),
#           legend.text = element_text(size = legend_title_size,family="Helvetica"),
#           strip.text.x = element_text(size = legend_text_size,family="Helvetica")) +
#     ylab('\u0394 Pain Ratings\n[LI - HI Exercise Pain Rating]')+
#     xlab(expression(paste('Training Volume [h/week]'))) +
#     ggtitle('')
# )
# 
# 
# # --------------------------------------------------------------------
# # Mixed Model with PEEP and Hotspin: Training Volume x Sex on Diff Pain Ratings
# #----------------------------------------------------------------------
# library(robustbase)
# 
# # statistics
# lme_model_main_df_combined_diff <- lmrob(diff_ints ~   hours_week*gender + treatment_order, data = df_combined_diff)
# summary(lme_model_main_df_combined_diff)
# confint(lme_model_main_df_combined_diff)
# 
# # Post hoc pairwise comparisons
# emm <- emmeans(lme_model_main_df_combined_diff, ~  hours_week*gender)
# 
# # Calculate correlations and p-values
# cor1 <- cor.test(df_combined_diff$diff_ints[df_combined_diff$gender==0], df_combined_diff$hours_week[df_combined_diff$gender==0])
# cor2 <- cor.test(df_combined_diff$diff_ints[df_combined_diff$gender==1], df_combined_diff$hours_week[df_combined_diff$gender==1])
# 
# # Create text labels for annotations
# label1 <- sprintf("Males (PEEP): r = %.2f, p = %.3f", cor1$estimate, cor1$p.value)
# label2 <- sprintf("Females (PEEP + HotSpin): r = %.2f, p = %.3f", cor2$estimate, cor2$p.value)
# 
# (corr_fitness_gender <- 
#     ggplot(df_combined_diff,aes(x = hours_week,y=diff_ints,fill = gender, colour = gender,group = gender))+
#     geom_point(size = 1,shape = 21,alpha = 0.9,show.legend = F, colour = 'black')+
#     geom_smooth(method = 'lm',alpha = 0.1,size = 1,se = T,show.legend = F, fullrange = F)+
#     theme_classic()+
#     scale_color_manual(labels = c("Female (PEEP + HotSpin)", "Male (PEEP)"),
#                        values = c("#DC8166", "#00c4b4")) +
#     scale_fill_manual(labels = c("Female (PEEP + HotSpin)", "Male (PEEP)"),
#                       values = c("#DC8166", "#00c4b4")) +
#     geom_hline(yintercept = 0,colour = 'black',size = 0.5)+
#     guides(fill=guide_legend(title=""))+
#     theme(legend.key.size = unit(0.25, 'cm'),axis.title = element_text(size = axis_title_size,family="Helvetica"),
#           axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),
#           plot.title = element_text(size = plot_title_size,family="Helvetica"),
#           legend.title = element_text(size = legend_title_size),
#           legend.text = element_text(size = legend_title_size,family="Helvetica"),
#           strip.text.x = element_text(size = legend_text_size,family="Helvetica")) +
#     ylab('\u0394 Pain Ratings\n[LI - HI Exercise Pain Rating]')+
#     xlab(expression(paste('Training Volume [h/week]'))) +
#     annotate("text", x = Inf, y = Inf, hjust = 6.45, vjust = 4.0,
#              label = label1, colour = "#00c4b4", size = 3.5, parse = FALSE) +
#     annotate("text", x = Inf, y = Inf, hjust = 4.6, vjust = 6.0,
#              label = label2, colour = "#DC8166", size = 3.5, parse = FALSE)+
#     ggtitle('')
# )
# 
# # --------------------------------------------------------------------
# # Mixed Model with PEEP and Hotspin: Fitness Level x Group on Diff Pain Ratings
# #----------------------------------------------------------------------
# 
# lme_model_main <- lmrob(diff_ints ~   group*hours_week+treatment_order, data = df_combined_diff)
# summary(lme_model_main)
# 
# 
# # Calculate correlations and p-values
# cor1 <- cor.test(df_combined_diff$diff_ints[df_combined_diff$group==1], df_combined_diff$hours_week[df_combined_diff$group==1])
# cor2 <- cor.test(df_combined_diff$diff_ints[df_combined_diff$group==2], df_combined_diff$hours_week[df_combined_diff$group==2])
# cor3 <- cor.test(df_combined_diff$diff_ints[df_combined_diff$group==3], df_combined_diff$hours_week[df_combined_diff$group==3])
# 
# # Create text labels for annotations
# label1 <- sprintf("Females (PEEP): r = %.2f, p = %.3f", cor1$estimate, cor1$p.value)
# label2 <- sprintf("Males (PEEP): r = %.2f, p = %.3f", cor2$estimate, cor2$p.value)
# label3 <- sprintf("Females (HotSpin): r = %.2f, p = %.3f", cor3$estimate, cor3$p.value)
# 
# # Add annotations in the respective group colors
# (corr_fitness_group <- 
#     ggplot(df_combined_diff, aes(x = hours_week, y = diff_ints, fill = group, colour = group, group = group)) +
#     geom_point(size = 1, alpha = .9, shape = 21, show.legend = TRUE, colour = 'black') +
#     geom_smooth(method = 'lm', alpha = 0.1, size = 1, se = TRUE, show.legend = FALSE, fullrange = FALSE) +
#     theme_classic() +
#     scale_color_manual(labels = c("Female (PEEP)", "Male (PEEP)", "Female (HotSpin)"),
#                        values = c("#00c4b4", "#be00c4", "#fcac00")) +
#     scale_fill_manual(labels = c("Female (PEEP)", "Male (PEEP)", "Female (HotSpin)"),
#                       values = c("#00c4b4", "#be00c4", "#fcac00")) +
#     geom_hline(yintercept = 0, colour = 'black', size = 0.5) +
#     #guides(fill = guide_legend(title = "")) +
#     theme(
#       legend.key.size = unit(0.25, 'cm'),
#       axis.title = element_text(size = axis_title_size, family = "Helvetica"),
#       axis.text = element_text(size = axis_text_size, colour = 'black', family = "Helvetica"),
#       plot.title = element_text(size = plot_title_size, family = "Helvetica"),
#       legend.title = element_text(size = legend_title_size),
#       legend.text = element_text(size = legend_title_size, family = "Helvetica"),
#       strip.text.x = element_text(size = legend_text_size, family = "Helvetica")
#     ) +
#     ylab('\u0394 Pain Ratings\n[LI - HI Exercise Pain Rating]') +
#     xlab(expression(paste('Training Volume [h/week]'))) +
#     ggtitle('') +
#     annotate("text", x = Inf, y = Inf, hjust = 6.0, vjust = 2.0,
#              label = label1, colour = "#be00c4", size = 3.5, parse = FALSE) +
#     annotate("text", x = Inf, y = Inf, hjust = 6.5, vjust = 4.0,
#              label = label2, colour = "#00c4b4", size = 3.5, parse = FALSE) +
#     annotate("text", x = Inf, y = Inf, hjust = 5.6, vjust = 6.0,
#              label = label3, colour = "#fcac00", size = 3.5, parse = FALSE)
# )
# 
# ggarrange(corr_fitness_overall+rremove('legend'), corr_fitness_group  +rremove('legend'),corr_fitness_gender +rremove('legend'),ncol = 3, nrow = 1, labels = c("A", "B","C"), font.label = list(size = 11))
# 
# #ggsave(paste(save_path,'correlation_training_volume_samples.svg'), width = 17, height = 7, units = "cm")
# #ggsave(paste(save_path,'correlation_training_volume_samples.png'), width = 17, height = 7, units = "cm")
# 
# # -- Correlation Fitenss and Training volumne
# summary_lmer_models_main<- df_combined_heat %>%
#   group_by(subject,pwc,hours_week,gender,group)%>%
#   summarise_at(c('pain_rating'),mean,na.rm = T)
# 
# cor.test(summary_lmer_models_main$pwc,summary_lmer_models_main$hours_week)
# 
# (corr_fitness_travol <- ggplot(summary_lmer_models_main,aes(x = pwc,y=hours_week,colour = "#024873",fill = "#024873"))+
#     geom_point(size = 0.5,aes(colour = exercise_intensity),shape = 21,colour = 'black',alpha = 0.6,show.legend = F)+
#     geom_smooth(method = 'lm',alpha = 0.1,size = 1,se = T,show.legend = F)+
#     theme_classic()+
#     scale_color_manual(labels = c("SAL"),
#                        values = c("#024873")) +
#     scale_fill_manual(labels = c("SAL"),
#                       values = c("#024873"))+
#     geom_hline(yintercept = 0,colour = 'black',size = 0.5)+
#     theme(axis.title = element_text(size = axis_title_size,family="Helvetica"),
#           axis.text = element_text(size = axis_text_size,colour = 'black',family="Helvetica"),
#           plot.title = element_text(size = plot_title_size,family="Helvetica"),legend.title = element_text(size = axis_title_size,family="Helvetica"),
#           legend.text = element_text(size = legend_title_size,family="Helvetica"),
#           strip.text.x = element_text(size = legend_text_size,family="Helvetica")) +
#   #guides(fill=guide_legend(title=""))+
#   ylab('Training Volume (h/week)')+xlab(expression(paste('FTP [Watt * kg'^-1,']')))+
#   #stat_cor(method = "pearson",alternative = 'two.sided',label.sep = "\n", size = 2,show.legend = F,colour = 'black')+
#   ggtitle('FTP and Training Volume'))
# 
# ggarrange(corr_fitness_travol+rremove('legend'),ncol = 1, nrow = 1, labels = c("A"), font.label = list(size = 11))
# 
# #ggsave(paste(save_path,'correlation_fitness_training_volume_samples.svg'), width = 5, height = 5, units = "cm")
# #ggsave(paste(save_path,'correlation_fitness_training_volume_samples.png'), width = 5, height = 5, units = "cm")
# 

##################################################################################
##       Arrange Graphs
##################################################################################


# -------------- Sample comparison (MAIN)
ggarrange(
  raincloud_pwc_group,                 
  raincloud_training_h_w_group+ rremove("y.text") + rremove("y.ticks") + rremove("legend"),  # plot 2
  labels = c("A", "B"),
  font.label = list(size = 11),align = 'h', widths = c(1.4, 0.8) ,
  ncol = 2
)

#ggsave(paste(save_path,'raincloud_group_comparison.png'), width = 15, height = 7, units = "cm")
#ggsave(paste(save_path,'raincloud_group_comparison.svg'), width = 15, height = 7, units = "cm")



# -------------- Exercise Intervention
ggarrange(power_cycling,rel_power_cycling,hr_cycling, borg_rating,ncol = 4, nrow = 1,labels = c("A","B","C","D"),font.label = list(size = 12.5),widths = c(0.75,0.75,0.75,0.75))

#ggsave(paste(save_path,'exercise_params.svg'), width = 20, height = 7, units = "cm")
#ggsave(paste(save_path,'exercise_params.png'), width = 20, height = 7, units = "cm")




# ---------------Max ratings 
  ggarrange(max_ratings_whole_across_VAS+rremove('legend'), max_ratings_whole+rremove('legend'), deltas_param_max, labels = c("A", "B","C"),
            font.label = list(size = 11), align = 'hv', ncol = 3, widths = c(0.9, 1.1, 0.9))


#ggsave(paste(save_path,'overview_results_pain_ratings.svg'), width = 17, height = 12, units = "cm")
#ggsave(paste(save_path,'overview_results_pain_ratings.png'), width = 17, height = 12, units = "cm")

  
  #---------- Comparison PEEP Sample 
  
  # Heat
  ggarrange(
    ggarrange(f_vs_f_vs_m_plot+rremove('legend'), deltas_param_samples, labels = c("A", "B"), 
              font.label = list(size = 11), align = 'hv', ncol = 2, widths = c(1.1, 1)),
    ggarrange(f_vs_f_vs_m_plot_vas70+rremove('legend'), deltas_param_samples_vas70, labels = c("C", "D"), 
              font.label = list(size = 11), align = 'hv', ncol = 2, widths = c(1.1, 1)),nrow = 2
  )
  
#ggsave(paste(save_path,'peep_hotspin_overview_results_heat_pain_ratings.svg'), width = 13, height = 13, units = "cm", device = "svg")
#ggsave(paste(save_path,'peep_hotspin_overview_results_heat_pain_ratings.png'), width = 13, height = 13, units = "cm", device = "png")
  
  
  
  # Correlation FTP , Sex and Pain (HEAT)
  ggarrange(peep_hotspin_heat_saline_fitness,peep_hotspin_heat_saline_fitness_gender,peep_hotspin_heat_saline_fitness_group,ncol = 3,labels = c("A", "B", "C"), font.label = list(size = 11))
#  ggsave(paste(save_path,'comp_peep_hotspin_fitness_sex.svg'), width = 15, height = 6, units = "cm")
#  ggsave(paste(save_path,'comp_peep_hotspin_fitness_sex.png'), width = 15, height = 6, units = "cm")
  
#=======================================================

# SUPPLEMENTS

#=====================================================

# -------------- Calibration and online ratings
ggarrange(ggarrange(
  raincloud_pressure +rremove('legend'),raincloud_heat+rremove('legend'),labels = c("A","B"),font.label = list(size = 11),
  ncol = 2, nrow = 1),
  ggarrange(ratings_pressure,ratings_heat+rremove('ylab'),labels = c("C","D"),font.label = list(size = 11),widths = c(1,1)),
  ncol = 1, nrow = 2,align = 'hv')

#ggsave(paste(save_path,'heat_pressure_calib_online_ratings_supplements.svg'), width = 10, height = 10, units = "cm")
#ggsave(paste(save_path,'heat_pressure_calib_online_ratings_supplements.png'), width = 10, height = 10, units = "cm")

#---------------- Menstural cycle pahse
ggarrange(cycle_phase_distribution,pain_exp, mood_plot,labels = c("A","B","C"), font.label = list(size = 11), ncol = 3, align = 'hv',widths = c(0.8, 1,1))

#ggsave(paste(save_path,'cyclePhase_expectation_mood_supplements.png'), width = 15, height = 7, units = "cm")
#ggsave(paste(save_path,'cyclePhase_expectation_mood_supplements.svg'), width = 15, height = 7, units = "cm")

# -------------- Sample comparison (SUPPLEMENTS)
ggarrange(
  raincloud_weight_group+ rremove("legend"), 
  raincloud_height_group+ rremove("legend") + rremove("y.text") + rremove("y.ticks") , # plot 3
  raincloud_age_group + rremove("y.text") + rremove("y.ticks") + rremove("legend"),           # plot 4
  labels = c("A", "B", "C" ),
  font.label = list(size = 12),align = 'h', widths = c(1.2, 0.8, 0.8) ,
  ncol = 3
)

#ggsave(paste(save_path,'raincloud_group_comparison_supplements.png'), width = 20, height = 7, units = "cm")
#ggsave(paste(save_path,'raincloud_group_comparison_supplements.svg'), width = 20, height = 7, units = "cm")



#-------------------- Online Ratings Exercise
ggarrange(pressure_online_rating_across_VAS+rremove('legend'),heat_online_rating_across_VAS+rremove('legend'),labels = c("A","B"),font.label = list(size = 11),widths = c(1,1))

#ggsave(paste(save_path,'overview_results_online_ratings_supplements.svg'),  width = 10, height = 7, units = "cm")
#ggsave(paste(save_path,'overview_results_online_ratings_supplements.png'), width = 10, height = 7, units = "cm")



#-------------------- Online Ratings Exercise (VAS)

ggarrange(online_rating_VAS,ncol = 1,common.legend = T,legend = 'none',align = 'hv',font.label = list(size = 11))

#ggsave(paste(save_path,'VAS_results_online_ratings.svg'), width = 11, height = 7, units = "cm")
#ggsave(paste(save_path,'VAS_results_online_ratings.png'), width = 11, height = 7, units = "cm")


#--------------------
ggarrange(
  ggarrange(pressure_online_rating_across_VAS+rremove('legend'),heat_online_rating_across_VAS+rremove('legend'),labels = c("A","B"),font.label = list(size = 11),widths = c(1,1)),
  ggarrange(max_ratings_whole+rremove('legend'), deltas_param_max, labels = c("C", "D"),
          font.label = list(size = 11), align = 'hv', ncol = 2, widths = c(1.1, 0.9)),nrow = 2
)
#ggsave(paste(save_path,'overview_results_online_max_ratings.svg'), width = 15, height = 12, units = "cm")




#---------------- Menstural cycle pahse
ggarrange(cycle_phase_distribution,pain_exp, mood_plot,labels = c("A","B","C"), font.label = list(size = 11), ncol = 3, align = 'hv',widths = c(0.8, 1,1))

#ggsave(paste(save_path,'raincloud_pwc_group.svg'), width = 8, height = 10, units = "cm")
#ggsave(paste(save_path,'expectation_cycle_phase_mood.png'), width = 24, height = 8, units = "cm")


# ------------ FTP and Difference Pain ratings
#ggarrange(correlation_plot+rremove('legend'), correlation_plot_plateau+rremove('legend'), correlation_plot_max+rremove('legend'), ncol = 1, nrow = 3, labels = c("A", "B", "C"), font.label = list(size = 11))

#ggsave(paste(save_path,'correlation_FTP_ratings.svg'), width = 13, height = 20, units = "cm")
#ggsave(paste(save_path,'correlation_FTP_ratings.png'), width = 13, height = 20, units = "cm")

# ------------ FTP and Difference Pain ratings at VAS 70
ggarrange(correlation_plot_VAS70+rremove('legend'), correlation_plot_plateau_VAS70  +rremove('legend'), correlation_plot_max_VAS70+rremove('legend'), ncol = 1, nrow = 3, labels = c("A", "B", "C"), font.label = list(size = 11))

#ggsave(paste(save_path,'correlation_FTP_ratings_VAS70.svg'), width = 13, height = 20, units = "cm")
#ggsave(paste(save_path,'correlation_FTP_ratings_VAS70.png'), width = 13, height = 20, units = "cm")


# Pressure

ggarrange(
  ggarrange(f_vs_f_vs_m_plot_pressure+rremove('legend'), deltas_param_samples_pressure, labels = c("A", "B"), 
            font.label = list(size = 11), align = 'hv', ncol = 2, widths = c(1.1, 1)),
  ggarrange(f_vs_f_vs_m_plot_vas70_pressure+rremove('legend'), deltas_param_samples_vas70_pressure, labels = c("C", "D"), 
            font.label = list(size = 11), align = 'hv', ncol = 2, widths = c(1.1, 1)),nrow = 2
)

#ggsave(paste(save_path,'peep_hotspin_overview_results_pressure_pain_ratings.svg'), width = 13, height = 13, units = "cm")
#ggsave(paste(save_path,'peep_hotspin_overview_results_pressure_pain_ratings.png'), width = 13, height = 13, units = "cm")


# Correlation Pressure
ggarrange(peep_hotspin_pressure_saline_fitness,peep_hotspin_pressure_saline_fitness_gender,peep_hotspin_pressure_saline_fitness_group,ncol = 3,labels = c("A", "B", "C"), font.label = list(size = 11))
#ggsave(paste(save_path,'comp_peep_hotspin_pressure_fitness_sex.svg'), width = 15, height = 6, units = "cm")
#ggsave(paste(save_path,'comp_peep_hotspin_pressure_fitness_sex.png'), width = 15, height = 6, units = "cm")


# Correlation Training Volume , Sex and Pain Heat
ggarrange(peep_hotspin_heat_saline_fitness_travol,peep_hotspin_heat_saline_fitness_gender_travol,peep_hotspin_heat_saline_fitness_group_travol,ncol = 3,labels = c("A", "B", "C"), font.label = list(size = 11))
#ggsave(paste(save_path,'comp_peep_hotspin_travol_sex:heat.svg'), width = 15, height = 6, units = "cm")
#ggsave(paste(save_path,'comp_peep_hotspin_travol_sex_heat.png'), width = 15, height = 6, units = "cm")

# Correlation Training Volume , Sex and Pain Pressure
ggarrange(peep_hotspin_pressure_saline_fitness_travol,peep_hotspin_pressure_saline_fitness_gender_travol,peep_hotspin_pressure_saline_fitness_group_travol,ncol = 3,labels = c("A", "B", "C"), font.label = list(size = 11))
#ggsave(paste(save_path,'comp_peep_hotspin_travol_sex_pressure.svg'), width = 15, height = 6, units = "cm")
#ggsave(paste(save_path,'comp_peep_hotspin_travol_sex_pressure.png'), width = 15, height = 6, units = "cm")


# Correlation Training Volume , Sex and Pain Heat + Pressure
ggarrange(
  ggarrange(peep_hotspin_pressure_saline_fitness_travol,peep_hotspin_pressure_saline_fitness_gender_travol,peep_hotspin_pressure_saline_fitness_group_travol,ncol = 3,labels = c("A", "B", "C"), font.label = list(size = 11)),
            ggarrange(peep_hotspin_heat_saline_fitness_travol,peep_hotspin_heat_saline_fitness_gender_travol,peep_hotspin_heat_saline_fitness_group_travol,ncol = 3,labels = c("D", "E", "F"), font.label = list(size = 11)), 
  nrow  = 2
)
            

  
 ggsave(paste(save_path,'comp_peep_hotspin_travol_sex.svg'), width = 17, height = 12, units = "cm")
ggsave(paste(save_path,'comp_peep_hotspin_travol_sex.png'), width = 17, height = 12, units = "cm")
