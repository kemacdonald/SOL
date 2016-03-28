#### Bayesian Regression using BayesFactor Package ######

#### Setup 
rm(list=ls())

# Packages
library(reshape)
library(R2jags)
library(rjags)
library(BayesFactor)
library(dplyr)
library(ggplot2)
library(magrittr)
library(tidyr)
theme_set(theme_bw())

##### Load and clean up the data
data <- read.csv("../../analysis/eye_movements/sol_ss_all.csv")

data %<>%  
    filter(age_group_collapsed == "Kids", value_cat == "Target") %>% 
    mutate(C_D_count = ifelse(is.na(C_D_count), 0, C_D_count),
           total_trials_shifting = C_T_count + C_D_count,
           Sub.Num = as.character(Sub.Num)) %>% 
    dplyr::select(Sub.Num, age_group, age_peek_months, signs_produced, 
                  C_T_count,
                  total_trials_shifting, 
                  mean_correct_rt, 
                  median_rt = median_ct_rt,
                  C_T_prop, 
                  prop_looking, 
                  mean_prop_looking_TD)

data_voc <- data %>% filter(is.na(signs_produced) == F)

#### Fit Acc model
bf_acc <- regressionBF(mean_prop_looking_TD ~ age_peek_months + signs_produced, data = data_voc)

# what's the best model based on bayes factor?
which.max(bf_acc)
best_comp = head(bf_acc) / max(bf_acc)
bf_best_acc <- lmBF(prop_looking ~ age_peek_months, data = data_voc)
lm_acc <- lm(prop_looking ~ age_peek_months, data = data_voc)

# sample from the posterior
chains_acc = posterior(bf_best_acc, iterations = 10000)
summary(chains_acc)
# plot
plot(chains_acc)

#### Fit RT model
bf_rt = regressionBF(median_rt ~ age_peek_months + signs_produced, data = data_voc)
bf_rt

which.max(bf_rt)
bf_best_rt <- lmBF(median_rt ~ signs_produced, data = data_voc)

# sample from posterior 
chains_rt = posterior(bf_best_rt, iterations = 10000)
summary(chains_rt)

# plot
plot(chains_rt)
