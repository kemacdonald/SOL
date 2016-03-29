library(polspline)

##### Load and clean up the data
d <- read.csv("../../analysis/eye_movements/sol_ss_all.csv")

d_all <- d %>% 
    filter(value_cat == "Target" | value_cat == "Distractor") %>% 
    mutate(C_D_count = ifelse(is.na(C_D_count), 0, C_D_count),
           total_trials_shifting = C_T_count + C_D_count,
           Sub.Num = as.character(Sub.Num)) %>% 
    dplyr::select(Sub.Num, age_peek_months, age_group, signs_produced, C_T_count, total_trials_shifting, 
                  mean_correct_rt, median_rt = median_ct_rt,C_T_prop, prop_looking, mean_prop_looking_TD,
                  value_cat, age_group_collapsed)

d %<>%  
    filter(age_group_collapsed == "Kids", value_cat == "Target") %>% 
    mutate(C_D_count = ifelse(is.na(C_D_count), 0, C_D_count),
           total_trials_shifting = C_T_count + C_D_count,
           Sub.Num = as.character(Sub.Num)) %>% 
    dplyr::select(Sub.Num, age_peek_months, signs_produced, C_T_count, total_trials_shifting, 
                  mean_correct_rt, median_rt = median_ct_rt,C_T_prop, prop_looking, mean_prop_looking_TD,
                  age_group_collapsed)

d_voc <- d %>% filter(is.na(signs_produced) == F)

####### Create, initialize, and adapt the model:

parameters = c( "beta0" ,  "beta1" , "sigma", 
                "zbeta0" , "zbeta1" , "zsigma" )
adaptSteps = 500  # Number of steps to "tune" the samplers
burnInSteps = 1000
nChains = 1 
thinSteps = 1
numSavedSteps = 20000
nIter = ceiling( ( numSavedSteps * thinSteps ) / nChains )

dataList_age = list(
    x = d$age_peek_months,
    y = d$mean_prop_looking_TD
)

parameters = c( "beta0" ,  "beta1" , "sigma", 
                "zbeta0" , "zbeta1" , "zsigma" )

# Get samples
samples <- jags(data = dataList_age, parameters.to.save = parameters,
                model.file="accuracy_model_null.txt", n.chains=1, n.iter=nIter, n.burnin = burnInSteps,
                n.thin=1, DIC=T)

samples_lm <- jags(data = dataList_age, parameters.to.save = parameters,
                model.file="accuracy_model.txt", n.chains=1, n.iter=nIter, n.burnin = burnInSteps, 
                n.thin=1, DIC=T)

df_null <- data.frame(beta0 = samples$BUGSoutput$sims.list$beta0, beta1 = samples$BUGSoutput$sims.list$beta1)
df_lm <- data.frame(beta0 = samples_lm$BUGSoutput$sims.list$beta0, beta1 = samples_lm$BUGSoutput$sims.list$beta1)

##### Get relevant parameter values for null model
age.seq <- seq(from = 10, to = 60, by = 1)
mu.link <- function(x) {df_null$beta0 + df_null$beta1*x}
mu <- sapply(age.seq, mu.link)
mu.mean <- apply(mu, 2, mean)
mu.HPDI <- apply(mu, 2, HPDI, prob = 0.95)
mu.HPDI_tidy <- data.frame(age_peek_months = age.seq, 
                           hpdi_lower = mu.HPDI[1,], hpdi_upper = mu.HPDI[2,],
                           mu.mean = mu.mean)

##### Get relevant parameter values for linear model
mu.link.lm <- function(x) {df_lm$beta0 + df_lm$beta1*x}
mu.lm <- sapply(age.seq, mu.link.lm)
mu.mean.lm <- apply(mu.lm, 2, mean)
mu.HPDI.lm <- apply(mu.lm, 2, HPDI, prob = 0.95)
mu.HPDI_tidy.lm <- data.frame(age_peek_months = age.seq, 
                           hpdi_lower = mu.HPDI.lm[1,], hpdi_upper = mu.HPDI.lm[2,],
                           mu.mean = mu.mean.lm)

## Visualize parameters

# Intercept-only model
ggplot(aes(x = beta1), data = df_null) + 
    geom_histogram() 

# Intercept and Slope model
ggplot(aes(x = beta1), data = df_lm) + 
    geom_histogram() 

## now plot
ggplot(data = d) +
    geom_line(aes(x = age_peek_months, y = mu.mean), data = mu.HPDI_tidy, size = 2) +
    geom_line(aes(x = age_peek_months, y = mu.mean.lm), data = mu.HPDI_tidy.lm, size = 2) +
    geom_point(aes(age_peek_months, mean_prop_looking_TD), color = "black", size = 5.5) +
    geom_point(aes(age_peek_months, mean_prop_looking_TD), color = "grey50", size = 4) +
    geom_ribbon(aes(x = age_peek_months, ymin = hpdi_lower, ymax = hpdi_upper), data = mu.HPDI_tidy,
                alpha = 0.2) +
    geom_ribbon(aes(x = age_peek_months, ymin = hpdi_lower, ymax = hpdi_upper), data = mu.HPDI_tidy.lm,
                alpha = 0.2) +
    ylab("Accuracy") +
    xlab("Child's Age (months)") +
    coord_cartesian(xlim=c(15, 55), ylim=c(0.25, 0.95)) +
    theme(axis.title.x = element_text(colour="grey30",size=22,
                                      angle=0,hjust=0.5,vjust=0,face="plain"),
          axis.title.y = element_text(colour="grey30",size=22,
                                      hjust=0.5,vjust=0.5,face="plain"),
          axis.text.x = element_text(colour="grey20",size=18,
                                     angle=0,hjust=0.5,vjust=0,face="plain"),
          axis.text.y = element_text(colour="grey20",size=18,
                                     angle=0,hjust=0.5,vjust=0,face="plain"),
          panel.grid.major=element_blank())

# Compute p(D) for posterior and prior: Savage-Dickey Method to get Bayes Factor
# Fits a density using spliens to approx. log-density
# uses 1997 knot and deletion algorithm
fit.posterior <- logspline(df_lm$beta1, lbound = 0)
posterior <- dlogspline(0, fit.posterior) # pdf @ beta=0
prior <- 2*dnorm(0) # height of order--restricted prior
posterior/prior # bayes factor

