######### Bayesian Analysis of Data on Development of Real-Time ASL Procesing Skills ###############

# Clear workspace
rm(list=ls())

# set wd
setwd("~/Documents/Projects/SOL/SOL-GIT/analysis/bayesian-data-analysis/")

# Packages
library(reshape)
library(polspline)
library(BayesFactor)
library(rethinking)
library(R2jags)
library(rjags)
library(dplyr)
library(ggplot2)
library(magrittr)
library(tidyr)
theme_set(theme_bw())

# set seed so we get the same simulation results
set.seed(seed = 10)

# Helper functions
find_mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
}
HDIofMCMC = function( sampleVec , credMass=0.95 ) {
    # Computes highest density interval from a sample of representative values,
    #   estimated as shortest credible interval.
    # Arguments:
    #   sampleVec
    #     is a vector of representative values from a probability distribution.
    #   credMass
    #     is a scalar between 0 and 1, indicating the mass within the credible
    #     interval that is to be estimated.
    # Value:
    #   HDIlim is a vector containing the limits of the HDI
    sortedPts = sort( sampleVec )
    ciIdxInc = ceiling( credMass * length( sortedPts ) )
    nCIs = length( sortedPts ) - ciIdxInc
    ciWidth = rep( 0 , nCIs )
    for ( i in 1:nCIs ) {
        ciWidth[ i ] = sortedPts[ i + ciIdxInc ] - sortedPts[ i ]
    }
    HDImin = sortedPts[ which.min( ciWidth ) ]
    HDImax = sortedPts[ which.min( ciWidth ) + ciIdxInc ]
    HDIlim = c( HDImin , HDImax )
    return( HDIlim )
}

##### Load and clean up the data
d <- read.csv("../../analysis/eye_movements/sol_ss_all.csv")

d %<>%  
    filter(age_group_collapsed == "Kids", value_cat == "Target") %>% 
    mutate(C_D_count = ifelse(is.na(C_D_count), 0, C_D_count),
           total_trials_shifting = C_T_count + C_D_count,
           Sub.Num = as.character(Sub.Num)) %>% 
    select(Sub.Num, age_peek_months, signs_produced, C_T_count, total_trials_shifting, 
                  median_rt = median_ct_rt, mean_prop_looking_TD, age_group_collapsed)

####### Standardize the data
d %<>% 
    mutate(age.s = (age_peek_months - mean(age_peek_months)) / sd(age_peek_months),
           acc.s = (mean_prop_looking_TD - mean(mean_prop_looking_TD)) / sd(mean_prop_looking_TD),
           rt.s = (median_rt - mean(median_rt)) / sd(median_rt))

d_voc <- d %>% 
    filter(is.na(signs_produced) == F) %>% 
    mutate(voc.s = (signs_produced - mean(signs_produced)) / sd(signs_produced))

####### Set MCMC parameters (same across all models)
parameters = c( "zbeta0" , "zbeta1" , "zsigma" )
adaptSteps = 500  # Number of steps to "tune" the samplers
burnInSteps = 1000
nChains = 1 
thinSteps = 1
numSavedSteps = 50000
nIter = ceiling( ( numSavedSteps * thinSteps ) / nChains )

######## Accuracy-Age model
dataList_age = list(
    x = d$age.s,
    y = d$acc.s,
    Ntotal = length(d$age.s)
)

dataList_age_prior = list(
    #y = d$mean_prop_looking_TD,
    x = d$age.s,
    Ntotal = length(d$age.s)
)

# Get samples conditioning on data
samples <- jags(data = dataList_age, parameters.to.save = parameters,
                model.file="accuracy_model.txt", n.chains=nChains, 
                n.iter=nIter, n.burnin = burnInSteps,
                n.thin=1, DIC=T)

df <- data.frame(zbeta0 = samples$BUGSoutput$sims.list$zbeta0, 
                 zbeta1 = samples$BUGSoutput$sims.list$zbeta1)

# transform standardized prior to measurement scale
df %<>% mutate(
    beta1 = zbeta1 * sd(d$mean_prop_looking_TD) / sd(d$age_peek_months),
    beta0 = zbeta0 * sd(d$mean_prop_looking_TD) + mean(d$mean_prop_looking_TD) - zbeta1 * mean(d$age_peek_months) * sd(d$mean_prop_looking_TD) / sd(d$age_peek_months)
    )
                 
# Get samples from prior
samples_prior <- jags(data = dataList_age_prior, parameters.to.save = parameters,
                model.file="accuracy_model.txt", n.chains=nChains, 
                n.iter=nIter, n.burnin = burnInSteps,
                n.thin=1, DIC=F)

df_prior <- data.frame(zbeta0 = samples_prior$BUGSoutput$sims.list$zbeta0, 
                       zbeta1 = samples_prior$BUGSoutput$sims.list$zbeta1)

# transform standardized prior to measurement scale
df_prior %<>% 
    mutate( beta1 = zbeta1 * sd(d$mean_prop_looking_TD) / sd(d$age_peek_months),
            beta0 = zbeta0 * sd(d$mean_prop_looking_TD) + mean(d$mean_prop_looking_TD) - (zbeta1 * mean(d$age_peek_months) * sd(d$mean_prop_looking_TD)) / sd(d$age_peek_months)
            )

### Plot prior predictive on measurement scale to make sure we get reasonable values
ggplot() + geom_density(aes(df_prior$beta0)) 
ggplot() + geom_density(aes(df_prior$beta1))

df_prior_samples <- sample_n(select(df_prior, beta0, beta1), size = 100)

ggplot(data = d) +
    geom_point(aes(age_peek_months, mean_prop_looking_TD), color = "black", size = 5.5) +
    geom_point(aes(age_peek_months, mean_prop_looking_TD), color = "grey50", size = 4) +
    geom_abline(aes(intercept = 0.32, slope = beta1), data = df_prior_samples, alpha = 0.3) +
    xlim(0,55) +
    ylim(0, 1)

# Compute p(D) for posterior and prior: Savage-Dickey Method to get Bayes Factor
# Fits a density using splines to approx. log-density
# uses 1997 knot and deletion algorithm
fit.posterior <- logspline(df$zbeta1, lbound = 0)
posterior <- dlogspline(0, fit.posterior) # pdf @ beta=0

# get prior density using logspline
fit.prior <- logspline(df_prior$zbeta1, lbound = 0)
prior <- dlogspline(0, fit.prior) # pdf @ beta=0
prior_analytic <- 2*dnorm(0, 1) # pdf @ beta=0
bf_acc_age <- posterior/prior # bayes factor
1/bf_acc_age

##### Get relevant parameter values from posterior distribution
post_mode_b1 <- round(find_mode(df$beta1), 3)     # MAP of slope
HDI_b1 <- HDIofMCMC( df$beta1 , credMass = 0.95 ) # HDI of slope
post_mode_b0 <- round(find_mode(df$beta0), 3)     # Map of intercept

## Visualize posterior distribution on slope of regression line
ggplot(aes(x = beta1), data = df) + 
    geom_density() +
    geom_vline(xintercept = HDI_b1[1], color = 'red') +
    geom_vline(xintercept = HDI_b1[2], color = 'red')

## simulate posteior distribution of mu for each value of age to get HDI for MAP regression line
mu.link <- function(x) {df$beta0 + df$beta1*x}
age.seq <- seq(from = 10, to = 60, by = 1)
mu <- sapply(age.seq, mu.link)
mu.mean <- apply(mu, 2, mean)
mu.HPDI <- apply(mu, 2, HPDI, prob = 0.95)
mu.HPDI_tidy <- data.frame(age_peek_months = age.seq, 
                           hpdi_lower = mu.HPDI[1,], hpdi_upper = mu.HPDI[2,],
                           mu.mean = mu.mean)

ggplot(data = d) +
    geom_line(aes(x = age_peek_months, y = mu.mean), data = mu.HPDI_tidy, size = 2) +
    geom_point(aes(age_peek_months, mean_prop_looking_TD), color = "black", size = 5.5) +
    geom_point(aes(age_peek_months, mean_prop_looking_TD), color = "grey50", size = 4) +
    geom_ribbon(aes(x = age_peek_months, ymin = hpdi_lower, ymax = hpdi_upper), data = mu.HPDI_tidy,
                alpha = 0.2) +
    ylab("Accuracy") +
    xlab("Child's Age (months)") +
    coord_cartesian(xlim=c(15, 55), ylim=c(0.25, 0.95)) +
    ggtitle(bquote(list(alpha==.(as.character(round(find_mode(df$beta0), 2))), 
                        beta==.(as.character(round(find_mode(df$beta1), 3))),
                        HDI==.(paste("95% HDI [", round(HDI_b1[1],3) , "," , round(HDI_b1[2], 3), "]"))))) +
    theme(axis.title.x = element_text(colour="grey30",size=22,
                                      angle=0,hjust=0.5,vjust=0,face="plain"),
          axis.title.y = element_text(colour="grey30",size=22,
                                      hjust=0.5,vjust=0.5,face="plain"),
          axis.text.x = element_text(colour="grey20",size=18,
                                     angle=0,hjust=0.5,vjust=0,face="plain"),
          axis.text.y = element_text(colour="grey20",size=18,
                                     angle=0,hjust=0.5,vjust=0,face="plain"),
          panel.grid.major=element_blank())

##################### Accuracy-Vocab Model #####################

dataList_voc = list(
    x = d_voc$voc.s,
    y = d_voc$acc.s,
    Ntotal = length(d_voc$age.s)
)

dataList_voc_prior = list(
    #y = d_voc$mean_prop_looking_TD,
    x = d_voc$voc.s,
    Ntotal = length(d_voc$age.s)
)

# Get samples
samples <- jags(data = dataList_voc, parameters.to.save = parameters,
                model.file="accuracy_model.txt", n.chains=nChains, n.iter=nIter, 
                n.thin=1, DIC=T, n.burnin=burnInSteps)

df <- data.frame(zbeta0 = samples$BUGSoutput$sims.list$zbeta0, 
                 zbeta1 = samples$BUGSoutput$sims.list$zbeta1)

# transform standardized prior to measurement scale
df %<>%
    mutate( beta1 = zbeta1 * sd(d_voc$mean_prop_looking_TD) / sd(d_voc$signs_produced),
            beta0 = zbeta0 * sd(d_voc$mean_prop_looking_TD) + mean(d_voc$signs_produced) - (zbeta1 * mean(d_voc$signs_produced) * sd(d_voc$mean_prop_looking_TD)) / sd(d_voc$signs_produced)
    )

# Get samples from prior
samples_prior <- jags(data = dataList_voc_prior, parameters.to.save = parameters,
                      model.file="accuracy_model.txt", n.chains=nChains, 
                      n.iter=nIter, n.burnin = burnInSteps,
                      n.thin=1, DIC=F)

df_prior <- data.frame(zbeta0 = samples_prior$BUGSoutput$sims.list$zbeta0, 
                       zbeta1 = samples_prior$BUGSoutput$sims.list$zbeta1)

# transform standardized prior to measurement scale
df_prior %<>% 
    mutate( beta1 = zbeta1 * sd(d_voc$mean_prop_looking_TD) / sd(d_voc$signs_produced),
            beta0 = zbeta0 * sd(d_voc$mean_prop_looking_TD) + mean(d_voc$signs_produced) - (zbeta1 * mean(d_voc$signs_produced) * sd(d_voc$mean_prop_looking_TD)) / sd(d_voc$signs_produced)
    )


### Plot prior predictive on measruement scale to make sure we get reasonable slope and intercept values
ggplot() + geom_density(aes(df_prior$beta0)) 
ggplot() + geom_density(aes(df_prior$beta1)) 

# Compute p(D) for posterior and prior: Savage-Dickey Method to get Bayes Factor
# Fits a density using splines to approx. log-density
# uses 1997 knot and deletion algorithm
fit.posterior <- logspline(df$beta1, lbound = 0)
posterior <- dlogspline(0, fit.posterior) # pdf @ beta=0

# get prior density using logspline
fit.prior <- logspline(df_prior$beta1, lbound = 0)
prior <- dlogspline(0, fit.prior) # pdf @ beta=0
bf_acc_voc <- posterior/prior # bayes factor
1/bf_acc_voc

####### Compute p(D) for posterior and prior: Savage-Dickey Method to get Bayes Factor

# Fits a density using spliens to approx. log-density
# uses 1997 knot and deletion algorithm
fit.posterior <- logspline(df_acc_voc$beta1, lbound = 0)
posterior <- dlogspline(0, fit.posterior) # pdf @ beta=0
prior <- 2*dnorm(0) # height of order--restricted prior
bf_voc_acc <- posterior/prior # bayes factor

##### Visualize parameter values

## get MAP and HDI
post_mode_b1_accvoc <- round(find_mode(df_acc_voc$beta1), 3)
post_mode_b0_accvoc <- round(find_mode(df_acc_voc$beta0), 3)
HDI_b1_accvoc <- HDIofMCMC(df_acc_voc$beta1 , credMass = 0.95 )

## Visualize posterior distribution on slope of regression line
ggplot(aes(x = beta1), data = df_acc_voc) + 
    geom_density() +
    geom_vline(xintercept = HDI_b1_accvoc[1], color = 'red') +
    geom_vline(xintercept = HDI_b1_accvoc[2], color = 'red')

## simulate to get HDI for regression line
vocab.seq <- seq(from = 15, to = 85, by = 1)
mu.link <- function(x) {df_acc_voc$beta0 + df_acc_voc$beta1*x}
mu <- sapply(vocab.seq, mu.link)
mu.mean <- apply(mu, 2, mean)
mu.HPDI <- apply(mu, 2, HPDI, prob = 0.95)
mu.HPDI_tidy <- data.frame(signs_produced = vocab.seq, hpdi_lower = mu.HPDI[1,], hpdi_upper = mu.HPDI[2,])

ggplot(data = d_voc) +
    geom_abline(intercept = post_mode_b0_accvoc, slope = post_mode_b1_accvoc, size = 2) +
    geom_point(aes(signs_produced, mean_prop_looking_TD), color = "black", size = 5.5) +
    geom_point(aes(signs_produced, mean_prop_looking_TD), color = "grey50", size = 4) +
    geom_ribbon(aes(x = signs_produced, ymin = hpdi_lower, ymax = hpdi_upper), data = mu.HPDI_tidy,
                alpha = 0.2) +
    ylab("Accuracy") +
    xlab("Signs Produced") +
    coord_cartesian(xlim=c(18, 82), ylim=c(0.25, 0.95)) +
    ggtitle(bquote(list(alpha==.(as.character(round(post_mode_b0_accvoc, 2))), 
                        beta==.(as.character(round(post_mode_b1_accvoc, 3))),
                        HDI=.(paste("95% HDI [", round(HDI_b1_accvoc[1],3) , "," , round(HDI_b1_accvoc[2], 3), "]"))))) +
    theme(axis.title.x = element_text(colour="grey30",size=22,
                                      angle=0,hjust=0.5,vjust=0,face="plain"),
          axis.title.y = element_text(colour="grey30",size=22,
                                      hjust=0.5,vjust=0.5,face="plain"),
          axis.text.x = element_text(colour="grey20",size=18,
                                     angle=0,hjust=0.5,vjust=0,face="plain"),
          axis.text.y = element_text(colour="grey20",size=18,
                                     angle=0,hjust=0.5,vjust=0,face="plain"),
          panel.grid.major=element_blank())

######## RT-Vocab MODEL: Latent mixture + Linear Regression ###########

dataList = list(
    x = d_voc$signs_produced,
    y = d_voc$median_rt,
    n_correct = d_voc$C_T_count,
    n_trials = d_voc$total_trials_shifting
)

# Specify initial values for phi and indicator var z
myinits <-  list(list(phi = 0.75, z = round(runif(length(d_voc$total_trials_shifting)))))

# Specify parameters to monitor (same for both RT models)
parameters = c( "beta0" ,  "beta1" , "sigma", 
                "zbeta0" , "zbeta1" , "zsigma",
                "phi" , "z")

# Get samples
samples <- jags(data = dataList, parameters.to.save = parameters,
                model.file="rt_model.txt", n.chains=nChains, n.iter=nIter, n.burnin = burnInSteps,
                n.thin=1, DIC=T)

df_rt_voc <- data.frame(beta0 = samples$BUGSoutput$sims.list$beta0,
                 beta1 = samples$BUGSoutput$sims.list$beta1,
                 phi = samples$BUGSoutput$sims.list$phi,
                 z = samples$BUGSoutput$sims.list$z)

####### Compute p(D) for posterior and prior: Savage-Dickey Method to get Bayes Factor

# Fits a density using spliens to approx. log-density
# uses 1997 knot and deletion algorithm
fit.posterior <- logspline(df_rt_voc$beta1, ubound = 0)
posterior <- dlogspline(0, fit.posterior) # pdf @ beta=0
prior <- 2*dnorm(0) # height of order--restricted prior
bf_rt_voc <- posterior/prior # bayes factor

##### Visualize parameter values

## get HDI and MAP
HDI_rt_voc_b1 <- HDIofMCMC( df_rt_voc$beta1 , credMass = 0.95 )
post_mode_b1_rt_voc <- round(find_mode(df_rt_voc$beta1), 3)
post_mode_b0_rt_voc <- round(find_mode(df_rt_voc$beta0), 3)

# plot
ggplot(aes(x = beta1), data = df_rt_voc) + 
    geom_density() +
    geom_vline(xintercept = HDI_rt_voc_b1[1], color = 'red') +
    geom_vline(xintercept = HDI_rt_voc_b1[2], color = 'red')

## simulate posterior means for each value of vocab to get HDI for regression line
vocab.seq <- seq(from = 15, to = 85, by = 1)
mu.link <- function(x) {df_rt_voc$beta0 + df_rt_voc$beta1*x}
mu <- sapply(vocab.seq, mu.link)
mu.mean <- apply(mu, 2, mean)
mu.HPDI <- apply(mu, 2, HPDI, prob = 0.95)
mu.HPDI_tidy <- data.frame(signs_produced = vocab.seq, hpdi_lower = mu.HPDI[1,], hpdi_upper = mu.HPDI[2,],
                           mu.rt = mu.mean)

# plot
ggplot(data = d_voc) +
    geom_line(aes(x = signs_produced, y = mu.mean), data = mu.HPDI_tidy, size = 2) +
    geom_point(aes(signs_produced, median_rt), color = "black", size = 5.5) +
    geom_point(aes(signs_produced, median_rt), color = "grey50", size = 4) +
    geom_ribbon(aes(x = signs_produced, ymin = hpdi_lower, ymax = hpdi_upper), data = mu.HPDI_tidy,
                alpha = 0.2) +
    ylab("Reaction Time (ms)") +
    xlab("Signs Produced") +
    coord_cartesian(xlim=c(18, 82), ylim=c(700, 2000)) +
    ggtitle(bquote(list(alpha==.(as.character(round(post_mode_b0_rt_voc, 2))), 
                        beta==.(as.character(round(post_mode_b1_rt_voc, 2))),
                        HDI=.(paste("95% HDI [", round(HDI_rt_voc_b1[1],2) , "," , round(HDI_rt_voc_b1[2], 2), "]"))))) +
    theme(axis.title.x = element_text(colour="grey30",size=22,
                                      angle=0,hjust=0.5,vjust=0,face="plain"),
          axis.title.y = element_text(colour="grey30",size=22,
                                      hjust=0.5,vjust=0.5,face="plain"),
          axis.text.x = element_text(colour="grey20",size=18,
                                     angle=0,hjust=0.5,vjust=0,face="plain"),
          axis.text.y = element_text(colour="grey20",size=18,
                                     angle=0,hjust=0.5,vjust=0,face="plain"),
          panel.grid.major=element_blank())

######## RT-Age MODEL: Latent mixture + Linear Regression ############

# Specify the data in a list, for later shipment to JAGS:
dataList = list(
    x = d$age_peek_months,
    y = d$median_rt,
    n_correct = d$C_T_count,
    n_trials = d$total_trials_shifting
)

# Specify initial values
myinits <-  list(list(phi = 0.75, z = round(runif(length(d$total_trials_shifting)))))

####### Create, initialize, and adapt the model:
samples <- jags(data = dataList, parameters.to.save = parameters,
                model.file="rt_model.txt", n.chains=1, n.iter=nIter, 
                n.burnin = burnInSteps, n.thin=1, DIC=T)

df_rt_age <- data.frame(beta0 = samples$BUGSoutput$sims.list$beta0,
                        beta1 = samples$BUGSoutput$sims.list$beta1,
                        phi = samples$BUGSoutput$sims.list$phi,
                        z = samples$BUGSoutput$sims.list$z)

####### Compute p(D) for posterior and prior: Savage-Dickey Method to get Bayes Factor

# Fits a density using spliens to approx. log-density
# uses 1997 knot and deletion algorithm
fit.posterior <- logspline(df_rt_age$beta1, ubound = 0)
posterior <- dlogspline(0, fit.posterior) # pdf @ beta=0
prior <- 2*dnorm(0) # height of order--restricted prior
bf_rt_age <- posterior/prior # bayes factor

##### Visualize parameter values

## get HDI and MAP
HDI_rt_age_b1 <- HDIofMCMC(df_rt_age$beta1 , credMass = 0.95 )
post_rt_age_mode_b1 <- round(find_mode(df_rt_age$beta1), 3)
post_rt_age_mode_b0 <- round(find_mode(df_rt_age$beta0), 3)

## plot
ggplot(aes(x = beta1), data = df_rt_age) + 
    geom_density() +
    geom_vline(xintercept = HDI_rt_age_b1[1], color = 'red') +
    geom_vline(xintercept = HDI_rt_age_b1[2], color = 'red')

# simulate to get HDI around MAP regression line
mu.link <- function(x) {df_rt_age$beta0 + df_rt_age$beta1*x}
age.seq <- seq(from = 10, to = 60, by = 1)
mu <- sapply(age.seq, mu.link)
mu.mean <- apply(mu, 2, mean)
mu.HPDI <- apply(mu, 2, HPDI, prob = 0.95)
mu.HPDI_tidy <- data.frame(age_peek_months = age.seq, hpdi_lower = mu.HPDI[1,], hpdi_upper = mu.HPDI[2,],
                           mu.rt = mu.mean)

# now plot
ggplot(data = d) +
    geom_line(aes(x = age_peek_months, y = mu.rt), data = mu.HPDI_tidy, size = 2) +
    geom_point(aes(age_peek_months, median_rt), color = "black", size = 5.5) +
    geom_point(aes(age_peek_months, median_rt), color = "grey50", size = 4) +
    geom_ribbon(aes(x = age_peek_months, ymin = hpdi_lower, ymax = hpdi_upper), data = mu.HPDI_tidy,
                alpha = 0.2) +
    ylab("Reaction Time (ms)") +
    xlab("Child's Age (months)") +
    coord_cartesian(xlim=c(15, 55), ylim=c(700, 1800)) +
    ggtitle(bquote(list(alpha==.(as.character(round(post_rt_age_mode_b0, 2))), 
                        beta==.(as.character(round(post_rt_age_mode_b1, 2))),
                        HDI=.(paste("95% HDI [", round(HDI_rt_age_b1[1],3) , "," , round(HDI_rt_age_b1[2], 3), "]"))))) +
    theme(axis.title.x = element_text(colour="grey30",size=22,
                                      angle=0,hjust=0.5,vjust=0,face="plain"),
          axis.title.y = element_text(colour="grey30",size=22,
                                      hjust=0.5,vjust=0.5,face="plain"),
          axis.text.x = element_text(colour="grey20",size=18,
                                     angle=0,hjust=0.5,vjust=0,face="plain"),
          axis.text.y = element_text(colour="grey20",size=18,
                                     angle=0,hjust=0.5,vjust=0,face="plain"),
          panel.grid.major=element_blank())

######## RT-Acc MODEL: Latent mixture + Linear Regression ############
# Specify the data in a list, for later shipment to JAGS:
dataList = list(
    x = d$mean_prop_looking_TD,
    y = d$median_rt,
    n_correct = d$C_T_count,
    n_trials = d$total_trials_shifting
)

# Specify initial values
myinits <-  list(list(phi = 0.75, z = round(runif(length(d$total_trials_shifting)))))

####### Create, initialize, and adapt the model:
jagsModel = jags.model( "rt_model.txt" , data=dataList , inits=myinits , 
                        n.chains=nChains , n.adapt=adaptSteps )

# Burn-in:
update( jagsModel , n.iter=burnInSteps )

# Get samples
samples <- jags(data = dataList, parameters.to.save = parameters,
                model.file="rt_model.txt", n.chains=1, n.iter=10000, 
                n.thin=1, DIC=T)

df_acc_rt <- data.frame(beta0 = samples$BUGSoutput$sims.list$beta0,
                        beta1 = samples$BUGSoutput$sims.list$beta1,
                        phi = samples$BUGSoutput$sims.list$phi,
                        z = samples$BUGSoutput$sims.list$z)

####### Compute p(D) for posterior and prior: Savage-Dickey Method to get Bayes Factor

# Fits a density using spliens to approx. log-density
# uses 1997 knot and deletion algorithm
fit.posterior <- logspline(df_acc_rt$beta1, ubound = 0)
posterior <- dlogspline(0, fit.posterior) # pdf @ beta=0
prior <- 2*dnorm(0) # height of order--restricted prior
bf_acc_rt <- posterior/prior # bayes factor

##### Visualize parameter values

## get HDI and MAP
HDI_acc_rt_b1 <- HDIofMCMC(df_acc_rt$beta1 , credMass = 0.95 )
post_acc_rt_mode_b1 <- round(find_mode(df_acc_rt$beta1), 3)
post_acc_rt_mode_b0 <- round(find_mode(df_acc_rt$beta0), 3)

## plot
ggplot(aes(x = beta1), data = df_acc_rt) + 
    geom_density() +
    geom_vline(xintercept = HDI_acc_rt_b1[1], color = 'red') +
    geom_vline(xintercept = HDI_acc_rt_b1[2], color = 'red')

# simulate to get HDI around MAP regression line
mu.link <- function(x) {df_acc_rt$beta0 + df_acc_rt$beta1*x}
acc.seq <- seq(from = 0.2, to = 0.9, by = 0.1)
mu <- sapply(acc.seq, mu.link)
mu.mean <- apply(mu, 2, mean)
mu.HPDI <- apply(mu, 2, HPDI, prob = 0.95)
mu.HPDI_tidy <- data.frame(mean_prop_looking_TD = acc.seq, hpdi_lower = mu.HPDI[1,], hpdi_upper = mu.HPDI[2,],
                           mu.rt = mu.mean)

ggplot(data = d) +
    geom_line(aes(mean_prop_looking_TD, mu.rt), data = mu.HPDI_tidy, size = 2) +
    geom_point(aes(mean_prop_looking_TD, median_rt), color = "black", size = 5.5) +
    geom_point(aes(mean_prop_looking_TD, median_rt), color = "grey50", size = 4) +
    geom_ribbon(aes(x = mean_prop_looking_TD, ymin = hpdi_lower, ymax = hpdi_upper), data = mu.HPDI_tidy,
                alpha = 0.2) +
    ylab("Reaction Time (ms)") +
    xlab("Mean Accuracy") +
    coord_cartesian(xlim=c(0.2, 0.9), ylim=c(700, 1800)) +
    ggtitle(bquote(list(alpha==.(as.character(round(post_acc_rt_mode_b0, 2))), 
                        beta==.(as.character(round(post_acc_rt_mode_b1, 2))),
                        HDI=.(paste("95% HDI [", round(HDI_acc_rt_b1[1],2) , "," , round(HDI_acc_rt_b1[2], 2), "]"))))) +
    theme(axis.title.x = element_text(colour="grey30",size=22,
                                      angle=0,hjust=0.5,vjust=0,face="plain"),
          axis.title.y = element_text(colour="grey30",size=22,
                                      hjust=0.5,vjust=0.5,face="plain"),
          axis.text.x = element_text(colour="grey20",size=18,
                                     angle=0,hjust=0.5,vjust=0,face="plain"),
          axis.text.y = element_text(colour="grey20",size=18,
                                     angle=0,hjust=0.5,vjust=0,face="plain"),
          panel.grid.major=element_blank())

####### Visualize posterior distribution on group membership and group level success

# grab subject numbers and add to data frame
colnames(df_rt_age)[which(names(df_rt_age)=="z.1"):which(names(df_rt_age)=="z.29")] <- as.character(d$Sub.Num)

# melt data frame
df_mixture <- df_rt_age %>% dplyr::select(-beta0, -beta1)
df_melt <- reshape::melt(df_mixture[,2:length(df_mixture)], variable.name = "Sub.Num", value.name = "group_membership")

# change factor label
df_melt %<>% mutate(group_membership_factor = factor(df_melt$value, labels = c("G", "K"))) 

a <- qplot(data=melt(df_rt_age$phi),x=value,geom='histogram',binwidth=0.008)+
    theme_bw()+
    xlab('phi')+
    xlim(0,1)

b <- ggplot(data=df_melt,aes(x=group_membership_factor, 
                             fill=group_membership_factor)) +
    geom_bar(stat="count") +
    facet_wrap(~variable,scales='fixed')+
    guides(fill=F) +
    langcog::scale_fill_solarized() +
    theme_bw()+
    xlab('Group Membership') +
    ylab("Count")

gridExtra::grid.arrange(a,b,nrow=1, top='Group success rate and individual group membership', 
             widths = c(2, 4))
