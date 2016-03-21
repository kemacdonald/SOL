# Clear workspace
rm(list=ls())
# Packages
library(reshape)
library(R2jags)
library(rjags)
library(dplyr)
library(ggplot2)
library(magrittr)
library(tidyr)
theme_set(theme_bw())

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
data <- read.csv("../../analysis/eye_movements/sol_ss_all.csv")

data %<>%  
    filter(age_group_collapsed == "Kids", value_cat == "Target") %>% 
    mutate(C_D_count = ifelse(is.na(C_D_count), 0, C_D_count),
           total_trials_shifting = C_T_count + C_D_count,
           Sub.Num = as.character(Sub.Num)) %>% 
    dplyr::select(Sub.Num, age_peek_months, 
                  signs_produced, 
                  C_T_count,
                  total_trials_shifting, 
                  mean_correct_rt, 
                  median_rt = median_ct_rt,
                  C_T_prop, 
                  prop_looking, 
                  mean_prop_looking_TD)

data_voc <- data %>% filter(is.na(signs_produced) == F)

##### Structure the data
rt <-data$mean_correct_rt
rt_median <- data$median_rt
n_correct <- data$C_T_count
n_trials <- data$total_trials_shifting  
acc <- data$mean_prop_looking_TD
age <- data$age_peek_months
vocab <- data$signs_produced

# For Vocab models we need to remove one participant for whom we didn't collect CDI data
rt_voc <- data_voc$mean_correct_rt
rt_median_voc <- data_voc$median_rt
n_correct_voc <- data_voc$C_T_count
n_trials_voc <- data_voc$total_trials_shifting  
acc_voc <- data_voc$mean_prop_looking_TD
age_voc <- data_voc$age_peek_months
vocab_voc <- data_voc$signs_produced

######## Accuracy-Age model

dataList_age = list(
    x = age,
    y = acc
)


modelString = "
        # Standardize the data:
        data {
        Ntotal <- length(y)
        xm <- mean(x)
        ym <- mean(y)
        xsd <- sd(x)
        ysd <- sd(y)
        for ( i in 1:length(y) ) {
        zx[i] <- ( x[i] - xm ) / xsd
        zy[i] <- ( y[i] - ym ) / ysd
        }
        }
        # Specify the model for standardized data:
        model {
        for ( i in 1:Ntotal ) {
             zy[i] ~ dt( zbeta0 + zbeta1 * zx[i] , 1/(zsigma)^2 , nu )
        }
        # Priors vague on standardized scale:
        zbeta0 ~ dnorm( 0 , 1/(10)^2 )  
        zbeta1 ~ dnorm( 0 , 1/(10)^2 )
        zsigma ~ dunif( 1.0E-3 , 1.0E+3 )
        nu <- nuMinusOne+1
        nuMinusOne ~ dexp(1/29.0)
        # Transform to original scale:
        beta1 <- zbeta1 * ysd / xsd  
        beta0 <- zbeta0 * ysd  + ym - zbeta1 * xm * ysd / xsd 
        sigma <- zsigma * ysd
        }
        " 
writeLines( modelString , con="acc_age_model.txt" )

###### RUN THE CHAINS
parameters = c( "beta0" ,  "beta1" , "sigma", 
                "zbeta0" , "zbeta1" , "zsigma", "nu")
adaptSteps = 500  # Number of steps to "tune" the samplers
burnInSteps = 1000
nChains = 1 
thinSteps = 1
numSavedSteps = 20000
nIter = ceiling( ( numSavedSteps * thinSteps ) / nChains )

####### Create, initialize, and adapt the model:
jagsModel <- jags.model( "acc_age_model.txt" , data=dataList_age , inits=NULL , 
                        n.chains=nChains , n.adapt=adaptSteps )

# Burn-in:
update( jagsModel , n.iter=burnInSteps )

# Get samples
samples <- jags(data = dataList_age, parameters.to.save = parameters,
                model.file="acc_age_model.txt", n.chains=1, n.iter=10000, 
                n.thin=1, DIC=T)

df <- data.frame(beta0 = samples$BUGSoutput$sims.list$beta0,
                 beta1 = samples$BUGSoutput$sims.list$beta1)

##### Visualize parameter values

## get MAP and HDI
post_mode <- round(find_mode(df$beta1), 3)
HDI <- HDIofMCMC( df$beta1 , credMass = 0.95 )

## Visualize posterior distribution on slope of regression line
ggplot(aes(x = beta1), data = df) + 
    geom_density() +
    geom_vline(xintercept = HDI[1], color = 'red') +
    geom_vline(xintercept = HDI[2], color = 'red')

#### Plot MAP regression line with scatter plot
slope_acc_age <- samples$BUGSoutput$median$beta1
intercept_acc_age <- samples$BUGSoutput$median$beta0

# grab 20 samples from posterior 
post_20 <- dplyr::sample_n(df, size = 20)

ggplot(data = data, aes(age, acc)) +
    geom_abline(intercept = intercept_acc_age, slope = slope_acc_age, size = 2) +
    geom_abline(data = post_20, aes(intercept = beta0, slope = beta1), size = 0.2, alpha = 0.5) +
    geom_point(color = "black", size = 5.5) +
    geom_point(color = "grey50", size = 4) +
    # geom_hline(yintercept = 0.33, linetype = "dashed") +
    ylab("Accuracy") +
    xlab("Child's Age (months)") +
    # annotate("text", x = 45, y = 0.4, label = paste("MAP =", post_mode), size=8) +
    ggtitle(bquote(list(alpha==.(as.character(round(intercept_acc_age, 2))), 
                        beta==.(as.character(post_mode)),
                        HDI=.(paste("95% HDI [", round(HDI[1],3) , "," , round(HDI[2], 3), "]"))))) +
    theme(axis.title.x = element_text(colour="grey30",size=22,
                                      angle=0,hjust=0.5,vjust=0,face="plain"),
          axis.title.y = element_text(colour="grey30",size=22,
                                      hjust=0.5,vjust=0.5,face="plain"),
          axis.text.x = element_text(colour="grey20",size=18,
                                     angle=0,hjust=0.5,vjust=0,face="plain"),
          axis.text.y = element_text(colour="grey20",size=18,
                                     angle=0,hjust=0.5,vjust=0,face="plain"),
          panel.grid.major=element_blank())

####### Accuracy-Vocab Model ###########

dataList_voc = list(
    x = vocab_voc,
    y = acc_voc
)

modelString = "
# Standardize the data:
data {
Ntotal <- length(y)
xm <- mean(x)
ym <- mean(y)
xsd <- sd(x)
ysd <- sd(y)
for ( i in 1:length(y) ) {
zx[i] <- ( x[i] - xm ) / xsd
zy[i] <- ( y[i] - ym ) / ysd
}
}
# Specify the model for standardized data:
model {
for ( i in 1:Ntotal ) {
zy[i] ~ dt( zbeta0 + zbeta1 * zx[i] , 1/(zsigma)^2 , nu )
}
# Priors vague on standardized scale:
zbeta0 ~ dnorm( 0 , 1/(10)^2 )  
zbeta1 ~ dnorm( 0 , 1/(10)^2 )
zsigma ~ dunif( 1.0E-3 , 1.0E+3 )
nu <- nuMinusOne+1
nuMinusOne ~ dexp(1/29.0)
# Transform to original scale:
beta1 <- zbeta1 * ysd / xsd  
beta0 <- zbeta0 * ysd  + ym - zbeta1 * xm * ysd / xsd 
sigma <- zsigma * ysd
}
" 
writeLines( modelString , con="acc_vocab_model.txt" )

###### RUN THE CHAINS
parameters = c( "beta0" ,  "beta1" , "sigma", 
                "zbeta0" , "zbeta1" , "zsigma", "nu")
adaptSteps = 500  # Number of steps to "tune" the samplers
burnInSteps = 1000
nChains = 1 
thinSteps = 1
numSavedSteps = 20000
nIter = ceiling( ( numSavedSteps * thinSteps ) / nChains )

####### Create, initialize, and adapt the model:
jagsModel <- jags.model( "acc_vocab_model.txt" , data=dataList_voc , inits=NULL , 
                         n.chains=nChains , n.adapt=adaptSteps )

# Burn-in:
update( jagsModel , n.iter=burnInSteps )

# Get samples
samples <- jags(data = dataList_voc, parameters.to.save = parameters,
                model.file="acc_vocab_model.txt", n.chains=1, n.iter=10000, 
                n.thin=1, DIC=T)

df_acc_voc <- data.frame(beta0 = samples$BUGSoutput$sims.list$beta0,
                 beta1 = samples$BUGSoutput$sims.list$beta1)

##### Visualize parameter values

## get MAP and HDI
post_mode_acc_voc <- round(find_mode(df_acc_voc$beta1), 3)
HDI_acc_voc <- HDIofMCMC( df_acc_voc$beta1 , credMass = 0.95 )

## Visualize posterior distribution on slope of regression line
ggplot(aes(x = beta1), data = df_acc_voc) + 
    geom_density() +
    geom_vline(xintercept = HDI_acc_voc[1], color = 'red') +
    geom_vline(xintercept = HDI_acc_voc[2], color = 'red')

#### Plot MAP regression line with scatter plot
slope_acc_voc <- samples$BUGSoutput$median$beta1
intercept_acc_voc <- samples$BUGSoutput$median$beta0

# grab 20 samples from posterior 
post_20_acc_voc <- dplyr::sample_n(df_acc_voc, size = 20)

ggplot(data = data_voc, aes(signs_produced, mean_prop_looking_TD)) +
    geom_abline(intercept = intercept_acc_voc, slope = slope_acc_voc, size = 2) +
    geom_abline(data = post_20_acc_voc, aes(intercept = beta0, slope = beta1), size = 0.2, alpha = 0.5) +
    geom_point(color = "black", size = 5.5) +
    geom_point(color = "grey50", size = 4) +
    #geom_hline(yintercept = 0.33, linetype = "dashed") +
    ylab("Accuracy") +
    xlab("Signs Produced") +
   # annotate("text", x = 65, y = 0.37, label = paste("MAP =", post_mode), size=8) +
    ggtitle(bquote(list(alpha==.(as.character(round(intercept_acc_voc, 2))), 
                        beta==.(as.character(round(slope_acc_voc, 3))),
                        HDI=.(paste("95% HDI [", round(HDI_acc_voc[1],3) , "," , round(HDI_acc_voc[2], 3), "]"))))) +
    theme(axis.title.x = element_text(colour="grey30",size=22,
                                      angle=0,hjust=0.5,vjust=0,face="plain"),
          axis.title.y = element_text(colour="grey30",size=22,
                                      hjust=0.5,vjust=0.5,face="plain"),
          axis.text.x = element_text(colour="grey20",size=18,
                                     angle=0,hjust=0.5,vjust=0,face="plain"),
          axis.text.y = element_text(colour="grey20",size=18,
                                     angle=0,hjust=0.5,vjust=0,face="plain"),
          panel.grid.major=element_blank())

######## RT-Vocab MODEL: Latent mixture + Linear Regression

# Specify the data in a list, for later shipment to JAGS:
dataList = list(
    x = vocab_voc,
    y = rt_median_voc,
    n_correct = n_correct_voc,
    n_trials = n_trials_voc
)

# Specify initial values
myinits <-  list(list(phi = 0.75, z = round(runif(length(n_trials_voc)))))

modelString = "
    # Standardize the data:
    data {
        Ntotal <- length(y)
        xm <- mean(x)
        ym <- mean(y)
        xsd <- sd(x)
        ysd <- sd(y)
        for ( i in 1:length(y) ) {
            zx[i] <- ( x[i] - xm ) / xsd
            zy[i] <- ( y[i] - ym ) / ysd
        }
    }
    # Specify the model for standardized data:
    model {
        # Latent Mixture Model
        for (i in 1:Ntotal) {
            z[i] ~ dbern(0.75)
        }
        # First Group Guesses
        psi <- 0.5
        # Second Group Has Some Unknown Greater Rate Of Success
        phi ~ dunif(0.5,1)
        # Data Follow Binomial With Rate Given By Each Persons Group Assignment
        for (i in 1:Ntotal){
            theta[i] <- equals(z[i],0)*psi+equals(z[i],1)*phi
            n_correct[i] ~ dbin(theta[i],n_trials[i])
        }
        # Bayesian Linear Regression
        for ( i in 1:Ntotal ) {
            zbeta0[i] <- equals(z[i],0)*nuisance_beta0+equals(z[i],1)*true_beta0
            zbeta1[i] <- equals(z[i],0)*nuisance_beta1+equals(z[i],1)*true_beta1
            zsigma[i] <- equals(z[i],0)*nuisance_sigma+equals(z[i],1)*true_sigma
           # zy[i] ~ dt( zbeta0[i] + zbeta1[i] * zx[i] , 1/(zsigma[i])^2 , nu )
           zy[i] ~ dnorm( zbeta0[i] + zbeta1[i] * zx[i] , 1/(zsigma[i])^2 )
        }
        # Priors vague on standardized scale set for both linear regressions
        true_beta0 ~ dnorm( 0 , 1/(10)^2 )
        nuisance_beta0 ~ dnorm( 0 , 1/(10)^2 )  
        true_beta1 ~ dnorm( 0 , 1/(10)^2 )
        nuisance_beta1 ~ dnorm( 0 , 1/(10)^2 )
        true_sigma ~ dunif( 1.0E-3 , 1.0E+3 )
        nuisance_sigma ~ dunif( 1.0E-3 , 1.0E+3 )
        nu <- nuMinusOne+1
        nuMinusOne ~ dexp(1/29.0)
        # Transform to original scale:
        beta1 <- true_beta1 * ysd / xsd  
        beta0 <- true_beta0 * ysd  + ym - true_beta1 * xm * ysd / xsd 
        sigma <- zsigma * ysd
    }"

writeLines( modelString , con="rt-voc.txt" )

###### RUN THE CHAINS
parameters = c( "beta0" ,  "beta1" , "sigma", 
                "zbeta0" , "zbeta1" , "zsigma", "nu" ,
                "phi" , "z")
adaptSteps = 500  # Number of steps to "tune" the samplers
burnInSteps = 1000
nChains = 1 
thinSteps = 1
numSavedSteps = 20000
nIter = ceiling( ( numSavedSteps * thinSteps ) / nChains )

####### Create, initialize, and adapt the model:
jagsModel = jags.model( "rt-voc.txt" , data=dataList , inits=myinits , 
                        n.chains=nChains , n.adapt=adaptSteps )

# Burn-in:
update( jagsModel , n.iter=burnInSteps )

# Get samples
samples <- jags(data = dataList, parameters.to.save = parameters,
                model.file="rt-voc.txt", n.chains=1, n.iter=10000, 
                n.thin=1, DIC=T)

df_rt_voc <- data.frame(beta0 = samples$BUGSoutput$sims.list$beta0,
                 beta1 = samples$BUGSoutput$sims.list$beta1,
                 phi = samples$BUGSoutput$sims.list$phi,
                 z = samples$BUGSoutput$sims.list$z)

##### Visualize parameter values

## get HDI
HDI_rt_voc <- HDIofMCMC( df_rt_voc$beta1 , credMass = 0.95 )

## Visualize posterior distribution on slope of regression line
post_mode <- round(find_mode(df_rt_voc$beta1), 3)

ggplot(aes(x = beta1), data = df_rt_voc) + 
    #geom_histogram() +
    geom_density() +
    #annotate("text", x = 0, y = 200, label = paste("Mode =", post_mode), size=5) +
    geom_vline(xintercept = HDI_rt_voc[1], color = 'red') +
    geom_vline(xintercept = HDI_rt_voc[2], color = 'red')

#### Plot MAP regression line with scatter plot
slope_rt_voc <- samples$BUGSoutput$median$beta1
intercept_rt_voc <- samples$BUGSoutput$median$beta0

# grab 20 samples from posterior 
post_20_rt_voc <- df_rt_voc %>% 
    select(beta0, beta1) %>% 
    sample_n(., size = 20)

ggplot(data = data_voc, aes(signs_produced, median_rt)) +
    geom_abline(intercept = intercept_rt_voc, slope = slope_rt_voc, size = 2) +
    geom_abline(data = post_20_rt_voc, aes(intercept = beta0, slope = beta1), size = 0.2, alpha = 0.5) +
    geom_point(color = "black", size = 5.5) +
    geom_point(color = "grey50", size = 4) +
    ylab("Reaction Time (ms)") +
    xlab("Signs Produced") +
    # annotate("text", x = 65, y = 1650, label = paste("MAP =", post_mode), size=8) +
    ggtitle(bquote(list(alpha==.(as.character(round(intercept_rt_voc, 2))), 
                        beta==.(as.character(round(slope_rt_voc, 2))),
                        HDI=.(paste("95% HDI [", round(HDI_rt_voc[1],2) , "," , round(HDI_rt_voc[2], 2), "]"))))) +
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
    x = age,
    y = rt_median,
    n_correct = n_correct,
    n_trials = n_trials
)

# Specify initial values
myinits <-  list(list(phi = 0.75, z = round(runif(length(n_trials)))))

modelString = "
# Standardize the data:
data {
Ntotal <- length(y)
xm <- mean(x)
ym <- mean(y)
xsd <- sd(x)
ysd <- sd(y)
for ( i in 1:length(y) ) {
zx[i] <- ( x[i] - xm ) / xsd
zy[i] <- ( y[i] - ym ) / ysd
}
}
# Specify the model for standardized data:
model {
# Latent Mixture Model
for (i in 1:Ntotal) {
z[i] ~ dbern(0.5)
}
# First Group Guesses
psi <- 0.5
# Second Group Has Some Unknown Greater Rate Of Success
phi ~ dunif(0.5,1)
# Data Follow Binomial With Rate Given By Each Persons Group Assignment
for (i in 1:Ntotal){
theta[i] <- equals(z[i],0)*psi+equals(z[i],1)*phi
n_correct[i] ~ dbin(theta[i],n_trials[i])
}
# Bayesian Linear Regression
for ( i in 1:Ntotal ) {
zbeta0[i] <- equals(z[i],0)*nuisance_beta0+equals(z[i],1)*true_beta0
zbeta1[i] <- equals(z[i],0)*nuisance_beta1+equals(z[i],1)*true_beta1
zsigma[i] <- equals(z[i],0)*nuisance_sigma+equals(z[i],1)*true_sigma
# zy[i] ~ dt( zbeta0[i] + zbeta1[i] * zx[i] , 1/(zsigma[i])^2 , nu )
zy[i] ~ dnorm( zbeta0[i] + zbeta1[i] * zx[i] , 1/(zsigma[i])^2 )
}
# Priors vague on standardized scale set for both linear regressions
true_beta0 ~ dnorm( 0 , 1/(10)^2 )
nuisance_beta0 ~ dnorm( 0 , 1/(10)^2 )  
true_beta1 ~ dnorm( 0 , 1/(10)^2 )
nuisance_beta1 ~ dnorm( 0 , 1/(10)^2 )
true_sigma ~ dunif( 1.0E-3 , 1.0E+3 )
nuisance_sigma ~ dunif( 1.0E-3 , 1.0E+3 )
nu <- nuMinusOne+1
nuMinusOne ~ dexp(1/29.0)
# Transform to original scale:
beta1 <- true_beta1 * ysd / xsd  
beta0 <- true_beta0 * ysd  + ym - true_beta1 * xm * ysd / xsd 
sigma <- zsigma * ysd
}"

writeLines( modelString , con="rt-age.txt" )

###### RUN THE CHAINS
parameters = c( "beta0" ,  "beta1" , "sigma", 
                "zbeta0" , "zbeta1" , "zsigma", "nu" ,
                "phi" , "z")
adaptSteps = 500  # Number of steps to "tune" the samplers
burnInSteps = 1000
nChains = 1 
thinSteps = 1
numSavedSteps = 20000
nIter = ceiling( ( numSavedSteps * thinSteps ) / nChains )

####### Create, initialize, and adapt the model:
jagsModel = jags.model( "rt-age.txt" , data=dataList , inits=myinits , 
                        n.chains=nChains , n.adapt=adaptSteps )

# Burn-in:
update( jagsModel , n.iter=burnInSteps )

# Get samples
samples <- jags(data = dataList, parameters.to.save = parameters,
                model.file="rt-voc.txt", n.chains=1, n.iter=10000, 
                n.thin=1, DIC=T)

df_rt_age <- data.frame(beta0 = samples$BUGSoutput$sims.list$beta0,
                        beta1 = samples$BUGSoutput$sims.list$beta1,
                        phi = samples$BUGSoutput$sims.list$phi,
                        z = samples$BUGSoutput$sims.list$z)

##### Visualize parameter values

## get HDI
HDI_rt_age <- HDIofMCMC( df_rt_age$beta1 , credMass = 0.95 )

## Visualize posterior distribution on slope of regression line
post_mode <- round(find_mode(df_rt_age$beta1), 3)

ggplot(aes(x = beta1), data = df_rt_age) + 
    #geom_histogram() +
    geom_density() +
    #annotate("text", x = 0, y = 200, label = paste("Mode =", post_mode), size=5) +
    geom_vline(xintercept = HDI_rt_age[1], color = 'red') +
    geom_vline(xintercept = HDI_rt_age[2], color = 'red')

#### Plot MAP regression line with scatter plot
slope_rt_age <- samples$BUGSoutput$median$beta1
intercept_rt_age <- samples$BUGSoutput$median$beta0

# grab 20 samples from posterior 
post_20_rt_age <- df_rt_age %>% 
    select(beta0, beta1) %>% 
    sample_n(., size = 20)

ggplot(data = data, aes(age_peek_months, median_rt)) +
    geom_abline(intercept = intercept_rt_age, slope = slope_rt_age, size = 2) +
    geom_abline(data = post_20_rt_age, aes(intercept = beta0, slope = beta1), size = 0.2, alpha = 0.5) +
    geom_point(color = "black", size = 5.5) +
    geom_point(color = "grey50", size = 4) +
    ylab("Reaction Time (ms)") +
    xlab("Child's Age (months)") +
    # annotate("text", x = 65, y = 1650, label = paste("MAP =", post_mode), size=8) +
    ggtitle(bquote(list(alpha==.(as.character(round(intercept_rt_voc, 2))), 
                        beta==.(as.character(round(slope_rt_voc, 2))),
                        HDI=.(paste("95% HDI [", round(HDI_rt_age[1],2) , "," , round(HDI_rt_age[2], 2), "]"))))) +
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
colnames(df)[which(names(df)=="z.1"):which(names(df)=="z.28")] <- as.character(data_voc$Sub.Num)

# melt data frame
df_mixture <- df %>% dplyr::select(-beta0, -beta1)
df_melt <- reshape::melt(df_mixture[,2:length(df_mixture)], variable.name = "Sub.Num", value.name = "group_membership")

# change factor label
df_melt %<>% mutate(group_membership_factor = factor(df_melt$value, labels = c("G", "K"))) 

a <- qplot(data=melt(df$phi),x=value,geom='histogram',binwidth=0.008)+
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