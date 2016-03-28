######## Sample from prior predictive ########

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
    dplyr::select(Sub.Num, age_peek_months, signs_produced, C_T_count,
                  total_trials_shifting, mean_correct_rt, median_rt = median_ct_rt,
                  C_T_prop, prop_looking, mean_prop_looking_TD)

#### Model 
# Here we don't send the data to jags, so we can sample from the prior

dataList_age = list(
    # y = acc ,
    x = (data$age_peek_months - mean(data$age_peek_months)) / sd(data$age_peek_months) ,
    Ntotal = length(data$age_peek_months)
)

modelString = "
# Specify the model for standardized data:
model {
for ( i in 1:Ntotal ) {
zy[i] ~ dnorm( zbeta0 + zbeta1 * x[i] , 1/(zsigma)^2 )
}
# Priors vague on standardized scale:
zbeta0 ~ dnorm( 0 , 1/(10)^2 )  
zbeta1 ~ dnorm(0, 1.0E-2) I(0, )
zsigma ~ dunif( 1.0E-3 , 1.0E+3 )
}
" 
writeLines( modelString , con="acc_age_prior_predictive.txt" )

###### RUN THE CHAINS
parameters = c("zbeta0" , "zbeta1" , "zsigma" )
adaptSteps = 500  # Number of steps to "tune" the samplers
burnInSteps = 1000
nChains = 1 
thinSteps = 1
numSavedSteps = 20000
nIter = ceiling( ( numSavedSteps * thinSteps ) / nChains )

####### Create, initialize, and adapt the model:
jagsModel <- jags.model( "acc_age_prior_predictive.txt" , data=dataList_age , inits=NULL , 
                         n.chains=nChains , n.adapt=adaptSteps )

# Burn-in:
update( jagsModel , n.iter=burnInSteps )

# Get samples
samples <- jags(data = dataList_age, parameters.to.save = parameters,
                model.file="acc_age_prior_predictive.txt", n.chains=1, n.iter=10000, 
                n.thin=1, DIC=F)

df <- data.frame(beta0 = samples$BUGSoutput$sims.list$zbeta0,
                 beta1 = samples$BUGSoutput$sims.list$zbeta1)

##### Visualize parameter values
qplot(df$beta0, df$beta1)

## get MAP and HDI
post_mode <- round(find_mode(df$beta1_d), 3)
HDI <- HDIofMCMC( df$beta1 , credMass = 0.95 )

## Visualize posterior distribution on slope of regression line
ggplot(aes(x = beta1), data = df) + 
    geom_density() +
    geom_vline(xintercept = HDI[1], color = 'red') +
    geom_vline(xintercept = HDI[2], color = 'red')

#### Plot MAP regression line with scatter plot
slope_acc_age <- samples$BUGSoutput$median$zbeta1
intercept_acc_age <- samples$BUGSoutput$median$zbeta0

# grab 20 samples from posterior 
post_20 <- dplyr::sample_n(df, size = 100)

ggplot(data = df, aes(age_peek_months, mean_prop_looking_TD)) +
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