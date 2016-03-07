# Clear workspace
rm(list=ls())
# Packages
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

##### Load data

data <- read.csv("../../analysis/eye_movements/sol_ss_all.csv")

data %<>%  
    filter(age_group_collapsed == "Kids", value_cat == "Target") %>% 
    mutate(C_D_count = ifelse(is.na(C_D_count), 0, C_D_count),
           total_trials_shifting = C_T_count + C_D_count,
           Sub.Num = as.character(Sub.Num)) %>% 
    dplyr::select(Sub.Num, age_peek_months, 
                  signs_produced, 
                  C_T_count,
                  total_trials_shifting, mean_correct_rt, 
                  log_median_rt = median_ct_rt.x, 
                  median_rt = median_ct_rt.y,
                  C_T_prop, prop_looking)

data_voc <- data %>% filter(is.na(signs_produced) == F)

##### Structure the data

rt <-data$mean_correct_rt
rt_median <- data$median_rt
n_correct <- data$C_T_count
n_trials <- data$total_trials_shifting  
acc <- data$prop_looking
age <- data$age_peek_months
vocab <- data$signs_produced

# For Vocab models we need to remove one participant for whom we didn't get CDI
rt_voc <-data_voc$mean_correct_rt
rt_median_voc <- data_voc$median_rt
n_correct_voc <- data_voc$C_T_count
n_trials_voc <- data_voc$total_trials_shifting  
acc_voc <- data_voc$prop_looking
age_voc <- data_voc$age_peek_months
vocab_voc <- data_voc$signs_produced

# Create matrix for multiple regression: each column is a predictor and each row is data point 

x_matrix <- cbind(age_voc, vocab_voc)

# Specify the data in a list, for later shipment to JAGS:

dataList = list(
    x = x_matrix,
    y = rt_voc,
    n_correct = n_correct,
    n_trials = n_trials,
    Nx = dim(x_matrix)[2] ,
    Ntotal = dim(x_matrix)[1]
)

# Specify initial values

myinits <-  list(list(phi = 0.75, z = round(runif(length(n_correct)))))

######## THE MODEL

modelString = "
    # Standardize the data:
      data {
        ym <- mean(y)
        ysd <- sd(y)
        for ( i in 1:Ntotal ) {
          zy[i] <- ( y[i] - ym ) / ysd
        }
        for ( j in 1:Nx ) {
          xm[j]  <- mean(x[,j])
          xsd[j] <-   sd(x[,j])
          for ( i in 1:Ntotal ) {
            zx[i,j] <- ( x[i,j] - xm[j] ) / xsd[j]
          }
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
        zsigma[i] <- equals(z[i],0)*nuisance_sigma+equals(z[i],1)*true_sigma
        # need to draw values for both predictors
        for ( j in 1:Nx ) {
            zbeta[i,] <- equals(z[i],0)*nuisance_beta[i,]+equals(z[i],1)*true_beta[,]
        }

        # zy[i] ~ dt( zbeta0[i] + zbeta1[i] * zx[i] , 1/(zsigma[i])^2 , nu )
        zy[i] ~ dt( zbeta0[i] + sum( zbeta[1:Nx] * zx[i,1:Nx] ) , 1/zsigma[i]^2 , nu )
    }
    # Priors vague on standardized scale; set for both linear regressions
    true_beta0 ~ dnorm( 0 , 1/(10)^2 )
    nuisance_beta0 ~ dnorm( 0 , 1/(10)^2 )  
    for ( j in 1:Nx ) {
        true_beta[j] ~ dnorm( 0 , 1/(10)^2 )
        nuisance_beta[j] ~ dnorm( 0 , 1/(10)^2 )
    }
    true_sigma ~ dunif( 1.0E-3 , 1.0E+3 )
    nuisance_sigma ~ dunif( 1.0E-3 , 1.0E+3 )
    nu <- nuMinusOne+1
    nuMinusOne ~ dexp(1/29.0)
    # Transform to original scale:
    beta[1:Nx] <- ( zbeta[1:Nx] / xsd[1:Nx] )*ysd
    beta0 <- zbeta0*ysd  + ym - sum( zbeta[1:Nx] * xm[1:Nx] / xsd[1:Nx] )*ysd
    sigma <- zsigma * ysd
}"

writeLines( modelString , con="TEMPmodel.txt" )

###### RUN THE CHAINS

parameters = c( "beta0" ,  "beta" ,  "sigma", 
                "zbeta0" , "zbeta1" , "zsigma", "nu" ,
                "phi" , "z")
adaptSteps = 500  # Number of steps to "tune" the samplers
burnInSteps = 1000
nChains = 1 
thinSteps = 1
numSavedSteps = 20000
nIter = ceiling( ( numSavedSteps * thinSteps ) / nChains )

####### Create, initialize, and adapt the model:

jagsModel = jags.model( "TEMPmodel.txt" , data=dataList , inits=myinits , 
                        n.chains=nChains , n.adapt=adaptSteps )

# Burn-in:
update( jagsModel , n.iter=burnInSteps )

# Get samples
samples <- jags(data = dataList, parameters.to.save = parameters,
                model.file="TEMPmodel.txt", n.chains=1, n.iter=5000, 
                n.thin=1, DIC=T)


df <- data.frame(beta0 = samples$BUGSoutput$sims.list$true_beta0,
                 beta1 = samples$BUGSoutput$sims.list$true_beta1,
                 phi = samples$BUGSoutput$sims.list$phi,
                 z = samples$BUGSoutput$sims.list$z)

##### Visualize parameter values

## get HDI
HDI <- HDIofMCMC( df$beta1 , credMass = 0.95 )

## Visualize posterior distribution on slope of regression line
post_mode <- round(find_mode(df$beta1), 3)

ggplot(aes(x = beta1), data = df) + 
    geom_histogram() +
    annotate("text", x = 0, y = 90, label = paste("Mode =", post_mode), size=5) +
    geom_vline(xintercept = HDI[1], color = 'red') +
    geom_vline(xintercept = HDI[2], color = 'red')

####### Visualize posterior distribution on group membership and group level success

# grab subject numbers and add to data frame
colnames(df)[which(names(df)=="z.1"):which(names(df)=="z.29")] <- as.character(data$Sub.Num)

# melt data frame
df_mixture <- df %>% dplyr::select(-beta0, -beta1)
df_melt <- melt(df_mixture[,2:length(df_mixture)], variable.name = "Sub.Num", value.name = "group_membership")

# change factor label
df_melt %<>% mutate(group_membership_factor = factor(df_melt$group_membership, 
                                                     labels = c("G", "K"))) 

a <- qplot(data=melt(df$phi),x=value,geom='histogram',binwidth=0.008)+
    theme_bw()+
    xlab('phi')+
    xlim(0,1)

b <- ggplot(data=df_melt,aes(x=group_membership_factor, 
                             fill=group_membership_factor)) +
    geom_bar(stat="count") +
    facet_wrap(~Sub.Num,scales='fixed')+
    guides(fill=F) +
    scale_fill_solarized() +
    theme_bw()+
    xlab('Group Membership') +
    ylab("Count")

grid.arrange(a,b,nrow=1, top='Group success rate and individual group membership', 
             widths = c(2, 4))