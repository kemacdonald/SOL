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

HDIofMCMC = function( sampleVec , credMass=0.90 ) {
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
                  total_trials_shifting, 
                  mean_correct_rt, 
                  median_rt = median_ct_rt,
                  C_T_prop, 
                  prop_looking)

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
    y = rt_median_voc,
    Nx = dim(x_matrix)[2] ,
    Ntotal = dim(x_matrix)[1]
)


######## THE ACC MODEL

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
    for ( i in 1:Ntotal ) {
    zy[i] ~ dt( zbeta0 + sum( zbeta[1:Nx] * zx[i,1:Nx] ) , 1/(zsigma)^2 , nu )
    }
    # Priors vague on standardized scale:
    zbeta0 ~ dnorm( 0 , 1/2^2 )  
    for ( j in 1:Nx ) {
    zbeta[j] ~ dnorm( 0 , 1/2^2 )
    }
    zsigma ~ dunif( 1.0E-5 , 1.0E+1 )
    nu <- nuMinusOne+1
    nuMinusOne ~ dexp(1/29.0)
    # Transform to original scale:
    beta[1:Nx] <- ( zbeta[1:Nx] / xsd[1:Nx] )*ysd
    beta0 <- zbeta0*ysd  + ym - sum( zbeta[1:Nx] * xm[1:Nx] / xsd[1:Nx] )*ysd
    sigma <- zsigma*ysd
    }"
writeLines( modelString , con="TEMPmodel.txt" )

###### RUN THE CHAINS

parameters = c( "beta0" ,  "beta" ,  "sigma", 
                "zbeta0" , "zsigma", "nu")
adaptSteps = 500  # Number of steps to "tune" the samplers
burnInSteps = 1000
nChains = 1 
thinSteps = 1
numSavedSteps = 20000
nIter = ceiling( ( numSavedSteps * thinSteps ) / nChains )

####### Create, initialize, and adapt the model:

jagsModel = jags.model( "TEMPmodel.txt" , data=dataList , inits=NULL , 
                        n.chains=nChains , n.adapt=adaptSteps )

# Burn-in:
update( jagsModel , n.iter=burnInSteps )

# Get samples
samples <- jags(data = dataList, parameters.to.save = parameters,
                model.file="TEMPmodel.txt", n.chains=1, n.iter=10000, 
                n.thin=1, DIC=T)

df <- data.frame(beta0 = samples$BUGSoutput$sims.list$beta0,
                 ageBeta = samples$BUGSoutput$sims.list$beta[,1],
                 vocabBeta = samples$BUGSoutput$sims.list$beta[,2])

###### Visualize posterior distribution on slope of regression line

### Age Coef
post_mode <- round(find_mode(df$ageBeta), 3)
HDI <- HDIofMCMC( df$ageBeta , credMass = 0.95 )

ggplot(aes(x = ageBeta), data = df) + 
    geom_density() +
    geom_vline(xintercept = HDI[1], color = 'red') +
    geom_vline(xintercept = HDI[2], color = 'red')

## Vocab coef
post_mode <- round(find_mode(df$vocabBeta), 3)
HDI <- HDIofMCMC( df$vocabBeta , credMass = 0.95 )

ggplot(aes(x = vocabBeta), data = df) + 
    geom_density() +
    geom_vline(xintercept = HDI[1], color = 'red') +
    geom_vline(xintercept = HDI[2], color = 'red')

############# RT Model: Guessing Model + Multiple Regression ##################

dataList = list(
    x = x_matrix,
    y = rt_median_voc,
    Nx = dim(x_matrix)[2] ,
    Ntotal = dim(x_matrix)[1],
    n_correct = n_correct_voc,
    n_trials = n_trials_voc
)

# Specify initial values
myinits <-  list(list(phi = 0.75, z = round(runif(length(n_trials_voc)))))

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
         # Specify the Latent Mixture Model:
        model {
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
         # Specify the linear regression model for standardized data:
            for ( i in 1:Ntotal ) {
                # choose parameters based on output of latent mixture
                zbeta0[i] <- equals(z[i],0)*nuisance_beta0+equals(z[i],1)*true_beta0
                zbeta[i, 1:Nx] <- equals(z[i],0)*nuisance_beta+equals(z[i],1)*true_beta
                # multiple regression 
                zy[i] ~ dt( zbeta0[i] + sum( zbeta[i,1:Nx] * zx[i,1:Nx] ) , 1/zsigma^2 , nu )
            }
            # Priors vague on standardized scale:
            true_beta0 ~ dnorm( 0 , 1/2^2 )
            nuisance_beta0 ~ dnorm( 0 , 1/2^2 )
            for ( j in 1:Nx ) {
                true_beta[j] ~ dnorm( 0 , 1/2^2 )
                nuisance_beta[j] ~ dnorm( 0 , 1/2^2 )
            }
            zsigma ~ dunif( 1.0E-5 , 1.0E+1 )
            nu <- nuMinusOne+1
            nuMinusOne ~ dexp(1/29.0)
            # Transform to original scale:
            beta[1:Nx] <- ( true_beta[1:Nx] / xsd[1:Nx] )*ysd
            beta0 <- true_beta0*ysd  + ym - sum( true_beta[1:Nx] * xm[1:Nx] / xsd[1:Nx] )*ysd
            sigma <- zsigma*ysd
        }"

writeLines( modelString , con="TEMPmodel.txt" )

###### RUN THE CHAINS

parameters = c( "beta0" ,  "beta" ,  "sigma", 
                "zbeta0" , "zbeta" , "zsigma", "nu", 
                "z", "phi")

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

df <- data.frame(beta0 = samples$BUGSoutput$sims.list$beta0,
                 ageBeta = samples$BUGSoutput$sims.list$beta[,1],
                 vocabBeta = samples$BUGSoutput$sims.list$beta[,2],
                 phi = samples$BUGSoutput$sims.list$phi,
                 z = samples$BUGSoutput$sims.list$z)

###### Visualize posterior distribution on slope of regression line

### Age Coef
post_mode <- round(find_mode(df$ageBeta), 3)
HDI <- HDIofMCMC( df$ageBeta , credMass = 0.95 )

ggplot(aes(x = ageBeta), data = df) + 
    geom_density() +
    geom_vline(xintercept = HDI[1], color = 'red') +
    geom_vline(xintercept = HDI[2], color = 'red')

### Vocab coef
post_mode <- round(find_mode(df$vocabBeta), 3)
HDI <- HDIofMCMC( df$vocabBeta , credMass = 0.95 )

ggplot(aes(x = vocabBeta), data = df) + 
    geom_density() +
    geom_vline(xintercept = HDI[1], color = 'red') +
    geom_vline(xintercept = HDI[2], color = 'red')

####### Visualize posterior distribution on group membership and group level success

# grab subject numbers and add to data frame
colnames(df)[which(names(df)=="z.1"):which(names(df)=="z.28")] <- as.character(data_voc$Sub.Num)

# melt data frame
df_mixture <- df %>% dplyr::select(-ageBeta, -vocabBeta, -beta0)
df_melt <- reshape2::melt(df_mixture[,2:length(df_mixture)], 
                          variable.name = "Sub.Num", 
                          value.name = "group_membership")

# change factor label
df_melt %<>% mutate(group_membership_factor = factor(df_melt$group_membership, 
                                                     labels = c("G", "K"))) 

a <- qplot(data= reshape2::melt(df$phi),x=value,geom='histogram',binwidth=0.008)+
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

gridExtra::grid.arrange(a,b,nrow=1, top='Group success rate and individual group membership', 
             widths = c(2, 4))