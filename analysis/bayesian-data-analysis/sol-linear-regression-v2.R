# Packages
library(rjags)
library(dplyr)
library(ggplot2)
library(magrittr)
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
    dplyr::select(Sub.Num, age_peek_months, signs_produced, C_T_count,
                  total_trials_shifting, mean_correct_rt, median_ct_rt,
                  C_T_prop, prop_looking)

##### Structure the data

y = data$prop_looking
x = data$age_peek_months

# Specify the data in a list, for later shipment to JAGS:

dataList = list(
    x = x ,
    y = y
)

######## THE MODEL

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
    }"

writeLines( modelString , con="TEMPmodel.txt" )

###### RUN THE CHAINS

parameters = c( "beta0" ,  "beta1" ,  "sigma", 
                "zbeta0" , "zbeta1" , "zsigma", "nu" )
adaptSteps = 500  # Number of steps to "tune" the samplers
burnInSteps = 1000
nChains = 4 
thinSteps = 1
numSavedSteps = 20000
nIter = ceiling( ( numSavedSteps * thinSteps ) / nChains )

####### Create, initialize, and adapt the model:

jagsModel = jags.model( "TEMPmodel.txt" , data=dataList , #inits=initsList , 
                        n.chains=nChains , n.adapt=adaptSteps )

# Burn-in:
update( jagsModel , n.iter=burnInSteps )


codaSamples = coda.samples( jagsModel , variable.names=parameters , 
                            n.iter=nIter , thin=thinSteps )

samples <- jags(data = dataList, parameters.to.save = parameters,
                model.file="TEMPmodel.txt", n.chains=1, n.iter=2000, 
                n.thin=1, DIC=T)


df <- data.frame(beta0 = samples$BUGSoutput$sims.list$beta0,
                 beta1 = samples$BUGSoutput$sims.list$beta1)

##### Visualize parameter values

## get HDI

HDI <- HDIofMCMC( df$beta1 , credMass = 0.95 )

## Visualize posterior distribution on slope of regression line

post_mode <- round(find_mode(df$beta1), 3)

ggplot(aes(x = beta1), data = df) + 
    geom_histogram(binwidth = .0005) +
    annotate("text", x = 0, y = 90, label = paste("Mode =", post_mode), size=5) +
    geom_vline(xintercept = HDI[1], color = 'red') +
    geom_vline(xintercept = HDI[2], color = 'red')
