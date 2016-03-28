############# Fit model using map2stan ###############
library(rethinking)
library(dplyr)
library(ggplot2)
library(rjags)

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


### Fit model
m_acc_age <- map2stan(
    alist(
        mean_prop_looking_TD ~ dnorm(mu , sigma),
        mu <- Intercept + bAge*age_peek_months , 
        Intercept ~ dnorm(0,10),
        bAge ~ dnorm(0, 10),
        sigma ~ dcauchy(0,2)
    ),
    data = select(data_voc, age_peek_months, mean_prop_looking_TD)
)

### Model output
precis(m_acc_age)

### Visualize
post <- as.data.frame(extract.samples(m_acc_age_s))
pairs(post)
ggplot(aes(b_age_peek_months), data = post) + geom_density() 
    
## DIC and WAIC
show(m_acc_age_s)

## Trace plot
plot(m_acc_age_s)

######## NULL Model comparison

m_acc_age_null <- map2stan(
    alist(
        mean_prop_looking_TD ~ dnorm(mu , sigma),
        mu <- Intercept , 
        Intercept ~ dnorm(0,10),
        sigma ~ dcauchy(0,2)
    ),
    data = select(data_voc, mean_prop_looking_TD)
)

m_acc_age_voc <- map2stan(
    alist(
        mean_prop_looking_TD ~ dnorm(mu , sigma),
        mu <- Intercept + bAge*age_peek_months + bVoc*signs_produced, 
        Intercept ~ dnorm(0,10),
        bAge ~ dnorm(0, 10),
        bVoc ~ dnorm(0, 10),
        sigma ~ dcauchy(0,2)
    ),
    data = select(data_voc, mean_prop_looking_TD, age_peek_months, signs_produced)
)

acc_age_models <- rethinking::compare(m_acc_age , m_acc_age_null, m_acc_age_voc)

plot(acc_age_models, SE=T, dSE=T)

plot(coeftab(m_acc_age_null, m_acc_age, m_acc_age_voc))

### visualize parameters multiple regression
post_multi <- as.data.frame(extract.samples(m_acc_age_voc))
plot(m_acc_age_voc)
plot(post_multi)
