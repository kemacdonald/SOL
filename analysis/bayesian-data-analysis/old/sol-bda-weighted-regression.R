
############ Regression using posterior guessing probability as a weight ##############

## summarize guessing posterior probability for each participant
ss_post <- df_melt %>%
    dplyr::rename(Sub.Num = variable) %>% 
    group_by(Sub.Num) %>% 
    summarise(guessing_count = n() - sum(value),
              post_prob_guessing = 1 - sum(value) / n()) 

## join with data and standardize
data <- left_join(data, ss_post, by = "Sub.Num") %>% 
    mutate(post_prob_guessing.s = (post_prob_guessing - mean(post_prob_guessing)) / sd(post_prob_guessing),
           guessing_div_mean = post_prob_guessing / mean(post_prob_guessing))

# grab weights
weights <- data %>% 
    filter(is.na(signs_produced) == F) %>% # include or remove kid who doesn't have cdi 
    select(post_prob_guessing.s) 

weights <- weights[['post_prob_guessing.s']] # convert weights from df to vector
weights <- weights + .000001

dataList = list(
    x = vocab_voc,
    y = rt_median_voc,
    w = weights
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
zw[i] <- w[i]
}
}
# Specify the model for standardized data:
model {
for ( i in 1:Ntotal ) {
zy[i] ~ dt( zbeta0 + zbeta1 * zx[i] , 1/(zw[i] * zsigma)^2 , nu )
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
writeLines( modelString , con="acc_vocab_weighted_model.txt" )

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
jagsModel <- jags.model( "acc_vocab_weighted_model.txt" , data=dataList , inits=NULL , 
                         n.chains=nChains , n.adapt=adaptSteps )

# Burn-in:
update( jagsModel , n.iter=burnInSteps )

# Get samples
samples <- jags(data = dataList, parameters.to.save = parameters,
                model.file="acc_vocab_weighted_model.txt", n.chains=1, n.iter=10000, 
                n.thin=1, DIC=T)

df_rt_voc <- data.frame(beta0 = samples$BUGSoutput$sims.list$beta0,
                        beta1 = samples$BUGSoutput$sims.list$beta1)

##### Visualize parameter values

## get MAP and HDI
post_mode <- round(find_mode(df_rt_voc$beta1), 3)
HDI <- HDIofMCMC( df_rt_voc$beta1 , credMass = 0.95 )

## Visualize posterior distribution on slope of regression line
ggplot(aes(x = beta1), data = df_rt_voc) + 
    geom_density() +
    geom_vline(xintercept = HDI[1], color = 'red') +
    geom_vline(xintercept = HDI[2], color = 'red')