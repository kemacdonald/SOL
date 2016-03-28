######### Bayesian Analysis of Data on Development of Real-Time ASL Procesing Skills ###############

######### Categorical Analyses ##########

# Clear workspace
rm(list=ls())

# Packages
library(reshape)
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

d_all <- d %>% 
    filter(value_cat == "Target" | value_cat == "Distractor") %>% 
    mutate(C_D_count = ifelse(is.na(C_D_count), 0, C_D_count),
           total_trials_shifting = C_T_count + C_D_count,
           Sub.Num = as.character(Sub.Num)) %>% 
    dplyr::select(Sub.Num, age_peek_months, age_group, signs_produced, C_T_count, total_trials_shifting, 
                  mean_correct_rt, median_rt = median_ct_rt,C_T_prop, prop_looking, mean_prop_looking_TD,
                  value_cat, age_group_collapsed, hearing_status_participant)

d %<>%  
    filter(age_group_collapsed == "Kids", value_cat == "Target") %>% 
    mutate(C_D_count = ifelse(is.na(C_D_count), 0, C_D_count),
           total_trials_shifting = C_T_count + C_D_count,
           Sub.Num = as.character(Sub.Num)) %>% 
    dplyr::select(Sub.Num, age_peek_months, signs_produced, C_T_count, total_trials_shifting, 
                  mean_correct_rt, median_rt = median_ct_rt,C_T_prop, prop_looking, mean_prop_looking_TD,
                  age_group_collapsed, hearing_status_participant)

d_voc <- d %>% filter(is.na(signs_produced) == F)

##### Fit models with age group (younger, older, adults) as categorial predictor of Target/Distractor looking
d_model <- d_all %>% 
    mutate(age_group_older = ifelse(age_group == ">= 27 Months", 1, 0),
           age_group_adults = ifelse(age_group == "Adults", 1, 0),
           age_group_clean = ifelse(age_group_older == 1, "Older", 
                                    ifelse(age_group_adults == 1, "Adult", 
                                           "Younger")),
           age_group_adults_collapsed = ifelse(age_group == "Adults", 1, 0),
           hearing_group = ifelse(hearing_status_participant == "hearing", 1, 0))

## Accuracy by hearing status
acc_hearing_stat_cat_targ.m <- map(
    alist(
        mean_prop_looking_TD ~ dnorm(mu, sigma),
        mu <- Intercept + b.Hearing*hearing_group,
        Intercept ~ dnorm(0.6, 10), 
        b.Hearing ~ dnorm(0, 1),
        sigma ~ dunif(0, 10)
    ),
    data = filter(d_model, value_cat == "Target", age_group_collapsed == "Kids")
)

precis(acc_hearing_stat_cat_targ.m)

## Accuracy by age group
acc_age_cat_targ.m <- map(
    alist(
        mean_prop_looking_TD ~ dnorm(mu, sigma),
        mu <- Intercept + b.Older*age_group_older + b.Adults*age_group_adults,
        Intercept ~ dnorm(0.6, 10), 
        b.Older ~ dnorm(0, 1),
        b.Adults ~ dnorm(0, 1),
        sigma ~ dunif(0, 10)
    ),
    data = filter(d_model, value_cat == "Target")
)

## Accuracy by age group
acc_age_cat_targ.m.collapsed <- map(
    alist(
        mean_prop_looking_TD ~ dnorm(mu, sigma),
        mu <- Intercept + b.Adults*age_group_adults,
        Intercept ~ dnorm(0.6, 10), 
        b.Adults ~ dnorm(0, 1),
        sigma ~ dunif(0, 10)
    ),
    data = filter(d_model, value_cat == "Target")
)

## Distractor looking by age group
acc_age_cat_dist.m <- map(
    alist(
        mean_prop_looking_TD ~ dnorm(mu, sigma),
        mu <- Intercept + b.Older*age_group_older + b.Adults*age_group_adults,
        Intercept ~ dnorm(0.6, 10), 
        b.Older ~ dnorm(0, 1),
        b.Adults ~ dnorm(0, 1),
        sigma ~ dunif(0, 10)
    ),
    data = filter(d_model, value_cat == "Distractor")
)

## RT by age group
rt_age_cat.m <- map2stan(
    alist(
        median_rt ~ dnorm(mu, sigma),
        mu <- Intercept + b_Older*age_group_older + b_Adults*age_group_adults,
        Intercept ~ dnorm(1300, 100), 
        b_Older ~ dnorm(0, 500),
        b_Adults ~ dnorm(0, 500),
        sigma ~ dunif(0, 500)
    ),
    data = select(filter(d_model, value_cat == "Target"), median_rt, age_group_older, age_group_adults)
)

## RT by hearing status
rt_hearing_stat.m <- map2stan(
    alist(
        median_rt ~ dnorm(mu, sigma),
        mu <- Intercept + b_Hearing*hearing_group,
        Intercept ~ dnorm(1300, 100), 
        b_Hearing ~ dnorm(0, 500),
        sigma ~ dunif(0, 500)
    ),
    data = select(filter(d_model, value_cat == "Target", age_group_collapsed == "Kids"), 
                  median_rt, hearing_group)
)

## RT by age group collapsed
rt_age_cat_collapsed.m <- map2stan(
    alist(
        median_rt ~ dnorm(mu, sigma),
        mu <- Intercept + b_Adults*age_group_adults_collapsed,
        Intercept ~ dnorm(1300, 100), 
        b_Adults ~ dnorm(0, 500),
        sigma ~ dunif(0, 500)
    ),
    data = select(filter(d_model, value_cat == "Target"), median_rt, age_group_adults_collapsed)
)

#### Get samples so we can explore posterior distributions ########

# target and distractor looking
post_ages.acc <- extract.samples(acc_age_cat_targ.m)
post_ages.acc_collapsed <- extract.samples(acc_age_cat_targ.m.collapsed)
post_ages.dist <- extract.samples(acc_age_cat_dist.m)

# rt models
post_ages.rt <- extract.samples(rt_age_cat.m)
post_ages.rt_collapsed <- extract.samples(rt_age_cat_collapsed.m)

# hearing status
post_hearing_stat.rt <- extract.samples(rt_hearing_stat.m)
post_hearing_stat.acc <- extract.samples(acc_hearing_stat_cat_targ.m)


######### Store values from posterior

# Accuracy: hearing status
mu.hearing <- post_hearing_stat.acc$Intercept
mu.deaf <- post_hearing_stat.acc$b.Hearing + post_hearing_stat.acc$Intercept
post_df_hearing_stat.acc <- data.frame(mu.hearing, mu.deaf)
precis(post_df_hearing_stat.acc)

# Accuracy: compute averages for each group (intercept is youngest kids)
mu.younger <- post_ages.acc$Intercept
mu.older <- post_ages.acc$b.Older + post_ages.acc$Intercept
mu.adults <- post_ages.acc$b.Adults + post_ages.acc$Intercept
post_df.acc <- data.frame(mu.younger, mu.older, mu.adults)
precis(post_df.acc)

# Accuracy: compute averages for each group (intercept is youngest kids) collapsing
mu.kids.collapsed <- post_ages.acc_collapsed$Intercept
mu.adults.collapsed <- post_ages.acc_collapsed$b.Adults + post_ages.acc_collapsed$Intercept
post_df.acc.collapsed <- data.frame(mu.kids.collapsed, mu.adults.collapsed)
precis(post_df.acc.collapsed)

# Distractor: compute averages for each group (intercept is youngest kids)
mu.younger.dist <- post_ages.dist$Intercept
mu.older.dist <- post_ages.dist$b.Older + post_ages.dist$Intercept
mu.adults.dist <- post_ages.dist$b.Adults + post_ages.dist$Intercept
post_df.dist <- data.frame(mu.younger.dist, mu.older.dist, mu.adults.dist)
precis(post_df.dist)

# Merge Target and Distractor looking posterior samples
post_df.dist$looking_type <- "Distractor"
post_df.acc$looking_type <- "Target"
post_df.acc_all <- rbind(post_df.acc, post_df.dist)

# RT: all age groups
mu.younger.rt <- as.numeric(post_ages.rt$Intercept)
mu.older.rt <- as.numeric(post_ages.rt$b_Older + post_ages.rt$Intercept)
mu.adults.rt <- as.numeric(post_ages.rt$b_Adults + post_ages.rt$Intercept)
post_df.rt <- data.frame(mu.younger.rt, mu.older.rt, mu.adults.rt) 
precis(post_df.rt)

# RT: age group collapsing across kids/adults
mu.kids.rt.collapsed <- as.numeric(post_ages.rt_collapsed$Intercept)
mu.adults.rt.collapsed <- as.numeric(post_ages.rt_collapsed$b_Adults + post_ages.rt_collapsed$Intercept)
post_df.rt.collapsed <- data.frame(mu.kids.rt.collapsed, mu.adults.rt.collapsed) 
precis(post_df.rt.collapsed)

# RT: hearing status
mu.hearing.rt <- as.numeric(post_hearing_stat.rt$Intercept)
mu.deaf.rt <- as.numeric(post_hearing_stat.rt$b_Hearing + post_hearing_stat.rt$Intercept)
post_hearing_stat.rt <- data.frame(mu.hearing.rt, mu.deaf.rt) 
precis(post_hearing_stat.rt)

################# Test contrasts #################

# hearing status acc
diff.deaf.hearing <- mu.hearing - mu.deaf
quantile(diff.deaf.hearing, probs = c(0.025, 0.5, 0.975))

# older vs. younger acc
diff.old.young <- mu.younger - mu.older
quantile(diff.old.young, probs = c(0.025, 0.5, 0.975))

# older vs. younger rt
diff.old.young.rt <- mu.younger.rt - mu.older.rt
quantile(diff.old.young.rt, probs = c(0.025, 0.5, 0.975))

# adults vs. kids rt
diff.old.young.rt.collapsed <- mu.kids.rt.collapsed - mu.adults.rt.collapsed
quantile(diff.old.young.rt.collapsed, probs = c(0.025, 0.5, 0.975))

# hearing status rt
diff.deaf.hearing.rt <- mu.hearing.rt - mu.deaf.rt
quantile(diff.deaf.hearing.rt, probs = c(0.025, 0.5, 0.975))

# adults vs. kids acc
diff.old.young.acc.collapsed <- mu.kids.collapsed - mu.adults.collapsed
quantile(diff.old.young.acc.collapsed, probs = c(0.025, 0.5, 0.975))

# target vs. distractor looking at each time point
diff.targ.dist.young <- mu.younger - mu.younger.dist
quantile(diff.targ.dist.young, probs = c(0.025, 0.5, 0.975))

diff.targ.dist.older <- mu.older - mu.older.dist
quantile(diff.targ.dist.older, probs = c(0.025, 0.5, 0.975))

diff.targ.dist.adults <- mu.adults - mu.adults.dist
quantile(diff.targ.dist.adults, probs = c(0.025, 0.5, 0.975))

################# Posterior means and Plots #################

# Accuracy: MAP and HDI for each age group
post_df.acc_all %<>% gather(key = age_group, value = post_sample, mu.younger:mu.adults) 

ms_all <- post_df.acc_all %>%
    group_by(age_group, looking_type) %>% 
    summarise(mean_post = mean(post_sample),
              HDI_lower = HDIofMCMC(post_sample, credMass = 0.95)[1],
              HDI_upper = HDIofMCMC(post_sample, credMass = 0.95)[2]) %>% 
    mutate(HDI_lower = ifelse(HDI_lower <= 0, 0, HDI_lower))


######### Now plot bar graph of MAP and HDI for target/distractor looking ############

# RT: MAP and HDI for each age group
post_df.rt %<>% gather(key = age_group, value = post_sample)

ms_all_rt <- post_df.rt %>%
    group_by(age_group) %>% 
    summarise(mean_post = mean(post_sample),
              HDI_lower = HDIofMCMC(post_sample, credMass = 0.95)[1],
              HDI_upper = HDIofMCMC(post_sample, credMass = 0.95)[2]) 

# Bar Plot
all_acc <- ggplot(aes(x = age_group, y = mean_post, fill = looking_type), data = ms_all) +
    geom_bar(stat = "identity", position="dodge") +
    geom_linerange(aes(ymin = HDI_lower,
                       ymax = HDI_upper),
                   position = position_dodge(width = 0.9)) +
    xlab("Age Group") +
    ylab("Proportion Looking") +
    coord_cartesian(ylim=c(0, 0.9)) +
    scale_fill_brewer(type = "qual", palette = "Set1") +
    guides(fill=guide_legend(title=NULL)) +
    scale_x_discrete(labels=c("Younger", "Older", "Adults")) +
    theme_bw() +
    theme(plot.title = element_text(face = "bold", size = 20),
          axis.title.x = element_text(colour="grey40",size=18,
                                      angle=0,hjust=0.5,vjust=0,face="plain"),
          axis.title.y = element_text(colour="grey40",size=18,
                                      angle=90,hjust=0.5,vjust=0.5,face="plain"),
          axis.text.x = element_text(colour="grey20",size=18,
                                     angle=0,hjust=0.5,vjust=0,face="plain"),
          axis.text.y = element_text(colour="grey20",size=18,
                                     angle=0,hjust=0.5,vjust=0,face="plain"),
          plot.margin = unit(c(0.5,2,1,1), "cm"),
          legend.position = c(0.2, 0.85),
          legend.text = element_text(size = 16)) +
    ggtitle("Accuracy") 

all_rt <- ggplot(aes(x= age_group, y = mean_post, fill = age_group), data = ms_all_rt) +
    geom_bar(stat = "identity") +
    geom_linerange(aes(ymax = HDI_upper, ymin = HDI_lower)) +
    scale_fill_brewer(type = "qual", palette = "Set1") +
    guides(fill=F) +
    xlab("Age Group") +
    ylab("Mean RT (ms)") +
    scale_x_discrete(labels=c("Younger", "Older", "Adults")) +
    theme_bw() +
    theme(plot.title = element_text(face = "bold", size = 20),
          axis.title.x = element_text(colour="grey40",size=18,
                                      angle=0,hjust=0.5,vjust=0,face="plain"),
          axis.title.y = element_text(colour="grey40",size=18,
                                      angle=90,hjust=0.5,vjust=0.5,face="plain"),
          axis.text.x = element_text(colour="grey20",size=18,
                                     angle=0,hjust=0.5,vjust=0,face="plain"),
          axis.text.y = element_text(colour="grey20",size=18,
                                     angle=0,hjust=0.5,vjust=0,face="plain"),
          plot.margin = unit(c(0.5,2,1,1), "cm")) +
    ggtitle("Reaction Time") 

gridExtra::grid.arrange(all_acc, all_rt, ncol = 2)
