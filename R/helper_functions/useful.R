library(grid)
library(ggplot2)
library(bootstrap)
library(lme4)
library(stringr)
library(plotrix)
library(reshape2)
library(plyr)
library(car)
library(ggm)
library(psych)
library(tidyr)
library(dplyr)

## add some style elements for ggplot2
theme_set(theme_classic())

## standard error of the mean
sem <- function (x) {
  sd(x,na.rm=TRUE) / sqrt(length(x))
}

## NA functions
na.mean <- function(x) {mean(x,na.rm=T)}
na.median <- function(x) {median(x,na.rm=T)}
na.sum <- function(x) {sum(x,na.rm=T)}
na.sd <- function(x) {sd(x,na.rm=T)}

## convert to number
to.n <- function(x) {
  as.numeric(as.character(x))
}

## inverse logistic
inv.logit <- function (x) {
  exp(x) / (1 + exp(x)) 
}

## number of unique subs
n.unique <- function (x) {
  length(unique(x))
}

## for bootstrapping 95% confidence intervals
theta <- function(x,xdata,na.rm=T) {mean(xdata[x],na.rm=na.rm)}
ci.low <- function(x,na.rm=T) {
  mean(x,na.rm=na.rm) - quantile(bootstrap(1:length(x),1000,theta,x,na.rm=na.rm)$thetastar,.025,na.rm=na.rm)}
ci.high <- function(x,na.rm=T) {
  quantile(bootstrap(1:length(x),1000,theta,x,na.rm=na.rm)$thetastar,.975,na.rm=na.rm) - mean(x,na.rm=na.rm)}

## for basic plots, add linear models with correlations
lm.txt <- function (p1,p2,x=7.5,yoff=.05,lt=2,c="black",data=data)
{
  l <- lm(p2 ~ p1)
  regLine(l,lty=lt,col=c)
  cl <- coef(l)
  text(x,cl[1] + cl[2] * x + yoff,
       paste("r = ",sprintf("%2.2f",sqrt(summary(l)$r.squared)),
            getstars(anova(l)$"Pr(>F)"[1]),sep=""),
       xpd="n")
}

## get stars for significance testing
getstars <- function(x) {
  if (x > .1) {return("")}
  if (x < .001) {return("***")}
  if (x < .01) {return("**")}
  if (x < .05) {return("*")}
}

## Multiple plot function
# note, from internet. 
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

# anonymize subject ids by giving them a value 1:num_subjects
anonymize.sids <- function(df, subject_column_label) {
  subj_col = which(names(df) == subject_column_label) # get workerid column index
  temp <- data.frame(workerid = unique(df[,subj_col])) # make new df of unique workerids
  temp$subid <- 1:length(unique(df[,subj_col])) # make list of subids
  index <- match(df[,subj_col], temp$workerid) 
  df$subids <- temp$subid[index]
  df[,subj_col] <- NULL 
  df$subids  = as.factor(df$subids)
  return(df)
}

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

HDIofMCMC = function( sampleVec , credMass=0.95 ) {
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

find_mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


# function to run simulations over multiple windows
# inputs are a list of the data frames (with the ss data), predictor variable name,
# outcome variable name, and the name of the text file with the JAGs model code
# returns a single, tidy data frame with simulations for each window
window_simulations_fun <- function(list_of_dfs, predictor, outcome, model_type, prior_list) {
  df_final <- data.frame()
  counter <- 0
  
  # check what predictors are to get which variables to unstandardize
  if ( predictor == "age.s") {
    unstandardized_pred <- "age_peek_months"
  } else if ( predictor == "voc.s" ) {
    unstandardized_pred <- "signs_produced"
  } else if ( predictor == "acc.s") {
    unstandardized_pred <- "mean_prop_looking_TD" 
  } else {
    unstandardized_pred <- "median_rt"
  }
  
  if ( !(model_type %in% c("accuracy", "rt", "rt_no_mix")) ) {
    stop("Incorrect model type specification")
  }
  
  for (df in list_of_dfs) {
    for (prior in prior_list) {
      counter <- counter + 1
      window_name <- unique(df$window)
      prior_name <- prior
      model_path <- paste("models/", model_type, "_model.txt", sep = "")
      
      dataList <- list(
        y = unlist(df[outcome]),
        x = unlist(df[predictor]),
        Ntotal = length(unlist(df[predictor])),
        n_correct = unlist(df["C_T_count"]),
        n_trials = unlist(df["total_trials_shifting"]),
        prior = prior_name
      )
      
      ## Draw samples from the prior and posterior for given window and a given prior
      samples <- jags(data = dataList, parameters.to.save = parameters,
                      model.file = model_path, n.chains=nChains, 
                      n.iter=nIter, n.burnin = burnInSteps,
                      n.thin=thinSteps, DIC=F)
      
      ## Grab just the parameters we care about
      if ( model_type == "accuracy") {
        
        df_tmp <- data.frame(window = window_name,
                             slope_prior = prior_name,
                             zbeta0 = samples$BUGSoutput$sims.list$zbeta0, 
                             zbeta1 = samples$BUGSoutput$sims.list$zbeta1,
                             zbeta1_prior = samples$BUGSoutput$sims.list$zbeta1_prior,
                             zbeta0_prior = samples$BUGSoutput$sims.list$zbeta0_prior,
                             stringsAsFactors = F)
        
        ## Unstandardized variable
        unstandardized_outcome <- "mean_prop_looking_TD"
      } else if ( model_type == "rt" ) {
        
        df_tmp <- data.frame(window = window_name,
                             slope_prior = prior_name,
                             zbeta0 = samples$BUGSoutput$sims.list$true_beta0,
                             zbeta1 = samples$BUGSoutput$sims.list$true_beta1,
                             zbeta0_prior = samples$BUGSoutput$sims.list$true_beta0_prior,
                             zbeta1_prior = samples$BUGSoutput$sims.list$true_beta1_prior,
                             z = samples$BUGSoutput$sims.list$z,
                             stringsAsFactors = F)
        ## Unstandardized variable
        unstandardized_outcome <- "median_rt"
        
      } else if (model_type == "rt_no_mix") {
        
        df_tmp <- data.frame(window = window_name,
                             slope_prior = prior_name,
                             zbeta0 = samples$BUGSoutput$sims.list$zbeta0, 
                             zbeta1 = samples$BUGSoutput$sims.list$zbeta1,
                             zbeta1_prior = samples$BUGSoutput$sims.list$zbeta1_prior,
                             zbeta0_prior = samples$BUGSoutput$sims.list$zbeta0_prior,
                             stringsAsFactors = F)
        ## Unstandardized variable
        unstandardized_outcome <- "median_rt"
      }
      
      # standardize the data
      sd_unstand_out <- sd(unlist(df[unstandardized_outcome]))
      mean_unstand_out <- mean(unlist(df[unstandardized_outcome]))
      sd_unstand_pred <- sd(unlist(df[unstandardized_pred]))
      mean_unstand_pred <- mean(unlist(df[unstandardized_pred]))
      
      df_tmp %<>%  
        mutate(beta1 = zbeta1 * sd_unstand_out / sd_unstand_pred,
               beta1_prior = zbeta1_prior * sd_unstand_out / sd_unstand_pred,
               beta0 = zbeta0 * sd_unstand_out + mean_unstand_out - 
                 (zbeta1 * mean_unstand_pred * sd_unstand_out) / sd_unstand_pred,
               beta0_prior = zbeta0_prior * sd_unstand_out + mean_unstand_out - 
                 (zbeta1_prior * mean_unstand_pred * sd_unstand_out) / sd_unstand_pred,
               iteration = seq(1, nrow(.)))
      # store in final data frame
      df_final <- rbind(df_final, df_tmp)
    }
  }
  return(df_final)
} 

# Compute p(D) for posterior and prior: Savage-Dickey Method to get Bayes Factor
# Fits a density using splines to approx. log-density
# uses 1997 knot and deletion algorithm
compute_bf <- function (df, posterior_column_name, prior_column_name, bound) {
  if(bound == "upper") {
    fit.posterior <- logspline(df[posterior_column_name], ubound = 0)
    posterior <- dlogspline(0, fit.posterior) # pdf @ beta=0
    # get prior density using logspline
    fit.prior <- logspline(df[prior_column_name], ubound = 0)
    prior <- dlogspline(0, fit.prior) # pdf @ beta=0
    bf_rt_age <- posterior / prior
    1/bf_rt_age
  } else {
    fit.posterior <- logspline(df[posterior_column_name], lbound = 0)
    posterior <- dlogspline(0, fit.posterior) # pdf @ beta=0
    # get prior density using logspline
    fit.prior <- logspline(df[prior_column_name], lbound = 0)
    prior <- dlogspline(0, fit.prior) # pdf @ beta=0
    bf_rt_age <- posterior / prior
    1/bf_rt_age
  }
}

#### This is a slightly different version of the previous window simulation function
#### I made this version to do the linear regression with and without the gussing model 
#### Bad programming practices here in that I'm copying a lot of code, but it works

window_simulations_fun_no_mix <- function(list_of_dfs, predictor, outcome, 
                                          model_type, prior_list, path) {
  df_final <- data.frame()
  counter <- 0
  
  # check what predictors are to get which variables to unstandardize
  if ( predictor == "age.s") {
    unstandardized_pred <- "age_peek_months"
  } else if ( predictor == "voc.s" ) {
    unstandardized_pred <- "signs_produced"
  } else if ( predictor == "acc.s") {
    unstandardized_pred <- "mean_prop_looking_TD" 
  } else {
    unstandardized_pred <- "median_rt"
  }
  
  if ( !(model_type %in% c("accuracy", "rt", "rt_no_mix")) ) {
    stop("Incorrect model type specification")
  }
  
  for (df in list_of_dfs) {
    for (prior in prior_list) {
      counter <- counter + 1
      window_name <- unique(df$window)
      prior_name <- prior
      model_path <- paste(path, model_type, "_model.txt", sep = "")
      
      dataList <- list(
        y = unlist(df[outcome]),
        x = unlist(df[predictor]),
        Ntotal = length(unlist(df[predictor])),
        n_correct = unlist(df["C_T_count"]),
        n_trials = unlist(df["total_trials_shifting"]),
        prior = prior_name
      )
      
      ## Draw samples from the prior and posterior for given window and a given prior
      samples <- jags(data = dataList, parameters.to.save = parameters,
                      model.file = model_path, n.chains=nChains, 
                      n.iter=nIter, n.burnin = burnInSteps,
                      n.thin=thinSteps, DIC=F)
      
      ## Grab just the parameters we care about
      if ( model_type == "accuracy") {
        
        df_tmp <- data.frame(window = window_name,
                             slope_prior = prior_name,
                             zbeta0 = samples$BUGSoutput$sims.list$zbeta0, 
                             zbeta1 = samples$BUGSoutput$sims.list$zbeta1,
                             zbeta1_prior = samples$BUGSoutput$sims.list$zbeta1_prior,
                             zbeta0_prior = samples$BUGSoutput$sims.list$zbeta0_prior,
                             stringsAsFactors = F)
        
        ## Unstandardized variable
        unstandardized_outcome <- "mean_prop_looking_TD"
      } else if ( model_type == "rt" ) {
        
        df_tmp <- data.frame(window = window_name,
                             slope_prior = prior_name,
                             zbeta0 = samples$BUGSoutput$sims.list$true_beta0,
                             zbeta1 = samples$BUGSoutput$sims.list$true_beta1,
                             zbeta0_prior = samples$BUGSoutput$sims.list$true_beta0_prior,
                             zbeta1_prior = samples$BUGSoutput$sims.list$true_beta1_prior,
                             z = samples$BUGSoutput$sims.list$z,
                             stringsAsFactors = F)
        ## Unstandardized variable
        unstandardized_outcome <- "median_rt"
        
      } else if (model_type == "rt_no_mix") {
        
        df_tmp <- data.frame(window = window_name,
                             slope_prior = prior_name,
                             zbeta0 = samples$BUGSoutput$sims.list$zbeta0, 
                             zbeta1 = samples$BUGSoutput$sims.list$zbeta1,
                             zbeta1_prior = samples$BUGSoutput$sims.list$zbeta1_prior,
                             zbeta0_prior = samples$BUGSoutput$sims.list$zbeta0_prior,
                             stringsAsFactors = F)
        ## Unstandardized variable
        unstandardized_outcome <- "median_rt"
      }
      
      # standardize the data
      sd_unstand_out <- sd(unlist(df[unstandardized_outcome]))
      mean_unstand_out <- mean(unlist(df[unstandardized_outcome]))
      sd_unstand_pred <- sd(unlist(df[unstandardized_pred]))
      mean_unstand_pred <- mean(unlist(df[unstandardized_pred]))
      
      df_tmp %<>%  
        mutate(beta1 = zbeta1 * sd_unstand_out / sd_unstand_pred,
               beta1_prior = zbeta1_prior * sd_unstand_out / sd_unstand_pred,
               beta0 = zbeta0 * sd_unstand_out + mean_unstand_out - 
                 (zbeta1 * mean_unstand_pred * sd_unstand_out) / sd_unstand_pred,
               beta0_prior = zbeta0_prior * sd_unstand_out + mean_unstand_out - 
                 (zbeta1_prior * mean_unstand_pred * sd_unstand_out) / sd_unstand_pred,
               iteration = seq(1, nrow(.)))
      # store in final data frame
      df_final <- rbind(df_final, df_tmp)
    }
  }
  return(df_final)
} 