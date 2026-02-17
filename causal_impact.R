#Causal Impact models

#Libraries----
library(CausalImpact)
library(dplyr)
library(zoo)
library(coda)

#Import Data----

causalImpact_data <- read.csv("causalImpact_data.csv")
rownames(causalImpact_data) <- causalImpact_data$date


#Define pre and post intervention period----

pre.period <- as.Date(c("2010-01-01", "2012-01-15"))
post.period <- as.Date(c("2012-01-16", "2017-12-31"))

#Data Interpolation----

causalImpact_data$date <- as.Date(causalImpact_data$date)
causalImpact_data_interpolated <- causalImpact_data %>% mutate(across(where(is.numeric),~ na.approx(.x, rule = 2)))

#basic Causal Impact model----

impact <- CausalImpact(causalImpact_data_interpolated, pre.period, post.period)
causalImpact_data_interpolated$date <- causalImpact_data$date
plot(impact)
summary(impact)
summary(impact, "report")

#model with seasonality and higher number of iterations----
impact2 <- CausalImpact(causalImpact_data_interpolated, pre.period, post.period,
                        model.args = list(niter = 5000, nseasons = 52,
                          season.duration = 1,
                          prior.level.sd = 0.01))
plot(impact2)
summary(impact2)
#> summary(impact2)
#Posterior inference {CausalImpact}
#
#Average         Cumulative     
#Actual                   37              11407          
#Prediction (s.d.)        43 (3.1)        13289 (950.5)  
#95% CI                   [37, 49]        [11357, 15150] 
#
#Absolute effect (s.d.)   -6.1 (3.1)      -1882.5 (950.5)
#95% CI                   [-12, 0.16]     [-3743, 49.82] 
#
#Relative effect (s.d.)   -14% (6.3%)     -14% (6.3%)    
#95% CI                   [-25%, 0.44%]   [-25%, 0.44%]  
#
#Posterior tail-area probability p:   0.02764
#Posterior probability of an effect:  97.236%
#
#For more details, type: summary(impact, "report")


bsts_model <- impact2$model$bsts.model
coef_samples <- bsts_model$coefficients

# Observation noise samples
sigma_samples <- bsts_model$sigma.obs

# State contributions (list of matrices)
state_samples <- bsts_model$state.contributions

# Counterfactual prediction samples
cf_samples <- impact2$impact$pred  
coef_mcmc  <- as.mcmc(coef_samples)
sigma_mcmc <- as.mcmc(sigma_samples)

x11()
plot(coef_mcmc)

#Trace plot for observation noise
x11()
plot(sigma_mcmc)

#Posterior density plots for coefficients
x11()
densplot(coef_mcmc)

#Plot state components
x11()
plot(bsts_model, "state")

#Inclusion Probability
PlotBstsCoefficients(bsts_model)

#Save interpolated dataset---- 
write.csv(causalImpact_data_interpolated, file = "causalImpact_data_interpolated.csv")
