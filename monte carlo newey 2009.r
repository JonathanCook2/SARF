######################################################################
# File name: monte carlo newey 2009.R
##
# Description: This file compares SARF with a procedure based on 
#				Newey's (2009) procedure for recovering regression
#				coefficients. Estimating using naive linear regression
#				and random forest are also provided for reference.
#
#				These results are presented in the appendix of 
#				"Sample-selection-adjsuted random forests"
##
# Contents: 1. Load libraries 
# 2. Set simulation parameters
# 3. Create matrixes for storing simulaiton results
# 4. Run simulation with a linear relationship
# 5. Run simulation with a non-linear relationship
##
# Date created: 12/11/2019
# Last edited: 5/12/2022
# Author: Jonathan Cook
######################################################################
# 1 Load libraries 
library(grf)
library(MASS)
library(svMisc)

######################################################################
# 2 Set simulation parameters
min_node <- 1
num_trees <- 1000

sim <- 1000 # Number of simulations
obs <- 10000 # Observations used for each simulation

range_x <- 3
var_epsilon <- 2
v1 <- matrix(c(var_epsilon, (0.9*var_epsilon^.5), (0.9*var_epsilon^.5),1), 2, 2)

######################################################################
# 3 Create matrixes for storing simulaiton results
newey_reg_coef <- heckman_intercept <- 
naive_reg_coef <- naive_intercept <- rep(NA, sim)

naive_ols_MSE <- newey_MSE <- naive_RF_MSE <- SARF_MSE <- rep(NA, sim)
naive_ols_bias <- newey_bias <- naive_RF_bias <- SARF_bias <- rep(NA, sim)

######################################################################
# 4 Run simulation with a linear relationship

set.seed(1234)

for(i in 1:sim){
# Generate sample
x <- runif(obs,-range_x,range_x)
z <- runif(obs,-range_x,range_x)

draw <- mvrnorm(n = obs, mu = rep(0, 2), Sigma = v1)
latent_y <- x + draw[,1]
latent_s <- x + z + draw[,2]
s <- 1*(latent_s > 0)
y <- ifelse(s == 1, latent_y, NA)

# Naive linear regression 
ols <- lm(y ~ x)
naive_intercept[i] <- coef(ols)[1]
naive_reg_coef[i] <- coef(ols)[2]

# Naive RF 
naive_RF <- regression_forest(as.matrix(x[which(s == 1)]), as.matrix(y[which(s == 1)]), min.node.size = min_node, num.trees = num_trees, honesty = TRUE)

# Find probabilties of selection
selection_RF <- regression_forest(as.matrix(cbind(x, z)), as.matrix(s), min.node.size = min_node, num.trees = num_trees, honesty = TRUE)
p_hat_RF <- predict(selection_RF, as.matrix(cbind(x, z)))$predictions

# Newey (2009) with intercept identified at infinity
p_hat_RF_sq <- p_hat_RF^2
p_hat_RF_cu <- p_hat_RF^3
newey_reg_RF <- lm(y ~ x + p_hat_RF + p_hat_RF_sq + p_hat_RF_cu)
newey_reg_coef[i] <- coef(newey_reg_RF)[2]
heckman_intercept[i] <- mean(y[which(p_hat_RF > .95 & s == 1)] - coef(newey_reg_RF)[2]*x[which(p_hat_RF > .95 & s == 1)])

# SARF
SARF <- regression_forest(as.matrix(cbind(x[which(s == 1)], p_hat_RF[which(s == 1)])), as.matrix(y[which(s == 1)]), min.node.size = min_node, num.trees = num_trees, honesty = TRUE)


# Random sample for finding MSE and bias
runif(obs,-range_x,range_x)
runif(obs,-range_x,range_x)
draw <- mvrnorm(n = obs, mu = rep(0, 2), Sigma = v1)
latent_y <- x + draw[,1]

# Make predictions
naive_ols_pred <- naive_intercept[i] + naive_reg_coef[i]*x
naive_RF_pred <- predict(naive_RF, as.matrix(x))$predictions
newey_pred <- heckman_intercept[i] + newey_reg_coef[i]*x
SARF_pred <- predict(SARF, as.matrix(cbind(x,rep(1,obs))))$predictions

# Find MSEs
naive_ols_MSE[i] <- mean((latent_y - naive_ols_pred)^2)
naive_RF_MSE[i] <- mean((latent_y - naive_RF_pred)^2)
newey_MSE[i] <- mean((latent_y - newey_pred)^2) 
SARF_MSE[i] <- mean((latent_y - SARF_pred)^2)

# Find bias
naive_ols_bias[i] <- mean(naive_ols_pred - latent_y)
naive_RF_bias[i] <- mean(naive_RF_pred - latent_y)
newey_bias[i] <- mean(newey_pred - latent_y)
SARF_bias[i] <- mean(SARF_pred - latent_y)
progress(100*(i / sim))
}

# Naive OLS
mean(naive_intercept)
sd(naive_intercept)
mean(naive_reg_coef)
sd(naive_reg_coef)
mean(naive_ols_MSE)
sd(naive_ols_MSE)
mean(naive_ols_bias)
sd(naive_ols_bias)

# Corrected OLS
mean(heckman_intercept)
sd(heckman_intercept)
mean(newey_reg_coef)
sd(newey_reg_coef)
mean(newey_MSE)
sd(newey_MSE)
mean(newey_bias)
sd(newey_bias)

# Naive RF
mean(naive_RF_MSE)
sd(naive_RF_MSE)
mean(naive_RF_bias)
sd(naive_RF_bias)

# SARF
mean(SARF_MSE)
sd(SARF_MSE)
mean(SARF_bias)
sd(SARF_bias)


######################################################################
# 5 Run simulation with a non-linear relationshi

set.seed(1234)

for(i in 1:sim){
# Make x and z correlated
x <- runif(obs,-range_x,range_x)
z <- runif(obs,-range_x,range_x)
# Make u and epsilon bivar norm
draw <- mvrnorm(n = obs, mu = rep(0, 2), Sigma = v1)

#latent_y <- pnorm(x) + draw[,1]
latent_y <- x^2 + draw[,1]
latent_s <- x - z + draw[,2]
s <- 1*(latent_s > 0)
y <- ifelse(s == 1, latent_y, NA)


# Naive linear regression 
ols <- lm(y ~ x)
naive_intercept[i] <- coef(ols)[1]
naive_reg_coef[i] <- coef(ols)[2]

# Naive RF 
naive_RF <- regression_forest(as.matrix(x[which(s == 1)]), as.matrix(y[which(s == 1)]), min.node.size = min_node, num.trees = num_trees, honesty = TRUE)

# Find probabilties of selection
selection_RF <- regression_forest(as.matrix(cbind(x, z)), as.matrix(s), min.node.size = min_node, num.trees = num_trees, honesty = TRUE)
p_hat_RF <- predict(selection_RF, as.matrix(cbind(x, z)))$predictions

# Newey (2009) with intercept identified at infinity
p_hat_RF_sq <- p_hat_RF^2
p_hat_RF_cu <- p_hat_RF^3
newey_reg_RF <- lm(y ~ x + p_hat_RF + p_hat_RF_sq + p_hat_RF_cu)
newey_reg_coef[i] <- coef(newey_reg_RF)[2]
heckman_intercept[i] <- mean(y[which(p_hat_RF > .95 & s == 1)] - coef(newey_reg_RF)[2]*x[which(p_hat_RF > .95 & s == 1)])

# SARF
SARF <- regression_forest(as.matrix(cbind(x[which(s == 1)], p_hat_RF[which(s == 1)])), as.matrix(y[which(s == 1)]), min.node.size = min_node, num.trees = num_trees, honesty = TRUE)

# Random sample for finding MSE and bias
x <- runif(obs,-range_x,range_x)
z <- runif(obs,-range_x,range_x)
draw <- mvrnorm(n = obs, mu = rep(0, 2), Sigma = v1)
latent_y <- x^2 + draw[,1]

# Make predictions
naive_ols_pred <- naive_intercept[i] + naive_reg_coef[i]*x
naive_RF_pred <- predict(naive_RF, as.matrix(x))$predictions
newey_pred <- heckman_intercept[i] + newey_reg_coef[i]*x
SARF_pred <- predict(SARF, as.matrix(cbind(x,rep(1,obs))))$predictions

# Find MSEs
naive_ols_MSE[i] <- mean((latent_y - naive_ols_pred)^2)
naive_RF_MSE[i] <- mean((latent_y - naive_RF_pred)^2)
newey_MSE[i] <- mean((latent_y - newey_pred)^2) 
SARF_MSE[i] <- mean((latent_y - SARF_pred)^2)

# Find bias
naive_ols_bias[i] <- mean(naive_ols_pred - latent_y)
naive_RF_bias[i] <- mean(naive_RF_pred - latent_y)
newey_bias[i] <- mean(newey_pred - latent_y)
SARF_bias[i] <- mean(SARF_pred - latent_y)
progress(100*(i / sim))
}

# Naive OLS
mean(naive_intercept)
sd(naive_intercept)
mean(naive_reg_coef)
sd(naive_reg_coef)
mean(naive_ols_MSE)
sd(naive_ols_MSE)
mean(naive_ols_bias)
sd(naive_ols_bias)

# Corrected OLS
mean(heckman_intercept)
sd(heckman_intercept)
mean(newey_reg_coef)
sd(newey_reg_coef)
mean(newey_MSE)
sd(newey_MSE)
mean(newey_bias)
sd(newey_bias)

# Naive RF
mean(naive_RF_MSE)
sd(naive_RF_MSE)
mean(naive_RF_bias)
sd(naive_RF_bias)

# SARF
mean(SARF_MSE)
sd(SARF_MSE)
mean(SARF_bias)
sd(SARF_bias)




