######################################################################
# File name: Figures 1 and 2.R
##
# Description: This file creates simulated data subject to sample
# selection, then estimates a standard random forest
# and a sample-selection-adjusted random forest.
#
# The plots created by this code are Figures 1 and 2 in "Sample-
#	selection-adjusted random forests"
##
# Contents: 1. Load libraries 
# 2. Set simulation parameters
# 3. Create matrixes for storing simulaiton results
# 4. Run simulation
# 5. Examine results
# 5.a MSE and bias for predictions
# 5.b Partial dependence plot for p_hat
##
# Date created: 12/11/2019
# Last edited: 5/12/2022
# Author: Jonathan Cook
######################################################################
# 1 Load libraries 
library(MASS)
library(mvtnorm)
library(grf)
library(svMisc)

######################################################################
# 2 Set simulation parameters
rho <- .9	# Correlation in unobservables
n <- 10000	# Sample size
sigma_sq <- 1 # Variance of error term for outcomes
range_of_x <- 2 # The max and min values of x

varcovar <- matrix(c(sigma_sq, (sigma_sq^.5)*rho, (sigma_sq^.5)*rho, 1), nrow = 2, ncol = 2)

# Function to generate latent variables
gen_latent_var <- function(outcome_vars, selection_vars){
	error <- mvrnorm(n, rep(0, 2), varcovar)
	latent_y <- sin(outcome_vars) + error[, 1]
	latent_s <- selection_vars[, 1]^2 - selection_vars[, 2]^2 + error[, 2]

	return(cbind(latent_y, latent_s))
}
# Function to calculate MSE
MSE <- function(predictions, actual){
	return(mean((predictions - actual)^2, na.rm=TRUE))
}
# Function to calculate bias
bias <- function(predictions, actual){
	return(mean(predictions - actual, na.rm=TRUE))
}
######################################################################
# 3 Create matrixes for storing simulaiton results
actual_outcomes <-
	observed_outcomes <-
		naive_HF_predictions <-
			adjusted_HF_predictions <- matrix(NA, nrow=n, ncol=1)
######################################################################
# 4 Run simulation 

set.seed(1)

# Draw data (training sample)
x <- runif(n, -range_of_x, range_of_x)
z <- runif(n, -range_of_x, range_of_x)
latent_vars <- gen_latent_var(outcome_vars = x, selection_vars = cbind(x, z))
selection <- 1*(latent_vars[, 2] > 0)
y <- ifelse(selection == 1, latent_vars[,1], NA) 
latent_y <- latent_vars[, 1]

# RF for prob selection
full_data <- as.data.frame(cbind(y, selection, x, z))
HF <- regression_forest(as.matrix(full_data[, 3:4]), as.matrix(full_data[, 2]), min.node.size=1, num.trees=1000, honesty=TRUE)
p_hat_HF <- predict(HF,as.matrix(full_data[, 3:4]))$predictions

# Plot p_hat versus actual probability of selection 
true_prob <- pnorm(x^2 - z^2)
plot(true_prob, p_hat_HF, xlab="Actual", ylab="Prediction")

# RF for outcome
full_data <- as.data.frame(cbind(y, x, p_hat_HF))
full_data <- full_data[complete.cases(full_data), ]
naive_HF <- regression_forest(as.matrix(full_data[, 2]), as.matrix(full_data[, 1]), min.node.size=1, num.trees=1000, honesty=TRUE)

# Adjusted RF for outcome
outcome_HF <- regression_forest(as.matrix(full_data[, 2:3]), as.matrix(full_data[,1]), min.node.size=1, num.trees=1000, honesty=TRUE)

# Draw data (testing sample)
x <- runif(n, -range_of_x, range_of_x)
z <- runif(n, -range_of_x, range_of_x)
latent_vars <- gen_latent_var(outcome_vars = x, selection_vars = cbind(x, z))
selection <- 1*(latent_vars[, 2] > 0)
y <- ifelse(selection == 1, latent_vars[, 1], NA) 
actual_outcomes[, 1] <- latent_vars[, 1] 	
observed_outcomes[, 1] <- y
	
###################################################################### 
# 5. Examine results
# 5.a. MSE and bias for predictions

# Predictions
full_data <- as.data.frame(cbind(y, x))
naive_HF_predictions <- predict(naive_HF, as.matrix(full_data[, 2]))$predictions
naive_MSE <- MSE(naive_HF_predictions, actual_outcomes)

all_MSEs <- all_p <- rep(NA, 101)

for(i in 0:100){
	p_tilde <- i / 100
	p_hat_HF <- runif(n, p_tilde, 1)
	full_data <- as.data.frame(cbind(y, x, p_hat_HF))
	adjusted_HF_predictions <- predict(outcome_HF, as.matrix(full_data[, 2:3]))$predictions
	all_MSEs[(i + 1)] <- MSE(adjusted_HF_predictions, actual_outcomes)
	all_p[(i + 1)] <- p_tilde
	progress(100*(i / 101))
}

plot(all_p, all_MSEs, type="l", lwd=2, xlab=expression(tilde(p)), ylab="MSE")
abline(h=naive_MSE, lty=3, lwd=2)  


# 5.b. Partial dependence plot for p_hat 
set.seed(1)
x <- runif(n, -range_of_x, range_of_x)
z <- runif(n, -range_of_x, range_of_x)
latent_vars <- gen_latent_var(outcome_vars = x, selection_vars = cbind(x, z))
selection <- 1*(latent_vars[, 2] > 0)
y <- ifelse(selection == 1, latent_vars[, 1], NA)

# RF for prob selection
full_data <- as.data.frame(cbind(y, selection, x, z))
HF <- regression_forest(as.matrix(full_data[, 3:4]),as.matrix(full_data[, 2]), min.node.size=1, num.trees=1000, honesty=TRUE)
p_hat_HF <- predict(HF, as.matrix(full_data[, 3:4]))$predictions

# RF for outcome
full_data <- as.data.frame(cbind(y, x, p_hat_HF))
full_data <- full_data[complete.cases(full_data), ]

# Corrected RF for outcome
outcome_HF <- regression_forest(as.matrix(full_data[, 2:3]), as.matrix(full_data[, 1]), min.node.size=1, num.trees=1000, honesty=TRUE)
orig.data <- full_data

pdp_p_hat <- rep(NA, nrow=dim(full_data)[1])
size <- dim(full_data)[1]
for(i in 1:size){
	full_data[, 3] <- orig.data[i, 3]
	preds <- predict(outcome_HF, full_data[, 2:3], estimate.variance=FALSE)
	pdp_p_hat[i] <- mean(preds$predictions)
	full_data <- orig.data
	progress(100*(i / size))
}

mean_f <- mean(full_data[, 2])
correct_pdp <- mean_f + rho*(sigma_sq^.5)*dnorm(qnorm(orig.data[, 3])) / pnorm(qnorm(orig.data[, 3])) 

for_pdp_plot <- cbind(orig.data[,3], pdp_p_hat, correct_pdp)
for_pdp_plot <- for_pdp_plot[order(orig.data[,3]), ]
plot(for_pdp_plot[, 1], for_pdp_plot[, 2], type= "l", lwd=2, xlab="Probability of selection", ylab="PDP")
lines(for_pdp_plot[, 1], for_pdp_plot[, 3], col="red", lwd =2, lty=2)

# Note: Need to check the range of x for observed and full population, as RF does not extrapolate well.
# Range of x
c(min(x), max(x))
# Range of x for observed data
c(min(x[which(selection == 1)]), max(x[which(selection == 1)]))
