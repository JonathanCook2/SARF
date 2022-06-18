# Sample-selection-adjusted random forests (SARF)

## Description

Code used in the paper 

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Cook, Jonathan (forthcoming): "Sample-selection-adjusted random forests," _International Journal of Data Science and Analytics_. [http://doi.org//10.1007/s41060-022-00337-w](http://doi.org//10.1007/s41060-022-00337-w) 

is available here. Please see this paper for more details and cite this paper when using the sample-selection-adjusted random forests (SARF) procedure.

A predictive model that is trained with non-randomly selected samples can offer biased predictions for the population. Cook (forthcoming) discusses when non-random selection is a problem. For the applications in which it is a problem, Cook presents a procedure for adjusting the predictions of random forest to account for non-random sampling of the training data. This adjustment results in more accurate predictions for the population.

The file "Figures 1 and 2.R" simulates data subject to sample selection then performs the SARF procedure on the simulated data. A naive random forest is also fitted using the selected sample for comparison. After forming predictions, the partial depenedence plot for predicted probabilies is compared to the correct value. The plots of partial dependence plot and the MSE as a function of "p tilde" are Figures 1 and 2 in Cook (forthcoming).)

The file "Monte Carlo Newey 2009.r" provides a Monte Carlo experiment that compares SARF with a procedure based on Newey (2009), a naive linear regression, and a naive random forest regression. These results are presented in the appendix of Cook (forthcoming).

## References

Cook, Jonathan (forthcoming): "Sample-selection-adjusted random forests," _International Journal of Data Science and Analytics_. [http://doi.org//10.1007/s41060-022-00337-w](http://doi.org//10.1007/s41060-022-00337-w) 

Newey, Whitney K. (2009): "Two‐step series estimation of sample selection models," _The Econometrics Journal_, 12, S217-S229. [https://doi.org/10.1111/j.1368-423X.2008.00263.x](https://doi.org/10.1111/j.1368-423X.2008.00263.x)
