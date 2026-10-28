# datasets_tosend/eva/

* Datasets to send:
	- processed skipTrack dataset with clean_clue function + the following variance in cycle length measures:
		- raw variance
		- $\tau_i$ summary statistics (fit with 10,000 iterations to MCMC, no covariates)
		- $c_{ij}$s
		- MCMC convergence diagnostics (see previous writeups)
	- $\tau_i$ kernel-smoothed posterior densities
	- skipTrack "as designed" ((fit with 10,000 iterations to MCMC) -- *supplemental*
		- stress as a covariate
		- stress + age as covariates
		- MCMC convergence diagnostics (see previous writeups)
	- Custom script/function for sampling efficiently from these posteriors ^ for MIME-type sampling (propagation of $\tau_i$ uncertainty to regression coefficient estimates)

* `mu_cyclelen_dists1.Rmd` -- Generates mean cycle length distributions -- a) raw, b) c_ij-corrected, c) median from SkipTrack, d) mean derived from SkipTrack median (under distributional assumptions) 


