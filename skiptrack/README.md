# skiptrack

* `clean_clue.R` : Function to convert Clue data to skipTrack-style cycle tracking data. In particular:
	- Clue data are logged as "days". A person may or may not enter a day as a "period" or "spotting" day (or enter birth control data). This function converts those entries in cycles according to the following rules...
		* Each "period" day is part of of a period.
		* A "spotting" day is only part of a period if it is immediately preceded or succeeded by a "period" day ('buffer' = 1 day.)
		* Each "period" must be >7 days apart, otherwise entries are considered to be part of the same period.
		* Cycle length = days between 1st day of current cycle & 1st day of next cycle (therefore, final period recorded is dropped; it has no cycle length, since the first day of the next period is not observed.)
		* Cycle number = index j for cycle j observed for individual i. That is, each cycle can be indexed as $cycle_{ij}$: cycle j for individual i. And each $cycle_{ij}$ has a corresponding cycle length.
		* Cycle lengths > 365 days are excluded.
		* If individuals report < 2 cycles (i.e. only 1), they are excluded from the dataset. Since the outcome of interest is individual-level cycle variability, having a single cycle length for any given individual is uninformative.

