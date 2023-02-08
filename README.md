# MarketsMatterSim

This repository contains code and input data sets that run a simulation assessing RMSE of within- and out-of-market comparison groups in a difference-in-differences analysis.  The goal of the repository is to share the code and data sets for other researchers to use, either to validate the results of a (hopefully forthcoming) research methods article or to use as a template for their own simulation studies.

This repository is intended as an archive rather than a living or working project and thus is *read-only*.

Files are organized into data and code folders.

## Code

1. 2021-09-20_hierarchical_modeling.Rmd:  R script that processes background data to produce simulation parameters and prepares a data set of parameters for each scenario to be tested.  The script then takes as input the output of the simulate function and summarizes it.

2. 2021-09-20_simulate_function.R: R script that runs the simulation, given the input parameters produced in (1). 

3. generate_report_figures.Rmd: Script that produced the summary figures included in the publication.

## Data

* 2021-09-19_matching_scenarios_params.xlsx: Simulation input parameters specifying the difference in slopes between within-market and out-of-market comparisons.  These values are user-specified but grounded in analysis of the Comprehensive Primary Care Plus demonstration.
* 2021-09-20_simulation_summaries.csv: Summarized simulation results, produced by script (2).
* CPC_plus_counts.csv: Simulation input describing the number of primary care practices per market, for selected vs. non-selected markets, in the Comprehensive Primary Care Plus demonstration.
* cvd_totexp_analytic_data_060321.xlsx: Simulation input describing the mean and standard deviation of Medicare spending in each hospital referral region (HRR) per month in data from the Comprehensive Primary Care Plus demonstration.
* slope_sds_asof08042021.csv: Simulation input describing standard deviations across practice- and market-level intercepts as well as practice- and market-level slopes.  These estimates were calculated from a Bayesian hierarchical model fit to data from the Comprehensive Primary Care Plus demonstration.
