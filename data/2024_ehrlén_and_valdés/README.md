# Selection favors high spread and asymmetry of flower opening dates within plant individuals

[https://doi.org/10.5061/dryad.2jm63xsz1](https://doi.org/10.5061/dryad.2jm63xsz1)

This data includes data on individuals (data_ids_clean.csv) and on individual flowers (data_id_flowers_clean.csv) coming from recordings of 5287 individual flowers of *Lathyrus vernus* belonging to 495 flowering events (i.e. one individual flowering in one year) over three years. It also includes weather data from a nearby meteorological station (weather_clean.csv). The code used in the manuscript is on the file 5_code_analyses_dryad.Rmd.

## Description of the data and file structure

The file data_ids_clean.csv contains the following columns:

*   id: identifier of the plant individual
*   fr_set: fruit set (proportion of flowers that developed to mature fruits)
*   n_seed: total number of seeds
*   fitness: total number of intact seeds
*   prop_seed_preyed: proportion of preyed seeds
*   year: the year of the flowering event
*   var: variance of opening dates
*   skew: skewness of opening dates
*   kurt: kurtosis of opening dates
*   avFD_v: mean opening date (expressed as the number of days after the vernal equinox in each year)
*   n_fl_std: number of flowers (standardized)
*   avFD_std: mean opening date (standardized)
*   var_std: variance of opening dates (standardized)
*   skew_std: skewness of opening dates (standardized)
*   kurt_std: kurtosis of opening dates (standardized)
*   fitness_rel: relativized fitness
*   n_fl: number of flowers
*   n_mat_intact_fr: number of mature intact fruits
*   prop_seeds_escpred: proportion of seeds escaping predation
*   n_fl_log: number of flowers (log)
*   n_fl_log_std: number of flowers (log), standardized
*   npr_seeds: number of non-predated seeds
*   pr_seeds: number of predated seeds

The file data_id_flowers_clean.csv contains the following columns:

*   id: identifier of the plant individual
*   year: the year of the flowering event
*   opening_date_c: opening date of each individual flower (expressed as a calendar date)
*   opening_date_v: opening date of each individual flower (expressed as the number of days after the vernal equinox in each year)

The file weather_clean.csv contains the following columns:

*   date: the date of the recording
*   year: the year of the recording
*   month: the month of the recording
*   day: the day of the recording
*   mean: the daily mean temperature recorded (ºC)
*   min: the daily minimum temperature recorded (ºC)
*   max: the daily maximum temperature recorded (ºC)

The file 5_code_analyses_dryad.Rmd is an R notebook all the code needed to run the analyses on the manuscript. This code is uploaded as software to Zenodo linked in this publication. 
