RARtrials: Response Adaptive Randomization in Clinical Trials
============================================================

**RARtrials** is designed for simulating some popular response-adaptive randomization methods in the literature with comparisons of each treatment group to a control group under no delay and delayed (time between treatment and outcome availability) scenarios. All the designs are based on one-sided tests with a choice from values of 'upper' and 'lower'. The general assumption is that binary outcomes follow Binomial distributions, while continuous outcomes follow normal distributions. Additionally, the number of patients accrued in the population follows a Poisson process and users can specify the enrollment rate of patients enrolled in the trial. 

Install RAR from CRAN with:

```r
install.packages('RAR')
```

Alternatively, install the RAR package from github with:


```r
#install.packages('devtools')
devtools::install_github("yayayaoyaoyao/RARtrials")
```

# Usage
There are two main groups of functions: those for simulating trials, which begin with `sim_`, and other functions that constitute the code for `sim_` with varying names. Functions included in this R package are as follows:

- `sim_RPTW` for the Randomized Play-the-Winner rule with binary outcomes in two-armed trials (Wei and Durham, 1978);

- `sim_dabcd_small_var` for the doubly adaptive biased coin design targeting Neyman allocation and RSIHR allocation using minimal variance strategy with binary outcomes in trials with up to five arms (Biswas and Mandal, 2004; Atkinson and Biswas, 2013) and `dabcd_small_var` calculates the allocation probabilities with available data using this method;

- `sim_dabcd_max_power` for the doubly adaptive biased coin design targeting Neyman allocation and RSIHR allocation using maximal power strategy with binary outcomes in trials with up to five arms and up to three arms respectively (Tymofyeyev, Rosenberger, and Hu, 2007; Jeon and Hu, 2010; Bello and Sabo, 2016) and `dabcd_max_power` calculates the allocation probabilities with available data using this method;

- `sim_A_optimal_known_var`, `sim_A_optimal_unknown_var`, `sim_Aa_optimal_known_var`, `sim_Aa_optimal_unknown_var`, `sim_RSIHR_optimal_known_var` and `sim_RSIHR_optimal_unknown_var` for Neyman allocation ($A_a$-optimal allocation and $A$-optimal allocation) and generalized RSIHR allocation subject to constraints for continuous outcomes with known and unknown variances in trials with up to five arms (Sverdlov and Rosenberger, 2013; Biswas and Mandal, 2004; Atkinson and Biswas, 2013);

- `sim_brar_binary`, `sim_brar_known_var` and `sim_brar_unknown_var` for Bayesian response-adaptive randomization using the Thall & Wathen method for binary outcomes, continuous outcomes with known and unknown variances in trials with up to five arms (Thall and Wathen, 2007); `brar_select_au_binary`, `brar_select_au_known_var` and `brar_select_au_unknown_var` can select appropriate $a_U$ using this method under null hypotheses; Functions start with `pgreater_` calculate the posterior probability of stopping a treatment group due to futility around $1\%$; Functions start with `pmax_` calculate the posterior probability that a particular arm is the best in a trial; `convert_gamma_to_chisq`, `convert_chisq_to_gamma` and `update_par_nichisq` are particular set-up for continuous outcomes with unknown variances;


- `sim_flgi_binary`, `sim_flgi_known_var` and `sim_flgi_unknown_var` for the forward-looking Gittins index rule and the controlled forward-looking Gittins index rule for binary outcomes and continuous outcomes with known and unknown variances in trials with up to five arms (Villar, Wason, and Bowden, 2015; Williamson and Villar, 2019); `flgi_cut_off_binary`, `flgi_cut_off_flgi_known_var` and `flgi_cut_off_flgi_unknown_var` can select cut-off values at the final stage for statistical inference; `Gittins` provides Gittins indices for binary reward processes and normal reward processes with known and unknown variance for certain discount factors.


