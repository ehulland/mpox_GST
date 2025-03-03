# Code and data compendium for the manuscript: '"What’s in a name?": Using mpox as a case study to understand the importance of communication, advocacy, and information accuracy in disease nomenclature'

## Erin N. Hulland, Marie-Laure Charpignon, Ghinwa Y. Hayek, Angel N. Desai, and Maimuna S. Majumder

______

This document describes the analytic step through with a focus on the elements of the R code included in this repository and details not articulated in the main manuscript or appendix. 

### Data acquisition and preprocessing step

GST data were downloaded from https://trends.google.com/trends/explore for `monkeypox` and `mpox` simultaneously between the period of November 1st, 2022 and December 1st, 2023. This was done for each of the 184 countries included in the analysis. For countries with one or more of the five UN / WHO languages (Appendix Table 1), this search was conducted again using the term `mpox` and the foreign-language term (Appendix Table 2) simultaneously. Covariate data were also acquired at this time.

### R codes for the main analyses 

These are the steps for running each piece of code in the main pipeline and the general details about each step.

Start with code **00_data_import_cleaning.R**, making sure to specify your own working directory. This code can be run top-to-bottom.

This import and cleaning code brings in all the GST data and covariate data, standardizes country names and spelling, and creates our two outcome variables: our **binary** variable of being an adopter country or not and our **continuous** variable defined as the fraction of all searches that were for `mpox`. Additional covariate processing is tackled here, as described in appendix section 3.1. At the end of this code, we have the following outputs:

 * **full_data**: a dataset with 192 countries included and all covariate data. This dataset is used in imputation Approach 2 and 3 where variable transformations are done per imputation.
* **full_dat**: a subset of full_data with 188 countries included (excluding those missing GDP), and with all transformations done to the covariates prior to building a model. This dataset is used only in imputation Approach 1 (complete case analysis).
* **listwise_del_ctrys**: This is a list of the 154 countries included in the complete case analysis (Approach 1, also referred to as listwise deletion)
* **final_166_ctrys**: This is a list of the 166 countries included in the GAI-imputation only method (Approach 2)
* **M**: This is a matrix of all the correlations between covariates in the analysis
* **get_betas**: This is a function created to pull n draws from a multivariate random normal distribution (using the MVRNORM function from the *MASS* package) with specified coefficient estimates and variance-covariance matrix




Next run **01_mpox_models_with_all_imputation.R**. This code calls in **00_data_import_cleaning.R** to be run automatically if you do not want to run it manually. 

This code is broken down into two major sections for each the binary and continuous outcome, and then three subsections (one for each imputation strategy outlined in the paper: Approach 1, complete case analysis; Approach 2, GAI-only imputation; Approach 3, full imputation). Each imputation strategy is replicated across both outcomes.

The specific steps taken in each subsection are described below. 

2.1) For the first imputation strategy of complete case analysis (i.e. no imputation at all), we build all 18 univariable models separately. We then take ten draws from a multivariate random normal distribution for each model. For each model, we consider whether the variable was significant at a threshold of 0.1 and if it was, we include it in our multivariable model, otherwise we exclude it; highly collinear variables are also excluded (see Table 1 of the manuscript). After constructing our multivariable model, we again take 10 draws from a multivariate random normal distribution, accounting for the variance-covariance matrix of the multivariable model. From these sets of 10 draws, we then calculate the median, 2.5th and 97.5th percentiles to generate our 95% uncertainty intervals. The results of this section include two forest plots per outcome: one for the univariable models (Appendix Figure S8a binary outcome and S9a continuous outcome) and one for the multivariable model with inclusion at 0.1 (Appendix Figure S8b binary outcome and S9b continuous outcome).

2.2) For the second imputation strategy, we remove all rows of data with missing data other than GAI and only impute GAI and perform complete case analysis. We perform 10 imputations using predictive mean matching to impute GAI. Within each imputation, we construct our 18 univariable models and for each one we pull one draw of the multivariate random normal distribution using MVRNORM and the **get_betas** function defined above. Still within a singular imputation, we then determine which variables are significant at the 0.1 threshold; these variables are then used to construct the multivariable model. The multivariable model is then evaluated within the imputation by taking one multivariate random normal draw. All ten imputations (and MVRNORM draws) are then summarized jointly, using the median, 2.5th, and 97.5th percentiles to generate the 95% uncertainty interval. The results of this section include two forest plots per outcome: one for the univariable models (Appendix Figure S10a binary outcome and S11a continuous outcome) and one for the multivariable model with inclusion at 0.1 (Appendix Figure S10b binary outcome and S11b continuous outcome)


2.3) For the final imputation strategy, we impute all variables with any missing data. Outside of the imputation strategy, this method is identical to step 2.2. We perform 10 imputations using predictive mean matching to impute GAI. Within each imputation, we construct our 18 univariable models and for each one we pull one draw of the multivariate random normal distribution using MVRNORM and the **get_betas** function defined above. Still within a singular imputation, we then determine which variables are significant at the 0.1 threshold; these variables are then used to construct the multivariable model. The multivariable model is then evaluated within the imputation by taking one multivariate random normal draw. All ten imputations (and MVRNORM draws) are then summarized jointly, using the median, 2.5th, and 97.5th percentiles to generate the 95% uncertainty interval. We additionally consider a threshold of 0.05 for multivariable model inclusion as a sensitivity analysis; these results can be seen in Appendix Figures 13 a) binary and b) continuous). The results of this section include three forest plots per outcome: one for the univariable models (Main Manuscript, Figure 1a binary and Figure2a continuous), one for the multivariable model with inclusion at 0.1 (Main Manuscript, Figure 1b binary and Figure 2b continuous), and one for the multivariable model with inclusion at 0.05 (Appendix Figure S14a and S14b)

### R codes for the sensitivity analyses

In order to run the numerous sensitivity analyses included in this manuscript, we include the following codes:

**02_mpox_models_foreign_language.R**: This code performs the analysis for the foreign language analysis results shown in Appendix Figures S6 and S7. It only includes imputation Approach 3 (full imputation). This code uses the same imputation strategy as step 2.3 above, and produces two forest plots per outcome: one for the univariable models (Appendix Figure S6a binary outcome and S7a continuous outcome) and one for the multivariable model with inclusion at 0.1 (Appendix Figure S6b binary outcome and S7b continuous outcome). Importantly, to get the accurate results for this analysis, a *different* preprocessing and importation code needs to be run, as specified at the top of the code: 

**00_data_import_cleaning_foreign_language.R**: This code is very similar to the code in step 1 above, but swaps in the foreign language `monkeypox` search intensity for english search intensity. It also does not produce correlation matrix **M**.

Lastly, at the end of this script there is code to reproduce the counterfactual scenerii for Ecuador and Madagascar where the second foreign language was unavailable. These results are visible in Appendix Figure S20. 

**03_mpox_models_mpox_case_countries.R**: This code performs the analysis for the set of countries that had every recorded a human or mammalian mpox case since its discovery in humans in the 1970s, results shown in Appendix Figures S12 and S13. It only includes imputation Approach 3 (full imputation). This code uses the same imputation strategy as step 2.3 above, and produces two forest plots per outcome: one for the univariable models (Appendix Figure S12a binary outcome and S13a continuous outcome) and one for the multivariable model with inclusion at 0.1 (Appendix Figure S12b binary outcome and S13b continuous outcome). This code uses the original data import and cleaning code. 

**04_mpox_models_tobit.R**: This code performs the analysis using a Tobit model for the zero-inflated continuous outcome, with results shown in Appendix Figure S15. It only includes imputation Approach 3 (full imputation). This code uses the same imputation strategy as step 2.3 above, and produces two forest plots, only for the continuous outcome: one for the univariable models (Appendix Figure S15a) and one for the multivariable model with inclusion at 0.1 (Appendix Figure S15b). This code uses the original data import and cleaning code. 

**05_mpox_models_bonferroni_correction.R**: This code performs the analysis for the Bonferroni corrected p-values for inclusion into the model and multivariable model significance It only includes imputation Approach 3 (full imputation). Because model selection is done per-imputation, this code summarizes the results of the univariable models and multivariable models based on significance threshold per imputation. 

**06a_mpox_models_rockchalk.R**: This code performs the analysis using the *rockchalk* package in R in lieu of the *MASS* package, which allows for seeding of the model results, and thus exact replication. Results are shown in Appendix Figures S16 and S17. It only includes imputation Approach 3 (full imputation). This code uses the same imputation strategy as step 2.3 above, and produces two forest plots per outcome: one for the univariable models (Appendix Figure S16a binary outcome and S17a continuous outcome) and one for the multivarible model with inclusion at 0.1 (Appendix Figure S16b binary outcome and S17b continuous outcome). This code uses the original data import and cleaning code. 

**06a_mpox_models_mvtnorm.R**: This code performs the analysis using the *mvtnorm* package in R in lieu of the *MASS* package, which allows for seeding of the model results, and thus exact replication. Results are shown in Appendix Figures S18 and S19. It only includes imputation Approach 3 (full imputation). This code uses the same imputation strategy as step 2.3 above, and produces two forest plots per outcome: one for the univariable models (Appendix Figure S18a binary outcome and S19a continuous outcome) and one for the multivarible model with inclusion at 0.1 (Appendix Figure S18b binary outcome and S19b continuous outcome). This code uses the original data import and cleaning code. 
