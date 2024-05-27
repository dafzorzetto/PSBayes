## Code files:

### General functions (in **`src`** folder):
 - **`functions_simulations.R`**:
   
    functions to simulate different settings
 - **`CASBAH.R`**:

   Gibbs sampler to estimate our proposed model: CASBAH
 - **`model_with_SUN.R`**:

   functions to estimate CASBAH with the conjugate prior for the probit (SUN distribution)

   suggestion: use CASBAH.R to avoid matrices dimensional problems
- **`competitor_SLM.R`**:

  functions to estimate the competitor model, proposed by <a href=https://www.tandfonline.com/doi/abs/10.1198/jasa.2011.ap10425> _"A Bayesian Semiparametric Approach to Intermediate Variables in Causal Inference"_ </a> by Schwartz, Li, and Mealli
 - **`plots_simulations.R`**:

   functions to visualize the simulated data

 - **`visualization_results_functions.R`**:

   function to visualize the estimated principal causal effects

### Reproduce the results in the Simulation Study Section (in **`simulation study`** folder):
 - **`1_simulation_data.R`**:

   generate the different settings
 - **`2_estimation_models.R`**:

   estimate CASBAH and competitor model for the different settings
 - **`3_bias_MSE_ARI.R`**:

   estimate bias and MSE for all the models and settings, and the adjusted ARI for CASBAH
 - **`4_visualizing_results.R`**:

   functions to visualize the results: bias and MSE

### Reproduce the results in the Application Section (in **`application`** folder):
 - **`1_merge_dataset.R`**:
   
    steps to merge the dataset
- **`2_model.R`**:
   
    estimation of CASBAH for the real dataset
- **`3_plot_paper.R`**:
   
    code to reproduce the results plots
