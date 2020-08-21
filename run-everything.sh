#!/bin/bash

# RECEIPE TO FOLLOW TO RUN EVERYTHING
# remember that these are roughly in order,
#  as some scripts require outputs of others

# 1: make sure all the input data are set up
# jelly biomass data
./prepare-process_inputs_GZbiomass.R

# 2: set up physical and biological data inputs
./prepare-process_woa_temp.R
./prepare-process_sst.R
./prepare-process_inputs_phytoplankton.R
./prepare-process_inputs_zooplankton.R
./prepare-process_inputs_export.R

# 3: model structure prep
Rscript prepare-process_global_sinking_rates.R
Rscript prepare-process_global_sinking_rates_2.R
Rscript prepare-compute_biomes.R
Rscript prepare-tunicate_allometry.R
Rscript prepare-MCparameters.R
Rscript prepare-fill_NAs_gz_model_inputs.R

# 4: setup and run model
./run-setup_biomass_mod_cases.R
./run-carbon_flux_model.R 1
./run-carbon_flux_model.R 2
./run-carbon_flux_model.R 3
./run-setup_and_execute_MC_ensemble_runs.R 1 1 # this calls run-export_flux.R
./run-setup_and_execute_MC_ensemble_runs.R 2 1
./run-setup_and_execute_MC_ensemble_runs.R 3 1

# 5: extract results
./run-extract_MC_ensemble_results.R 1
./run-extract_MC_ensemble_results.R 2
./run-extract_MC_ensemble_results.R 3

# 6: collating, plotting results
./analysis-results_by_biome.R 1
./analysis-results_by_biome.R 2
./analysis-results_by_biome.R 3
./analysis-MC_accepted_sims_stats.R
./analysis-results_by_biome_combResults.R
./analysis-generate_ensemble_means.R
./analysis-1deg_results.R
./analysis-martin_curve_calc.R 
./analysis-sinking_flux.R 1
./analysis-sinking_flux.R 2
./analysis-sinking_flux.R 3
./analysis-sinking_flux_combBiomes.R
./analysis-extract_point_comparisons.R 1
./analysis-extract_point_comparisons.R 2
./analysis-extract_point_comparisons.R 3
./analysis-plot_point_comparisons.R

