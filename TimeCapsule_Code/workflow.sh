


####################
# To reproduce the real MCSA data results run:

R CMD INSTALL hmm_1.1-6.tar.gz

for seed in {1..10}
do
Rscript RunFile.r $seed RealData Bayesian no_scale 30 Real_data_output/
done

# The following line produces trace plots and histograms of the MCMC samples.  See the 
# supplementary materials.
Rscript OutFile.r Real_data_output/ 10 no_scale
# The following line produces figures 2, 3, 4, 8, 9, 10, 11, and Table 3.
Rscript OutFile_figures.r Real_data_output/ 10 no_scale

# This code uses 30 threads in parallel.  The for loop above can be run in an embarrassingly 
# parallel fashion.  The code for each $seed may take about 3 days to run.
###################



###################
# To reproduce the sensitivity analysis for the real MCSA data results run:

R CMD INSTALL hmm_1.1-6.tar.gz

for seed in {1..10}
do
Rscript RunFile.r $seed RealData Bayesian 10x_scale 30 Sensitivity_test/
done

# The following line produces trace plots and histograms of the MCMC samples.  See the 
# supplementary material.
Rscript OutFile.r Sensitivity_test/ 10 10x_scale
# The following line produces alternative Figures 2, 3, 4, 8, 9, 10, 11, and Table 3 for the
# MCMC samples based on diffuse prior specifications.
Rscript OutFile_figures.r Sensitivity_test/ 10 10x_scale

# This code uses 30 threads in parallel.  The for loop above can be run in an embarrassingly 
# parallel fashion.  The code for each $seed may take about 3 days to run.
###################



###################
# To reproduce the synthetic MCSA simulation results run:

R CMD INSTALL hmm_1.1-6.tar.gz

for seed in {1..50}
do
Rscript Simulate.r $seed Simulation_output/
Rscript RunFile.r $seed Simulation Bayesian no_scale 30 Simulation_output/
done

# The following line produces Figure 7 in the manuscript, as well as trace plots, 
# histograms, and box plots for the MCSA simulation results.  See the supplementary material.
Rscript OutFile_sim.r Simulation_output/ 50

# This code uses 30 threads in parallel.  The for loop above can be run in an embarrassingly 
# parallel fashion.  The code for each $seed may take about 2 days to run.
###################



###################
# To reproduce the CAV data simulation results run:

R CMD INSTALL hmm_1.1-6.tar.gz

for seed in {1..100}
do
Rscript Simulate_cav.r $seed FALSE CAV_sim/
Rscript RunFile_cav.r $seed 8 CAV_sim/
done

# The following line produces Figure 5 in the manuscript, as well as box plots for all 
# parameters.  See the supplementary material.
Rscript OutFile_boxplot_cav.r CAV_sim/ 100 FALSE


for seed in {1..100}
do
Rscript Simulate_cav.r $seed TRUE CAV_sim_population/
Rscript RunFile_cav.r $seed 8 CAV_sim_population/
done

# The following line produces Figure 6 in the manuscript, as well as box plots for all 
# parameters.  See the supplementary material.
Rscript OutFile_boxplot_cav.r CAV_sim_population/ 100 TRUE

# This code uses 8 threads in parallel.  The for loops above can each be run in an 
# embarrassingly parallel fashion.  The code for each $seed may take about 1 day to run.
###################
