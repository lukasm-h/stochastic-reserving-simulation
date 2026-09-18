# stochastic-reserving-simulation
R code for the simulation and comparison of actuarial reserving practices for different model specifications

This repository contains the R code used for the simulation study in my master's thesis on the comparison of stochastic reserving methods for different (correct/incorrect) model specifications. 

Some parts of the simulation and scripts use functions from the comprehensive ChainLadder R package (Gesmann et al., 2026). 

## Simulation design

We focus on the Mack-, ODP- and Lognormal-model of stochastic reserving. 

The aim of the study is to get information on the risk of a misspecified model for the stochastic reserving process. The underlying question is: What is the expected difference of results, if the true model is X but the actuary (wrongly) specifies model Y?
The 'results' include the estimated ultimate/reserve, MSEP and multiple reserve quantiles. 

To achieve that, the following (replicable) process is followed: 

0. Specify a parameter set (including parameters for Mack, ODP and Lognormal model) that is somewhat 'common'/'prototypical' for your branch of insurance. You can do that by calibrating on a 'prototypical' triangle. All results are conditional on this parameter set. 

1. Simulate for each generating model (Mack, ODP, Lognormal) 10,000 triangles, using the parameters specified in (0) and following exactly the distributional properties of the models.

2. For each triangle, fit each fitting model (Mack, ODP, Lognormal) and record the corresponding estimated ultimate/reserve, MSEP and reserve quantiles (using analytical and Bootstrap methods).

3. For each triangle, compare the estimated quantities when fitting the 'correct' underlying (generating) model and when fitting the other two 'misspecified' models. 

4. Average or plot the observed differences in (3) for each combination of generating and fitting model. 

## Repository structure

- 'Simulation_Routine.R'
	
	Main simulation script. It is intended to be executed section by section, to give the opportunity to observe the process and check the diagnostics. Given a parameter set, the triangles for each generating model are simulated, for each triangle the fitting models are fitted and the differences of results are obtained and exported as CSV.  

- 'Parameter_calibration.R'

	Contains the calibration function. For a given 'prototypical' triangle, all the needed parameters for the simulation are obtained and saved. 

- 'Functions.R'

	All the analytical and Bootstrap-based functions used in the "Simulation_Routine.R"
	This is the mathematical core of the simuation.

- 'Export_Statistics.R'

	Function for exporting the simulation results to CSV. 

- 'Taylor_Ashe_1983_parameters.rds'

	Parameters calibrated on the triangle from Taylor and Ashe (1983, 59-60), used as the standard parameter set for the simulation, enabling testing of the script without supplying own parameter set or calibration triangle. 

## Usage

Important prerequisite: The simulation is based on a 'prototypical' parameter set supplied by the user. It needs to include (modulo model) equivalent parameters for Mack-, ODP- and Lognormal-model. For easy usage, the function in 'Parameter_calibration.R' can be used to obtain the exactly right format of parameters for a given 'prototypical' triangle, supplied in an Excel file. For testing purposes, the parameters calibrated on the triangle from Taylor and Ashe (1983, pp. 59-60) are included in the repository and set as the standard parameters in the main 'Simulation_Routine.R'.

The main process is contained in 'Simulation_Routine.R'. It is intended to be executed section by section, to give the opportunity to observe the process and check the diagnostics. Please check the comments and (small) manual input required in that script. 
Given a parameter set, the triangles for each generating model are simulated, for each triangle the fitting models are fitted and the differences of results are obtained and exported as CSV.

Please note: The simulation can be computationally intensive, so you can adjust the number of generated triangles and Bootstrap samples used for the fitting process in section 0 of that script. 

## R environment / required packages
Used packages:

- ChainLadder
- readxl
- cli

The dependencies are saved using 'renv', which can be used to restore the environment as well. 

## References

Gesmann, M., Murphy, D., Zhang, Y., Carrato, A., Wüthrich, M., Concina, F., & Dal Moro, E. 2026. *ChainLadder: Statistical Methods and Models for Claims Reserving in General Insurance*. R package version 0.2.21. https://doi.org/10.32614/CRAN.package.ChainLadder

Taylor, G.C. & Ashe, F.R. 1983. Second moments of estimates of outstanding claims. *Journal of Econometrics 23*(1), 37–61. https://doi.org/10.1016/0304-4076(83)90074-X

For a general introduction to the relevant models and methods, see e.g.

Wüthrich, Mario V. & Merz, Michael. 2008. *Stochastic Claims Reserving Methods in
Insurance*. John Wiley & Sons Ltd.
