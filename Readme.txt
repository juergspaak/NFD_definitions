Pyhton code to compute niche and fitness differences (N and F) as defined by
Intuitive and broadly applicable definitions of niche and fitness differences", J.W.Spaak, F. deLaender
DOI: https://doi.org/10.1101/482703 

To compute NFD for a mathematical model:

Encode the differential equations into a function `f` that returns the per capita growth rate of the species. i.e.

dN/dt = N*f(N)

To compute the parameters simply call (from numerical_NFD):

pars = NFD_model(f)

For more than two species one has to specifiy `n_spec`. The code automatically computes equilibrium densities, checks for stability and feasibility of said equilibria and computes the invasion growth rates.
Further information must be provided if automatic solver can't find stable equilibria. Examples can be found in "Example,compute NFD.py" and "Complicated examples for NFD.py".

To compute the parameters for experimental data use "NFD_for_experiments.py". An example can be found in "plot_Figure4.py".

The code is available in Python, can however be used in R aswell by using the package "reticulate" (Note that python must be installed on the computer to run reticulate). For this see "Example,compute NFD.R".

Available files are:

Pyhton code:
	numerical_NFD.py
		Contains the actual algorithm to compute niche and fitness differences.
		Only the function "NFD_model" should be used from this file, the other functions are
		called automatically by NFD_model.
		Examples on how to use NFD_model can be found in the files "Example,compute NFD.py" and "Complicated exampled for NFD.py"
	
	NFD_for_experiments.py
		Contains code to compute niche and fitness differences for an experiment, as designed in the manuscript.
		Only the function "NFD_experiment" should be used from this file, the other functions are
		called automatically by NFD_model.
		Examples on how to use NFD_experiment can be found in "plot_Figure4.py"

	Example,compute NFD.py
		Contains simple examples on how to compute NFD.
		For most models these examples are sufficient.

	Complicated examples for NFD.py
		Contains more advanced examples, that include non-monotonic per capita growth rates.
		The code is written such that it raises various errors from the NFD_model function for illustration
		To run the code without errors set "no_error = True" at the beginning of the file

	plot_Figure*
		plots the corresponding figure from the manuscript
		To plot all figures at once use plot_all_figures.py

R code:
	Example,compute NFD.R
		Shows how to use reticulate to compute NFD with R

Data files (csv)
	Exp_NFD,densities_BS*.csv
		Contains the experimentally measured densities for the cyanobacteria strains BS4 and BS5
		These are used in plot_Figure4.py

	Lights_and_absorbtions.csv
		Contains the incoming light spectrum used during the experiment and the
		absorption spectra from the two cyanobacteria.
		These are used in plot_FigureA3.py

	