Pyhton code to compute niche and fitness differences (N and F) as defined by
Intuitive and broadly applicable definitions of niche and fitness differences", J.W.Spaak, F. deLaender
DOI: https://doi.org/10.1101/482703 
and by "Building modern coexistence theory from the ground up: the role of community assembly"
By J.W. Spaak and S. J. Schreiber

R code Author: Sebastian J. Schreiber
Python code Author: Jurg W. Spaak

This repository combines automated code to compute the invasion graphs defined in Spaak and Schreiber as well as code to compute niche and fitness differneces as defined in Spaak and De Laender. 

################################################################################################################
To compute niche and fitness differences according to Spaak and De Laender:

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

##################################################################################################################
To compute teh invasion graph and determine permanence:

We provide automated code for step 1 for Lotka-Volterra community models, steps 2-4 for any model type, and steps 5-6 for niche and fitness differences for Lotka-Volterra community models (\ref{ap:NFD}).
We provide two sets of codes, one for any community model and one for Lotka-Volterra community models only.

For any community model, given the invasion scheme (invasion growth rates for all possible sub-communities), we compute the invasion graph, check for cycles and find an assembly end state.
Additionally, we provide a visualisation of the invasion graph as well as the assembly end state.
The relevant functions for this are found in scheme_plot.py for Python and xxx

For a Lotka-Volterra community model, given the intrinsic growth rates \mu and the interaction matrix, we compute the invasion scheme, invasion graph, check for cycles and find an assembly end state.
Additionally, we compute the niche and fitness differences of all species (presnt and absent in the assembly end state).
The relevant functions for this are found in scheme_plot.py for Python and xxx

###########################################################################
invasion_graph_main_functions.R contains functions for producing an invasion scheme for a Lotka-Volterra community model, producing the invasion graph from an invasion scheme, assessing properties of the invasion graph (e.g. acyclic or not, permanent communities, -i communities), and plotting the invasion graph. These functions were orginally developed for the article "Permanence via invasion graphs: Incorporating community assembly into Modern Coexistence Theory" by Josef Hofbauer and Sebastian J. Schreiber in the Journal of Mathematical Biology (2022) and made available at 10.5281/zenodo.7111753

invasion_graph_aux_functions.R contains functions for an alternative plot of the invasion graph, finding assembly end states and the associated n-1 communities, and finding cycles in a cyclic graph. 

invasion_graph_example.R illustrates using the functions from  invasion_graph_main_functions.R and invasion_graph_aux_functions.R. The illustration corresponds to the communities from Figures 3C and D in the manuscript. 

non-equilibrium-figure.R produces the plot in Box 1 of the manuscript. 

###########################################################################
Python code:
scheme_plot.py
	Contains automated code to compute invasion graphs and invasion schemes
	The most important functions are "do_all" and "do_all_LV", which perform the above mentioned code.
	For examples see plot_NFD_example.py

############################################################################
Files used for the figures in the paper
plot_base_schemes.py
	Plots figure 2 from the manuscript

plot_NFD_example.py
	Plots figure 3 from the manuscript
	Uses Vandermeer_1969.csv and Geijzendorffer 2011.csv as input

plot_decomposition_example.py
	Plots figure 4 from the manucript

letten_2018_functions.py
	Contains model functions for Letten et al. 2018 paper
	These functions are used in "plot_decomposition_example.py"

Data files
Letten_2018.csv
	Contains the species traits needed for letten_2018_functions.py
	Data is from Letten et al. 2018, PNAS

Vandermeer_1969.csv
	Contains the species interaction matrix for the Vandermeer example.
	Columns correspond to the effect of species i on species j

Geijzendorffer 2011.csv
	Contains the species interaction matrix for the Geijzendorffer example.
	Columns correspond to the effect of species i on species j



	