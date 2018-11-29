Pyhton code to compute niche and fitness differences (N and F) as defined by
The unified Niche and Fitness definition, J.W.Spaak, F. deLaender
DOI:

To compute NFD for an mathematical model:

Encode the differential equations into a function `f` that returns the per capita growth rate of the species. i.e.

dN/dt = N*f(N)

To compute the parameters simply call (from numerical_NFD):

pars = NFD_model(f)

For more than two species one has to specifiy `n_spec`. The code automatically computes equilibrium densities, checks for stability and feasibility of said equilibria and computes the invasion growth rates. Further information must be provided if automatic solver can't find stable equilibria. Examples can be found in "Example,compute NFD.py" and "Complicated exampled for NFD.py".

To compute the parameters for experimental data use "NFD_for_experiments.py". An example can be found in "Example,NFD for experiments.py".

The code is available in Python, can however be used in R aswell by using the package "reticulate" (Note that python must be installed on the computer to run reticulate). For this see "Example,compute NFD.R".

Furthermore the files to generate the figures of the paper can be found in
plot_annual_plant_definitions.py and plots.py