# Sea otter effects on butter clams in Glacier Bay, Alaska

Code for analyses conducted in:

Leach, Clinton B., *et al.*, Revealing the extent of sea otter impacts on bivalve prey through multi-trophic monitoring and mechanistic models, In Review, *Journal of Animal Ecology*

## Abstract

1.  Sea otters are apex predators that can exert considerable influence over the nearshore communities they occupy. Since facing near extinction in the early 1900s, sea otters are making a remarkable recovery in Southeast Alaska, particularly in Glacier Bay, the largest protected tidewater glacier fjord in the world. The expansion of sea otters across Glacier Bay offers both a challenge to monitoring and stewardship and an unprecedented opportunity to study the top-down effect of a novel apex predator across a diverse and productive ecosystem.
2.  Our goal was to integrate monitoring data across trophic levels, space, and time to quantify and map the predator-prey interaction between sea otters and butter clams (*Saxidomus gigantea*), one of the dominant large bivalves in Glacier Bay and a favored prey of sea otters.
3.  We developed a spatially-referenced mechanistic differential equation model of butter clam dynamics that combined both environmental drivers of local population growth and estimates of otter abundance from aerial survey data. We embedded this model in a Bayesian statistical framework and fit it to clam survey data from 43 intertidal and subtidal sites across Glacier Bay.
4.  Prior to substantial sea otter expansion, we found that butter clam density was structured by an environmental gradient driven by distance from glacier (represented by latitude) and a quadratic effect of current speed. Estimates of sea otter attack rate revealed spatial heterogeneity in sea otter impacts and a negative relationship with local shoreline complexity.
5.  Sea otter exploitation of productive butter clam habitat substantially reduced the abundance and altered the distribution of butter clams across Glacier Bay, with potential cascading consequences for nearshore community structure and function. Spatial variation in estimated sea otter predation processes further suggests that community context and local environmental conditions mediate the top-down influence of sea otters on a given prey. Overall, our framework provides high-resolution insights about the interaction among components of this food web and could be applied to a variety of other systems involving invasive species, epidemiology, or migration.

## Repository

This repository contains the code to process the data and align with sea otter abundance estimates and environmental covariates, run the MCMC chain to fit the model, process the MCMC output, and generate publication figures.

All code is contained in the `src` directory. The code expects to find the following datasets in a `data` directory:

-   `data/prey`: Intertidal and subtidal sea otter prey sampling data from <https://doi.org/10.5066/P9LODH0Z>
-   `data/otters`: output from sea otter ecological diffusion model (Lu *et al.* 2019) from ...
-   `data/current`: raster of root-mean-squared estimates of current speed from tidal circulation (Drew *et al*. 2013) from ...
-   `data/PH6502`: shapefile of Glacier Bay shoreline, downloaded from the NOAA shoreline data explorer at <https://nsde.ngs.noaa.gov/> (Project ID PH6502)

Files in the `src` directory are as follows:

-   `align_data.r`: generates all of the required inputs for the model, including the clam abundance observations, the sea otter abundance time series for each prey sampling site, and the environmental covariates for each site (latitude, current speed, and shoreline complexity)

-   `compute_lambda.r`: Computes the posterior mean sea otter abundance raster from the Lu *et al.* 2019 model output; requires `glb_specs.r` and `lambda_helpers.r`

-   `plots.r`: Includes all of the code to generate figures and tables from the manuscript, as well as additional exploratory figures; expects an `output` directory from which to read in MCMC output produced by `run.jl`, and an `output/figures` directory to write figures to

-   `model.jl`: defines structs to contain all of the model constants and parameters

-   `process.jl`: defines the process model ODE

-   `run.jl`: sets up and runs the MCMC

-   `sample.jl`: defines all of the machinery for the MCMC including the log-likelihood function, functions to sample each of the parameters, and function to execute the full MCMC
