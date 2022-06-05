# SamplingRespCurves
This repo contains the R code used to analyse the effect of data sampling on the response curves estimated by species distribution models.

Scripts and content:
- Wide_species: R code to gather climatic and topographic data for the analyses, and for generating _Dianthus sperandii_ (species with wide distribution). Also, the script includes the code to: 1) simulate the sampling of species occurrence data; 2) get and compare results (bias, variance and mean squared error); and create diffetent plots (maps of sampling effort, plot showing correlation between true _vs_ estimated occurrence probabilities, plot comparing true _vs_ simulated species response curves).
- Rare_species: basically, same content as in 'Wide_species', but the code is to generate _Dianthus tundrae_ (species with narrow distribution), and then get and compare results as done for _D. sperandii_.
- Prop_stratified, Env_uniform, Systematic, Proximity, Topograhic: R scripts to perform the different sampling strategies in R. 
