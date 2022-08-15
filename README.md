# SamplingRespCurves
This repo contains the R code used to analyse the effect of data sampling on the response curves estimated by species distribution models.

Scripts and content:
- Bio_Elev_data: script reporting code for getting and processing elevation and bioclimatic (CHELSA) data.
- Sampling_strategies: script reporting code for generating data and functions used to simulate the tested sampling approaches (random, proportional stratified, road proximity, uniform within the environmental space, systematic within the geographic space, topographic).
- Dianthus_sperandii, Dianthus_tundrae: scripts reporting code used for simulations related to the virtual specis with wide (D. sperandii) and narrow (D. tundrae) geographic distribution.
- Figures: script reporting code for reproducing some of the figures presented in the manuscript (e.g., simulated species' response curves and geographic distributions). 

The order in which the code should be run is the following: Bio_Elev_data -> Sampling_strategies -> Dianthus_sperandii -> Dianthus_tundrae -> Figures.

R version used: R version 4.1.0 (2021-05-18)

Version of the main R packages used:
- raster ‘3.5.9’
- sf ‘1.0.0’
- ggplot2 ‘3.3.5’
- elevatr ‘0.4.2’
- car ‘3.0.11’

