# LEFM with Firn Layers
If this code is used, please cite "T. Clayton, R. Duddu, T. Hageman, E. Martinez-Paneda, The influence of firn-layer material properties on surface crevasse propagation in glaciers and ice shelves. The Cryosphere (2024), [10.5194/egusphere-2024-660](https://doi.org/10.5194/egusphere-2024-660) "

Matlab code performing Linear-elastic fracture propagation through icesheets due to the presence of meltwater. Notably, the stress models used include the effect of surface firn layers, and their impact on the predicted crevasse depth is compared to that of solid ice.

# Usage
The file main.m performs all procedures needed to study the crevasse propagation depth for both the cases with firn included, and without. At the top of this file, the material parameters used throughout this analysis can be changed, after which a parametric sweep is performed over the different meltwater height ratios. 

The other files included are implementations of LEFM for the different stress profiles resulting from the inclusion of firn. These files are called via main.m, and do not need to be run by themselves.

 
