# Markov-process-analysis
## About the project
---
This is one of my graduate research projects. I derived analytical expressions of Markov process time series using linear algebra to construct complex enzyme kinetic models, and implemented the expressions using MATLAB scripts to analyze experimental data, which leads to new insights in enzymatic mechanisms. The derivation part is detailed in a manuscript, which is about to be submitted, and will be updated here once it's published. At this point, I only include the MATLAB code that implements the analytical expressions to fit experimental data.\
The experimental data were taken from Hendry, G. and Wydrzynski, T. Biochemistry 2003, 42, 6209.,
Hillier, W. and Wydrzynski, T. Phys. Chem. Chem. Phys. 2004, 6, 4882.,
and de Lichtenberg, C. and Messinger, J. Phys. Chem. Chem. Phys. 2020, 22, 12894.


## Getting started
---
I used [MATLAB](https://www.mathworks.com/products/get-matlab.html?s_tid=gn_getml) for this project. Download the whole repository including the `.csv` files containing experimental data and run the `Kinetic_modeling.m` file. Several input values can be changed according to the need of the researcher, including `k1`, `k2`, `btrp`, `save`, and the experimental data for the fit.
