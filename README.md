# M_-M_K-

This code is meant to provide realistic masses and uncertainties of M dwarfs given a user-provided 2MASS K and distance (and uncertainties). The code will then read in the included posterior from Mann et al. (in prep) to estimate the error arising from scatter in the relation itself. The output is a posterior on stellar mass (Solar units). You can also request 1D (mean and sigma) errors if requested.

You can include [Fe/H] if it is known. In this scenario the code will use a Mk-Mass-[Fe/H] relation. Note that the [Fe/H] term is technically not statistically significant. So this is mostly for illustrative purposes or tests on stars with extreme metallicities (e.g., subdwarfs). 

The code requires scalar values for K, distance, and corresponding uncertainties. If you instead of a set of values (e.g., an asymmetric or otherwise non-Gaussian posterior on distance) the suggested solution is to read in the posterior and run the code on each point separately, and then combine the resulting posteriors on mass.

The success of the associated paper depends on the ease of use of this software, so please feel free to send questions/suggestions to mann.andrew.w [at] gmail.com (or open an issue). 

The code is available in IDL, but I promise to make a python version soon (under some sufficiently vague definition of soon). 
 
