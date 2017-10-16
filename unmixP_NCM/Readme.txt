Sampled NCM (Sampled Normal Compositional Model) ReadMe

------------------------------------------------------------------------------------

NOTE: If the Sampled NCM Algorithm is used in any publication or presentation, the following reference must be cited: A. Zare, P. Gader, and G. Casella, “Sampling piecewise convex unmixing and endmember extraction,” IEEE Trans. Geosci. Remote Sens., vol. 51, iss. 3, pp. 1655-1665, 2013. 

------------------------------------------------------------------------------------

The command to run the Sampled NCM code:

[P] = unmixP_NCM(X,E,cov,Parameters)

The input data is a dxN matrix of N input data points each of d spectral bands (dimensionality).

The parameters input is a struct with the following fields:

Parameters.NumberIterations = 1000; % number of iterations.

The parameters structure can be generated using the unmixP_NCM_Parameters.m function.

--------------------------------------------------------------------------------------

Real Hyperspectral Data Demo:

To Run the Sampled NCM Algorithm, with the example real hyperspectral data set run the scipt -

demo_unmixP_NCM.m

This script will estimate proportions for the simulated hyperspectral data set and plot the results.

--------------------------------------------------------------------------------------

If you have any questions, please contact:

Alina Zare
Electrical and Computer Engineering
University of Missouri
327 Engineering Building West
Columbia, MO 65211
zarea@missouri.edu