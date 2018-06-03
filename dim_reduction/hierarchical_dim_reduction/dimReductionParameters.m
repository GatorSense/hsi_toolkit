function [Parameters] = dimReductionParameters()

%function [Parameters] = dimReductionParameters()
%
% Parameters for Hierarchical Dimensionality Reduction
% Computes KL-divergence between every pair of bands in image and merges
% bands that are similar based on KL-divergence.
%
% Author: Alina Zare
% Email Address: azare@ufl.edu
% Created: September 12, 2008
% Latest Revision: June 3, 2018
%%

%Parameters for Hierarchical Dimensionality Reduction
Parameters.numBands = 7; %Reduced dimensionality size
Parameters.type = 'complete'; %Type of Hierarchy used in clustering
Parameters.showH = 0; %Set to 1 to show clustering, 0 otherwise
Parameters.NumCenters = 255; %Number of centers used in computing KL-divergence
