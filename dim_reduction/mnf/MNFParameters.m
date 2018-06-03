function params = MNFParameters()

%function params = MNFParameters()
%
% Parameters for Maximum Noise Fraction code
%
% Author: Alina Zare
% Email Address: azare@ufl.edu
% Created: 2008
% Latest Revision: June 3, 2018
%%

params = struct();
params.dim_reduce = true;  % logical flag on whether to do a dimensionality reduction on final PCA step
params.en_pct = 0.98;       % if dim_reduce, percentage of eigenvalues to retain