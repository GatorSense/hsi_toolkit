function im_reduced = dim_reduction_demo(img, wavelengths)

%function im_reduced = dim_reduction_demo(img, wavelengths)
%
% Dimensionality reduction demo
%
% Inputs:  
%  img: hyperspectral data cube (n_row x n_col x n_bands)
%   wavelengths: vector containing wavelengths of each HSI band (n_bands x  1)
%
% Outputs:
%   im_reduced: cell array containing reduced data results
%
% Author: Alina Zare
% Email Address: azare@ufl.edu
% Created: June 3, 2018

i = 1;
addpath('hierarchical_dim_reduction');
im_reduced{i}.result = dimReduction(img, dimReductionParameters());
im_reduced{i}.method = 'hierarchical dimensionality reduction';

i = i+1;
addpath('mnf');
im_reduced{i}.result = MNF(img, MNFParameters());
im_reduced{i}.result = im_reduced{i}.result - min(min(im_reduced{i}.result(:)));
im_reduced{i}.result = im_reduced{i}.result / max(max(im_reduced{i}.result(:)));
im_reduced{i}.method = 'mnf';

%%
figure;
numR = 1;
numC = ceil((i+1)/numR); 
subplot(numR,numC,1); imagesc(get_RGB(img, wavelengths)); title('RGB image');
for i = 1:length(im_reduced)
	subplot(numR, numC, i+1); imagesc(im_reduced{i}.result(:,:,1:3)); title(im_reduced{i}.method);
end