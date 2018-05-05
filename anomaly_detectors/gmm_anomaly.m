function [gmm_out] = gmm_anomaly(hsi_img,mask,n_comp)
%
%[gmm_out] = gmm_anomaly(hsi_img,n_comp)
%
% Gaussian Mixture Model Anomaly Detector
%  fits GMM assuming entire image is background
%  computes negative log likelihood of each pixel in the fit model
%
% inputs:
%  hsi_image - n_row x n_col x n_band hyperspectral image
%  mask - binary image limiting detector operation to pixels where mask is true
%         if not present or empty, no mask restrictions are used
%  n_comp - number of Gaussian components to use
%
% outputs:
%   gmm_out - detector output image
%
% 8/7/2012 - Taylor C. Glenn
% 5/5/2018 - Edited by Alina Zare
	
	addpath(fullfile('..','util'));
	gmm_out = img_det(@gmm_det,hsi_img,[],mask,n_comp);
end

function gmm_data = gmm_det(hsi_data,~,n_comp)
	% fit the model
	gmm = gmdistribution.fit(hsi_data',n_comp,'Replicates',1,'Regularize',1e-6);

	% compute mixture likelihood of each pixel
	gmm_data = -log(pdf(gmm,hsi_data'));
end