function [gmrx_out] = gmrx_anomaly(hsi_img,mask,n_comp)
%
%[gmrx_out] = gmrx_anomaly(hsi_img,mask,n_comp)
%
% Gaussian Mixture RX Anomaly Detector
%  fits GMM assuming entire image is background
%  assigns pixels to highest posterior probability mixture component
%  computes pixel Mahlanobis distance to component mean
%
% inputs:
%  hsi_image - n_row x n_col x n_band hyperspectral image
%  mask - binary image limiting detector operation to pixels where mask is true
%         if not present or empty, no mask restrictions are used
%  n_comp - number of Gaussian components to use
%
% outputs:
%   gmrx_out - detector output image
%
% 8/7/2012 - Taylor C. Glenn - tcg@cise.ufl.edu
% 5/5/2018 - Edited by Alina Zare

	addpath(fullfile('..','util'));
	gmrx_out = img_det(@gmrx_det,hsi_img,[],mask,n_comp);
end

function gmrx_data = gmrx_det(hsi_data,~,n_comp)
	n_pix = size(hsi_data,2);

	% fit the model
	gmm = gmdistribution.fit(hsi_data',n_comp,'Replicates',1,'Regularize',1e-6);

	% cluster/assign mixture component to each pixel
	%  get mahalanobis distance to components
	[idx,nlogl,P,logpdf,M] = gmm.cluster(hsi_data');

	gmrx_data = zeros(1,n_pix);
	for i=1:n_pix
	    gmrx_data(i) = M(i,idx(i));
	end
end