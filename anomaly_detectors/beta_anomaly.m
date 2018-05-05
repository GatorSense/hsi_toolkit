function [beta_out] = beta_anomaly(hsi_img,mask)
%
%[beta_out] = beta_anomaly(hsi_img)
%
% Beta Distribution Anomaly Detector
%  fits beta distribution to each band assuming entire image is background
%  computes negative log likelihood of each pixel in the model
%
% inputs:
%  hsi_image - n_row x n_col x n_band hyperspectral image
%  mask - binary image limiting detector operation to pixels where mask is true
%         if not present or empty, no mask restrictions are used
%
% outputs:
%   beta_out - detector output image
%
% 8/24/2012 - Taylor C. Glenn 
% 5/5/2018 - Edited by Alina Zare

	if ~exist('mask','var'); mask = []; end

	addpath(fullfile('..','util'));
	beta_out = img_det(@beta_det,hsi_img,[],mask);
end

function beta_data = beta_det(hsi_data,~)
	[n_band,n_pt] = size(hsi_data);

	hsi_data(hsi_data <=0) = 1e-6;
	hsi_data(hsi_data >=1) = 1 - 1e-6;

	% fit the model
	alpha = zeros(n_band,1);
	beta = zeros(n_band,1);

	for i=1:n_band
	    params = betafit(hsi_data(i,:));
	    alpha(i) = params(1);
	    beta(i) = params(2);
	end

	% compute likelihood of each pixel
	like = zeros(n_band,n_pt);
	for i=1:n_band
	    like(i,:) = log(betapdf(hsi_data(i,:),alpha(i),beta(i)));
	end

	beta_data = -sum(like,1);
end
