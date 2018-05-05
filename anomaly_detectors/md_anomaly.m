function dist_img = md_anomaly(hsi_img,mask)
%
%function dist_img = md_anomaly(hsi_img,mask)
%
% Mahalanobis Distance anomaly detector
% uses global image mean and covariance as background estimates
%
% inputs:
%  hsi_image - n_row x n_col x n_band hyperspectral image
%  mask - binary image limiting detector operation to pixels where mask is true
%         if not present or empty, no mask restrictions are used
%
% outputs:
%   dist_img - detector output image
%
% 8/7/2012 - Taylor C. Glenn 
% 5/5/2018 - Edited by Alina Zare

	if ~exist('mask','var'); mask = []; end
	addpath(fullfile('..','util'));
	dist_img = img_det(@md_det,hsi_img,[],mask);
end

function dist_data = md_det(hsi_data,~)
	n_pix = size(hsi_data,2);

	mu = mean(hsi_data,2);
	sigma = cov(hsi_data');

	z = hsi_data - repmat(mu,[1,n_pix]);
	siginv = pinv(sigma);

	dist_data = zeros(1,n_pix);
	for i=1:n_pix
	    dist_data(i) = z(:,i)'*siginv*z(:,i);
	end
end
