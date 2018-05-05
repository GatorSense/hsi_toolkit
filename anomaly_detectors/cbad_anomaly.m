function [cbad_out,cluster_img] = cbad_anomaly(hsi_img,mask,n_cluster)
%
% cbad_out = cbad_anomaly(hsi_img,mask,n_cluster)
%
% Cluster Based Anomaly Detection (CBAD)
% Ref: Carlotto, Mark J. "A cluster-based approach for detecting man-made objects and changes in imagery." IEEE Transactions on Geoscience and Remote Sensing 43.2 (2005): 374-387.
%
% inputs:
%  hsi_image - n_row x n_col x n_band hyperspectral image
%  mask - binary image limiting detector operation to pixels where mask is true
%         if not present or empty, no mask restrictions are used
%  n_cluster - number of clusters to use
%
% outputs:
%  cbad_out - detector output image
%  cluster_img - cluster label image
%
% 8/7/2012 - Taylor C. Glenn 
% 5/5/2018 - Edited by Alina Zare

addpath(fullfile('..','util'))
[cbad_out,cluster_img] = img_det(@cbad_det,hsi_img,[],mask,n_cluster);

end

function [cbad_data,idx] = cbad_det(hsi_data,~,n_cluster)
	n_pix = size(hsi_data,2);

	% cluster the data
	[idx,C] = kmeans(hsi_data',n_cluster,'emptyaction','singleton');

	% get cluster statistics
	mu = cell(1,n_cluster);
	siginv = cell(1,n_cluster);
	for i=1:n_cluster
	    
	    z = hsi_data(:,idx == i);
	    
	    mu{i} = mean(z,2);
	    
	    siginv{i} = pinv(cov(z'));
	    
	end

	% compute Mahalanobis distance of each point to its cluster
	cbad_data = zeros(1,n_pix);

	for i=1:n_pix
	    z = hsi_data(:,i) - mu{idx(i)};
	    
	    cbad_data(i) = z'*siginv{idx(i)}*z;
	end
end

