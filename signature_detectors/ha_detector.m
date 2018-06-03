function [ha_out] = ha_detector(hsi_img,tgt_sig,mask,ems,n_comp)
%
%function [ha_out] = ha_detector(hsi_img,tgt_sig,mask,ems,n_comp)
%
% Hybrid Abundance Detector
%
% inputs:
%  hsi_image - n_row x n_col x n_band hyperspectral image
%  tgt_sig - target signature (n_band x 1 - column vector)
%  mask - binary image limiting detector operation to pixels where mask is true
%         if not present or empty, no mask restrictions are used
%  ems - background endmembers
%  n_comp - number of mixture components for abundance mixtures
%
% outputs:
%  ha_out - detector image
%
% 8/27/2012 - Taylor C. Glenn
% 6/2/2018 - Edited by Alina Zare

addpath(fullfile('..','util'));
ha_out = img_det(@ha_det,hsi_img,tgt_sig,mask,ems,n_comp);

end

function ha_data = ha_det(hsi_data,tgt_sig,ems,n_comp)

[~,n_pix] = size(hsi_data);

% unmix data with only background
P = unmix(hsi_data,ems); 

% unmix data with target signature as well
targ_P = unmix(hsi_data,[tgt_sig ems]); 

% model abundances with and without target
gmm_bg = gmdistribution.fit(P,n_comp,'Replicates',1,'Regularize',1e-6);

% compute mixture likelihood ratio of each pixel
n_em = size(ems,2);
ll_bg = log(pdf(gmm_bg,P))'/(n_em * n_comp);

ll_tgt = targ_P(:,1)' > 0.05;
hs_data = zeros(1,n_pix);

for i=1:n_pix
    z = hsi_data(:,i) - ems*P(i,:)';
    w = hsi_data(:,i) - [tgt_sig ems]*targ_P(i,:)';
    
    hs_data(i) = (z'*z) / (w'*w);
end

ha_data = hs_data + ll_tgt - ll_bg;


end
