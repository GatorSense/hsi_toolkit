function [kbcmf_out] = kbcmf_detector(hsi_img,tgt_sig,mask,k)
%
%function [kbcmf_out] = kbcmf_detector(hsi_img,tgt_sig,mask,k)
%
% K_Betas Conditional Matched Filters
%
% inputs:
%  hsi_image - n_row x n_col x n_band hyperspectral image
%  tgt_sig - target signature (n_band x 1 - column vector)
%  mask - binary image limiting detector operation to pixels where mask is true
%         if not present or empty, no mask restrictions are used
%  k - number of clusters
%
% outputs:
%  kbcmf_out - detector image
%
% 8/24/2012 - Taylor C. Glenn - tcg@cise.ufl.edu
%

if ~exist('mask','var'); mask = []; end
if ~exist('k','var'); k = 5; end

kbcmf_out = img_det(@kbcmf_det,hsi_img,tgt_sig,mask,k);

end

function kbcmf_data = kbcmf_det(hsi_data,tgt_sig,k)

[n_band,n_pix] = size(hsi_data);

% fit the model
[idx,ll,alpha,beta] = kbetas_cluster(hsi_data,k);

% make a matched filter for each background class
mu = cell(1,k);
siginv = cell(1,k);
filt = cell(1,k);

for i=1:k
    Z = hsi_data(:,idx==i);
    mu{i} = mean(Z,2);
    siginv{i} = pinv(cov(Z'));
   
    s = tgt_sig - mu{i};
    filt{i} = s'*siginv{i} / sqrt(s'*siginv{i}*s);
end

% run appropriate filter for cluster of each pixel
kbcmf_data = zeros(1,n_pix);

for i=1:n_pix
    kbcmf_data(i) = filt{idx(i)}*(hsi_data(:,i) - mu{idx(i)});
end

end