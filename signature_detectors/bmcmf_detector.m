function [bmcmf_out] = bmcmf_detector(hsi_img,tgt_sig,mask,n_comp)
%
%function [bmcmf_out] = bmcmf_detector(hsi_img,tgt_sig,mask,n_comp)
%
% Beta Mixture Conditional Matched Filters
%
% inputs:
%  hsi_image - n_row x n_col x n_band hyperspectral image
%  tgt_sig - target signature (n_band x 1 - column vector)
%  mask - binary image limiting detector operation to pixels where mask is true
%         if not present or empty, no mask restrictions are used
%  n_comp - number of mixture components
%
% outputs:
%  bmcmf_out - detector image
%
% 8/26/2012 - Taylor C. Glenn - tcg@cise.ufl.edu
%

if ~exist('mask','var'); mask = []; end
if ~exist('n_comp','var'); n_comp = 5; end

bmcmf_out = img_det(@bmcmf_det,hsi_img,tgt_sig,mask,n_comp);

end

function bmcmf_data = bmcmf_det(hsi_data,tgt_sig,n_comp)

[n_band,n_pix] = size(hsi_data);

% fit the model
[pi,alpha,beta,idx] = bmm_fit(hsi_data,n_comp,200);

% make a matched filter for each background class
mu = cell(1,n_comp);
siginv = cell(1,n_comp);
filt = cell(1,n_comp);

for i=1:n_comp
    Z = hsi_data(:,idx==i);
    if ~isempty(Z)
        mu{i} = mean(Z,2);
        siginv{i} = pinv(cov(Z'));
        
        s = tgt_sig - mu{i};
        filt{i} = s'*siginv{i} / sqrt(s'*siginv{i}*s);
    end
end

% run appropriate filter for cluster of each pixel
bmcmf_data = zeros(1,n_pix);

for i=1:n_pix
    bmcmf_data(i) = filt{idx(i)}*(hsi_data(:,i) - mu{idx(i)});
end

end